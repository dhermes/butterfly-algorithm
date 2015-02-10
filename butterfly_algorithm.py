from matplotlib import pyplot as plt
import numpy as np
from scipy.misc import factorial
import time

from _fortran_utils import speedup
intermediate_coeffs = speedup.intermediate_coeffs
coeff_new_level = speedup.coeff_new_level


DEBUG = True


if DEBUG:
    def info(msg):
        print(msg)
else:
    def info(msg):
        pass


def dft_kernel(t, s):
    return np.exp(- 1.0j * t * s)


def dft_data(N):
    N_vals = np.arange(N, dtype=np.float64)
    t = 2 * np.pi * N_vals
    s = N_vals / N
    return t, s


def create_sigma_bins(s, data, num_bins):
    """Puts source points from s into bins.

    Also puts the corresponding data to s in the bin so that
    the s -> D(s) mapping is preserved to use when creating
    a coefficient.

    Also returns a list of the endpoints of the bins.

    Assumes s is sorted.
    """
    # Due to ordering, max and min at beginning.
    bin_endpoints = np.linspace(s[0], s[-1], num_bins + 1)
    # This corresponds to iteration -1 of the loop.
    left_val, right_val = None, bin_endpoints[0]

    pairs = zip(s, data)  # Assumes same length.
    num_pairs = len(pairs)

    result = []
    sigma_vals = []
    s_index = 0
    # NOTE: Consider using numpy.digitize.
    for curr_bin_index in xrange(num_bins):
        left_val, right_val = right_val, bin_endpoints[curr_bin_index + 1]
        sigma_vals.append(0.5 * (left_val + right_val))
        curr_values = []
        while s_index < num_pairs and pairs[s_index][0] <= right_val:
            curr_values.append(pairs[s_index])
            s_index += 1
        result.append(tuple(curr_values))

    return sigma_vals, result


def bin_coefficient(bin_pairs, tau, sigma, alpha):
    """Computes the coefficient for a given bin.

    Assumes bin_pairs is a tuple of pairs (s, D(s)) and that these
    values represent all s in B(sigma).

    Computes D(tau, sigma, alpha) under these assumptions.
    """
    result = 0.0
    # NOTE: We could speed this up by vectorizing but when we call it
    #       `bin_pairs` will typically have just 1 value.
    for s_val, data_val in bin_pairs:
        result += dft_kernel(tau, s_val) * (s_val - sigma)**alpha * data_val
    return (- 1.0j)**alpha / factorial(alpha) * result


def initial_coefficients(s, data, t, L, R=8):
    """Computes the set of D(tau, sigma, alpha) coefficients when ell = 0.

    Returns two values:
    - A list of lists of the initial coefficients for tau the only center at
      level ell = 0 and all possible sigma values. Uses (R) as the truncation
      of the Taylor series
    - A dictionary of the s values in each bin; we no longer need s -> D(s)
      mapping after using the D(s) values to compute D(tau, sigma, alpha)
    """
    max_num_bins = 2**L

    sigma_vals, sigma_bins = create_sigma_bins(s, data, max_num_bins)
    # NOTE: We could find the maximum and minimum in a single
    #       pass through t (instead of separate calls).
    min_t = np.min(t)
    max_t = np.max(t)
    tau = 0.5 * (min_t + max_t)

    coefficients_list = []
    for index in xrange(max_num_bins):
        sigma = sigma_vals[index]
        bin_pairs = sigma_bins[index]
        coefficients = tuple(bin_coefficient(bin_pairs, tau, sigma, alpha)
                             for alpha in xrange(R))
        coefficients_list.append(coefficients)

    tau_endpoints = [(min_t, max_t)]
    return coefficients_list, sigma_vals, tau_endpoints


def coarsen_sigma_vals(sigma_vals):
    # Assumes N is even.
    N = len(sigma_vals)
    halfN, _ = divmod(N, 2)
    result = []
    for i in xrange(halfN):
        left = sigma_vals[2 * i]
        right = sigma_vals[2 * i + 1]
        result.append(0.5 * (left + right))
    return result


def refine_tau_endpoints(tau_endpoints):
    result = []
    for left, right in tau_endpoints:
        mid = 0.5 * (left + right)
        result.extend([
            (left, mid),
            (mid, right),
        ])
    return result


def update_coefficients(coefficients_list, sigma_vals, tau_endpoints,
                        factorial_values):
    new_sigma_vals = coarsen_sigma_vals(sigma_vals)
    new_tau_endpoints = refine_tau_endpoints(tau_endpoints)

    num_sigma = len(sigma_vals)
    num_new_sigma = len(new_sigma_vals)
    R = len(coefficients_list[0])

    # We order the values of `coefficients_list` according to the
    # tau index and then the sigma index.
    # [ (tau(0), sigma(0)), (tau(0), sigma(1)), ... (tau(1), sigma(0)), ...]

    new_coefficients_list = []
    for tau_plus_index, tau_pair in enumerate(new_tau_endpoints):
        tau_plus = 0.5 * (tau_pair[0] + tau_pair[1])
        tau_index, _ = divmod(tau_plus_index, 2)
        tau = 0.5 * (tau_endpoints[tau_index][0] + tau_endpoints[tau_index][1])

        for sigma_minus_index, sigma_minus in enumerate(new_sigma_vals):
            # We want the floor of (INDEX / 2).
            sigma_index = 2 * sigma_minus_index
            sigma_prime_index = 2 * sigma_minus_index + 1

            sigma = sigma_vals[sigma_index]
            sigma_prime = sigma_vals[sigma_prime_index]

            D_tau_sigma = coefficients_list[tau_index * num_sigma + sigma_index]
            D_tau_sigma_prime = coefficients_list[
                tau_index * num_sigma + sigma_prime_index]

            new_coeffs = tuple(
                coeff_new_level(D_tau_sigma, D_tau_sigma_prime, tau, tau_plus,
                                sigma, sigma_prime, sigma_minus, alpha,
                                factorial_values, R)
                for alpha in xrange(R)
            )
            new_coefficients_list.append(new_coeffs)

    return new_coefficients_list, new_sigma_vals, new_tau_endpoints


def match_with_tau(t_vals, tau_endpoints):
    sorted_indices = np.argsort(t_vals)
    num_bins = len(tau_endpoints)
    num_t = len(t_vals)

    result = {}
    t_index = 0
    # NOTE: Consider using numpy.digitize.
    for curr_bin_index in xrange(num_bins):
        left_val, right_val = tau_endpoints[curr_bin_index]
        tau = 0.5 * (left_val + right_val)
        curr_values = []
        while t_index < num_t:
            sorted_index = sorted_indices[t_index]
            curr_val = t_vals[sorted_index]
            if curr_val > right_val:
                break
            curr_values.append((sorted_index, curr_val))
            result[sorted_index] = (tau, curr_bin_index)
            t_index += 1

    return result


def make_f_hat(t, tau_map, sigma, R, coefficients_list):
    f_hat = []
    for t_index, t_val in enumerate(t):
        tau, tau_index = tau_map[t_index]
        # K_prime = dft_kernel(t_val - tau, sigma) * (t_val - tau)**alpha
        f_hat.append(dft_kernel(t_val - tau, sigma) * np.dot(
            (t_val - tau)**np.arange(R),
            coefficients_list[tau_index],
        ))
    return f_hat


def approximate_f_hat(t, s, data, R=8):
    L = int(np.floor(np.log2(len(s))))

    start = time.time()
    coefficients_list, sigma_vals, tau_endpoints = initial_coefficients(
        s, data, t, L, R=R)
    duration = time.time() - start
    info('initial_coefficients time: %g' % (duration,))

    factorial_values = factorial(range(R))
    for _ in xrange(L):
        start = time.time()
        coefficients_list, sigma_vals, tau_endpoints = update_coefficients(
            coefficients_list, sigma_vals, tau_endpoints, factorial_values)
        duration = time.time() - start
        info('loop %d time: %g' % (_, duration))

    # There should be exactly one left.
    sigma, = sigma_vals

    start = time.time()
    tau_map = match_with_tau(t, tau_endpoints)
    duration = time.time() - start
    info('match_with_tau time: %g' % (duration,))

    start = time.time()
    f_hat = make_f_hat(t, tau_map, sigma, R, coefficients_list)
    duration = time.time() - start
    info('make_f_hat time: %g' % (duration,))
    return np.array(f_hat)


def simple_correctness_test(L=5, R=11):
    N = 2**L
    t, s = dft_data(N)
    data = np.random.random(t.shape)
    f_hat = approximate_f_hat(t, s, data, R=R)
    fft_f_hat = np.fft.fft(data, n=N)
    print('2-norm: %g' % (np.linalg.norm(f_hat - fft_f_hat, ord=2),))
    print('sup-norm: %g' % (np.linalg.norm(f_hat - fft_f_hat, ord=np.inf),))
