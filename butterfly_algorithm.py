from matplotlib import pyplot as plt
import numpy as np
from scipy.misc import comb as combinations
from scipy.misc import factorial
import time


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
    left_val = 0
    right_val = bin_endpoints[1]

    pairs = zip(s, data)  # Assumes same length.
    num_pairs = len(pairs)

    result = []
    sigma_vals = []
    s_index = 0
    for curr_bin_index in xrange(num_bins):
        sigma_vals.append(0.5 * (left_val + right_val))
        curr_values = []
        while s_index < num_pairs and pairs[s_index][0] <= right_val:
            curr_values.append(pairs[s_index])
            s_index += 1
        result.append(tuple(curr_values))
        if curr_bin_index < num_bins - 1:
            left_val, right_val = right_val, bin_endpoints[curr_bin_index + 2]

    return sigma_vals, result


def bin_coefficient(bin_pairs, tau, sigma, alpha):
    """Computes the coefficient for a given bin.

    Assumes bin_pairs is a tuple of pairs (s, D(s)) and that these
    values represent all s in B(sigma).

    Computes D(tau, sigma, alpha) under these assumptions.
    """
    result = 0.0
    for s_val, data_val in bin_pairs:
        result += dft_kernel(tau, s_val) * (s_val - sigma)**alpha * data_val
    return (- 1.0j)**alpha / factorial(alpha) * result


def initial_coefficients(s, data, t, R=8, L=None):
    """Computes the set of D(tau, sigma, alpha) coefficients when ell = 0.

    Returns two values:
    - A list of lists of the initial coefficients for tau the only center at
      level ell = 0 and all possible sigma values. Uses (R) as the truncation
      of the Taylor series
    - A dictionary of the s values in each bin; we no longer need s -> D(s)
      mapping after using the D(s) values to compute D(tau, sigma, alpha)
    """
    if L is None:
        L = int(np.floor(np.log2(len(s))))
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
    N = len(sigma_vals)
    if N == 1:
        raise ValueError('Cannot coarsen a singleton.')
    halfN, _ = divmod(N, 2)  # Assumes N is even
    result = []
    for i in xrange(halfN):
        left = sigma_vals[2 * i]
        right = sigma_vals[2 * i + 1]
        result.append(0.5 * (left + right))
    return result


def refine_tau_endpoints(tau_endpoints):
    N = len(tau_endpoints)
    result = []
    for i in xrange(N):
        left, right = tau_endpoints[i]
        mid = 0.5 * (left + right)
        result.extend([
            (left, mid),
            (mid, right),
        ])
    return result


def intermediate_coeffs(D_tau_sigma, D_tau_sigma_prime, tau, tau_plus,
                        sigma, sigma_prime, beta, R):
    sub_result = sub_result_prime = 0.0
    for gamma in xrange(beta, R):
        comb_val = (combinations(gamma, beta) *
                    (tau_plus - tau)**(gamma - beta))
        sub_result += comb_val * D_tau_sigma[gamma]
        sub_result_prime += comb_val * D_tau_sigma_prime[gamma]

    sub_result *= dft_kernel(tau_plus - tau, sigma)
    sub_result_prime *= dft_kernel(tau_plus - tau, sigma_prime)

    return sub_result, sub_result_prime


def coeff_new_level(D_tau_sigma, D_tau_sigma_prime, tau, tau_plus,
                    sigma, sigma_prime, sigma_minus, alpha, R):
    result = 0.0
    for beta in xrange(alpha + 1):
        sub_result, sub_result_prime = intermediate_coeffs(
            D_tau_sigma, D_tau_sigma_prime, tau, tau_plus,
            sigma, sigma_prime, beta, R)

        partial = (sigma - sigma_minus)**(alpha - beta) * sub_result
        partial_prime = ((sigma_prime - sigma_minus)**(alpha - beta) *
                         sub_result_prime)
        result += (
            (- 1.0j)**(alpha - beta) *
            (partial + partial_prime) / factorial(alpha - beta)
        )
    return result


def update_coefficients(coefficients_list, sigma_vals, tau_endpoints):
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
                                sigma, sigma_prime, sigma_minus, alpha, R)
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


def approximate_f_hat(t, s, data, R=8):
    L = int(np.floor(np.log2(len(s))))
    coefficients_list, sigma_vals, tau_endpoints = initial_coefficients(
        s, data, t, R=R, L=L)
    for _ in xrange(L):
        coefficients_list, sigma_vals, tau_endpoints = update_coefficients(
            coefficients_list, sigma_vals, tau_endpoints)

    sigma, = sigma_vals

    tau_map = match_with_tau(t, tau_endpoints)
    f_hat = []
    for t_index, t_val in enumerate(t):
        tau, tau_index = tau_map[t_index]
        # K_prime = dft_kernel(t_val - tau, sigma) * (t_val - tau)**alpha
        f_hat.append(dft_kernel(t_val - tau, sigma) * np.dot(
            (t_val - tau)**np.arange(R),
            coefficients_list[tau_index],
        ))
    return np.array(f_hat)


def simple_correctness_test(L=5, R=11):
    N = 2**L
    t, s = dft_data(N)
    data = np.random.random(t.shape)
    start = time.time()
    f_hat = approximate_f_hat(t, s, data, R=R)
    duration = time.time() - start
    print 'duration:', duration
    fft_f_hat = np.fft.fft(data, n=N)
    print '2-norm:', np.linalg.norm(f_hat - fft_f_hat, ord=2)
    print 'sup-norm:', np.linalg.norm(f_hat - fft_f_hat, ord=np.inf)
