from matplotlib import pyplot as plt
import numpy as np
from scipy.misc import factorial


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
    result = {}
    curr_bin_index = 0

    max_s = s[-1]  # Due to ordering.
    min_s = s[0]  # Due to ordering.
    bin_width = (max_s - min_s) / float(num_bins)
    right_val = bin_width

    curr_values = []
    for pair in zip(s, data):
        # pair[0] == s_value
        if pair[0] <= right_val:
            curr_values.append(pair)
        else:
            result[curr_bin_index] = tuple(curr_values)
            # Increment until s_value is contained.
            while pair[0] > right_val:
                curr_bin_index += 1
                right_val += bin_width
            curr_values = [pair]

    # Add the final set of values.
    result[curr_bin_index] = tuple(curr_values)

    bin_endpoints = np.linspace(min_s, max_s, num_bins + 1)
    return bin_endpoints, result


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


def initial_coefficients(s, data, t, R=8):
    """Computes the set of D(tau, sigma, alpha) coefficients when ell = 0.

    Returns two values:
    - A list of lists of the initial coefficients for tau the only center at
      level ell = 0 and all possible sigma values. Uses (R) as the truncation
      of the Taylor series
    - A dictionary of the s values in each bin; we no longer need s -> D(s)
      mapping after using the D(s) values to compute D(tau, sigma, alpha)
    """
    L = int(np.floor(np.log2(len(s))))
    max_num_bins = 2**L

    bin_endpoints, sigma_bins = create_sigma_bins(s, data, max_num_bins)
    # NOTE: We could find the maximum and minimum in a single
    #       pass through t (instead of separate calls).
    min_t = np.min(t)
    max_t = np.max(t)
    tau = 0.5 * (min_t + max_t)

    coefficients_list = []
    s_bin_values = {}
    for index in xrange(max_num_bins):
        left, right = bin_endpoints[index:index + 2]
        sigma = 0.5 * (left + right)
        bin_pairs = sigma_bins.get(index, ())
        s_bin_values[index] = (left, right,
                               tuple(pair[0] for pair in bin_pairs))
        coefficients = [bin_coefficient(bin_pairs, tau, sigma, alpha)
                        for alpha in xrange(R)]
        coefficients_list.append(coefficients)

    return coefficients_list, s_bin_values
