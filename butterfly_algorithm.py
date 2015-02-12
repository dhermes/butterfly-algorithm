import numpy as np
from itertools import izip


def compute_f_hat(f, t, s, kernel_func):
    f_hat = np.zeros(f.shape, dtype=np.complex128)
    for k in xrange(len(f)):
        # Vectorized update.
        f_hat[k] = np.sum(kernel_func(t[k], s) * f)

    return f_hat


def dft_kernel(t, s):
    return np.exp(- 1.0j * t * s)


def bin_s(s_values, L):
    min_s = np.min(s_values)
    max_s = np.max(s_values)

    num_bins = 2**L
    s_bins = np.linspace(min_s, max_s, num_bins + 1)

    bin_width = s_bins[1] - s_bins[0]
    curr_sigma = 0.5 * (s_bins[1] + s_bins[0])

    # Cheat so that max_s is in the last bin; since digitize checks
    # BINS[i - 1] <= VAL < BINS[i].
    s_bins[-1] += 0.1

    # Bin info is a list of pairs. Each pair will contain two key
    # bin information pieces: (sigma, (list of s in bin)).
    bin_info = [[curr_sigma, {}]]
    for _ in xrange(num_bins - 1):
        curr_sigma += bin_width
        bin_info.append([curr_sigma, {}])

    for i, bin_index in enumerate(np.digitize(s_values, s_bins)):
        bin_info[bin_index - 1][1][i] = s_values[i]

    return bin_info


def create_initial_data(s_values, tau, actual_data, L, M):
    # Compute D(tau, sigma, alpha)
    bin_info = bin_s(s_values, L)
    result = []
    for sigma, value_dict in bin_info:
        # alpha = 0 case
        s_with_initial = [(s, dft_kernel(tau, s) * actual_data[i])
                          for i, s in value_dict.items()]
        # We need to make sure s is in the same order as K(tau, s) D(s),
        # so we keep the pair since dict.items() is non-deterministic.
        s_data = np.array([pair[0] for pair in s_with_initial])
        s_data = (s_data - sigma) * (-1.0j)
        values = np.array([pair[1] for pair in s_with_initial])

        alpha_vals = [np.sum(values)]
        for alpha in xrange(1, M):
            values = values * s_data / alpha
            alpha_vals.append(np.sum(values))

        # Keep tau, sigma and the values.
        result.append((tau, sigma, np.array(alpha_vals).reshape(M, 1)))

    return result


def compute_t_by_bins(t_values, coeff_vals, L):
    min_t = np.min(t_values)
    max_t = np.max(t_values)
    M = len(coeff_vals[0][-1])

    num_bins = 2**L
    t_bins = np.linspace(min_t, max_t, num_bins + 1)
    # Row vector exponents; will be multiplied by the column
    # vector alpha values.
    exponents = np.arange(M, dtype=np.float64).reshape(1, M)

    # Cheat so that max_t is in the last bin; since digitize checks
    # BINS[i - 1] <= VAL < BINS[i].
    t_bins[-1] += 0.1

    result = []
    for t_val, bin_index in izip(t_values, np.digitize(t_values, t_bins)):
        tau, sigma, alpha_vals = coeff_vals[bin_index - 1]
        delta = t_val - tau
        value_dot_prod = (delta**exponents).dot(alpha_vals)
        # The result is a 1x1 array.
        result.append(dft_kernel(delta, sigma) * value_dot_prod[0, 0])

    return np.array(result)


def A1(M, delta, eye_func=np.eye):
    result = eye_func(M)
    result[0, 1] = delta

    for col in xrange(2, M):
        prev_val = result[0, col - 1]

        # Pascal's triangle does not apply at ends. The zero term
        # already set on the diagonal.
        result[0, col] = delta * prev_val
        for row in xrange(1, col):
            curr_val = result[row, col - 1]
            result[row, col] = prev_val + delta * curr_val
            prev_val = curr_val

    return result


def mult_diag(val, A, M, diag):
    first_row = max(0, -diag)
    last_row = min(M, M - diag)
    for row in xrange(first_row, last_row):
        A[row, row + diag] *= val


def A_update(A_val, scale_multiplier=0.5, upper_diags=True):
    M = A_val.shape[0]

    # If not `upper_diags` we want the lower diagonal.
    diag_mult = 1 if upper_diags else -1
    # We don't need to update the main diagonal since exponent=0.
    scale_factor = 1
    for diagonal in xrange(1, M):
        scale_factor *= scale_multiplier
        mult_diag(scale_factor, A_val, M, diagonal * diag_mult)


def complex_eye(M):
    return np.eye(M, dtype=np.complex128)


def set_diag(val, A, M, diag):
    first_row = max(0, -diag)
    last_row = min(M, M - diag)
    for row in xrange(first_row, last_row):
        A[row, row + diag] = val


def A2(M, delta, eye_func=complex_eye, imag=1.0j):
    result = eye_func(M)
    new_delta = -imag * delta

    diagonal_value = 1
    for sub_diagonal in xrange(1, M):
        diagonal_value = diagonal_value * new_delta / sub_diagonal
        set_diag(diagonal_value, result, M, -sub_diagonal)

    return result


def increase_tau_refinement(values, num_tau, num_sigma, update_func):
    result = []
    for tau_index in xrange(num_tau):
        begin = tau_index * num_sigma
        end = begin + num_sigma
        # We need to hold the right values until all the left values
        # have been added.
        right_updated = []

        # Assumes num_sigma is even.
        left_vals, right_vals = values[begin:end:2], values[1 + begin:end:2]
        for left_val, right_val in izip(left_vals, right_vals):
            new_left, new_right = update_func(left_val, right_val)
            result.append(new_left)
            right_updated.append(new_right)

        result.extend(right_updated)

    return result


def make_update_func(A1_minus, A1_plus, A2_minus, A2_plus, delta_T):

    def update_func(left_val, right_val):
        tau, sigma, alpha_vals_left = left_val
        tau_same, sigma_prime, alpha_vals_right = right_val
        # We expect the pair to share tau, and don't check to avoid slowdown.
        # if tau_same != tau:
        #     import ipdb
        #     ipdb.set_trace()

        sigma_minus = 0.5 * (sigma + sigma_prime)
        tau_left = tau - delta_T
        tau_right = tau + delta_T

        k_minus_sigma = dft_kernel(-delta_T, sigma)
        k_plus_sigma = 1.0 / k_minus_sigma
        k_minus_sigma_prime = dft_kernel(-delta_T, sigma_prime)
        k_plus_sigma_prime = 1.0 / k_minus_sigma_prime

        new_alpha_vals_left = (
            k_minus_sigma * A2_minus.dot(A1_minus.dot(alpha_vals_left)) +
            k_minus_sigma_prime * A2_plus.dot(A1_minus.dot(alpha_vals_right))
        )
        new_alpha_vals_right = (
            k_plus_sigma * A2_minus.dot(A1_plus.dot(alpha_vals_left)) +
            k_plus_sigma_prime * A2_plus.dot(A1_plus.dot(alpha_vals_right))
        )

        new_left_val = (tau_left, sigma_minus, new_alpha_vals_left)
        new_right_val = (tau_right, sigma_minus, new_alpha_vals_right)
        return new_left_val, new_right_val

    return update_func


def solve(s, t, data, L, M):
    min_t = np.min(t)
    max_t = np.max(t)
    tau = 0.5 * (min_t + max_t)
    coeff_vals = create_initial_data(s, tau, data, L, M)

    # ell = 0
    delta_S = (np.max(s) - np.min(s)) / 2**(L + 1)
    delta_T = (max_t - min_t) * 0.25
    num_tau, num_sigma = 1, 2**L

    A1_minus = A1(M, -delta_T)
    A1_plus = A1(M, delta_T)
    A2_minus = A2(M, -delta_S)
    A2_plus = A2(M, delta_S)

    for ell in xrange(1, L + 1):
        update_func = make_update_func(A1_minus, A1_plus, A2_minus,
                                       A2_plus, delta_T)
        coeff_vals = increase_tau_refinement(coeff_vals, num_tau,
                                             num_sigma, update_func)

        num_tau = num_tau * 2
        num_sigma = num_sigma / 2
        delta_T *= 0.5

        A_update(A1_plus, scale_multiplier=0.5, upper_diags=True)
        A_update(A1_minus, scale_multiplier=0.5, upper_diags=True)
        A_update(A2_plus, scale_multiplier=2.0, upper_diags=False)
        A_update(A2_minus, scale_multiplier=2.0, upper_diags=False)

    return compute_t_by_bins(t, coeff_vals, L)
