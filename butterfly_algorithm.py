from matplotlib import pyplot as plt
import numpy as np


def compute_f_hat(f, t, s, kernel_func):
    f_hat = np.zeros(f.shape, dtype=np.complex128)
    for k in xrange(len(f)):
        # Vectorized update.
        f_hat[k] = np.sum(kernel_func(t[k], s) * f)

    return f_hat


def dft_kernel(t, s):
    return np.exp(- 1.0j * t * s)


def dft_data(N):
    N_vals = np.arange(N, dtype=np.float64)
    t = 2 * np.pi * N_vals
    s = N_vals / N
    return t, s


def solve_problem_slides():
    # http://www.mathworks.com/help/matlab/math/fast-fourier-transform-fft.html
    data = np.load('resources/bluewhale.npz')
    X = data['X']
    sampling_rate = int(data['rate'])

    blue_whale_begin = 24500 - 1
    blue_whale_end = 31000 - 1

    blue_whale_call = X[blue_whale_begin:blue_whale_end + 1]
    size_call = len(blue_whale_call)

    N = int(2**np.ceil(np.log2(size_call)))

    blue_whale_call = np.hstack([
        blue_whale_call,
        np.zeros(N - len(blue_whale_call)),
    ])

    dft_whale_call = np.fft.fft(blue_whale_call, n=N)


class CoefficientsOwner(object):
    """General information about a set target and source data.

    Assumes:
    - s, t are real
    - s is in ascending order
    """

    def __init__(self, L, s, t):
        self.L = L
        self.s = s
        # Don't need t for the coefficients, just need tau.

        self.max_num_bins = 2**self.L
        # NOTE: We could find the maximum and minimum in a single
        #       pass through t (instead of separate calls).
        self.t_endpoints = tuple(
            np.linspace(np.min(t), np.max(t), self.max_num_bins + 1).tolist())
        self.s_endpoints = tuple(
            np.linspace(np.min(s), np.max(s), self.max_num_bins + 1).tolist())

        self.s_values_by_bin = self._s_values_by_bin()

    def _tau_val(self, level, tau_index):
        stride = 2**(self.L - level)
        left_val = self.t_endpoints[stride * tau_index]
        right_val = self.t_endpoints[stride * (tau_index + 1)]
        return 0.5 * (left_val + right_val)

    def _sigma_val(self, level, sigma_index):
        stride = 2**level
        left_val = self.s_endpoints[stride * sigma_index]
        right_val = self.s_endpoints[stride * (sigma_index + 1)]
        return 0.5 * (left_val + right_val)

    def _s_values_by_bin(self):
        result = {}
        curr_bin_index = 0
        right_val = self.s_endpoints[curr_bin_index + 1]

        # Assumes s is sorted.
        curr_values = []
        for s_value in self.s:
            if s_value <= right_val:
                curr_values.append(s_value)
            else:
                result[curr_bin_index] = tuple(curr_values)
                # Increment until s_value is contained.
                while s_value > right_val:
                    curr_bin_index += 1
                    right_val = self.s_endpoints[curr_bin_index + 1]
                curr_values = [s_value]

        # Add the final set of values.
        result[curr_bin_index] = tuple(curr_values)

        return result


class DataCoefficient(object):

    def __init__(self, owner, level, tau_index, sigma_index):
        self.owner = owner
        self.level = level

        self.tau_index = tau_index
        self.tau = self.owner._tau_val(self.level, self.tau_index)

        self.sigma_index = sigma_index
        self.sigma = self.owner._sigma_val(self.level, self.sigma_index)


def make_intervals(t, s, L=None):
    if L is None:
        L = int(np.floor(np.log2(len(s))))

    owner = CoefficientsOwner(L, s, t)
    level = tau_index = 0
    return [
        DataCoefficient(owner, level, tau_index, sigma_index)
        for sigma_index in xrange(2**L)
    ]



def test_make_intervals(L=4):
    t, s = dft_data(2**L)
    intervals = make_intervals(t, s, L=L)
    return intervals
