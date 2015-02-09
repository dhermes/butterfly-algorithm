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

    def __init__(self, level, L, s, t):
        self.level = level
        self.L = L
        self.s = s
        self.t = t

        self.max_t = np.max(t)
        self.min_t = np.min(t)
        self.max_s = np.max(s)
        self.min_s = np.min(s)

        self.num_tau = 2**self.level
        self.num_sigma = 2**(self.L - self.level)

        self.tau_width = (self.max_t - self.min_t) / float(self.num_tau)
        self.sigma_width = (self.max_s - self.min_s) / float(self.num_sigma)

        self._tau_values = {}
        self._sigma_values = {}
        self._s_values_by_bin = {}

    def _tau_val(self, tau_index):
        if tau_index not in self._tau_values:
            interval_begin = self.min_t + tau_index * self.tau_width
            self._tau_values[tau_index] = interval_begin + 0.5 * self.tau_width
        return self._tau_values[tau_index]

    def _sigma_val(self, sigma_index):
        if sigma_index not in self._sigma_values:
            interval_begin = self.min_s + sigma_index * self.sigma_width
            self._sigma_values[sigma_index] = (interval_begin +
                                               0.5 * self.sigma_width)
        return self._sigma_values[sigma_index]

    def _get_s_values(self, sigma_index):
        if sigma_index not in self._s_values_by_bin:
            interval_begin = self.min_s + sigma_index * self.sigma_width
            if sigma_index == self.num_sigma - 1:
                interval_end = self.max_s
            else:
                interval_end = interval_begin + self.sigma_width

            s_values = [value for value in self.s
                        if interval_begin < value <= interval_end]

            if sigma_index == 0 and self.min_s in self.s:
                s_values.insert(0, self.min_s)

            self._s_values_by_bin[sigma_index] = s_values

        return self._s_values_by_bin[sigma_index]


class DataCoefficient(object):

    def __init__(self, owner, tau_index, sigma_index):
        self.owner = owner

        self.tau_index = tau_index
        self.tau = self.owner._tau_val(self.tau_index)

        self.sigma_index = sigma_index
        self.sigma = self.owner._sigma_val(self.sigma_index)

        self.s_values = self.owner._get_s_values(self.sigma_index)


def make_intervals(t, s, L=None):
    if L is None:
        L = int(np.floor(np.log2(len(s))))

    owner = CoefficientsOwner(0, L, s, t)
    return [
        DataCoefficient(owner, 0, sigma_index)
        for sigma_index in xrange(2**L)
    ]



def test_make_intervals(L=4):
    t, s = dft_data(2**L)
    intervals = make_intervals(t, s, L=L)
    return intervals
