from itertools import izip
import numpy as np

from butterfly_algorithm import dft_data
from butterfly_algorithm import initial_coefficients
from new_butterfly_algorithm import create_initial_data
from new_butterfly_algorithm import load_whale


N, data = load_whale()
t, s = dft_data(N)
L = int(np.floor(np.log2(N)))
M = 5

coefficients_list, sigma_vals, tau_endpoints = initial_coefficients(
    s, data, t, L, M=M)
tau = np.mean(tau_endpoints[0])
coeff_vals = create_initial_data(s, tau, data, L, M)

if len(coeff_vals) != N or len(coefficients_list) != N:
    raise ValueError('Wrong Size')


for OLD, NEW in izip(coefficients_list, coeff_vals):
    OLD = np.array(OLD)
    NEW = NEW[2].flatten()
    if not np.allclose(OLD, NEW):
        raise ValueError((OLD, NEW))
