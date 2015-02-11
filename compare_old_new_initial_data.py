from itertools import izip
import numpy as np

from butterfly_algorithm import dft_data
from butterfly_algorithm import initial_coefficients
from new_butterfly_algorithm import create_initial_data
from new_butterfly_algorithm import load_whale


def get_coefficients_lists():
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

    the_rest = (L, M, s, t, sigma_vals, tau_endpoints)
    return coefficients_list, coeff_vals, the_rest


def main():
    coefficients_list, coeff_vals, _ = get_coefficients_lists()

    for OLD, NEW in izip(coefficients_list, coeff_vals):
        OLD = np.array(OLD)
        NEW = NEW[2].flatten()
        if not np.allclose(OLD, NEW):
            raise ValueError((OLD, NEW))


if __name__ == '__main__':
    main()
