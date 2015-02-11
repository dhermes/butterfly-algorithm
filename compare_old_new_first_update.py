import numpy as np
from scipy.misc import factorial

from butterfly_algorithm import update_coefficients
from compare_old_new_initial_data import get_coefficients_lists
from new_butterfly_algorithm import A1
from new_butterfly_algorithm import A2
from new_butterfly_algorithm import make_update_func
from new_butterfly_algorithm import increase_tau_refinement


def main():
    coefficients_list, coeff_vals, the_rest = get_coefficients_lists()
    L, M, s, t, sigma_vals, tau_endpoints = the_rest
    delta_S = (np.max(s) - np.min(s)) / 2**(L + 1)
    delta_T = (np.max(t) - np.min(t)) * 0.25

    num_tau, num_sigma = 1, 2**L

    A1_minus = A1(M, -delta_T)
    A1_plus = A1(M, delta_T)
    A2_minus = A2(M, -delta_S)
    A2_plus = A2(M, delta_S)

    update_func = make_update_func(A1_minus, A1_plus, A2_minus,
                                   A2_plus, delta_T)
    coeff_vals = increase_tau_refinement(coeff_vals, num_tau,
                                         num_sigma, update_func)

    factorial_values = factorial(range(M))
    coefficients_list, sigma_vals, tau_endpoints = update_coefficients(
        coefficients_list, sigma_vals, tau_endpoints, factorial_values)

    print '\n',
    print 'M:', M
    print 'tau old %g vs. %g:' % (np.mean(tau_endpoints[0]), coeff_vals[0][0])
    print 'sigma old: %g vs. %g' % (sigma_vals[0], coeff_vals[0][1])
    print 'coeffs:', np.array(coefficients_list[0])
    print 'coeffs new:', coeff_vals[0][2].flatten()


if __name__ == '__main__':
    main()
