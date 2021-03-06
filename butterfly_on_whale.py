#!/usr/bin/env python

import argparse
import numpy as np
import time


from butterfly_algorithm import load_whale
from butterfly_algorithm import solve


def dft_data(N):
    N_vals = np.arange(N, dtype=np.float64)
    t = 2 * np.pi * N_vals
    s = N_vals / N
    return t, s


def main(M=11, L=None):
    _, data = load_whale()
    N = len(data)
    if L is None:
        L = int(np.floor(np.log2(N)))
    print 'N = %d points, M = %d truncated terms, L = %d refinements' % (
        N, M, L)
    t, s = dft_data(N)

    start = time.time()
    f_hat = solve(t, s, data, L, M=M)
    duration = time.time() - start
    print('total computation time: %g' % (duration,))

    fft_f_hat = np.fft.fft(data, n=N)
    print('2-norm: %g' % (np.linalg.norm(f_hat - fft_f_hat, ord=2),))
    print('sup-norm: %g' % (np.linalg.norm(f_hat - fft_f_hat, ord=np.inf),))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run Butterfly on whale test data.')
    parser.add_argument('--M', dest='M', type=int, default=11,
                        help='Size of Taylor series truncation.')
    parser.add_argument('--L', dest='L', type=int,
                        help='Number of grid refinement levels.')
    args = parser.parse_args()
    main(M=args.M, L=args.L)
