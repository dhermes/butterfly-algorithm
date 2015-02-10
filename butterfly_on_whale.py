#!/usr/bin/env python

import argparse
import numpy as np
import time


from butterfly_algorithm import approximate_f_hat
from butterfly_algorithm import dft_data


def get_whale():
    data = np.load('resources/bluewhale.npz')
    X = data['X']

    blue_whale_begin = 24500 - 1
    blue_whale_end = 31000 - 1

    blue_whale_call = X[blue_whale_begin:blue_whale_end + 1]
    size_call = len(blue_whale_call)

    N = int(2**np.ceil(np.log2(size_call)))

    blue_whale_call = np.hstack([
        blue_whale_call,
        np.zeros(N - len(blue_whale_call)),
    ])

    return blue_whale_call, N


def main(R=11):
    """Expected performance:

    R= 8 -- 1.9 seconds per loop; 2-norm error = O(10^(-1))
    R= 9 -- 2.3 seconds per loop; 2-norm error = O(10^(-2))
    R=10 -- 2.8 seconds per loop; 2-norm error = O(10^(-3))
    R=11 -- 3.3 seconds per loop; 2-norm error = O(10^(-4))
    R=12 -- 3.9 seconds per loop; 2-norm error = O(10^(-5))
    R=13 -- 4.5 seconds per loop; 2-norm error = O(10^(-6))
    R=14 -- 5.3 seconds per loop; 2-norm error = O(10^(-7))
    R=15 -- 6.0 seconds per loop; 2-norm error = O(10^(-8))
    """
    data, N = get_whale()
    print 'N = %d points, R = %d truncated terms' % (N, R)
    t, s = dft_data(N)

    start = time.time()
    f_hat = approximate_f_hat(t, s, data, R=R)
    duration = time.time() - start
    print('total computation time: %g' % (duration,))

    fft_f_hat = np.fft.fft(data, n=N)
    print('2-norm: %g' % (np.linalg.norm(f_hat - fft_f_hat, ord=2),))
    print('sup-norm: %g' % (np.linalg.norm(f_hat - fft_f_hat, ord=np.inf),))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run Butterfly on whale test data.')
    parser.add_argument('--R', dest='R', type=int, default=11,
                        help='Size of Taylor series truncation.')
    args = parser.parse_args()
    main(R=args.R)
