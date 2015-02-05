from matplotlib import pyplot as plt
import numpy as np


N_exp = 4  # 8 boxes
N = 2**N_exp


def make_example_plot(ell):
    t_size = 2**ell
    s_size = 2**(N_exp - ell)
    all_s = np.linspace(0, N, s_size + 1)
    all_t = np.linspace(0, N, t_size + 1)

    # Make all horizontal lines.
    for t in all_t:
        plt.plot(all_s, t * np.ones(s_size + 1),
                 color='b', marker='o')

    for s in all_s:
        plt.plot(s * np.ones(t_size + 1), all_t,
                 color='b', marker='o')

    plt.title(r'$\ell = %d$' % (ell,))
    ax = plt.gca()
    ax.set_frame_on(False)
    ax.set_xticklabels([])
    # H/T: http://stackoverflow.com/a/20416681/1068170
    for tic in ax.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    ax.set_yticklabels([])
    for tic in ax.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    plt.xlabel('$S$')
    plt.ylabel('$T$', rotation=0)
    plt.axis('scaled')
    plt.xlim(-1, N + 1)
    plt.ylim(-1, N + 1)
    plt.show()
