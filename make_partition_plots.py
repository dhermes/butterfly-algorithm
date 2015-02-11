from matplotlib import pyplot as plt
import numpy as np


N_exp = 4  # 8 boxes
N = 2**N_exp


def remove_axis_frame(ax):
    ax.set_frame_on(False)
    ax.yaxis.set_label_position('right')
    ax.set_xticklabels([])
    # H/T: http://stackoverflow.com/a/20416681/1068170
    for tic in ax.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    ax.set_yticklabels([])
    for tic in ax.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False


def make_1D_ax_just_s(ax, s_size, N_exp=N_exp):
    N = 2**N_exp
    all_s = np.linspace(0, N, s_size + 1)
    ax.plot(all_s, 0.5 * N * np.ones(all_s.shape),
            color='b', marker='o')
    remove_axis_frame(ax)
    ax.axis('scaled')
    ax.set_xlim(-1, N + 1)
    ax.set_ylim(-1, N + 1)


def make_1D_ax_just_t(ax, t_size, N_exp=N_exp):
    N = 2**N_exp
    all_t = np.linspace(0, N, t_size + 1)
    ax.plot(0.5 * N * np.ones(all_t.shape), all_t,
            color='b', marker='o')
    remove_axis_frame(ax)
    ax.axis('scaled')
    ax.set_xlim(-1, N + 1)
    ax.set_ylim(-1, N + 1)


def make_1D_ax_both(ax, ell, t_size, s_size, N_exp=N_exp):
    ax.set_title(r'$L = %d, \, \ell = %d$' % (N_exp, ell))

    N = 2**N_exp
    all_s, all_t = np.linspace(0, N, s_size + 1), np.linspace(0, N, t_size + 1)
    # Make all horizontal lines.
    for t in all_t:
        ax.plot(all_s, t * np.ones(s_size + 1),
                color='b', marker='o')
    # Make all vertical lines.
    for s in all_s:
        ax.plot(s * np.ones(t_size + 1), all_t,
                color='b', marker='o')

    remove_axis_frame(ax)
    ax.set_xlabel('$S$')
    ax.set_ylabel('$T$', rotation=0)
    ax.axis('scaled')
    ax.set_xlim(-1, N + 1)
    ax.set_ylim(-1, N + 1)


def make_1D_example_plot(ell, N_exp=N_exp):
    rows, cols = 1, 3
    fig, (ax1, ax2, ax3) = plt.subplots(rows, cols)

    t_size, s_size = 2**ell, 2**(N_exp - ell)
    make_1D_ax_just_s(ax1, s_size, N_exp=N_exp)
    make_1D_ax_both(ax2, ell, t_size, s_size, N_exp=N_exp)
    make_1D_ax_just_t(ax3, t_size, N_exp=N_exp)

    width, height = fig.get_size_inches()
    fig.set_size_inches(2 * width, 2 * height, forward=True)

    plt.show()


def make_2D_ax(ax, ell, N_exp=N_exp, title=None):
    if title is not None:
        ax.set_title(title)

    N = 2**N_exp
    size = 2**ell

    all_pts = np.linspace(0, N, size + 1)
    for val in all_pts:
        # Make horizontal line.
        ax.plot(all_pts, val * np.ones(size + 1),
                color='b', marker='o')
        # Make vertical line.
        ax.plot(val * np.ones(size + 1), all_pts,
                color='b', marker='o')

    remove_axis_frame(ax)
    ax.axis('scaled')
    ax.set_xlim(-1, N + 1)
    ax.set_ylim(-1, N + 1)


def make_2D_arrow(ax, N_exp=N_exp):
    N = 2**N_exp
    ax.annotate('', (0, 0.5 * N), (N, 0.5 * N),
                arrowprops={'arrowstyle':'<->'})
    # ax.arrow(
    remove_axis_frame(ax)
    ax.axis('scaled')
    ax.set_xlim(-1, N + 1)
    ax.set_ylim(-1, N + 1)


def make_2D_example_plot(ell, N_exp=N_exp):
    rows, cols = 1, 3
    fig, (ax1, ax2, ax3) = plt.subplots(rows, cols)

    make_2D_ax(ax1, N_exp - ell, N_exp=N_exp, title='$S$')
    make_2D_arrow(ax2, N_exp=N_exp)
    make_2D_ax(ax3, ell, N_exp=N_exp, title='$T$')

    width, height = fig.get_size_inches()
    fig.set_size_inches(2 * width, 2 * height, forward=True)

    plt.show()


def make_1D_centers(L=N_exp):
    rows, cols = 1, 1
    fig, ax = plt.subplots(rows, cols)

    N = 2**L

    all_s = np.linspace(0, N, N + 1)
    ax.plot(all_s, np.zeros(all_s.shape),
            color='black', marker='|', markersize=20)

    center_pts = all_s[:-1] + 0.5
    ax.plot(center_pts, np.zeros(center_pts.shape),
            color='red', marker='o', linestyle='None')

    remove_axis_frame(ax)
    ax.axis('scaled')
    ax.set_xlim(-1, N + 1)
    ax.set_ylim(-0.5, 0.5)

    width, height = fig.get_size_inches()
    fig.set_size_inches(2 * width, 2 * height, forward=True)

    plt.show()


def make_interval_at_stages(center_label, left_label, right_label):
    rows, cols = 1, 1
    fig, ax = plt.subplots(rows, cols)

    ax.plot([-1, 1], [0, 0], color='black', marker='|', markersize=20)
    ax.plot([0], [0], color='blue', marker='o', linestyle='None')
    ax.plot([-0.5, 0.5], [0, 0], color='red', marker='o', linestyle='None')

    remove_axis_frame(ax)

    ax.axis('scaled')
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.1, 0.1)
    ax.annotate(center_label, xy=(0.0, 0.0), xytext=(-0.025, -0.125),
                fontsize=24)
    ax.annotate(left_label, xy=(0.0, 0.0), xytext=(-0.5 - 0.025, -0.125),
                fontsize=24)
    ax.annotate(right_label, xy=(0.0, 0.0), xytext=(0.5 - 0.025, -0.125),
                fontsize=24)

    width, height = fig.get_size_inches()
    fig.set_size_inches(2 * width, 2 * height, forward=True)

    plt.show()


def make_err_plots(L10_M_choices, L11_M_choices, L12_M_choices, L13_M_choices,
                   L10_time_pairs, L11_time_pairs, L12_time_pairs,
                   L13_time_pairs):
    rows = cols = 2
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(rows, cols)

    L10_log_errors = [np.log10(pair[0]) for pair in L10_time_pairs]
    L11_log_errors = [np.log10(pair[0]) for pair in L11_time_pairs]
    L12_log_errors = [np.log10(pair[0]) for pair in L12_time_pairs]
    L13_log_errors = [np.log10(pair[0]) for pair in L13_time_pairs]
    true_min = np.min(L10_log_errors + L11_log_errors +
                      L12_log_errors + L13_log_errors)
    true_min -= 1.0
    true_max = np.max(L10_log_errors + L11_log_errors +
                      L12_log_errors + L13_log_errors)
    true_max += 1.0
    true_left = np.min(L10_M_choices + L11_M_choices +
                       L12_M_choices + L13_M_choices)
    true_left -= 3
    true_right = np.max(L10_M_choices + L11_M_choices +
                        L12_M_choices + L13_M_choices)
    true_right += 3

    ax1.plot(L10_M_choices, L10_log_errors, marker='o')
    ax1.set_ylabel(r'$\log_{10} ||e||_2$', rotation=0, fontsize=20, labelpad=40)
    ax1.set_title(r'$L = 10$', fontsize=20)
    ax1.set_xlim((true_left, true_right))
    ax1.set_ylim((true_min, true_max))

    ax2.plot(L11_M_choices, L11_log_errors, marker='o')
    ax2.set_title(r'$L = 11$', fontsize=20)
    ax2.set_xlim((true_left, true_right))
    ax2.set_ylim((true_min, true_max))

    ax3.plot(L12_M_choices, L12_log_errors, marker='o')
    ax3.set_ylabel(r'$\log_{10} ||e||_2$', rotation=0, fontsize=20, labelpad=40)
    ax3.set_title(r'$L = 12$', fontsize=20)
    ax3.set_xlim((true_left, true_right))
    ax3.set_ylim((true_min, true_max))

    ax4.plot(L13_M_choices, L13_log_errors, marker='o')
    ax4.set_title(r'$L = 13$', fontsize=20)
    ax4.set_xlim((true_left, true_right))
    ax4.set_ylim((true_min, true_max))

    width, height = fig.get_size_inches()
    fig.set_size_inches(2 * width, 2 * height, forward=True)

    plt.show()


def make_time_plots(L10_M_choices, L11_M_choices, L12_M_choices, L13_M_choices,
                    L10_time_pairs, L11_time_pairs, L12_time_pairs,
                    L13_time_pairs):
    rows = cols = 2
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(rows, cols)

    L10_log_errors = [np.log10(pair[0]) for pair in L10_time_pairs]
    L11_log_errors = [np.log10(pair[0]) for pair in L11_time_pairs]
    L12_log_errors = [np.log10(pair[0]) for pair in L12_time_pairs]
    L13_log_errors = [np.log10(pair[0]) for pair in L13_time_pairs]

    L10_time_only = [pair[1] for pair in L10_time_pairs]
    L11_time_only = [pair[1] for pair in L11_time_pairs]
    L12_time_only = [pair[1] for pair in L12_time_pairs]
    L13_time_only = [pair[1] for pair in L13_time_pairs]

    true_min = np.min(L10_log_errors + L11_log_errors +
                      L12_log_errors + L13_log_errors)
    true_min -= 1.0
    true_max = np.max(L10_log_errors + L11_log_errors +
                      L12_log_errors + L13_log_errors)
    true_max += 1.0
    true_left = np.min(L10_time_only + L11_time_only +
                       L12_time_only + L13_time_only)
    true_left -= 0.1
    true_right = np.max(L10_time_only + L11_time_only +
                        L12_time_only + L13_time_only)
    true_right += 0.1

    ax1.plot(L10_time_only, L10_log_errors, marker='o')
    ax1.set_ylabel(r'$\log_{10} ||e||_2$', rotation=0, fontsize=20, labelpad=40)
    ax1.set_title(r'$L = 10$', fontsize=20)
    ax1.set_xlim((true_left, true_right))
    ax1.set_ylim((true_min, true_max))

    ax2.plot(L11_time_only, L11_log_errors, marker='o')
    ax2.set_title(r'$L = 11$', fontsize=20)
    ax2.set_xlim((true_left, true_right))
    ax2.set_ylim((true_min, true_max))

    ax3.plot(L12_time_only, L12_log_errors, marker='o')
    ax3.set_ylabel(r'$\log_{10} ||e||_2$', rotation=0, fontsize=20, labelpad=40)
    ax3.set_xlabel('runtime', fontsize=20)
    ax3.set_title(r'$L = 12$', fontsize=20)
    ax3.set_xlim((true_left, true_right))
    ax3.set_ylim((true_min, true_max))

    ax4.plot(L13_time_only, L13_log_errors, marker='o')
    ax4.set_title(r'$L = 13$', fontsize=20)
    ax4.set_xlabel('runtime', fontsize=20)
    ax4.set_xlim((true_left, true_right))
    ax4.set_ylim((true_min, true_max))

    width, height = fig.get_size_inches()
    fig.set_size_inches(2 * width, 2 * height, forward=True)

    plt.show()