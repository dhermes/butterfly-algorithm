from butterfly_algorithm import make_intervals


def test_make_intervals(L=4):
    t, s = dft_data(2**L)
    intervals = make_intervals(t, s, L=L)
    return intervals
