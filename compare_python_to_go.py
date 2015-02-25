import numpy as np
from butterfly_algorithm import create_initial_data


s = np.array([0.0, 1.0])
min_s, max_s = 0.0, 1.0
tau = 1.0
actual_data = np.array([3.0, 4.0], dtype=np.complex128)
num_bins = 4
M = 2
sum_parts = create_initial_data(s, min_s, max_s, tau,
                                actual_data, num_bins, M)
for row in sum_parts:
    print row[:2], row[2].flatten()
