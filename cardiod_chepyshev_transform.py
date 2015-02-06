from matplotlib import pyplot as plt
import numpy as np

# See http://en.wikipedia.org/wiki/Discrete_Chebyshev_transform


num_points = 2**5
x = np.linspace(0, 1, num_points)
u = 4 * x**3 - 7 * x + 1

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
ax1.plot(x, u, marker='o')
ax1.axis('scaled')

u_hat = np.zeros(u.shape)
for m in range(len(u)):
    pm = 1.0 if m == 0 else 2.0
    u_hat[m] = pm * u.dot(np.cos(m * np.arccos(x))) / num_points

ax2.plot(x, u_hat, marker='o', color='r')
plt.show()
