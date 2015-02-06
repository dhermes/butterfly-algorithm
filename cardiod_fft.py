from matplotlib import pyplot as plt
import numpy as np


a = 0.5
num_points = 2**5
t = np.linspace(0, 2 * np.pi, num_points + 1)[:-1]
x = a * (2 * np.cos(t) - np.cos(2 * t))
y = a * (2 * np.sin(t) - np.sin(2 * t))

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
ax1.plot(
    np.hstack([x, x[0]]),
    np.hstack([y, y[0]]),
    marker='o',
)
ax1.axis('scaled')

f = x + 1j * y
f_hat = np.fft.fft(f, n=num_points)
x_hat = np.real(f_hat)
y_hat = np.imag(f_hat)

ax2.plot(
    np.hstack([x_hat, x_hat[0]]),
    np.hstack([y_hat, y_hat[0]]),
    marker='o',
    color='r',
)
plt.show()

# Remove far left and far right
far_left = np.argmin(x_hat)
far_right = np.argmax(x_hat)

# Drop the far left and far right.
plt.plot(
    np.delete(x_hat, (far_left, far_right)),
    np.delete(y_hat, (far_left, far_right)),
    marker='o',
    linestyle='None',
)
plt.show()
