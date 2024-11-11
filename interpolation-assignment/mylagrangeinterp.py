import numpy as np
import matplotlib.pyplot as plt
import time

points = np.loadtxt('points.txt')

def lagrange_interpolation(x, y, z):
    n = len(x)
    p = 0.0
    for i in range(n):
        L = 1.0
        for j in range(n):
            if j != i:
                L *= (z - x[j]) / (x[i] - x[j])
        p += y[i] * L
    return p

x_points = [points[i][0] for i in range(len(points))]
y_points = [points[i][1] for i in range(len(points))]

z_values = np.linspace(-6, 6, 100)
start = time.time()
w_values = [lagrange_interpolation(x_points, y_points, z) for z in z_values]
end = time.time()

print('Interpolation took', end-start)

fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))

ax1.plot(x_points, y_points, 'o', label="Data Points")
ax1.plot(z_values, w_values, 'r-', label="Lagrange Polynomial")
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$p(x)$")
ax1.legend()
ax1.set_title("Lagrange Interpolation Polynomial")
plt.tight_layout()
plt.show()
