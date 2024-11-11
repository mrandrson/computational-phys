import numpy as np
import matplotlib.pyplot as plt
import time

def divided_difference(x, y):
    n = len(x)
    F = np.zeros((n, n))
    for i in range(n):
        F[i][0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            F[i][j] = (F[i + 1][j - 1] - F[i][j - 1]) / (x[i + j] - x[i])
    return F[0]

def newton_poly(x, y, z):
    a = divided_difference(x, y)
    n = len(a)
    p = a[0]
    for i in range(1, n):
        term = a[i]
        for j in range(i):
            term *= (z - x[j])
        p += term
    return p

points = np.loadtxt('points.txt')

x_points = [points[i][0] for i in range(len(points))]
y_points = [points[i][1] for i in range(len(points))]

zvals = np.linspace(-6, 6, 100)

start = time.time()
w = [newton_poly(x_points, y_points, z) for z in zvals]
end = time.time()

print('Interpolation took', end-start)

fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))

ax1.plot(x_points, y_points, 'o', label="Data Points")
ax1.plot(zvals, w, 'r-', label="Newton Polynomial")
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$p(x)$")
ax1.legend()
ax1.set_title("Newton Interpolation Polynomial")
plt.tight_layout()
plt.show()
