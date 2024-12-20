import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def rotate(v, w, theta):
    w = w / np.linalg.norm(w)
    return (
        v * np.cos(theta)
        + np.cross(w, v) * np.sin(theta)
        + w * np.dot(w, v) * (1 - np.cos(theta))
    )

a = np.array([2, 3, 4])
b = np.array([2, 0, 4])

rotated_vector = rotate(a, b, np.pi)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.quiver(0, 0, 0, a[0], a[1], a[2], color='r', label='Vector a', arrow_length_ratio=0.1)

ax.quiver(0, 0, 0, b[0], b[1], b[2], color='g', label='Axis b', arrow_length_ratio=0.1)

ax.quiver(0, 0, 0, rotated_vector[0], rotated_vector[1], rotated_vector[2], color='b', label='Rotated a', arrow_length_ratio=0.1)

max_range = np.array([a, b, rotated_vector]).max()
ax.set_xlim([-max_range, max_range])
ax.set_ylim([-max_range, max_range])
ax.set_zlim([-max_range, max_range])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

plt.show()

