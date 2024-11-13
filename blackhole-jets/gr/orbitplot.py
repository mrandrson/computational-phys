import numpy as np
import matplotlib.pyplot as plt

data = np.fromfile("orbital_plane_positions.bin", dtype=np.float64).reshape(-1, 2)

x, y = data[:, 0], data[:, 1]

print("Data shape:", data.shape)

plt.plot(x, y)
plt.xlabel("x position")
plt.ylabel("y position")
plt.title("Orbital Plane Trajectory")
plt.show()

