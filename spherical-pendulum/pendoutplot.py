import struct
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

with open('output.bin', 'rb') as f:
    data = f.read()

data_size = len(data)
steps = data_size // (3 * 8)

x, y, z = [], [], []

for i in range(steps):
    offset = i * 3 * 8
    x_val, y_val, z_val = struct.unpack('ddd', data[offset:offset + 24])
    x.append(x_val)
    y.append(y_val)
    z.append(z_val)

plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

line, = ax.plot([], [], [], linewidth=1, label='Connecting Line')
bob, = ax.plot([], [], [], 'ro', label='Bob')

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    bob.set_data([], [])
    bob.set_3d_properties([])
    return line, bob

def update(frame):
    line.set_data([0, x[frame]], [0, y[frame]])
    line.set_3d_properties([0, z[frame]])

    bob.set_data([x[frame]], [y[frame]])
    bob.set_3d_properties([z[frame]])
    return line, bob

ani = FuncAnimation(fig, update, frames=steps, init_func=init, blit=True)

ani.save('spinning_pendulum.gif', writer='imagemagick', fps=30)

print("GIF saved as spinning_pendulum.gif")

