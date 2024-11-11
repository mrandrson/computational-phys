import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

total_time = 1000
N = total_time * 500
dt = total_time / N
m = 1
l = 1
g = 9.81

def dHdtheta(theta):
    return m * g * l * np.sin(theta)

def dHdp(ptheta):
    return ptheta / (m * l**2)

def update_verlet(theta, ptheta):
    phalfstep = ptheta - (dt / 2) * dHdtheta(theta)
    thetanext = theta + dt * dHdp(phalfstep)
    pthetanext = phalfstep - (dt / 2) * dHdtheta(thetanext)
    return thetanext, pthetanext

theta_verlet = np.zeros(N)
ptheta_verlet = np.zeros(N)

theta_verlet[0] = np.pi / 3
ptheta_verlet[0] = 0

for i in range(N - 1):
    theta_verlet[i + 1], ptheta_verlet[i + 1] = update_verlet(theta_verlet[i], ptheta_verlet[i])

time = np.arange(N) * dt

##Phase Plot##
'''
plt.figure(figsize=(12, 6))
plt.plot(theta_verlet, ptheta_verlet, label="Störmer-Verlet", color="blue")

plt.xlabel(r"$\theta$")
plt.ylabel(r"$\dot{\theta}$")
plt.title("Phase Diagram: Störmer-Verlet")
plt.legend()
plt.show()
'''

##Angle vs Time Plot##
'''
plt.plot(time, theta_verlet, label="Störmer-Verlet (Angle)", color="blue")

plt.xlabel("Time (s)")
plt.ylabel("Angle (rad)")
plt.legend()
plt.show()
'''

##Animation##
'''
fig, ax = plt.subplots()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')

line_verlet, = ax.plot([], [], '-o', color="blue", label="Störmer-Verlet")
ax.legend()

def init():
    line_verlet.set_data([], [])
    return line_verlet,

def animate(i):
    x_verlet = l * np.sin(theta_verlet[i])
    y_verlet = -l * np.cos(theta_verlet[i])

    line_verlet.set_data([0, x_verlet], [0, y_verlet])
    return line_verlet,

ani = FuncAnimation(fig, animate, frames=N, init_func=init, blit=True, interval=1)
plt.show()
'''
