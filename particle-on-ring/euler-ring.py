import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


total_time = 100
N = total_time*500
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

def update_euler(theta, ptheta):
    thetanext = theta + dt * dHdp(ptheta)
    pthetanext = ptheta - dt * dHdtheta(theta)
    return thetanext, pthetanext

theta_verlet = np.zeros(N)
ptheta_verlet = np.zeros(N)
theta_euler = np.zeros(N)
ptheta_euler = np.zeros(N)

theta_verlet[0] = np.pi / 3   
ptheta_verlet[0] = 0          

theta_euler[0] = np.pi / 3    
ptheta_euler[0] = 0

for i in range(N - 1):
    theta_verlet[i + 1], ptheta_verlet[i + 1] = update_verlet(theta_verlet[i], ptheta_verlet[i])
    
    theta_euler[i + 1], ptheta_euler[i + 1] = update_euler(theta_euler[i], ptheta_euler[i])

time = np.arange(N) * dt

##Phase Plot##
'''
plt.figure(figsize=(12, 6))
plt.plot(theta_verlet, ptheta_verlet, label="Störmer-Verlet ", color="blue")

plt.plot(theta_euler, ptheta_euler, label="Explicit Euler", color="red", linestyle="--")

plt.xlabel(r"$\theta$")
plt.ylabel(r"$\dot{\theta}$")
plt.title("Phase Diagram: Störmer-Verlet vs. Explicit Euler")
plt.legend()
plt.show()
'''

##Angle vs Time Plot##
'''
plt.plot(time, theta_verlet, label="Störmer-Verlet (Angle)", color="blue")

plt.plot(time, theta_euler, label="Explicit Euler (Angle)", color="red", linestyle="--")

plt.xlabel("Time (s)")
plt.ylabel("Angle (rad)")
plt.legend()
plt.show()
'''

##Animation##
fig, ax = plt.subplots()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')

line_verlet, = ax.plot([], [], '-o', color="blue", label="Störmer-Verlet")
line_euler, = ax.plot([], [], '-o', color="red", linestyle="--", label="Explicit Euler")
ax.legend()

def init():
    line_verlet.set_data([], [])
    line_euler.set_data([], [])
    return line_verlet, line_euler

def animate(i):
    x_verlet = l * np.sin(theta_verlet[i])
    y_verlet = -l * np.cos(theta_verlet[i])
    x_euler = l * np.sin(theta_euler[i])
    y_euler = -l * np.cos(theta_euler[i])
    
    line_verlet.set_data([0, x_verlet], [0, y_verlet])
    line_euler.set_data([0, x_euler], [0, y_euler])
    return line_verlet, line_euler

ani = FuncAnimation(fig, animate, frames=N, init_func=init, blit=True, interval=1)
plt.show()


##Energy Error Plot##
'''
def total_energy(theta, ptheta):
    kinetic_energy = ptheta**2 / (2 * m * l**2)
    potential_energy = m * g * l * (1 - np.cos(theta))
    return kinetic_energy + potential_energy

energy_verlet = np.zeros(N)
energy_euler = np.zeros(N)

for i in range(N):
    energy_verlet[i] = total_energy(theta_verlet[i], ptheta_verlet[i])
    energy_euler[i] = total_energy(theta_euler[i], ptheta_euler[i])

initial_energy_verlet = energy_verlet[0]
initial_energy_euler = energy_euler[0]

energy_error_verlet = (energy_verlet - initial_energy_verlet) / initial_energy_verlet
energy_error_euler = (energy_euler - initial_energy_euler) / initial_energy_euler

plt.figure(figsize=(12, 6))
plt.plot(time, energy_error_verlet, label="Störmer-Verlet Energy Error", color="blue")
plt.plot(time, energy_error_euler, label="Explicit Euler Energy Error", color="red", linestyle="--")

plt.xlabel("Time (s)")
plt.ylabel("Relative Energy Error")
plt.title("Relative Energy Error: Störmer-Verlet vs. Explicit Euler")
plt.legend()
plt.show()
'''
