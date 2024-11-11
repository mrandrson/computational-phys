import rebound
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

G = 1           # Gravitational constant
M_bh = 1e6      # Black hole mass
N = 1000         # Number of particles in the accretion disk
r_min = 1       # Innermost stable circular orbit for Schwarzschild
r_max = 10      # Outer edge of the disk

sim = rebound.Simulation()
sim.G = G
sim.add(m=M_bh)  # Add the central black hole

def orbital_velocity(r):
    return np.sqrt(G * M_bh / r)

for _ in range(N):
    r = np.random.uniform(r_min, r_max)
    theta = np.random.uniform(0, 2 * np.pi)
    x, y = r * np.cos(theta), r * np.sin(theta)
    vx, vy = -orbital_velocity(r) * np.sin(theta), orbital_velocity(r) * np.cos(theta)
    sim.add(m=0, x=x, y=y, vx=vx, vy=vy)

sim.integrator = "whfast"
sim.dt = 1e-4
sim.move_to_com()

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(-r_max, r_max)
ax.set_ylim(-r_max, r_max)
ax.set_title("Accretion Disk Simulation around Schwarzschild Black Hole")
ax.set_xlabel("x (AU)")
ax.set_ylabel("y (AU)")

particles, = ax.plot([], [], 'o', markersize=1, color="blue")
black_hole, = ax.plot(0, 0, 'o', color="black", markersize=10)

def init():
    particles.set_data([], [])
    return particles, black_hole

def update(frame):
    sim.integrate(frame * sim.dt)
    x_vals = [p.x for p in sim.particles[1:]]
    y_vals = [p.y for p in sim.particles[1:]]
    particles.set_data(x_vals, y_vals)
    return particles, black_hole

ani = FuncAnimation(fig, update, frames=1000, init_func=init, blit=True, interval=20)
plt.show()

