import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
k_e = 8.987551787e9       # Coulomb's constant in N·m²/C²
e_charge = 1.602176634e-19  # Elementary charge in C
amu = 1.66053906660e-27   # Atomic mass unit in kg
k_B = 1.380649e-23        # Boltzmann constant in J/K
epsilon_0 = 8.854187817e-12  # Vacuum permittivity in F/m

# Simulation Parameters
num_protons = 100         # Number of protons
num_electrons = 100       # Number of electrons
num_particles = num_protons + num_electrons
box_length = 1e-9         # Simulation box length in meters (1 nm)
temperature = 15e6        # Temperature in Kelvin (Sun's core temperature)
time_step = 1e-19         # Time step in seconds
num_steps = 1000          # Number of simulation steps
softening = 1e-12         # Softening parameter in meters

# Particle Properties
positions = np.random.uniform(0, box_length, (num_particles, 3))
velocities = np.zeros_like(positions)
forces = np.zeros_like(positions)

# Assign masses and charges
masses = np.zeros(num_particles)
charges = np.zeros(num_particles)

# Protons
m_p = 1.00727647 * amu  # Proton mass in kg
positions[:num_protons] = np.random.uniform(0, box_length, (num_protons, 3))
masses[:num_protons] = m_p
charges[:num_protons] = e_charge

# Electrons
m_e = 5.48579909070e-4 * amu  # Electron mass in kg
positions[num_protons:] = np.random.uniform(0, box_length, (num_electrons, 3))
masses[num_protons:] = m_e
charges[num_protons:] = -e_charge

# Initialize velocities from Maxwell-Boltzmann distribution
def initialize_velocities(masses, temperature):
    stddev = np.sqrt(k_B * temperature / masses)
    velocities = np.random.normal(0, stddev[:, np.newaxis], (num_particles, 3))
    return velocities

velocities = initialize_velocities(masses, temperature)

# Apply periodic boundary conditions
def apply_pbc(positions):
    positions %= box_length
    return positions

# Compute forces using Coulomb's law with softening and Debye screening
def compute_forces(positions, charges):
    forces = np.zeros_like(positions)
    potential_energy = 0.0
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            rij = positions[i] - positions[j]
            # Minimum image convention
            rij -= box_length * np.round(rij / box_length)
            r2 = np.dot(rij, rij) + softening**2
            r = np.sqrt(r2)
            # Coulomb force with softening
            force_magnitude = (k_e * charges[i] * charges[j]) / r2
            force_vector = force_magnitude * rij / r
            forces[i] += force_vector
            forces[j] -= force_vector
            # Potential energy
            potential_energy += (k_e * charges[i] * charges[j]) / r
    return forces, potential_energy

# Time evolution using Velocity Verlet algorithm
positions_history = []
kinetic_energies = []
potential_energies = []
total_energies = []

for step in range(num_steps):
    positions_history.append(positions.copy())
    
    # Compute forces and potential energy
    forces, potential_energy = compute_forces(positions, charges)
    
    # Update velocities (half step)
    velocities += 0.5 * (forces / masses[:, np.newaxis]) * time_step
    
    # Update positions
    positions += velocities * time_step
    positions = apply_pbc(positions)
    
    # Compute new forces and potential energy
    new_forces, potential_energy = compute_forces(positions, charges)
    
    # Update velocities (second half step)
    velocities += 0.5 * (new_forces / masses[:, np.newaxis]) * time_step
    
    # Compute energies
    kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
    total_energy = kinetic_energy + potential_energy
    kinetic_energies.append(kinetic_energy)
    potential_energies.append(potential_energy)
    total_energies.append(total_energy)

positions_history = np.array(positions_history)

# Visualization setup (projected onto 2D plane)
fig, ax = plt.subplots(figsize=(6, 6))
scat_protons = ax.scatter([], [], s=20, c='red', label='Protons')
scat_electrons = ax.scatter([], [], s=20, c='blue', label='Electrons')
ax.set_xlim(0, box_length)
ax.set_ylim(0, box_length)
ax.set_title("Hydrogen Plasma Simulation")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_aspect('equal')
ax.legend()

def update(frame):
    current_positions = positions_history[frame]
    scat_protons.set_offsets(current_positions[:num_protons, :2])
    scat_electrons.set_offsets(current_positions[num_protons:, :2])
    return scat_protons, scat_electrons

ani = FuncAnimation(fig, update, frames=num_steps, interval=20, blit=True)
plt.show()

# Plot energies
plt.figure()
plt.plot(kinetic_energies, label='Kinetic Energy')
plt.plot(potential_energies, label='Potential Energy')
plt.plot(total_energies, label='Total Energy')
plt.legend()
plt.xlabel('Time Step')
plt.ylabel('Energy (J)')
plt.title('Energy vs. Time')
plt.show()

