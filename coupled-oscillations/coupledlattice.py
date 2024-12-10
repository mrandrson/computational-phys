import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Simulation Parameters
num_cells = 10            # Number of unit cells in each dimension
lattice_constant = 5.26   # Lattice constant in Ångströms (Å) for Argon
mass = 39.948             # Mass of Argon atom in atomic mass units (amu)
epsilon = 0.0103          # Depth of potential well in eV for Argon
sigma = 3.4               # Finite distance at which the inter-particle potential is zero in Å
dt = 1.0                  # Time step in femtoseconds (fs)
num_steps = 1000          # Number of simulation steps
cutoff = 2.5 * sigma      # Cutoff distance for Lennard-Jones potential
temperature = 50000          # Initial temperature in Kelvin

# Conversion factors
eV_to_amu_A2_fs2 = 1.036427e-4  # 1 eV = 1.036427e-4 amu·Å²·fs⁻²
kB_eV_per_K = 8.617333262145e-5  # Boltzmann constant in eV/K

# Convert epsilon from eV to amu·Å²·fs⁻²
epsilon *= eV_to_amu_A2_fs2
kB = kB_eV_per_K * eV_to_amu_A2_fs2  # Boltzmann constant in amu·Å²·fs⁻²·K⁻¹

# Create a 2D square lattice
def create_lattice(num_cells, lattice_constant):
    positions = []
    for i in range(num_cells):
        for j in range(num_cells):
            positions.append([i * lattice_constant, j * lattice_constant])
    return np.array(positions)

positions = create_lattice(num_cells, lattice_constant)
num_atoms = positions.shape[0]

# Initialize velocities from Maxwell-Boltzmann distribution
def initialize_velocities(num_atoms, mass, temperature):
    stddev = np.sqrt(kB * temperature / mass)
    velocities = np.random.normal(0, stddev, (num_atoms, 2))
    velocities -= np.mean(velocities, axis=0)  # Remove net momentum
    return velocities

velocities = initialize_velocities(num_atoms, mass, temperature)

# Define the simulation box
box_length = num_cells * lattice_constant

# Apply periodic boundary conditions
def apply_pbc(positions):
    positions %= box_length
    return positions

# Compute forces and potential energy using the Lennard-Jones potential
def compute_forces(positions):
    forces = np.zeros_like(positions)
    potential_energy = 0.0
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            rij = positions[i] - positions[j]
            # Minimum image convention
            rij -= box_length * np.round(rij / box_length)
            r2 = np.dot(rij, rij)
            if r2 < cutoff**2:
                r6 = (sigma**2 / r2)**3
                r12 = r6**2
                lj_scalar = 24 * epsilon * (2 * r12 - r6) / r2
                force = lj_scalar * rij
                forces[i] += force
                forces[j] -= force
                # Potential energy
                potential_energy += 4 * epsilon * (r12 - r6)
    return forces, potential_energy

# Time evolution using the Velocity Verlet algorithm
positions_history = []
kinetic_energies = []
potential_energies = []
total_energies = []

for step in range(num_steps):
    positions_history.append(positions.copy())

    # Compute forces and potential energy
    forces, potential_energy = compute_forces(positions)
    
    # Update velocities (half step)
    velocities += 0.5 * (forces / mass) * dt

    # Update positions
    positions += velocities * dt
    positions = apply_pbc(positions)

    # Compute new forces and potential energy
    new_forces, potential_energy = compute_forces(positions)

    # Update velocities (second half step)
    velocities += 0.5 * (new_forces / mass) * dt

    # Compute energies
    kinetic_energy = 0.5 * mass * np.sum(velocities**2)
    total_energy = kinetic_energy + potential_energy
    kinetic_energies.append(kinetic_energy)
    potential_energies.append(potential_energy)
    total_energies.append(total_energy)

positions_history = np.array(positions_history)

# Visualization setup
fig, ax = plt.subplots(figsize=(6, 6))
scat = ax.scatter(positions[:, 0], positions[:, 1], s=20, c='blue')
ax.set_xlim(0, box_length)
ax.set_ylim(0, box_length)
ax.set_title("Lennard-Jones Crystal Lattice Simulation")
ax.set_xlabel("x (Å)")
ax.set_ylabel("y (Å)")
ax.set_aspect('equal')

def update(frame):
    current_positions = positions_history[frame]
    scat.set_offsets(current_positions)
    return scat,

ani = FuncAnimation(fig, update, frames=num_steps, interval=20, blit=True)
plt.show()

# Plot energies
plt.figure()
plt.plot(kinetic_energies, label='Kinetic Energy')
plt.plot(potential_energies, label='Potential Energy')
plt.plot(total_energies, label='Total Energy')
plt.legend()
plt.xlabel('Time Step')
plt.ylabel('Energy (amu·Å²·fs⁻²)')
plt.title('Energy vs. Time')
plt.show()

