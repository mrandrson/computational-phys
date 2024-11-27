import numpy as np
from multipoleforce import QuadTree, initialize_particles, directsum_acceleration

class StormerVerlet:
    def __init__(self, particles, acceleration_function, dt, Nmax=10, theta_lim=0.1, initial_velocities=None):
        """
        Initialize the Stormer-Verlet integrator.

        Parameters:
            particles (list): List of particles, where each particle is (x, y, mass).
            acceleration_function (callable): Function to compute accelerations.
                Should take `particles`, `position`, and optional parameters like `theta_lim`.
            dt (float): Timestep for integration.
            Nmax (int): Maximum number of particles in a leaf cell for QuadTree subdivision.
            theta_lim (float): Maximum opening angle for multipole acceleration.
            initial_velocities (list or None): List of initial velocities [(vx, vy), ...] for particles.
                Defaults to zero velocity for all particles if None.
        """
        self.particles = particles
        self.acceleration_function = acceleration_function
        self.dt = dt
        self.Nmax = Nmax
        self.theta_lim = theta_lim

        self.positions = np.array([[p[0], p[1]] for p in particles])
        self.masses = np.array([p[2] for p in particles])

        if initial_velocities is None:
            self.velocities = np.zeros_like(self.positions)
        else:
            self.velocities = np.array(initial_velocities)

        self.quadtree = QuadTree(self._particle_list(), self.Nmax)

    def _particle_list(self):
        """Return an updated list of particles with current positions and masses."""
        return [(self.positions[i, 0], self.positions[i, 1], self.masses[i]) for i in range(len(self.particles))]

    def compute_accelerations(self):
        """
        Compute the accelerations for all particles using the provided acceleration function.
        """
        self.quadtree = QuadTree(self._particle_list(), self.Nmax)

        accelerations = np.zeros_like(self.positions)
        for i, pos in enumerate(self.positions):
            accelerations[i] = self.acceleration_function(self.particles, tuple(pos), self.quadtree, self.theta_lim)
        return accelerations

    def integrate(self, steps):
        """
        Integrate the motion of the particles using the Störmer-Verlet scheme.

        Parameters:
            steps (int): Number of integration steps to perform.

        Returns:
            list: A list of position arrays at each timestep.
        """
        trajectories = [self.positions.copy()]

        accelerations = -self.compute_accelerations()

        for _ in range(steps):
            self.velocities += 0.5 * accelerations * self.dt

            self.positions += self.velocities * self.dt

            accelerations_new = -self.compute_accelerations()

            self.velocities += 0.5 * accelerations_new * self.dt

            accelerations = accelerations_new
            print(accelerations, _)

            trajectories.append(self.positions.copy())

        return trajectories


if __name__ == "__main__":
    numparticles = 100
    particles = initialize_particles(1, center=(0, 0), rmin=1.5e11, rmax=1.55e11, mmin=1.1e25, mmax=1e25) + initialize_particles(1, center=(0, 0), rmin=1e-2, rmax=1, mmin=1.989e30, mmax=2e30) 

    def multipole_acceleration(particles, position, quadtree, theta_lim):
        if quadtree is None:
            raise ValueError("QuadTree must be initialized for multipole acceleration.")
        return quadtree.compute_acceleration(position, theta_lim)

    def direct_acceleration(particles, position, quadtree, theta_lim):
        return directsum_acceleration(particles, position)

    dt = 1e6
    steps = 1000
    initial_velocities = [[0, 1e4], [0, 0]]
    integrator = StormerVerlet(particles, multipole_acceleration, dt, Nmax=10, theta_lim=0.1, initial_velocities=initial_velocities)
    trajectories = integrator.integrate(steps)
    import matplotlib.pyplot as plt
    n_particles = trajectories[0].shape[0]  # Number of particles
    trajectories = np.array(trajectories)  # Convert to a single array for easier slicing

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_title("Particle Trajectories")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

    # Plot each particle's trajectory
    for i in range(n_particles):
        x_coords = trajectories[:, i, 0]
        y_coords = trajectories[:, i, 1]
        ax.plot(x_coords, y_coords, label=f"Particle {i+1}")

    ax.grid()
    ax.legend()
    plt.show()

