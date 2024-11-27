import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
G = 6.674e-11

class QuadTree:
    def __init__(self, particles, Nmax, bounds=None):
        """
        QuadTree constructor.
        Particles: List of (x, y, mass) tuples.
        """
        self.particles = particles
        self.Nmax = Nmax
        self.children = []
        self.bounds = bounds or (
            min(p[0] for p in particles),
            max(p[0] for p in particles),
            min(p[1] for p in particles),
            max(p[1] for p in particles),
        )

        self.total_mass = sum(p[2] for p in particles)
        if self.total_mass > 0:
            self.center_of_mass_x = sum(p[0] * p[2] for p in particles) / self.total_mass
            self.center_of_mass_y = sum(p[1] * p[2] for p in particles) / self.total_mass
        else:
            self.center_of_mass_x, self.center_of_mass_y = 0, 0

        if len(particles) > Nmax:
            quadrants = self.subdivide(particles)
            for quadrant, bounds in quadrants:
                if quadrant:
                    self.children.append(QuadTree(quadrant, Nmax, bounds))

    def subdivide(self, particles):
        dmin_x = min(p[0] for p in particles)
        dmax_x = max(p[0] for p in particles)
        dmin_y = min(p[1] for p in particles)
        dmax_y = max(p[1] for p in particles)
        mid_x = (dmin_x + dmax_x) / 2
        mid_y = (dmin_y + dmax_y) / 2

        ch1, ch2, ch3, ch4 = [], [], [], []

        for p in particles:
            x, y, m = p
            if x <= mid_x and y > mid_y:
                ch1.append(p)
            elif x > mid_x and y > mid_y:
                ch2.append(p)
            elif x <= mid_x and y <= mid_y:
                ch3.append(p)
            elif x > mid_x and y <= mid_y:
                ch4.append(p)
        return [
            (ch1, (dmin_x, mid_x, mid_y, dmax_y)),
            (ch2, (mid_x, dmax_x, mid_y, dmax_y)),
            (ch3, (dmin_x, mid_x, dmin_y, mid_y)),
            (ch4, (mid_x, dmax_x, dmin_y, mid_y)),
        ]

    def highlight_active_cells(self, newpos, theta_lim, ax, show_mean=True, connect_mean=False):
        """Highlight active cells based on theta_lim."""
        size = max(
            ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2)**0.5
            for p in self.particles
        ) if len(self.particles) > 0 else 0
        distance = ((newpos[0] - self.center_of_mass_x)**2 + (newpos[1] - self.center_of_mass_y)**2)**0.5

        ang = size / distance if distance != 0 else float('inf')

        if len(self.children) == 0 or ang < theta_lim:
            xmin, xmax, ymin, ymax = self.bounds
            width = xmax - xmin
            height = ymax - ymin
            rect = Rectangle((xmin, ymin), width, height, edgecolor='#9467bd', facecolor='none', lw=1)
            ax.add_patch(rect)

            if show_mean:
                ax.scatter(self.center_of_mass_x, self.center_of_mass_y, color='#d62728', zorder=1)
            if connect_mean:
                ax.plot([newpos[0], self.center_of_mass_x], [newpos[1], self.center_of_mass_y], '-', color='#2ca02c', lw=0.5, zorder=0)
        else:
            for child in self.children:
                child.highlight_active_cells(newpos, theta_lim, ax, show_mean=show_mean, connect_mean=connect_mean)

    def collect_leaf_cells(self):
        """Recursively collect only the leaf cells of the QuadTree."""
        if not self.children:
            size = max(
                ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2)**0.5
                for p in self.particles
            ) if len(self.particles) > 0 else 0
            return [(self.bounds, self.particles, size, self.total_mass)]
        else:
            leaf_cells = []
            for child in self.children:
                leaf_cells.extend(child.collect_leaf_cells())
            return leaf_cells
    
    def _compute_multipole_term(self, order, dx, dy, r2):
        """
        Compute a single term of the multipole expansion.

        Parameters:
            order (int): The order of the multipole term.
            dx (float): x-displacement from the cell's center of mass.
            dy (float): y-displacement from the cell's center of mass.
            r2 (float): Square of the distance to the cell's center of mass.

        Returns:
            float: Contribution of the multipole term to the potential.
        """
        if order == 0:
            return -G*self.total_mass / (r2**0.5)
        elif order == 1:
            return 0
        elif order == 2:
            r = r2**0.5
            nx = dx / r  
            ny = dy / r  

            Qxx = sum(
                p[2] * (3 * (p[0] - self.center_of_mass_x)**2 - r2)
                for p in self.particles
            )
            Qyy = sum(
                p[2] * (3 * (p[1] - self.center_of_mass_y)**2 - r2)
                for p in self.particles
            )
            Qxy = sum(
                p[2] * 3 * (p[0] - self.center_of_mass_x) * (p[1] - self.center_of_mass_y)
                for p in self.particles
            )

            quad_term = (
                0.5 * (Qxx * nx**2 + Qyy * ny**2 + 2 * Qxy * nx * ny) / r**3
            )
            return -G * quad_term
            '''
            return (
                self.total_mass * (dx**2 - dy**2) / (2 * r2**2.5)
            )  # Quadrupole term (simplified example)
            '''
        else:
            raise NotImplementedError("Only up to second-order multipole terms are implemented.")

    def compute_potential(self, newpos, theta_lim, multipole_order=2):
        """
        Compute the gravitational potential at a given point using the multipole expansion.

        Parameters:
            newpos (tuple): The position (x, y) where the potential is evaluated.
            theta_lim (float): Maximum opening angle for active cells.
            multipole_order (int): Maximum order of the multipole expansion (default: 2).

        Returns:
            float: The gravitational potential at `newpos`.
        """
        size = max(
            ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2)**0.5
            for p in self.particles
        ) if len(self.particles) > 0 else 0
        distance = ((newpos[0] - self.center_of_mass_x)**2 + (newpos[1] - self.center_of_mass_y)**2)**0.5

        ang = size / distance if distance != 0 else float('inf')

        if len(self.children) == 0 or ang < theta_lim:
            if len(self.children) == 0:
                potential = 0.0
                for p in self.particles:
                    dx, dy = newpos[0] - p[0], newpos[1] - p[1]
                    r = (dx**2 + dy**2)**0.5
                    if r > 0:
                        potential -= G * p[2] / r
                return potential
            else:
                potential = 0.0
                dx, dy = newpos[0] - self.center_of_mass_x, newpos[1] - self.center_of_mass_y
                r2 = dx**2 + dy**2
                for order in range(multipole_order + 1):
                    potential += self._compute_multipole_term(order, dx, dy, r2)
                return potential

        potential = 0.0
        for child in self.children:
            potential += child.compute_potential(newpos, theta_lim, multipole_order)
        return potential


    def compute_acceleration(self, newpos, theta_lim, multipole_order=2):
        """
        Compute the gravitational acceleration at a given point using the multipole expansion.

        Parameters:
            newpos (tuple): The position (x, y) where the acceleration is evaluated.
            theta_lim (float): Maximum opening angle for active cells.
            multipole_order (int): Maximum order of the multipole expansion (default: 2).

        Returns:
            tuple: Gravitational acceleration (ax, ay).
        """
        size = max(
            ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2)**0.5
            for p in self.particles
        ) if len(self.particles) > 0 else 0
        distance = ((newpos[0] - self.center_of_mass_x)**2 + (newpos[1] - self.center_of_mass_y)**2)**0.5

        ang = size / distance if distance != 0 else float('inf')

        if len(self.children) == 0 or ang < theta_lim:
            if len(self.children) == 0:
                ax, ay = 0.0, 0.0
                for p in self.particles:
                    dx = newpos[0] - p[0]
                    dy = newpos[1] - p[1]
                    r2 = dx**2 + dy**2
                    r = r2**0.5
                    if r > 0:
                        a = -G * p[2] / r2
                        ax -= a * dx / r
                        ay -= a * dy / r
                return ax, ay
            else:
                ax, ay = 0.0, 0.0
                dx, dy = newpos[0] - self.center_of_mass_x, newpos[1] - self.center_of_mass_y
                r2 = dx**2 + dy**2
                r = r2**0.5
                for order in range(multipole_order + 1):
                    if order == 0:
                        a = -G * self.total_mass / r2
                        ax += a * dx / r
                        ay += a * dy / r
                    elif order == 2:
                        nx = dx / r
                        ny = dy / r
                        Qxx = sum(
                            p[2] * (3 * (p[0] - self.center_of_mass_x)**2 - r2)
                            for p in self.particles
                        )
                        Qyy = sum(
                            p[2] * (3 * (p[1] - self.center_of_mass_y)**2 - r2)
                            for p in self.particles
                        )
                        Qxy = sum(
                            p[2] * 3 * (p[0] - self.center_of_mass_x) * (p[1] - self.center_of_mass_y)
                            for p in self.particles
                        )
                        quad_term_x = 1.5 * (Qxx * nx + Qxy * ny) / r2
                        quad_term_y = 1.5 * (Qxy * nx + Qyy * ny) / r2
                        ax += -G * quad_term_x
                        ay += -G * quad_term_y
                return ax, ay

        ax, ay = 0.0, 0.0
        for child in self.children:
            child_ax, child_ay = child.compute_acceleration(newpos, theta_lim, multipole_order)
            ax += child_ax
            ay += child_ay
        return ax, ay


def directsum_acceleration(particles, newpos):
    """
    Compute the gravitational acceleration at a given point using direct summation.

    Parameters:
        particles (list): List of particles, where each particle is (x, y, mass).
        newpos (tuple): The position (x, y) where the acceleration is evaluated.

    Returns:
        tuple: Gravitational acceleration (ax, ay) using direct summation.
    """
    ax, ay = 0.0, 0.0
    for p in particles:
        dx = newpos[0] - p[0]
        dy = newpos[1] - p[1]
        r2 = dx**2 + dy**2
        r = r2**0.5
        if r > 0:
            a = -G * p[2] / r2
            ax += a * dx / r
            ay += a * dy / r
    return -ax, -ay


def directsum(particles, newpos):
    """
    Compute the gravitational potential at a given point using direct summation.

    Parameters:
        particles (list): List of particles, where each particle is (x, y, mass).
        newpos (tuple): The position (x, y) where the potential is evaluated.

    Returns:
        float: The gravitational potential at `newpos` using direct summation.
    """
    potential = 0.0
    for p in particles:
        dx = newpos[0] - p[0]
        dy = newpos[1] - p[1]
        r = (dx**2 + dy**2)**0.5
        if r > 0:
            potential -= G*p[2] / r
    return potential

def initialize_particles(n, center, rmin, rmax, mmin, mmax):
    """
    Initialize particles around a center point with random positions and masses.
    """
    radii = np.random.uniform(rmin, rmax, n)
    angles = np.random.uniform(0, 2 * np.pi, n)
    masses = np.random.uniform(mmin, mmax, n)

    x = center[0] + radii * np.cos(angles)
    y = center[1] + radii * np.sin(angles)

    return [(x[i], y[i], masses[i]) for i in range(n)]

def plot_quadtree(tree, ax):
    """Plot the QuadTree."""
    if tree.bounds:
        xmin, xmax, ymin, ymax = tree.bounds
        width = xmax - xmin
        height = ymax - ymin
        rect = Rectangle((xmin, ymin), width, height, edgecolor='k', facecolor='none', lw=1)
        ax.add_patch(rect)

    if len(tree.particles) > 0:
        x, y, _ = zip(*tree.particles)
        ax.scatter(x, y, s=10)

    for child in tree.children:
        plot_quadtree(child, ax)

numparticles = int(1e5)
particles1 = initialize_particles(numparticles, center=(0, 0), rmin=1e10, rmax=1e11, mmin=1e20, mmax=1e25)
particles2 = initialize_particles(numparticles, center=(2e11, 2e11), rmin =1e10, rmax=1e11, mmin = 1e15, mmax = 1e24)
particles = particles1+particles2
Nmax = numparticles*1e-1
quadtree = QuadTree(particles, Nmax)
##Plot Quadtree##
'''
fig, ax = plt.subplots(figsize=(8, 8))
plot_quadtree(quadtree, ax)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
'''

##Plot Quadtree With Active Cells##
'''
query_position = (-1e10, 2e11)
theta_lim = 0.1

fig, ax = plt.subplots(figsize=(10, 10))
quadtree.highlight_active_cells(query_position, theta_lim, ax, show_mean=True, connect_mean=True)

plt.scatter(*query_position, color='blue', label='Query Position')

plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()
'''

##Multipole Expansion##
'''
theta_limit = 0.1

query_position = (-1e10, 2e11)

direct_potential = directsum(particles, query_position)
print(f"Direct summation potential at {query_position}: {direct_potential}")

multipole_potential = quadtree.compute_potential(query_position, theta_limit)
print(f"Multipole expansion potential at {query_position}: {multipole_potential}")

error = abs(direct_potential - multipole_potential) / abs(direct_potential)
print(f"Relative error between direct and multipole expansion: {error:.2e}")

direct_ax, direct_ay = directsum_acceleration(particles, query_position)
print(f"Direct summation acceleration at {query_position}: ({direct_ax:.3e}, {direct_ay:.3e})")

multipole_ax, multipole_ay = quadtree.compute_acceleration(query_position, theta_limit)
print(f"Multipole expansion acceleration at {query_position}: ({multipole_ax:.3e}, {multipole_ay:.3e})")

error_x = abs(direct_ax - multipole_ax) / abs(direct_ax) if direct_ax != 0 else 0
error_y = abs(direct_ay - multipole_ay) / abs(direct_ay) if direct_ay != 0 else 0
print(f"Relative error in acceleration (x): {error_x:.2e}")
print(f"Relative error in acceleration (y): {error_y:.2e}")
'''
