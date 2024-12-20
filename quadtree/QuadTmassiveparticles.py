import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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

particles1 = initialize_particles(1000, center=(0, 0), rmin=1e10, rmax=1e11, mmin=1e20, mmax=1e25)
particles2 = initialize_particles(1000, center=(2e11, 2e11), rmin =1e10, rmax=1e11, mmin = 1e15, mmax = 1e24)
particles = particles1+particles2
Nmax = 25
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
query_position = (-2e11, 2e11)
theta_lim = 0.1

fig, ax = plt.subplots(figsize=(10, 10))
quadtree.highlight_active_cells(query_position, theta_lim, ax, show_mean=True, connect_mean=True)

plt.scatter(*query_position, color='blue', label='Query Position')

plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()
'''
