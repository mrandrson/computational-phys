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
            ##Debugging Comment##
            '''
            print(f"Order: {order}, dx: {dx}, dy: {dy}, r2: {r2}, Term: {-G*self.total_mass / (r2**0.5)}")
            '''
            return -G*self.total_mass / (r2**0.5)
        elif order == 1:
            return 0
        elif order == 2:
            r = r2**0.5
            nx = dx / r
            ny = dy / r

            Qxx = sum(
                p[2] * (3 * (p[0] - self.center_of_mass_x)**2 - ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2))
                for p in self.particles
            )
            Qyy = sum(
                p[2] * (3 * (p[1] - self.center_of_mass_y)**2 - ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2))
                for p in self.particles
            )
            Qxy = sum(
                p[2] * 3 * (p[0] - self.center_of_mass_x) * (p[1] - self.center_of_mass_y)
                for p in self.particles
            )

            quad_term = (
                (Qxx * nx**2 + Qyy * ny**2 + 2 * Qxy * nx * ny) -
                (Qxx + Qyy)
            ) / (2 * r**3)

            ##Debugging Comment##
            '''
            print(f"Order: {order}, dx: {dx}, dy: {dy}, r2: {r2}, Term: {quad_term}")
            '''
            return -G * quad_term
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
                ##Debugging Print##
                '''
                print('Using Multipole Expansion')
                '''
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
                        ax += a * dx / r
                        ay += a * dy / r
                return ax, ay
            else:
                ax, ay = 0.0, 0.0
                dx, dy = newpos[0] - self.center_of_mass_x, newpos[1] - self.center_of_mass_y
                r2 = dx**2 + dy**2
                r = r2**0.5
                for order in range(multipole_order + 1):
                    if order == 0:
                        a = -G * self.total_mass / r2
                        ##Debugging Print##
                        '''
                        print(f'Monopole Acceleration: {a*dx/r}, {a*dy/r}')
                        '''
                        ax += a * dx / r
                        ay += a * dy / r
                    elif order == 2:
                        Qxx = sum(
                            p[2] * (3 * (p[0] - self.center_of_mass_x)**2 - ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2))
                            for p in self.particles
                        )
                        Qyy = sum(
                            p[2] * (3 * (p[1] - self.center_of_mass_y)**2 - ((p[0] - self.center_of_mass_x)**2 + (p[1] - self.center_of_mass_y)**2))
                            for p in self.particles
                        )
                        Qxy = sum(
                            p[2] * 3 * (p[0] - self.center_of_mass_x) * (p[1] - self.center_of_mass_y)
                            for p in self.particles
                        )
                        r2 = dx**2 + dy**2
                        r = r2**0.5

                        nx = dx / r
                        ny = dy / r

                        quad_term_x = 0
                        quad_term_y = 0

                        quad_term_x += 3 * (Qxx * nx**2 * nx + Qxy * nx**2 * ny + Qxy * ny**2 * nx + Qyy * ny**2 * nx)
                        quad_term_y += 3 * (Qyy * ny**2 * ny + Qxy * ny**2 * nx + Qxy * nx**2 * ny + Qxx * nx**2 * ny)

                        quad_term_x -= (Qxx + Qyy) * nx
                        quad_term_y -= (Qxx + Qyy) * ny

                        quad_term_x *= G / (r**5)
                        quad_term_y *= G / (r**5)
                        ##Debugging Comment##
                        '''
                        print(f"Quadrupole Moments: Qxx={Qxx}, Qyy={Qyy}, Qxy={Qxy}")
                        print(f"Quadrupole Acceleration: {quad_term_x}, {quad_term_y}")
                        '''
                        ax += -quad_term_x
                        ay += -quad_term_y
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
    return ax, ay


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
particles2 = initialize_particles(0, center=(5e11, 5e11), rmin =1e10, rmax=1e11, mmin = 1e20, mmax = 1e25)
centralmass = initialize_particles(1, center=(0, 0), rmin = 1e-1, rmax =1, mmin = 1.989e30, mmax=2e30)
particles = particles1+particles2+centralmass
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

positions=np.logspace(0, 15, 20)
for i in range(len(positions)):
    query_position = (-positions[i], positions[i])

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
'''
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

def compute_magnitude(point):
    """Compute acceleration magnitude for a single grid point."""
    i, j, x, y = point
    query_position = (x, y)
    ax, ay = quadtree.compute_acceleration(query_position, 0.1)
    return i, j, np.sqrt(ax**2 + ay**2)

if __name__ == '__main__':
    x_min, x_max = -1.5e10, 1.5e10
    y_min, y_max = -1.5e10, 1.5e10
    grid_points = 35

    x_vals = np.linspace(x_min, x_max, grid_points)
    y_vals = np.linspace(y_min, y_max, grid_points)
    X, Y = np.meshgrid(x_vals, y_vals)

    grid_points_list = [(i, j, X[i, j], Y[i, j]) for i in range(grid_points) for j in range(grid_points)]

    acceleration_magnitude = np.zeros_like(X)

    with ProcessPoolExecutor(max_workers=16) as executor:
        results = executor.map(compute_magnitude, grid_points_list)

    for i, j, magnitude in results:
        acceleration_magnitude[i, j] = magnitude

    acceleration_magnitude[acceleration_magnitude == 0] = np.min(acceleration_magnitude[acceleration_magnitude > 0])

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(X, Y, np.log10(acceleration_magnitude), levels=50, cmap='viridis')
    cbar = plt.colorbar(contour, label='Log10 Acceleration Magnitude (m/s²)')
    cbar.set_label('Log10 Acceleration Magnitude (m/s²)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Logarithmic Acceleration Magnitude Over Grid')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
'''
