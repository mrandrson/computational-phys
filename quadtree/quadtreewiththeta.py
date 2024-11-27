import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class QuadTree:
    def __init__(self, points, Nmax, bounds=None):
        self.points = points
        self.Nmax = Nmax
        self.children = []
        self.bounds = bounds or (
            min(point[0] for point in points),
            max(point[0] for point in points),
            min(point[1] for point in points),
            max(point[1] for point in points),
        )

        if len(points) > Nmax:
            quadrants = self.subdivide(points)
            for quadrant, bounds in quadrants:
                if quadrant:
                    self.children.append(QuadTree(quadrant, Nmax, bounds))

    def subdivide(self, points):
        dmin_x = min(point[0] for point in points)
        dmax_x = max(point[0] for point in points)
        dmin_y = min(point[1] for point in points)
        dmax_y = max(point[1] for point in points)
        mid_x = (dmin_x + dmax_x) / 2
        mid_y = (dmin_y + dmax_y) / 2

        ch1, ch2, ch3, ch4 = [], [], [], []

        for point in points:
            x, y = point
            if x <= mid_x and y > mid_y:
                ch1.append(point)
            elif x > mid_x and y > mid_y:
                ch2.append(point)
            elif x <= mid_x and y <= mid_y:
                ch3.append(point)
            elif x > mid_x and y <= mid_y:
                ch4.append(point)
        print(f"Subdivision: ch1={len(ch1)}, ch2={len(ch2)}, ch3={len(ch3)}, ch4={len(ch4)}") 
        xmin, xmax, ymin, ymax = dmin_x, dmax_x, dmin_y, dmax_y
        return [
            (ch1, (xmin, mid_x, mid_y, ymax)),  
            (ch2, (mid_x, xmax, mid_y, ymax)),  
            (ch3, (xmin, mid_x, ymin, mid_y)),  
            (ch4, (mid_x, xmax, ymin, mid_y)),
        ]
    
    def highlight_active_cells(self, newpos, theta_lim, ax, show_mean=True, connect_mean=False):
        """
        Highlight the active areas of the QuadTree for a given position.
        
        Parameters:
            newpos (tuple): The position at which to compute the potential.
            theta_lim (float): Maximum opening angle (threshold for activeness).
            ax (matplotlib.axes.Axes): The axis object to draw on.
            show_mean (bool): Whether to show the mean position of the cell.
            connect_mean (bool): Whether to connect the mean and the query position.
        """
        if len(self.points) > 0:
            center_of_mass_x = sum(p[0] for p in self.points) / len(self.points)
            center_of_mass_y = sum(p[1] for p in self.points) / len(self.points)
        else:
            center_of_mass_x, center_of_mass_y = 0, 0

        size = max(
            ((p[0] - center_of_mass_x)**2 + (p[1] - center_of_mass_y)**2)**0.5
            for p in self.points
        ) if len(self.points) > 0 else 0
        distance = ((newpos[0] - center_of_mass_x)**2 + (newpos[1] - center_of_mass_y)**2)**0.5

        ang = size / distance if distance != 0 else float('inf')

        if len(self.children) == 0 or ang < theta_lim:
            xmin, xmax, ymin, ymax = self.bounds
            width = xmax - xmin
            height = ymax - ymin
            rect = Rectangle((xmin, ymin), width, height, edgecolor='#9467bd', facecolor='none', lw=1)
            ax.add_patch(rect)

            if show_mean:
                ax.scatter(center_of_mass_x, center_of_mass_y, color='#d62728', zorder=1)
            if connect_mean:
                ax.plot([newpos[0], center_of_mass_x], [newpos[1], center_of_mass_y], '-', color='#2ca02c', lw=0.5, zorder=0)
        else:
            for child in self.children:
                child.highlight_active_cells(newpos, theta_lim, ax, show_mean=show_mean, connect_mean=connect_mean)

    def collect_leaf_cells(self):
        """Recursively collect only the leaf cells of the QuadTree."""
        if not self.children:
            if self.points:
                center_of_mass_x = sum(p[0] for p in self.points) / len(self.points)
                center_of_mass_y = sum(p[1] for p in self.points) / len(self.points)
            else:
                center_of_mass_x, center_of_mass_y = 0, 0
            size = max(
                ((p[0] - center_of_mass_x)**2 + (p[1] - center_of_mass_y)**2)**0.5
                for p in self.points
            ) if self.points else 0

            return [(self.bounds, self.points, size)]
        else:
            leaf_cells = []
            for child in self.children:
                leaf_cells.extend(child.collect_leaf_cells())
            return leaf_cells

def sample_yx3(npos, min_radius=1e10, max_radius=1e11):
    np.random.seed(42)

    radii = (np.random.random(npos) * (max_radius**-1 - min_radius**-1) + min_radius**-1)**-1

    angles = np.random.uniform(0, 2 * np.pi, npos)

    x = radii * np.cos(angles)
    y = radii * np.sin(angles)

    return np.column_stack((x, y))

def plot_quadtree(tree, ax):
    if tree.bounds:
        xmin, xmax, ymin, ymax = tree.bounds

        width = xmax - xmin
        height = ymax - ymin
        rect = Rectangle((xmin, ymin), width, height, edgecolor='k', facecolor='none', lw=1)
        ax.add_patch(rect)

    if len(tree.points) > 0:
        x, y = zip(*tree.points)
        ax.scatter(x, y, s=10)

    for child in tree.children:
        plot_quadtree(child, ax)

points = sample_yx3(int(1e4))
Nmax = 1e2
quadtree = QuadTree(points, Nmax)
all_cells = quadtree.collect_leaf_cells()

for bounds, points, size in all_cells:
    print(f"Bounds: {bounds}, Number of Points: {len(points)}, Size: {size:.2f}")

'''
fig, ax = plt.subplots(figsize=(8, 8))
plot_quadtree(quadtree, ax)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
'''
AU = 1.5e11
query_position = (5*AU/(2**(1/2)), 5*AU/(2**(1/2))) 
theta_lim = 0.1          

fig, ax = plt.subplots(figsize=(8, 8))
quadtree.highlight_active_cells(query_position, theta_lim, ax, show_mean=True, connect_mean=True)

plt.gca().set_aspect('equal', adjustable='box')
plt.scatter(*query_position, color='blue', label='Evaluation Point')
# plt.legend()
plt.show()

