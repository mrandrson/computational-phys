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

    def collect_leaf_cells(self):
        """Recursively collect only the leaf cells of the QuadTree."""
        if not self.children:
            return [(self.bounds, self.points)]
        else:
            leaf_cells = []
            for child in self.children:
                leaf_cells.extend(child.collect_leaf_cells())
            return leaf_cells
'''
def sample_yx3(npos):
    np.random.seed(42)
    pos = np.random.random((npos, 2)) * 2.0 - 1.0
    pos[:, 1] *= 0.2
    pos[:, 1] += pos[:, 0]**3
    return pos
'''

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

for bounds, points in all_cells:
    print(f"Bounds: {bounds}, Number of Points: {len(points)}")

fig, ax = plt.subplots(figsize=(8, 8))
plot_quadtree(quadtree, ax)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
