class QuadTree:

    def __init__(self, data, mins = None, maxs = None, maxiterations = 3):
        
        self.data = [tuple(point) for point in data]
        
        if mins is None:
            mins = [min(coord[i] for coord in self.data) for i in range(2)]
        if maxs is None:
            maxs = [max(coord[i] for coord in self.data) for i in range(2)]

        self.mins = mins
        self.maxs = maxs
        self.maxiterations = maxiterations
        self.children = []
        self.subdivide()

    def childnode(self, dat, xmin, ymin, xmax, ymax, iterations):
        
        if iterations >= 0 or len(dat) == 0:
            return None

        xmid = (xmax-xmin)/2
        ymid = (ymax-ymin)/2

        data_q1 = dat[(dat[:, 0] < xmid) & (dat[:, 1] < ymid)]
        data_q2 = dat[(dat[:, 0] < xmid) & (dat[:, 1] >= ymid)]
        data_q3 = dat[(dat[:, 0] >= xmid) & (dat[:, 1] < ymid)]
        data_q4 = dat[(dat[:, 0] >= xmid) & (dat[:, 1] >= ymid)]

        children = []

        if data_q1.shape[0] > 0:
            children.append(QuadTree(data_q1, [xmin, ymin], [xmid, ymid], iterations - 1))
        if data_q2.shape[0] > 0:
            children.append(QuadTree(data_q2, [xmin, ymid], [xmid, ymax], iterations - 1))
        if data_q3.shape[0] > 0:
            children.append(QuadTree(data_q3, [xmid, ymin], [xmax, ymid], iterations - 1))
        if data_q4.shape[0] > 0:
            children.append(QuadTree(data_q4, [xmid, ymid], [xmax, ymax], iterations - 1))

        return children

    def subdivide(self):
        """Subdivide the quadtree into four child nodes."""
        if self.maxiterations <= 0 or len(self.data) == 0:
            return
        
        xmid = (self.mins[0] + self.maxs[0]) / 2
        ymid = (self.mins[1] + self.maxs[1]) / 2
        
        data_q1 = [point for point in self.data if point[0] < xmid and point[1] < ymid]
        data_q2 = [point for point in self.data if point[0] < xmid and point[1] >= ymid]
        data_q3 = [point for point in self.data if point[0] >= xmid and point[1] < ymid]
        data_q4 = [point for point in self.data if point[0] >= xmid and point[1] >= ymid]
        
        if data_q1:
            self.children.append(QuadTree(data_q1, [self.mins[0], self.mins[1]], [xmid, ymid], self.maxiterations - 1))
        if data_q2:
            self.children.append(QuadTree(data_q2, [self.mins[0], ymid], [xmid, self.maxs[1]], self.maxiterations - 1))
        if data_q3:
            self.children.append(QuadTree(data_q3, [xmid, self.mins[1]], [self.maxs[0], ymid], self.maxiterations - 1))
        if data_q4:
            self.children.append(QuadTree(data_q4, [xmid, ymid], [self.maxs[0], self.maxs[1]], self.maxiterations - 1))

if __name__ == "__main__":
    data = [(1, 1), (2, 3), (4, 2), (3, 5), (6, 7), (7, 8), (8, 9)]

    qt = QuadTree(data)

    print("Root Node:")
    print(f"Data: {qt.data}")
    print(f"Mins: {qt.mins}, Maxs: {qt.maxs}")
    print(f"Number of Children: {len(qt.children)}")

    if qt.children:
        print("\nFirst Child Node:")
        print(f"Data: {qt.children[0].data}")
        print(f"Mins: {qt.children[0].mins}, Maxs: {qt.children[0].maxs}")
