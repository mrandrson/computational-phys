import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = []
y = []
z = []

with open('C60.xyz', 'r') as file:
    lines = file.readlines()

    for line in lines[2:]:
        parts = line.split()
        if len(parts) == 4:
            x.append(float(parts[1]))
            y.append(float(parts[2]))
            z.append(float(parts[3]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='r', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D Plot of Fullerene Molecule')

plt.show()

