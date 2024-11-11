import random

totalpoints = int(1e7)
circlepoints = 0
for i in range(1, totalpoints):
    x = random.uniform(0, 1)
    y = random.uniform(0, 1)
    if x**2+y**2 <= 1:
        circlepoints += 1

pival = 4*circlepoints/totalpoints

print(pival)
