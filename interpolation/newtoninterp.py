import numpy as np
from matplotlib import pyplot as plt
#
# Function to evaluate polynomial
#
def poly(z):
 return 1.7+4.0*(z-3.8)-0.303*(z-3.8)*(z-4.2)-0.158*(z-3.8)*(z-4.2)*(z-5.3)
#
# Read original points to plot
#
x=[]
y=[]
with open('newton_points.txt','r') as fp:
  for i in range(4):
    line=fp.readline()
    tokens=line.split()
    x.append(float(tokens[0]))
    y.append(float(tokens[1]))
#
# Evaluate polynomial to plot from a
# to b using n points
#
a=0.0
b=8.0
n=100
#
# Incements will be h
#
h=(b-a)/n
#
# Create array for x values
#
z=np.arange(a,b,h)
#
# Evaluate polynomail 
#
p=poly(z)
#
# Plot polynomial and original points
#
plt.plot(x,y,'o')
plt.plot(z,p,'red')
plt.show()
