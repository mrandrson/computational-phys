import numpy as np
from matplotlib import pyplot as plt
#
# Function to evaluate Lagrange polynomial
#
def lprod(x,i,z):
  q=len(x)
  result=1.0
  for m in range(q):
    if m!=i:
      result=result*(z-x[m])/(x[i]-x[m])
  return result
#
# Function to evaluate polynomial
#
def eval_l(x,y,z):
  n=len(x)
  result=0.0
  for i in range(n):
    result=result+y[i]*lprod(x,i,z)
  return result
#
# Main program
#
#
# Read points used to generate polynomial
#
x=[]
y=[]
with open('lagrange_points.txt','r') as fp:
  for i in range(5):
    line=fp.readline()
    tokens=line.split()
    x.append(float(tokens[0]))
    y.append(float(tokens[1]))
#
# Evaluate poynomial from a to b
# using n points
#
a=-5.0
b=5.0
n=100
#
# Increment for points will be dx
#
dx=(b-a)/n
#
# Write polynomial points for plotting
#
#fptr=open('lagrange.txt','w')
#for i in range(n):
#  z.append(a+i*dx)
#  w=eval_l(x,y,z)
#  fptr.write('{} {}\n'.format(z,w))
#fptr.close()
#
# Generate array for x value
#
z=np.arange(a,b,dx)
#
# Evaluate polynomail
#
w=eval_l(x,y,z)
#
# Plot polynomial
# and original points
#
plt.plot(x,y,'o')
plt.plot(z,w,'red')
plt.show()
