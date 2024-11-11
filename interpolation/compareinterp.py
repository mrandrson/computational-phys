import numpy as np
from matplotlib import pyplot as plt


def polynewton(x, y):
    n = len(x)
    divided_diffs = [y.copy()]
    for i in range(1, n):
        column = []
        for j in range(n - i):
            diff = (divided_diffs[i - 1][j + 1] - divided_diffs[i - 1][j]) / (x[j + i] - x[j])
            column.append(diff)
        divided_diffs.append(column)
    
    def newton_poly(x_val):
        result = divided_diffs[0][0]
        prod = 1
        for i in range(1, n):
            prod *= (x_val - x[i - 1])
            result += divided_diffs[i][0] * prod
        return result
    
    return newton_poly

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
with open('points.txt','r') as fp:
  for i in range(5):
    line=fp.readline()
    tokens=line.split()
    x.append(float(tokens[0]))
    y.append(float(tokens[1]))
#
# Evaluate poynomial from a to b
# using n points
#
a=-6.0
b=6.0
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

newtonpolynomial = polynewton(x, y) 
##Add the residual calculation here.
lagrange_residuals = [y[i] - eval_l(x, y, x[i]) for i in range(len(x))]
newton_residuals = [y[i] - newtonpolynomial(x[i]) for i in range(len(x))]
print(lagrange_residuals)
print(newton_residuals)
plt.plot(x,y,'o')
plt.plot(z,w,'red', label = 'Lagrange Interpolation')
plt.plot(z, newtonpolynomial(z), label='Newton Interpolation')
plt.legend()
plt.show()
