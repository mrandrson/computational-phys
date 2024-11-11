import numpy as np

def derivative(f, x):
    h = 1e-6
    return (f(x+h)-f(x-h))/(2*h)

def newtonroot(f, x0):
    xarr = [x0]
    while abs(f(xarr[-1])) > 1e-8:
        fp = derivative(f, xarr[-1])
        xnext = xarr[-1] - f(xarr[-1]) / fp
        xarr.append(xnext)
    return xarr[-1]

func = lambda x: x-np.cos(x)

print(newtonroot(func, 2))
