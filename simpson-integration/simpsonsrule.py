import numpy as np

def simpson(f, x):

    funcseq = f(x[0])
    for i in range(1, len(x)-1):
        if i % 2 == 1:
            funcseq += 4*f(x[i])
        else:
            funcseq += 2*f(x[i])
    deltax = abs(x[1]-x[0])
    funcseq += f(x[-1])

    return (deltax/3)*funcseq

func = lambda x: np.cos(x)

xvec = np.linspace(0, np.pi, int(1e6)+1)

print(simpson(func, xvec))
