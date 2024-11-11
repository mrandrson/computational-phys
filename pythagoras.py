import matplotlib.pyplot as plt

def get_c(a, b):
    return (a**2+b**2)**(1/2)

def mod(a, b):
    if b == 0:
        raise ValueError("The divisor 'b' cannot be zero.")
    while a >= b:
        a = a - b
    return a

pythpairs = []
reps = int(1e3)
for i in range(1, reps):
    for j in range(1, reps):
        c = get_c(i, j)
        if mod(c, int(c)) == 0:
            #print(i, j)
            pythpairs.append([i, j])

#print(pythpairs)
x = [pythpairs[i][0] for i in range(len(pythpairs))]
y = [pythpairs[i][1] for i in range(len(pythpairs))]
#plt.xscale('log')
#plt.yscale('log')
plt.scatter(x, y) 
plt.show()  
