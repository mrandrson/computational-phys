A = [[4, -1, -1], [-2, 6, 1], [-1, 1, 7]]
b = [[3], [9], [-6]]

def l2norm(b):
    b2 =  0
    for i in range(len(b)):
        b2 += b[i][0]**2
    return b2**(1/2)

def gausssiedel(A, b, tol=1e-8, max_iterations=5000):
    n = len(b)
    x = [[0] for i in range(n)]

    for i in range(max_iterations):
        x_old = [[x[j][0]] for j in range(n)]

        for j in range(n):
            sigma = sum(A[j][k] * x[k][0] if k != j else 0 for k in range(n))
            x[j][0] = (b[j][0] - sigma) / A[j][j]

        diff = l2norm([[x[j][0] - x_old[j][0]] for j in range(n)])
        if diff < tol:
            return x

print('Solution:', gausssiedel(A, b))
