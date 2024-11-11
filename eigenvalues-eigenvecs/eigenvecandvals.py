def power_method(A, x, tol=1e-6, max_iterations=1000):
    n = len(A)
    lambda_old = 0

    for iteration in range(1, max_iterations + 1):
        y = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
        lambda_new = max(abs(yi) for yi in y)
        x = [y[i] / lambda_new for i in range(n)]
        
        if abs(lambda_new - lambda_old) < tol:
            print(f"Converged after {iteration} iterations.")
            return lambda_new, x
        
        lambda_old = lambda_new


A = [[1, 2, 0],
    [-2, 1, 2],
    [1, 3, 1]]

x = [1, 1, 1]

eigenvalue, eigenvector = power_method(A, x)
print("\nLargest Eigenvalue:", eigenvalue)
print("Corresponding Eigenvector:", eigenvector)

