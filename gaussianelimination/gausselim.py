def gaussian_elimination(A, b):
    n = len(b)
    
    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > abs(A[max_row][i]):
                max_row = k
        
        A[i], A[max_row] = A[max_row], A[i]
        b[i], b[max_row] = b[max_row], b[i]
        
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]
    print('Row Reduced Matrix: ', A)
    x = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        sum_ax = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (b[i] - sum_ax) / A[i][i]
    
    return x

A = [[3, 2, 1],
    [2, 3, 1],
    [1, 1, 4]]

b = [11, 13, 12]

solution = gaussian_elimination(A, b)
print("Solution:", solution)
