import numpy as np


# matrix (N + 1) * (N + 2) , where (N+2)st column which has index (N+1) is a right hand-side of equations
# returns vector solution for input system
def tridAlgo(matrix: np.ndarray) -> np.ndarray:
    # N + 1 : number of equations in the system
    N = matrix.shape[0] - 1
    # we search for solution as matrix[i, i] = alpha[i+1]*matrix[i+1, i+1]+ beta[i+1], i from N-1 to 0
    sol = np.ndarray(N+1)

    alpha = np.zeros(N+1)
    beta = np.zeros(N+1)

    alpha[1] = - matrix[0, 1]/matrix[0, 0]
    beta[1] = matrix[0, N+1]/matrix[0, 0]

    for i in range(1, N):
        z = - matrix[i, i] - matrix[i, i-1]*alpha[i]
        alpha[i+1] = matrix[i, i+1]/z
        beta[i+1] = (-matrix[i, N+1] + matrix[i, i-1]*beta[i])/z

    sol[N] = (- matrix[N, N+1] + matrix[N, N-1]*beta[N])/(- matrix[N, N] - matrix[N, N-1]*alpha[N])
    # from N-1 to 0
    for i in reversed(range(N)):
        sol[i] = alpha[i+1]*sol[i+1] + beta[i+1]
    return sol

def main():
    A = np.array([[2, -1, 0], [5, 4, 2], [0, 1, -3]])
    b = np.array([3, 6, 2])
    sol = np.array(np.linalg.solve(A, b))
    print(sol)

    C = np.array([[2, -1, 0, 3], [5, 4, 2, 6], [0, 1, -3, 2]])
    print(tridAlgo(C))



main()
