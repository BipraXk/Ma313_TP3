import numpy as np
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def methodeI(n):
    A = np.zeros((n, n))
    b = np.zeros((n, 1))

    for i in range(n):
        b[i] = np.sin((i + 1) / 11)  # car on commence a 0 donc on rajoute +1
        for j in range(n):
            if i == j:
                A[i, j] = 1
            else:
                A[i, j] = 1 / (5 + (2 * (i) - 3 * (j)) ** 2)
    return (A, b), n


def methodeII(n):
    c = np.zeros((n, n))
    d = np.zeros((n, 1))
    for i in range(n):
        d[i] = np.cos((i + 1) / 5)
        for j in range(n):
            if i == j:
                c[i, i] = 2
            else:
                c[i, j] = 1 / (1 + 4 * abs(i - j))
    return (c, d), n


def methodeIII(n):
    A = np.zeros((n, n))
    b = np.zeros((n, 1))

    alpha = 11
    beta = -3
    gamma = -7

    for i in range(0, n):
        b[i] = np.sin(i / 11)
        A[i, i] = alpha
        if i < n - 1:
            A[i, i + 1] = beta
            A[i + 1, i] = gamma
    print(A)
    return (A, b), n


def MIJacobi(methode):
    (A, b), n = methode

    # N = np.zeros((n, n))
    M = np.zeros((n, n))
    x = np.zeros((n, 1))

    array_norm = []
    array_k = []

    for i in range(0, n):
        M[i, i] = A[i, i]
    N = (M - A)
    print(M, N)
    k = 0
    norm_r = 3
    condition = True
    while condition:
        x = np.linalg.solve(M, N @ x + b)
        r = A @ x - b
        norm_r = np.linalg.norm(r)
        array_norm.append(norm_r)
        k += 1
        array_k.append(k)
        condition = (k < n) and (norm_r > 10 ** -10)

    plt.xlabel("Nombre d'itérations de k")
    plt.ylabel("Valeur de norme r^(k)")
    plt.semilogy(array_k, array_norm, label="Jacobi")
    plt.title("La norme de r^k en fonction de k")
    plt.show()

    return x, norm_r, k


def MIGaussSeidel(methode):
    (A, b), n = methode
    N = np.zeros((n, n))
    M = np.zeros((n, n))
    x = np.zeros((n, 1))
    array_norm = []  # evolution de la norme de r^k
    array_k = []  # evolution de k
    for i in range(0, n):
        for j in range(0, n):
            if i >= j:
                M[i, j] = A[i, j]
    N = M - A
    k = 0
    normr = 3
    condition = True
    while condition:
        x = np.linalg.solve(M, N @ x + b)
        r = A @ x - b
        norm_r = np.linalg.norm(r)
        array_norm.append(normr)
        k += 1
        array_k.append(k)
        condition = (k < n) and (norm_r > 10 ** -10)
    plt.axes(xlabel="nombre d'itérations de k", ylabel="valeur de norme r^k)")
    plt.semilogy(array_k, array_norm, label="GaussSeidel")
    plt.title("évolution de la norme de r^k en fonction de k")
    plt.show()
    return array_norm, k


MIJacobi(methodeI(150))
MIGaussSeidel(methodeI(150))