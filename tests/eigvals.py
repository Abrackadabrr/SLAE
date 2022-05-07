import numpy as np

f = open("/media/evgen/Big_disc/MIPT/2nd level/Chapter 4/SLAE/cmake-build-debug/matrix.txt")
A = []
for line in f.readlines():
    row = line.split(' ')
    A.append(row)


A = np.array(A, dtype=np.float64)
# print(A)
n = np.int32(np.sqrt(A.shape[1])//1)
print(n)
L = np.diag(np.diag(A, -1), -1) + np.diag(np.diag(A, -n), -n)
U = np.diag(np.diag(A, 1), 1) + np.diag(np.diag(A, n), n)
D = np.diag(np.diag(A, 0), 0)

New = np.eye(n*n) - np.linalg.inv(D)@A
mu = max(np.abs(np.linalg.eigvals(New)))
print('For sor omega = ', 1 + (mu/(1 + np.sqrt(1 - mu*mu)))**2)

New = np.linalg.inv(D + U) @ L @ np.linalg.inv(L + D) @ U
print('For sym GS rho = ', max(np.linalg.eigvals(New)))