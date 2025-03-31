
# convert from gij coefficients (non-factorial) in f(z^i, zbar^k)
# to   complex Cxy coefficients (non-factorial) in f(x^i, y^k)
# transpose this matrix and use it in 'Chua_nonlinear_response_scatter.py'

import numpy as np

def expand (Ni, Nk):
    # expand the function z^i*zbar^k
    # as a polynomial in (x^i, y^k)

    Z = np.zeros((Ni + Nk + 1), complex)    # coeff of x^i*y^k
    Z[0] = 1
    for i in range (Ni):
        for j in range (Ni + Nk, 0, -1):
            Z[j] = Z[j] + complex (0, 1)*Z[j - 1]   # calc (x + iy)^Ni
            #print (i, j, Z)
    for i in range (Nk):
        for j in range (Ni + Nk, 0, -1):
            Z[j] = Z[j] + complex (0, -1)*Z[j - 1]  # calc (x - iy)^Nk
            #print (i, j, Z)
    return Z

# generate a uniform polynomial of degree N

N = 5                       # quintic
M = np.zeros((N + 1, N + 1), complex)
for r in range (N + 1):
    row = expand (N - r, r)
    M[r] = row
    print ('(', N - r, ',', r, ')', row)
print ('\nM    = \n', M)
print ('\nMinv = \n', 32*np.linalg.inv(M))  # use 2^N
