
# fit a sine/cos wave to Delayed Logistic: frequency 0, 1, 2
# see: https://docs.scipy.org/doc/scipy-1.16.2/reference/generated/scipy.optimize.fsolve.html

# see attachment 'Fourier Analysis of Delayed Logistic', Dec 8, 2025, p. 6
# original (complex) parameters are: b0, b1, b2 (complex), w (real)
# assume origin of time is arbitrary, set b1 real
# internal function  parameters are: x[0] - x[4]
# map to : [w, b0, b1, b2R, b2I, b3R, b3I, ...]

import numpy as np
from scipy.optimize import fsolve

def func(x):
    deg = int((len(x) - 1)/2)                                   # highest Fourier frequency
    #print ('x len =', len(x), ',', deg)
    rad = x[0]*np.pi/180
    fw  = complex (np.cos(rad), np.sin(rad))                    # exp(iw)
    fb = np.full(deg + 1, complex(0, 0))                        # array of coeff
    fb[0] = complex (x[1], 0)                                   # bo real (constant)
    fb[1] = complex (x[2], 0)                                   # b1 real (arbitrary t = 0)
    if deg > 1:
        for i in range (deg - 1):
            fb[i + 2] = complex(x[2*i + 3], x[2*i + 4])
    #print ('internal fw, fb[i]:', fw, ',', fb)
    wrap = np.full((2*deg + 1), complex(0, 0))                  # row of b[j] from j = -N to N
    for j in range (deg + 1):
        wrap[j]       = fb[deg - j].conjugate()
        wrap[j + deg] = fb[j]
    #print ('wrap           [i]:', fw, ',', wrap)
    M = np.full((deg + 1, 2*deg + 1), complex(0, 0))            # array of coeff of constraints
    for j in range (deg + 1):
        M[j][deg - j] += fb[j]
        M[j][deg]     -= a*fb[j]
        for k in range (j, 2*deg + 1):
            M[j][k] += a*wrap[k]*wrap[2*deg - k + j]
    #print ('M =', M, '\n')

    ret = np.full(deg + 1, complex(0, 0))                       # array of complex returns
    for j in range (deg + 1):
        for k in range (2*deg + 1):
            ret[j] += np.power(fw, deg - k)*M[j][k]
    #ret[0] = a*fb[2]*fb[2].conjugate()*np.power(fw, 2) + a*fb[1]*fb[1].conjugate()*fw + a*fb[0]*fb[0] + fb[0] - a*fb[0] + a*fb[1]*fb[1].conjugate()*fw.conjugate() + a*fb[2]*fb[2].conjugate()*np.power(fw, -2)
    #ret[1] = (a*fb[2]*fb[1].conjugate() + fb[1])*fw + a*fb[0]*fb[1] - a*fb[1] + a*fb[0]*fb[1]*fw.conjugate() + a*fb[2]*fb[1].conjugate()*np.power(fw, -2)
    #ret[2] = fb[2]*np.power(fw, 2) + a*fb[2]*fb[0] - a*fb[2] + a*fb[1]*fb[1]*fw.conjugate() + a*fb[0]*fb[2]*np.power(fw, -2)
    # LINEAR
    #ret[0] = a*fb[1]*fb[1].conjugate()*fw + a*fb[0]*fb[0] + fb[0] - a*fb[0] + a*fb[1]*fb[1].conjugate()*fw.conjugate()
    #ret[1] = fb[1]*fw + a*fb[0]*fb[1] - a*fb[1] + a*fb[0]*fb[1]*fw.conjugate()
    #print ('internal ret:', ret)
    arr = [ret[0].real, ret[1].real, ret[1].imag]
    if deg > 1:
        for i in range (deg - 1):
            arr.append(ret[i + 2].real)
            arr.append(ret[i + 2].imag)
    #print ('internal arr:', arr)
    return arr

N = 3
a = 2.249 # 3
w = np.arccos((a - 1)/2)
b0 = 1/a
b1 = np.sqrt((a - 2)/a/a/(a - 1))
data_in = [w*180/np.pi, b0, b1]
#data_in = [23, 0.20, 0.16]
if N > 1:
    for i in range (2*(N - 1)):
        data_in.append(0)
#data_in = [0, 0.14961085948676142 , 0.13659693701379186, 0.113546733269272, 0, 0.0844215497293165, 0, 0, 0, 0, 0]
data_in = [40.37351966920139 , 0.3466225233198259 , 0.21587455604003164, 0.07341036725540125, -0.04725885602123077, 0.0006885692623692628, -0.03319]
print ('DL init = ', N, ',', a, ',', data_in, ',', func(data_in))

root = fsolve(func, data_in, full_output = False)
#print ('final =', [root[0], root[1], root[2], root[3], root[4]], ',', func([root[0], root[1], root[2], root[3], root[4]]))
print ('final =', N, ',', a, ',', root, ',', func(root))

print ('summary    , a, w, b0, b1, bi_real, bi_imag')
print ('solve_b0_bN,', a, ',', root[0], ',', root[1], ',', root[2], end='')
if N > 1:
    for i in range (2*(N - 1)):
        print (',', root[i + 3], end='')
print ()
