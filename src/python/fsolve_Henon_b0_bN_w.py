
# fit a sine/cos wave to Gonchenko GHM Henon model: frequency 0, 1, 2, 3, 4, 5
# see: https://docs.scipy.org/doc/scipy-1.16.2/reference/generated/scipy.optimize.fsolve.html

# see attachment 'Fourier Analysis of Henon Gonchenko map', Jan 3, 2026
# original (complex) parameters are: b0, b1, b2 (complex), w (real)
# assume origin of time is arbitrary, set b1 real
# internal function  parameters are: x[0] - x[4]
# map to : [w, b0, b1, b2R, b2I, b3R, b3I, ...] same mapping as for Delayed Logistic

import numpy as np
from scipy.optimize import fsolve

def fn_w(temp_beta, temp_R):
    return np.arccos((temp_beta - 1)*(temp_R - 2)/2/temp_R)

def fn_b0(temp_beta, temp_R):
    return (temp_beta - 1)/temp_R

def fn_b1(temp_alpha, temp_beta, temp_R):
    return np.sqrt((temp_alpha - fn_b0(temp_beta, temp_R)*(1 + temp_beta) + fn_b0(temp_beta, temp_R)*fn_b0(temp_beta, temp_R)*(temp_R - 1))/(2*temp_beta - temp_R*(temp_beta - 1)))

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
    M[0][deg] = alpha
    for j in range (deg + 1):
        M[j][deg - j] -= fb[j]
        M[j][deg + j] -= beta*fb[j]
        for k in range (j, 2*deg + 1):
            M[j][k]   += R*wrap[k]*wrap[2*deg - k + j]
            M[j][deg] -= wrap[k]*wrap[2*deg - k + j]
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

N     = 6
R     = -0.1
alpha = -0.36 # -0.36 # -0.20
beta  =  1.0240
#beta = 1 - R + R*np.sqrt(1 + alpha)
'''
# scan beta to calculate w (fundamental freq only)
print ('fundamental frequency w,', alpha, ',', beta, ',', R)
print ('beta, w, b0, b1')
for i in range (61):
    scan = beta + 0.0002*i
    print (scan, ',', fn_w(scan, R)*180/np.pi, ',', fn_b0(scan, R), ',', fn_b1(alpha, scan, R))
# end of scan
'''
data_in = [fn_w(beta, R)*180/np.pi, fn_b0(beta, R), fn_b1(alpha, beta, R)]
if N > 1:
    for i in range (2*(N - 1)):
        data_in.append(0)
before = len(data_in)
#data_in[0] = -data_in[0]
#data_in = [75.35562750743725 , -0.27772645812224633 , 0.40336269569838623, 0.05933752093042137, -0.006625433326948167, -0.01795378700728204, 0.0005602612734979972, -0.10498147167105115, 0.021566541817140617, 0.06792773879064566, -0.014383422826146436]
after = len(data_in)
if (before != after):
    print ("\nINPUT ERROR: len(data_in) has changed from", before, 'to', after)
    exit()
print ('Henon init = ', N, ',', alpha, ',', beta, ',', R, ',', data_in, ',', func(data_in))

root = fsolve(func, data_in, full_output = False)
#print ('final =', [root[0], root[1], root[2], root[3], root[4]], ',', func([root[0], root[1], root[2], root[3], root[4]]))
print ('final =', N, ',', alpha, ',', beta, ',', R, ',', root, ',', func(root))

print ('summary    , alpha, beta, R, w, b0, b1, bi_real, bi_imag')
print ('solve_b0_bN,', alpha, ',', beta, ',', R, ',', root[0], ',', root[1], ',', root[2], end='')
if N > 1:
    for i in range (2*(N - 1)):
        print (',', root[i + 3], end='')
print ()
