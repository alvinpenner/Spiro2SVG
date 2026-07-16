
# fit a sine/cos wave to Henon homogeneous sine-cos model at frequency 0, 1, 2
# assuming a given b0, and b1 ... bN real
# see: https://docs.scipy.org/doc/scipy-1.16.2/reference/generated/scipy.optimize.fsolve.html

# see attachment 'Fourier Analysis of Henon_sine_cos_alpha map', Jan 19, 2026
# Henon model parameter: alpha
# adjustable parameters are: w, b1 ... bN (assumed real), initial codition b0
# internal function  parameters are: x[0] - x[N]
# map to : [w, b1, b2]

import numpy as np
from scipy.optimize import fsolve

def func(x):
    deg = len(x) - 1                                            # highest Fourier frequency
    rad = x[0]*np.pi/180
    #print ('x len =', len(x), ',', deg)
    fw = np.cos(rad)                                            # cos(w)
    fb = np.zeros(deg + 1)                                      # array of coeff bi
    fb[0] = b0                                                  # bo (known constant)
    for i in range (1, deg + 1):
        fb[i] = x[i]                                            # b1 (internal variable)
    #print ('internal fw, fb[i]:', fw, ',', fb)
    wrap = np.zeros(2*deg + 1)                                  # row of b[j] from j = -N to N
    for j in range (deg + 1):
        wrap[j]       = fb[deg - j]
        wrap[j + deg] = fb[j]
    #print ('wrap           [i]:', fw, ',', wrap)

    M = np.zeros((deg + 1, 2*deg + 1))                          # array of coeff of constraints
    for j in range (deg + 1):
        M[j][deg - j] += fb[j]
        M[j][deg + j] += fb[j]
        M[j][deg]     -= 2*np.cos(alpha)*fb[j]
        for k in range (j, 2*deg + 1):
            M[j][deg] -= np.sin(alpha)*wrap[k]*wrap[2*deg - k + j]
            #print (j, ',', k, ',', 2*deg - k + j)
    #print ('M =', M, '\n')

    ret = np.zeros(deg + 1)                                     # array of real returns
    for j in range (deg + 1):
        for k in range (2*deg + 1):
            #ret[j] += M[j][k]*np.cos(j*rad)
            ret[j] += M[j][k]*np.cos((deg - k)*rad)
            #print ('ret,', j, ',', k, ',', ret[j])
    #print ('internal init', alpha, ',', rad, ',', b0, ',', fb[1], ',', fb[2], ',', np.cos(alpha), ',', np.sin(alpha))
    #print ('internal ret:', ret)
    return ret

N = 5
b0 = 0.054                                              # arbitrary input
H_cos = 0.24
alpha = np.arccos(H_cos)
print (H_cos, ',', b0, ',', alpha*180/np.pi)
w1 = np.arccos(np.cos(alpha) + b0*np.sin(alpha))*180/np.pi      # first-order w
#w1 = 0
data_in = [w1, np.sqrt(-b0*b0/2 + b0*(1 - np.cos(alpha))/np.sin(alpha))]
if N > 1:
    for i in range (2, N + 1):
        data_in.append(0)
before = len(data_in)
#print ('data_in,', data_in)
data_in = [71.76889567698953 , 0.022301514293393433, -0.0003364977971672534, 0.003955685643386366, -0.2024705622422827, -0.006207]
after = len(data_in)
if (before != after):
    print ("\nINPUT ERROR: len(data_in) has changed", before, 'to', after)
    exit()
print ('H_sin_cos = ', N, ',', H_cos, ',', b0, ',', data_in, ',', func(data_in))

root = fsolve(func, data_in, full_output = False)
#print ('final =', [root[0], root[1], root[2], root[3], root[4]], ',', func([root[0], root[1], root[2], root[3], root[4]]))
print ('final     = ', N, ',', H_cos, ',', b0, ',', root, ',', func(root))

print ('summary    , H_cos, b0, w1, w, b1, b2')
print ('solve_b0_bN,', H_cos, ',', b0, ',', w1, ',', root[0], ',', root[1], end='')
if N > 1:
    for i in range (2, len(data_in)):
        print (',', root[i], end='')
print ()

#   double check analytically for N = 2

if N == 20:
    print ("\n......................................")
    print ('analytical check for N = 2, w =', root[0])
    wrad = root[0]*np.pi/180
    b2 = -b0 + (np.cos(wrad) - np.cos(alpha))/np.sin(alpha)
    b1_temp = (np.cos(2*wrad) - np.cos(alpha))/np.sin(alpha) - b0
    b1 = np.sqrt(2*b2*b1_temp)
    print ('b1_b2  ,', b1, ',', b2)
    print ('check w,', 2*b0*(1 - np.cos(alpha)) - np.sin(alpha)*(b0*b0 + 2*b1*b1 + 2*b2*b2))

    bsin = np.cos(alpha) + b0*np.sin(alpha)
    rhs = 2*np.cos(2*wrad)*np.cos(wrad) - 2*(np.cos(2*wrad) + np.cos(wrad))*bsin + 2*bsin*bsin
    b1 = np.sqrt(rhs)/np.sin(alpha)
    rhs = np.cos(wrad)*np.cos(wrad) - 2*bsin*np.cos(wrad) + bsin*bsin
    b2 = np.sqrt(rhs)/np.sin(alpha)
    print ('test b1 b2,', b1, ',', b2)
    C3 = 8                              # C3*cos^3 + C2*cos^2 + C1*cos + C0 = 0
    C2 = -8*bsin + 2
    C1 = -4 - 8*bsin
    C0 = -2*b0*np.sin(alpha)*(1 - np.cos(alpha)) + b0*b0*np.sin(alpha)*np.sin(alpha) + 4*bsin + 6*bsin*bsin
    print ('test cubic w,', C3*np.cos(wrad)*np.cos(wrad)*np.cos(wrad) + C2*np.cos(wrad)*np.cos(wrad) + C1*np.cos(wrad) + C0)
    print ()
    coef = [C3, C2, C1, C0]
    print ('coef = ', coef)
    roots = np.roots(coef)
    #print ('\nroots :', roots)
    print ('roots')
    for i in range(len(roots)):
        if roots[i].imag == 0 and abs(roots[i].real) < 1:
            print (i, ',', roots[i].real, ',', roots[i].imag, ', w =', np.arccos(roots[i].real)*180/np.pi, ', test =', C3*roots[i]*roots[i]*roots[i] + C2*roots[i]*roots[i] + C1*roots[i] + C0)
        else:
            print (i, ',', roots[i].real, ',', roots[i].imag)
    # calculate dw_db0
    # see: attachment 'Fourier Analysis of Henon, Mar 11, 2026', p.10
    dw_db0 = -(4*np.cos(alpha) + 1)/(4*np.cos(alpha) + 2)
    print ('slope dw_db0,', H_cos, ',', dw_db0*180/np.pi)
