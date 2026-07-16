
# calculate dw/da at a bifurc for delayed logistic map
# see Fourier Analysis of Henon Gonchenko map, Jan 3, 2026, pp 10,11
# see Fourier Analysis of Delayed Logistic map, Dec 8, 2025, pp 6, 14

import numpy as np
from numpy import linalg as LA
import math

a = 2
# ..................................
b0    = 1/a
w     = np.arccos(0.5)
C1    = np.cos(w)
C2    = np.cos(2*w)
S1    = np.sin(w)
S2    = np.sin(2*w)
print ('a,', a)
print ('w_b0_C1_S1,', w*180/np.pi, ',', b0, ',', C1, ',', S1)
print ()

M = np.array([[ 0           , 2*a*b0 + 1 - a, 2*a*C1, 0          , 0],
              [-S1 - a*b0*S1, a + a*C1      , 0     , a*C1 + a*C2, -a*S1 + a*S2],
              [ C1 - a*b0*C1, -a*S1         , 0     , a*S1 - a*S2,  a*C1 + a*C2],
              [ 0           , 0             , a*C1  , C2 + a*b0*C2 + a*b0 - a, -S2 + a*b0*S2],
              [ 0           , 0             ,-a*S1  , S2 - a*b0*S2           ,  C2 + a*b0*C2 + a*b0 - a]])
print ('M =', M)
bx = np.array([ -b0*b0 + b0, -b0 + 1 - b0*C1, b0*S1, 0, 0])
print ('b =', bx)
db = LA.solve(M, bx)
print ()
print ('sol =', db)
