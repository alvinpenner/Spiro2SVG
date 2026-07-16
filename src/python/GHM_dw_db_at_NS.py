
# calculate dw/dbeta at a NS bifurc for a Generalized Henon map
# see Fourier Analysis of Henon Gonchenko map, Jan 3, 2026, pp 10,11
# see Gonchenko p. 7
# N.B. note this is a new and improved version modified on May 30, 2026
# see 'Fourier Analysis of Henon Gonchenko Map' Jan 3, 2026, p.16

import numpy as np
from numpy import linalg as LA
import math

#beta = 1.08 # 0.9739453 # 0.95238
alpha = -0.20 # -0.36 # -0.20
R     = -0.1
#alpha = (beta - 1)*(beta - 1 + 2*R)/R/R
beta = 1 - R + R*np.sqrt(1 + alpha)        # Gonchenko Sec. 3.1.1
# ..................................
b0    = (beta - 1)/R
w     = np.arccos(b0*(R - 2)/2)
C1    = np.cos(w)
C2    = np.cos(2*w)
#C3    = np.cos(3*w)
S1    = np.sin(w)
S2    = np.sin(2*w)
#S3    = np.sin(3*w)
print ('a_b_R,', alpha, ',', beta, ',', R)
print ('w_b0_C1_S1_C2_S2,', w*180/np.pi, ',', b0, ',', C1, ',', S1, ',', C2, ',', S2)
print ()
M = np.array([[ 0   ,  2*C1 - 2      , 2*R*C1 - 2,  0               ,  0],
              [ 2*S1,  R - 2 + R*C1  , 0         ,  R*C1 + R*C2 - 2 , -R*S1 + R*S2],
              [ 0   , -R*S1          , 0         ,  R*S1 - R*S2     ,  R*C1 + R*C2 - 2],
              [ 0   ,  0             , R*C1 - 1  , -2*C2 + 2*C1     ,  0],
              [ 0   ,  0             ,-R*S1      ,  0               , -2*C2 + 2*C1]])
print ('M =', M)
rhs = np.array([ b0, C1, -S1, 0, 0])
print ('rhs =', rhs)
db = LA.solve(M, rhs)
print ()
print ('d/dbeta =,', alpha, ',', beta, ',', R, ',', db[0], ',', db[1], ',', db[2], ',', db[3], ',', db[4])
#print ('M*x =,', np.matmul(M, db))

# theory 'Fourier Analysis of Henon Gonchenko Map', p. 17
# NOTE that the sign of (num_b0, num_1) has been reversed because of an unfortunate error
# in the right hand side of the Maxima matrix in 'Henon_Gonchenko_w_rev_2.txt'
# the current Python matrix right hand side above is correct

den    = -4*S1*S1*(R*R - R)*(2*C1 + 1)*(C1 - 1)
num_b0 = -S1*(2*R*R*R*(2*C1 - 1)*(C1 + 1) + 4*R*C1 - 8*R*R*C1*C1)
num_1  =   2*R*R*S1*(C1 - 1)*(4*C1*C1 + 4*C1 - 1) \
         - 2*R*S1*(C1 - 1)*(8*C1*C1 + 6*C1 + 3) \
         + 4*S1*(C1 - 1)*(4*C1 + 1)
print ()
print ('theory den num_b0 num_1 final,', den, ',', num_b0, ',', num_1, ',', (b0*num_b0 + num_1)/den)

# theory revised (page 23 of Fourier Analysis of GHM)

den = -4*S1*R*(R - 1)*(R - 2)*(2*C1 + 1)*(C1 - 1)
num = R*R*R*( -4*C1*C1 - 6*C1 + 2) \
    +   R*R*(-16*C1*C1*C1 +  4*C1*C1 + 26*C1 + 2) \
    +     R*( 32*C1*C1*C1            - 24*C1 - 16) \
                          - 32*C1*C1 + 24*C1 + 8
print ('theory reduced form     final,', num, ',', den, ',', num/den)

#Maxima_den = R*R*(C1*C1*(4*S1*S2 + 8*S1*S1) + C1*(-4*S1*S2 - 4*C2*S1*S1) - 4*C2*S1*S1) \
#           + R*(C1*(-4*S1*S2 - 12*S1*S1) + 4*S1*S2 + (8*C2 + 4)*S1*S1)
#Maxima_b0  = R*R*R*(C1*C1*(S2 - S1) + S1*S1*S2 + C1*S2 - S1*S1*S1  + C2*S1) + 2*R*S2 + R*R*(-S2 + C1*(2*S1 - 3*S2) + (- C2 - 1)*S1)

#theory_1_R2 = - 2*R*R*S1*(C1 - 1)*(4*C1*C1 + 4*C1 - 1)
#Maxima_1_R2 = R*R*(C1*(-2*S1*S1*S2 + 2*S1*S1*S1  - 4*C2*S1) + C1*C1*(2*S2 + 2*S1) + 2*S1*S1*S2 + C1*C1*C1*(2*S1 - 2*S2) - 2*S1*S1*S1)
#print (theory_1_R2, ',', Maxima_1_R2)
#print ('Maxima den num_b0 num_1      ,', Maxima_den, ',', Maxima_b0, ',', Maxima_1)
#theory_1_R1 = 2*R*S1*(C1 - 1)*(8*C1*C1 + 6*C1 + 3)
#Maxima_1_R1 = R*(C1*C1*(2*S2 - 12*S1) + C1*(6*C2*S1 - 2*S2) + 6*C2*S1)
#print (theory_1_R1, ',', Maxima_1_R1)
#Maxima_1_1 = (-8*C2 - 4)*S1 + 12*C1*S1
#theory_1_1 = -4*S1*(C1 - 1)*(4*C1 + 1)
#print (theory_1_1, ',', Maxima_1_1)
