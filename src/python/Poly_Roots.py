
# find roots of a polynomial in R^2 from the map:
# z_n+1 = g10*z + g21*z*z*z_bar + g32*z*z*z*z_bar*z_bar
# calculate torus radius and phase shift for a 2D quintic map with cylindrical symmetry
# https://numpy.org/doc/2.2/reference/generated/numpy.polynomial.polynomial.Polynomial.html#numpy.polynomial.polynomial.Polynomial
# input data for gij: 'Normalize_Cubic_Map.py'
# see book Quintic Map, p. 27

import math
import numpy as np
from numpy.polynomial import Polynomial as Poly

#poly_roots = " 8 , 2.01 , NaN , NaN , 0.5 , 0.8717797887081346 , 0.9851727777823239 , -1.7178495427732 , 9.253024795056934 , -11.31079301357881 "
poly_roots = " 0 , -0.36 , 1.0202 , -0.1 , 0.20997342155426682 , 0.9778106634320061 , 0.4614630674651996 , -0.1169033251579869 , 0.3165186565249657 , -0.1617245583414001 "
poly_roots = " 1 , -0.36 , 1.021 , -0.1 , 0.20986718436432175 , 0.9782476240886823 , 0.4614230909783845 , -0.11683031978240993 , 0.3162837940793529 , -0.16094571165801336 "
poly_roots = " 2 , -0.36 , 1.022 , -0.1 , 0.20973455993273862 , 0.9787934891366136 , 0.46137327645451964 , -0.11673926950595609 , 0.3159904917241424 , -0.15997279541654347 "
poly_roots = " 3 , -0.36 , 1.023 , -0.1 , 0.20960212624244484 , 0.9793389818716814 , 0.46132363465305365 , -0.11664844772209401 , 0.3156975004166936 , -0.15900058748621143 "
#poly_roots = " 4 , -0.36 , 1.024 , -0.1 , 0.20946988283227083 , 0.9798841031140071 , 0.4612741648226246 , -0.11655785348371614 , 0.31540482520440655 , -0.15802908240757196 "

hdr = poly_roots.split(',')         # iT, parm a, parm b, parm c, gij etc
g10 = complex(float(hdr[4]), float(hdr[5]))
g21 = complex(float(hdr[6]), float(hdr[7]))
g32 = complex(float(hdr[8]), float(hdr[9]))
# coef style: C0 + C1*x + C2*x*x + C3*x*x*x + C4*x*x*x*x
coef = Poly([g10*g10.conjugate() - 1,
             g10*g21.conjugate() + g21*g10.conjugate(),
             g10*g32.conjugate() + g21*g21.conjugate() + g32*g10.conjugate(),
             g21*g32.conjugate() + g32*g21.conjugate(),
             g32*g32.conjugate()])
print ('coef = ', coef)
roots = coef.roots()
print ('\nroots')
for i in range(len(roots)):
    print ('test ', i, ',', roots[i], ',', coef.coef[0] + coef.coef[1]*roots[i] + coef.coef[2]*roots[i]*roots[i] + coef.coef[3]*roots[i]*roots[i]*roots[i] + coef.coef[4]*roots[i]*roots[i]*roots[i]*roots[i], end='')
    if roots[i].real < 0:
        print ('\t Bad Data: R^2 negative')
    elif abs(roots[i].imag) > 0.000001:
        print ('\t Bad Data: R^2 complex')
    else:
        print ('\t R = ', np.sqrt(roots[i]))

print ('\none pass starting at point (R, 0)')
R = np.sqrt(abs(roots[2]))      # assume root[2] is real
map10 = g10*R + g21*R*R*R + g32*R*R*R*R*R
print ('root 2,', hdr[0], ',', hdr[1], ',', hdr[2], ',', hdr[3], ',', R, ',', abs(map10), ',', map10.real, ',', map10.imag, ',,', abs(g10), ',', math.atan(g10.imag/g10.real)*180/np.pi, ',,', math.atan(map10.imag/map10.real)*180/np.pi)
R = np.sqrt(abs(roots[3]))      # assume root[3] is real
map10 = g10*R + g21*R*R*R + g32*R*R*R*R*R
print ('root 3,', hdr[0], ',', hdr[1], ',', hdr[2], ',', hdr[3], ',', R, ',', abs(map10), ',', map10.real, ',', map10.imag, ',,', abs(g10), ',', math.atan(g10.imag/g10.real)*180/np.pi, ',,', math.atan(map10.imag/map10.real)*180/np.pi)
