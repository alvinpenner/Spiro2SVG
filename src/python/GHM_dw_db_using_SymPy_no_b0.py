
# \Henon_Gonchenko\GHM_dw_db_using_SymPy.py
# https://docs.sympy.org/latest/guides/solving/solve-matrix-equation.html
# use SymPy to duplicate the functionality of 'GHM_dw_db_at_NS.py'
# solve the linear equation: M*sol = rhs
# this has been simplified to re-express C2, S2, b0 in terms of C1, S1, R

import sympy
from sympy.matrices import Matrix

print ('SymPy ver =', sympy.__version__)
print ()
R, C1, S1 = sympy.symbols('R, C1, S1')

M = Matrix([[ 0   ,  2*C1 - 2      , 2*R*C1 - 2,  0               ,  0],
            [ 2*S1,  R - 2 + R*C1  , 0         ,  R*C1 + R*(2*C1*C1 - 1) - 2 , -R*S1 + R*2*S1*C1],
            [ 0   , -R*S1          , 0         ,  R*S1 - R*2*S1*C1,  R*C1 + R*(2*C1*C1 - 1) - 2],
            [ 0   ,  0             , R*C1 - 1  , -2*(2*C1*C1 - 1) + 2*C1     ,  0],
            [ 0   ,  0             ,-R*S1      ,  0               , -2*(2*C1*C1 - 1) + 2*C1]])
print ('M =', M)
print ()
rhs = Matrix([2*C1/(R - 2), C1, -S1, 0, 0])
print ('rhs =', rhs)
print ()
#sol = M.solve(rhs)
#print ('sol[0] =', sol[0])
#print ()
#print ('sol[1] =', sol[1])

# attempt to factor output

sol_0 = (4*C1**4*R**2 + 2*C1**3*R**3 + 2*C1**3*R**2 - 16*C1**3*R + C1**2*R**3 + 4*C1**2*R**2*S1**2 - 4*C1**2*R**2 + 16*C1**2 + 2*C1*R**3*S1**2 + C1*R**3 - 6*C1*R**2*S1**2 - 7*C1*R**2 + 12*C1*R - 12*C1 - R**3*S1**2 + 2*R**2*S1**2 - 3*R**2 + 8*R - 4)/(4*C1**2*R**3*S1 - 12*C1**2*R**2*S1 + 8*C1**2*R*S1 - 2*C1*R**3*S1 + 6*C1*R**2*S1 - 4*C1*R*S1 - 2*R**3*S1 + 6*R**2*S1 - 4*R*S1)
print ('sol_0 =', sol_0)
print ()
print ('factor sol_0 =', sympy.factor(sol_0))
print ()

sol_1 = (-4*C1**2*R - C1*R**2 + 3*C1*R + 4*C1 + 2*R - 4)/(C1*R**3 - 3*C1*R**2 + 2*C1*R - R**3 + 3*R**2 - 2*R)
print ('sol_1 =', sol_1)
print ()
print ('factor sol_1 =', sympy.factor(sol_1))
print ()

# in sol_0 replace S1**2 = 1 - C1**2

sol_0 = (4*C1**4*R**2 + 2*C1**3*R**3 + 2*C1**3*R**2 - 16*C1**3*R + C1**2*R**3 + 4*C1**2*R**2*(1 - C1**2) - 4*C1**2*R**2 + 16*C1**2 + 2*C1*R**3*(1 - C1**2) + C1*R**3 - 6*C1*R**2*(1 - C1**2) - 7*C1*R**2 + 12*C1*R - 12*C1 - R**3*(1 - C1**2) + 2*R**2*(1 - C1**2) - 3*R**2 + 8*R - 4)/(4*C1**2*R**3*S1 - 12*C1**2*R**2*S1 + 8*C1**2*R*S1 - 2*C1*R**3*S1 + 6*C1*R**2*S1 - 4*C1*R*S1 - 2*R**3*S1 + 6*R**2*S1 - 4*R*S1)
print ('simplified factor sol_0 =', sympy.factor(sol_0))
print ()
