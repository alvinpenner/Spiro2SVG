
# consider either a Delayed Logistic Map (Aronson) with parameter 'a'
# or a Generalized Henon Map (Gonchenko_Kuznetsov) with parameters 'alpha, beta, R'
# transform real Cxy coeff into uniform linear form
# (this is a cleaned-up version of 'Delayed_Logistic_cubic_variable_c1.py')
# March 30, 2025

import numpy as np
from numpy import linalg as LA
import math

print ('cubic header: iT, alpha, beta, gamma, NaN, NaN, NaN, NaN,Cx[0],Cx[1],Cx[2],Cx[3],Cx[4],Cx[5],Cx[6],Cx[7],Cx[8],Cx[9],Cy[0],Cy[1],Cy[2],Cy[3],Cy[4],Cy[5],Cy[6],Cy[7],Cy[8],Cy[9]')
print ()

DL = not True
#data = [2.22, 2.17]             # Delayed Logistic
data = [1.022, 1.0249777]             # Henon

for i in range(len(data)):
    print ('______________________________________________________')
    if DL:
        # x_n+1 = y_n
        # y_n+1 = a*y_n*(1 - x_n)
        Log_a = data[i]
        x0 = (Log_a - 1)/Log_a          # stationary point
        y0 = x0
        print ('check Logistic x0:', Log_a, ',', x0, ',', Log_a*x0*(1 - x0))
        Jac = np.array([[0, 1], [-Log_a*y0, Log_a*(1 - x0)]])
    else:                               # assume Henon
        #x_n+1 = y_n
        #y_n+1 = Henon_a - Henon_b*x_n - y_n*y_n + Henon_R*x_n*y_n
        Henon_a = -0.36
        Henon_b = data[i]               # scan beta
        Henon_R = -0.1
        x0 = (Henon_b + 1 - np.sqrt((Henon_b + 1)*(Henon_b + 1) - 4*Henon_a*(Henon_R - 1)))/2/(Henon_R - 1) # stationary point
        y0 = x0
        print ('check Henon x0:', Henon_b, ',', x0, ',', Henon_a - Henon_b*x0 - y0*y0 + Henon_R*x0*y0)
        Jac = np.array([[0, 1], [-Henon_b + Henon_R*y0, -2*y0 + Henon_R*x0]])

    # calculate uniform linear response (see 'fit_linear_response()' in Chua_y_vs_x.java)
    # see Book Chaos III, p. 50 for skew transform

    w, v = LA.eig(Jac)
    #print ('eig: \n', w)
    #print ('vec: \n', v)
    Cx10 = w[0].real                        # assume the eigenvalues are complex conjugate pairs
    Cy10 = w[0].imag
    print ('Cxy10 = ', Cx10, ',', Cy10, ',', abs(w[0]), ',', math.atan(Cy10/Cx10)*180/np.pi)
    Re_V21 = (v[1][0]/v[0][0]).real         # define the skew-transform matrix
    Im_V21 = (v[1][0]/v[0][0]).imag
    print ("\nJac =\n", Jac, ',', Re_V21, ',', Im_V21)
    U = np.array([[1, 1], [Re_V21 + Im_V21, Re_V21 - Im_V21]])
    Uinv = LA.inv(U)
    #print ("\nU =\n", U)
    #print ("\nUinv =\n", Uinv)
    print ("transform Jac:\n", np.matmul(Uinv ,np.matmul(Jac, U)))
    print ()

    # calculate transformed quadratic response in a format suitable
    # for Chua_Simul_3().hdr (Java)
    # Cxy = np.array([complex(Cx20, Cy20), complex(Cx11, Cy11), complex(Cx02, Cy02)])

    if DL:
        Cx = -Log_a*np.array([U[0][0]*U[1][0], U[0][0]*U[1][1] + U[0][1]*U[1][0], U[0][1]*U[1][1]])
    else:                                               # assume Henon
        Cx = np.array([-U[1][0]*U[1][0]   + Henon_R*U[0][0]*U[1][0],
                       -2*U[1][0]*U[1][1] + Henon_R*(U[0][0]*U[1][1] + U[0][1]*U[1][0]),
                       -U[1][1]*U[1][1]   + Henon_R*U[0][1]*U[1][1]])
    Cy = Uinv[1][1]*Cx
    Cx = Uinv[0][1]*Cx

    print ('    private static String hdr = "', end='') # create a Java header
    if DL:
        print (i, ',', Log_a, ', NaN, NaN, NaN, NaN, NaN, NaN', end='') # Cx, Cy summary for 'Chua_Simul_3'
    else:                                               # assume Henon
        print (i, ',', Henon_a, ',', Henon_b, ',', Henon_R, ', NaN, NaN, NaN, NaN', end='')
    print (" ,", 0, ",", Cx10, ",", -Cy10, end='')      # insert first-order x response
    for k in range(len(Cx)):
        print (',', Cx[k], end='')
    print (', 0, 0, 0, 0', end='')                      # pad the header to be cubic
    print (', 0, 0, 0, 0, 0', end='')                   # pad the header to be quartic
    print (" ,", 0, ",", Cy10, ",", Cx10, end='')       # insert first-order y response
    for k in range(len(Cy)):
        print (',', Cy[k], end='')
    print (', 0, 0, 0, 0', end='')                      # pad the header to be cubic
    print (', 0, 0, 0, 0, 0', end='')                   # pad the header to be quartic
    print ('";')
    print ()
