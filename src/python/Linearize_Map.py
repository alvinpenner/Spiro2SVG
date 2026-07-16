
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

DL = True
#data = [2.22, 2.17]             # Delayed Logistic
data = [1.02, 1.0202, 1.021, 1.022, 1.023, 1.024, 1.025, 1.026, 1.027, 1.028, 1.029]    # Henon
data = [1.0201, 1.0202, 1.0203, 1.0204, 1.0205, 1.0206, 1.0207, 1.0208, 1.0209, 1.0210, 1.0211, 1.0212, 1.0213, 1.0214, 1.0215, 1.0216, 1.0217, 1.0218, 1.0219, 1.0220]
data = [1.0106, 1.0107, 1.0108, 1.0109, 1.0110, 1.0111, 1.0112, 1.0113, 1.0114, 1.0115, 1.0116, 1.0117, 1.0118, 1.0119, 1.0120, 1.0121, 1.0122, 1.0123, 1.0124, 1.0125, 1.0126, 1.0127, 1.0128, 1.0129, 1.0130]
#data = [1.0271, 1.0272, 1.0273, 1.0274, 1.0275, 1.0276, 1.0277, 1.0278, 1.0279]
#data = [1.0106, 1.011, 1.012, 1.013, 1.014, 1.015, 1.016, 1.017, 1.018, 1.019, 1.020, 1.021, 1.022, 1.023, 1.024, 1.025, 1.026]
#data = [1.0249777]
data = [2.00, 2.01, 2.02, 2.03, 2.04, 2.05]
#data = [1.021, 1.0201, 1.0202, 1.0203, 1.0204, 1.0205]
#data = [1.0105572809000085]


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
        #Henon_a = -0.36 # -0.36 # -0.20
        Henon_R = -0.1
        Henon_b = data[i]               # scan beta
        Henon_a = (Henon_b - 1)*(Henon_b - 1 + 2*Henon_R)/Henon_R/Henon_R
        x0 = (Henon_b + 1 - np.sqrt((Henon_b + 1)*(Henon_b + 1) - 4*Henon_a*(Henon_R - 1)))/2/(Henon_R - 1) # stationary point
        y0 = x0
        print ('check Henon  x0,', Henon_b, ',', x0, ',', Henon_a - Henon_b*x0 - y0*y0 + Henon_R*x0*y0)
        print ('check Re/Im_V21,', (Henon_b - 1)*(Henon_R - 2)/2/Henon_R, ',', \
                                  -np.sqrt((2*Henon_R + (Henon_b - 1)*(Henon_R - 2))*(2*Henon_R - (Henon_b - 1)*(Henon_R - 2)))/2/Henon_R)
        print ('dRe_V21_dbeta  ,', -(Henon_R - 2)*(Henon_b - 1)/2/(2*(Henon_b - 1) + (3 - Henon_b)*Henon_R))
        print ('dIm_V21_dbeta  ,', -0.5*((Henon_R - 2)*(Henon_R - 2)*(Henon_b - 1)*(Henon_b - 1) + 4*(Henon_b - 1 + Henon_R)*Henon_R)/ \
                                    (2*(Henon_b - 1) + (3 - Henon_b)*Henon_R)/ \
                                    np.sqrt((2*Henon_R + (Henon_b - 1)*(Henon_R - 2))* \
                                            (2*Henon_R - (Henon_b - 1)*(Henon_R - 2))))
        Jac = np.array([[0, 1], [-Henon_b + Henon_R*y0, -2*y0 + Henon_R*x0]])

    # calculate uniform linear response (see 'fit_linear_response()' in Chua_y_vs_x.java)
    # see Book Chaos III, p. 50 for skew transform

    w, v = LA.eig(Jac)
    #print ('eig: \n', w)
    #print ('vec: \n', v)
    Cx10 = w[0].real                        # assume the eigenvalues are complex conjugate pairs
    Cy10 = w[0].imag
    Re_V21 = (v[1][0]/v[0][0]).real         # define the skew-transform matrix
    Im_V21 = (v[1][0]/v[0][0]).imag
    print ("\nJac =\n", Jac)
    print ('Cxy10 = ', Cx10, ',', Cy10, ',', abs(w[0]), ',', math.atan(Cy10/Cx10)*180/np.pi, ',', Re_V21, ',', Im_V21)
    U = np.array([[1, 1], [Re_V21 + Im_V21, Re_V21 - Im_V21]])
    Uinv = LA.inv(U)
    #print ("\nU =\n", U)
    print ("\nUinv =\n", Uinv)
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

    if DL:
        print ('analytical Cx,', -Log_a*(Re_V21 + Im_V21)/2/Im_V21, ',', -Log_a*2*Re_V21/2/Im_V21, ',', -Log_a*(Re_V21 - Im_V21)/2/Im_V21)
        print ('analytical Cy,',  Log_a*(Re_V21 + Im_V21)/2/Im_V21, ',',  Log_a*2*Re_V21/2/Im_V21, ',',  Log_a*(Re_V21 - Im_V21)/2/Im_V21)
    else:
        print ('analytical Cx,', (  -(Re_V21 + Im_V21)*(Re_V21 + Im_V21) + Henon_R*(Re_V21 + Im_V21))/2/Im_V21,
                            ',', (-2*(Re_V21 + Im_V21)*(Re_V21 - Im_V21) + Henon_R*2*Re_V21)/2/Im_V21,
                            ',', (  -(Re_V21 - Im_V21)*(Re_V21 - Im_V21) + Henon_R*(Re_V21 - Im_V21))/2/Im_V21)
        print ('analytical Cy = -Cx')
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

    # add Kuznetsov code for cubic model (assuming original g21 = 0)
    # see Delayed_Logistic_cubic_variable_c1.py for original code
    # for original, see 'Chua_2D_cubic_variable_c1.py' and 'transform_quartic.py'

    g10 = w[0]
    print ("g10            ,", g10)
    Cxy = np.array([complex(Cx[0], Cy[0]), complex(Cx[1], Cy[1]), complex(Cx[2], Cy[2])])
    g20 = 2*np.matmul(Cxy, np.array([complex(0.5, 0), complex(0, -0.5), complex(-0.5, 0)]))/2
    g11 = np.matmul(Cxy, np.array([complex(1.0, 0), complex(0,  0.0), complex( 1.0, 0)]))/2
    g02 = 2*np.matmul(Cxy, np.array([complex(0.5, 0), complex(0,  0.5), complex(-0.5, 0)]))/2
    print ("g20 g11 g02    ,", g20, ",", g11, ",", g02)

    # optional aside for theoretical calc of gij

    g20_R = -complex (Re_V21 - Im_V21, Re_V21 + Im_V21)
    g11_R = Re_V21*complex (1, -1)
    g02_R =  complex (Re_V21 + Im_V21, Re_V21 - Im_V21)
    g20_1 =  complex ( Re_V21*Re_V21 - 2*Re_V21*Im_V21 - Im_V21*Im_V21,  Re_V21*Re_V21 + 2*Re_V21*Im_V21 - Im_V21*Im_V21)
    g11_1 =  complex (-Re_V21*Re_V21 - Im_V21*Im_V21, Re_V21*Re_V21 + Im_V21*Im_V21)
    g02_1 =  complex (-Re_V21*Re_V21 - 2*Re_V21*Im_V21 + Im_V21*Im_V21, -Re_V21*Re_V21 + 2*Re_V21*Im_V21 + Im_V21*Im_V21)
    if DL:
        print ("analytical gij ,", -Log_a*g20_R/2/Im_V21, ',', -Log_a*g11_R/2/Im_V21, ',', -Log_a*g02_R/2/Im_V21)
        print ("actual gij     ,", Log_a*g20_R/2/Im_V21*Log_a*g11_R/2/Im_V21, ',', \
                                   np.abs(Log_a*g11_R/2/Im_V21)*np.abs(Log_a*g11_R/2/Im_V21), ',', \
                                   np.abs(Log_a*g02_R/2/Im_V21)*np.abs(Log_a*g02_R/2/Im_V21))
        print ("test   gij     ,",-Log_a*Log_a*Re_V21*complex (Re_V21, Im_V21)/2/Im_V21/Im_V21, ',', \
                                   Log_a*Log_a*Re_V21*Re_V21/2/Im_V21/Im_V21, ',', \
                                   Log_a*Log_a*(Re_V21*Re_V21 + Im_V21*Im_V21)/2/Im_V21/Im_V21)
        print ("test g20*g11   ,", g20*g11/(Log_a*Log_a/2/Im_V21/Im_V21))
        print ("test g11*_g11  ,", g11*g11.conjugate()/(Log_a*Log_a/2/Im_V21/Im_V21))
        print ("test g02*_g02  ,", g02*g02.conjugate()/(Log_a*Log_a/2/Im_V21/Im_V21))
        print ("cof  g20*g11   ,", (g10.conjugate() - 3 + 2*g10)/2/(g10*g10 - g10)/(g10.conjugate() - 1))
        print ("cof  g11*_g11  ,", 1/g10.conjugate()/(g10 - 1))
        print ("cof  g02*_g02  ,", 1/2/(g10*g10 - g10.conjugate()))
    else:
        print ("analytical gij ,", (g20_1 + Henon_R*g20_R)/2/Im_V21, ',', (g11_1 + Henon_R*g11_R)/2/Im_V21, ',', (g02_1 + Henon_R*g02_R)/2/Im_V21)
        V21_Re = (Henon_R - 2)*(Henon_b - 1)/2/Henon_R         # only at bifurc
        theta = math.acos(V21_Re)               # first order response at bifurc
        V21_Im = math.sin(theta)                # only at bifurc
        #dtheta_dbeta = -V21_Re*(2*Henon_R + Henon_b - 1)/(Henon_R*(Henon_b - 3) - 2*(Henon_b - 1))/V21_Im
        #d_mu_dbeta = (Henon_b - 1 + Henon_R)/(2*(Henon_b - 1) - Henon_R*(Henon_b - 3))
        dtheta_dbeta = V21_Re*(Henon_R - 2 + V21_Re)/(Henon_R - 2)/(1 - V21_Re)/V21_Im
        d_mu_dbeta = (Henon_R - 2 + 2*V21_Re)/2/(Henon_R - 2)/(1 - V21_Re)
        #print ('dtheta/d_mu_dbeta ,', dtheta_dbeta, ',', d_mu_dbeta)
        theo_g20 = (g10 - Henon_R)*g10*complex (1, 1)/2/Im_V21
        theo_g11 = (-g10*g10.conjugate() + Henon_R*(g10 + g10.conjugate())/2)*complex (1, -1)/2/Im_V21
        theo_g02 = (-g10.conjugate() + Henon_R)*g10.conjugate()*complex (1, 1)/2/Im_V21
        print ('theo  gij      ,', theo_g20, ',', theo_g11, ',', theo_g02)

        theo_c1  = theo_g20*theo_g11*(g10.conjugate() - 3 + 2*g10)/2/g10/(g10 - 1)/(g10.conjugate() - 1)
        theo_c1 += theo_g11*theo_g11.conjugate()/g10.conjugate()/(g10 - 1)
        theo_c1 += theo_g02*theo_g02.conjugate()/2/(g10*g10 - g10.conjugate())
        theo_c1 *= g10.conjugate()/abs(g10)                              # compensate for exp(-i*theta)
        print ('theo  c1 (1)    ,', theo_c1)

        theo_c1 = theo_g20*theo_g11*g10.conjugate()*(2*g10 - 1)*(g10*g10 + g10 + 1) \
                - 2*g10*theo_g11*theo_g11.conjugate()*(g10*g10 + g10 + 1) \
                - g10*theo_g02*theo_g02.conjugate()
        theo_c1  = -theo_c1/2/(g10*g10*g10 - 1)
        theo_c1 *= g10.conjugate()/abs(g10)                              # compensate for exp(-i*theta)
        print ('theo  c1 (5)    ,', theo_c1)

        # analytical components of c1

        #print ('g2011_org ,', (Henon_R*V21_Re - 1)*(g10 - Henon_R)*(2*g10 - 1)*(g10*g10 + g10 + 1))
        g2011Re = (Henon_R*V21_Re - 1)*(2*V21_Re + 1)*(8*V21_Re*V21_Re*V21_Re - 6*V21_Re \
                                                     - 2*V21_Re*V21_Re + 1 \
                                                     - 4*Henon_R*V21_Re*V21_Re + 2*Henon_R + Henon_R*V21_Re)
        g2011Im = (Henon_R*V21_Re - 1)*(2*V21_Re + 1)*V21_Im*(8*V21_Re*V21_Re - 2 \
                                                              - 2*V21_Re \
                                                              - 4*Henon_R*V21_Re \
                                                              + Henon_R)
        print ('g2011_anal,', g2011Re, ',', g2011Im)
        #print ('g1111_org ,', -(Henon_R*V21_Re - 1)*(Henon_R*V21_Re - 1)*2*g10*(g10*g10 + g10 + 1))
        g1111Re = -2*(Henon_R*V21_Re - 1)*(Henon_R*V21_Re - 1)*(4*V21_Re*V21_Re*V21_Re + 2*V21_Re*V21_Re - 2*V21_Re - 1)
        g1111Im = -2*(Henon_R*V21_Re - 1)*(Henon_R*V21_Re - 1)*V21_Im*(4*V21_Re*V21_Re + 2*V21_Re)
        print ('g1111_anal,', g1111Re, ',', g1111Im)
        #print ('g0202_org ,', -(1 + Henon_R*Henon_R - 2*Henon_R*V21_Re)*g10)
        g0202Re = -(Henon_R*Henon_R - 2*Henon_R*V21_Re + 1)*V21_Re
        g0202Im = -(Henon_R*Henon_R - 2*Henon_R*V21_Re + 1)*V21_Im
        print ('g0202_anal,', g0202Re, ',', g0202Im)
        sumRe = g2011Re + g1111Re + g0202Re
        sumIm = g2011Im + g1111Im + g0202Im
        print ('sumRe     ,', sumRe, ',', sumIm)
        u3_1Re = 4*V21_Re*V21_Re*V21_Re - 3*V21_Re - 1
        u3_1Im = -V21_Im*(4*V21_Re*V21_Re - 1)
        print ('u3_1Re/Im ,', u3_1Re, ',', u3_1Im)
        c1_Re = u3_1Re*sumRe - u3_1Im*sumIm         # multiply u3_1*sumRe_Im
        c1_Im = u3_1Re*sumIm + u3_1Im*sumRe
        # calculate ubar*c1_Im_Re                   # compensate for exp(-i*theta)
        print ('ubar*c1_Im/Re ,', V21_Re*c1_Im - V21_Im*c1_Re, ',', V21_Re*c1_Re + V21_Im*c1_Im)
        print ('ubar*c1_Im/Re ,', (V21_Re*c1_Im - V21_Im*c1_Re)/(V21_Re*c1_Re + V21_Im*c1_Im))

        # rev 2 of sumRe_Im
        sumRe_2  = Henon_R*Henon_R*(-8*V21_Re*V21_Re*V21_Re*V21_Re*V21_Re - 12*V21_Re*V21_Re*V21_Re*V21_Re +  2*V21_Re*V21_Re*V21_Re +  7*V21_Re*V21_Re +   V21_Re)
        sumRe_2 +=         Henon_R*(16*V21_Re*V21_Re*V21_Re*V21_Re*V21_Re + 20*V21_Re*V21_Re*V21_Re*V21_Re +  2*V21_Re*V21_Re*V21_Re -  8*V21_Re*V21_Re - 8*V21_Re - 2)
        sumRe_2 +=                                                        - 16*V21_Re*V21_Re*V21_Re*V21_Re - 12*V21_Re*V21_Re*V21_Re + 10*V21_Re*V21_Re + 7*V21_Re + 1
        sumIm_2  = V21_Im*Henon_R*Henon_R*(-8*V21_Re*V21_Re*V21_Re*V21_Re - 12*V21_Re*V21_Re*V21_Re -  2*V21_Re*V21_Re +   V21_Re - 1)
        sumIm_2 +=         V21_Im*Henon_R*(16*V21_Re*V21_Re*V21_Re*V21_Re + 20*V21_Re*V21_Re*V21_Re + 10*V21_Re*V21_Re + 2*V21_Re - 1)
        sumIm_2 +=                 V21_Im*(                                -16*V21_Re*V21_Re*V21_Re - 12*V21_Re*V21_Re + 2*V21_Re + 1)
        print ('sumRe rev ,', sumRe_2, ',', sumIm_2)

    # calculate (resonant) g21new as per Kuznetsov (after eliminating hij)
    
    g21new = g20*g11*(g10.conjugate() - 3 + 2*g10)/2/g10/(g10 - 1)/(g10.conjugate() - 1)
    g21new += g11*g11.conjugate()/g10.conjugate()/(g10 - 1)
    g21new += g02*g02.conjugate()/2/(g10*g10 - g10.conjugate())
    g21new *= g10.conjugate()/abs(g10)                              # compensate for exp(-i*theta)
    if DL:
        print (data[i], ',', abs(g10), ',', '%f' % (np.arctan2(g10.imag, g10.real)*180/np.pi), ',', g21new.real, ',', g21new.imag, ',', g21new.imag/g21new.real)
    else:
        print ('Henon alpha, beta, R, abs(g10), theta, g21_Re, g21_Im, g21_Im/Re, dtheta/dbeta, d_mu/dbeta')
        print (Henon_a, ',', data[i], ',', Henon_R, ',', abs(g10), ',', '%f' % (np.arctan2(g10.imag, g10.real)*180/np.pi), ',', g21new.real, ',', g21new.imag, ',', g21new.imag/g21new.real, ',', dtheta_dbeta, ',', d_mu_dbeta)

