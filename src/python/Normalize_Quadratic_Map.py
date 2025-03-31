
# read quartic or quintic map as real coefficients Cxy
# (possibly from a curve-fit as in 'Chua_nonlinear_response_scatter.py')
# (or as output from an analytical map, as in 'Linearize_Map.py')
# assuming the linear coeff have already been made uniform
# convert real Cxy to complex coefficients (mu, gij quad, gij cubic, gij quart, gij quintic)
# analytically remove quadratic coeff, keeping overflow up to quintic
# see attached sheets: Normal Form of Quartic Map, Sep. 4, 2024
#                    : Final Quintic Map, Oct 30, 2024, p.10

# conversions from complex gij to real Cxy, and back, are consistent with 'expand_z.py'
# after transposing and multiplying by appropriate column and row scales

# (this is a cleaned-up version of 'transform_quartic.py')
# March 30, 2025

import numpy as np
import sys

hdr = "1 , -0.36 , 1.0249777 , -0.1 , NaN, NaN, NaN, NaN , 0 , 0.209340771971716 , -0.9804167097780551, -0.7825746941328198, 0.9143655752339029, -0.2638931501893879, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9804167097780551 , 0.209340771971716, 0.7825746941328199, -0.914365575233903, 0.26389315018938797, 0, 0, 0, 0, 0, 0, 0, 0, 0";
hdr = " 748319900, 99.98, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.9825013637459529 , -0.1864677328353816, 2.776092308471567, 5.466422132156249, -10.94611227200229, 199.69552085384856, 20.200712288857567, 420.82882397251774, 193.5610960875177, 46400.518972194746, -16062.84643487343, 26286.867824220077, -34342.71181116979, -12440.582713416536, 665452.9007178544, 1034143.7989926529, -3289396.260370887, 3854251.213474339, -1204350.7863173261, 827435.974478076 , 0 , 0.1864677328353816 , 0.9825013637459529, 10.424266879411505, 12.791127351240718, -3.256163611293784, 160.05156592217443, 532.7845130845899, 15.930539230456176, 62.2348457853685, -16035.267685450168, 16375.765827976264, -10553.648028867776, 7615.965376287072, 619.1370473148544, -259342.44893400747, -365616.21954569453, 1338351.6071672018, -1732867.981882307, 503087.95661363663, -437362.994824153";

if len(hdr.split(',')) == 38:       # generic quartic model (8 + 15 + 15)
    hdr_incr = 15
    print ("success ! quartic header")
elif len(hdr.split(',')) == 50:     # generic quintic model (8 + 21 + 21)
    hdr_incr = 21
    print ("success ! quintic header")
else:
    print ("Bad header length, neither quartic nor quintic !")
    sys.exit()

def prod_z_z (N1, p4, p3, p2, p1, p0, N2, q4, q3, q2, q1, q0):
    # multiply two (uniform) z polynomials to produce a new z polynomial (quintic)
    # first  : Sum_N1 (gik * z^i * zbar^k), degree N1
    # second : Sum_N2 (gik * z^i * zbar^k), degree N2

    ret = np.zeros((6,6), complex)    # final accumulated coeff of w^i*wbar^k
                                      # generated by z^i*zbar^k
    P = np.array([p4, p3, p2, p1, p0])
    Q = np.array([q4, q3, q2, q1, q0])
    for i in range (N1 + 1):
        for j in range (N2 + 1):
            #print (i, ',', j, ',', N1 + N2 - i - j, ',', i + j, ',', P[i], ',', Q[j])
            ret[N1 + N2 - i - j][i + j] += P[i]*Q[j]
    return ret

def prod_z_zbar (N1, p4, p3, p2, p1, p0, N2, q4, q3, q2, q1, q0):
    # multiply two (uniform) z polynomials to produce a new z polynomial (quintic)
    # first  : Sum_N1 (gik     * z^i * zbar^k) , degree N1
    # second : Sum_N2 (gik_bar * z^k * zbar^i) , degree N2

    ret = np.zeros((6,6), complex)    # final accumulated coeff of w^i*wbar^k
                                      # generated by z^i*zbar^k
    P = np.array([p4, p3, p2, p1, p0])
    Q = np.array([q4.conjugate(), q3.conjugate(), q2.conjugate(), q1.conjugate(), q0.conjugate()])
    for i in range (N1 + 1):
        for j in range (N2 + 1):
            #print (i, ',', j, ',', N1 + N2 - i - j, ',', i + j, ',', P[i], ',', Q[j])
            ret[N1 - i + j][N2 - j + i] += P[i]*Q[j]
    return ret

def prod_zbar_zbar (N1, p4, p3, p2, p1, p0, N2, q4, q3, q2, q1, q0):
    # multiply two (uniform) z polynomials to produce a new z polynomial (quintic)
    # first  : Sum_N1 (gik_bar * z^k * zbar^i) , degree N1
    # second : Sum_N2 (gik_bar * z^k * zbar^i) , degree N2

    ret = np.zeros((6,6), complex)    # final accumulated coeff of w^i*wbar^k
                                      # generated by z^i*zbar^k
    P = np.array([p4.conjugate(), p3.conjugate(), p2.conjugate(), p1.conjugate(), p0.conjugate()])
    Q = np.array([q4.conjugate(), q3.conjugate(), q2.conjugate(), q1.conjugate(), q0.conjugate()])
    for i in range (N1 + 1):
        for j in range (N2 + 1):
            #print (i, ',', j, ',', N1 + N2 - i - j, ',', i + j, ',', P[i], ',', Q[j])
            ret[i + j][N1 + N2 - i - j] += P[i]*Q[j]
    return ret

def concat (Ni, Nj):
    # expand function z^Ni*zbar^Nj in terms of w^l*wbar^m
    # calculate Xnew = Xold*Xin^Ni*Xin_bar*Nj
    Xold = np.zeros((6,6), complex)    # old accumulated coeff of w^l*wbar^m
    Xold[0][0] = 1
    Xnew = np.zeros((6,6), complex)    # new accumulated coeff of w^l*wbar^m
    #print ('\n', Ni, Nj, 'X0\n', X0.transpose().conjugate())
    for ni in range (Ni + Nj):
        if ni < Ni:
            Xin = X0                            # expand fxn z
        else:
            Xin = X0.transpose().conjugate()    # expand fxn zbar
        #print (ni, '\n', Xin)
        Xnew = np.zeros((6,6), complex)         # new accumulated coeff of w^l*wbar^m
        for oldi in range (6):
            for oldj in range (6 - oldi):
                for ini in range (6 - oldi - oldj):
                    for inj in range (6 - oldi - oldj - ini):
                        Xnew[oldi + ini][oldj + inj] += Xold[oldi][oldj]*Xin[ini][inj]
        #print (ni, '\n', Xnew)
        Xold = Xnew            
    return Xnew

def print_arr (arr):
    # print coeff of uniform polynomials
    print ('---------------')
    for i in range (len(arr)):
        for j in range (i + 1):
            print (arr[i - j][j], ',', end='')
        print ()
    print ('---------------')
    return

# parse a header given in the format of 'Chua_Simul_3.java'

g10 = complex(float(hdr.split(',')[9]), float(hdr.split(',')[9 + hdr_incr]))
print ("\ng10                , ", g10)

# for conversion matrices, use transpose of 'Minv' from 'expand_z.py'
# note gij coeff use factorial notation, Cxy do not

Cxy = np.array([complex(float(hdr.split(',')[11]), float(hdr.split(',')[11 + hdr_incr])), complex(float(hdr.split(',')[12]), float(hdr.split(',')[12 + hdr_incr])), \
                complex(float(hdr.split(',')[13]), float(hdr.split(',')[13 + hdr_incr]))])
print ("Cxy quadratic,", Cxy)
g20 = 2*np.matmul(Cxy, np.array([complex(1, 0), complex(0, -1), complex(-1, 0)]))/4 
g11 = 1*np.matmul(Cxy, np.array([complex(2, 0), complex(0,  0), complex( 2, 0)]))/4
g02 = 2*np.matmul(Cxy, np.array([complex(1, 0), complex(0,  1), complex(-1, 0)]))/4
print ("\ng20 g11 g02        , ", g20, ",", g11, ",", g02)

h20 = g20/(g10*g10 - g10)
h11 = g11/(g10*g10.conjugate() - g10)
h02 = g02/(g10.conjugate()*g10.conjugate() - g10)
print ("\nh20 h11 h02        , ", h20, ",", h11, ",", h02)

Cxy = np.array([complex(float(hdr.split(',')[14]), float(hdr.split(',')[14 + hdr_incr])), complex(float(hdr.split(',')[15]), float(hdr.split(',')[15 + hdr_incr])), \
                complex(float(hdr.split(',')[16]), float(hdr.split(',')[16 + hdr_incr])), complex(float(hdr.split(',')[17]), float(hdr.split(',')[17 + hdr_incr]))])
#print ("Cxy cubic,", Cxy) # see Chua_2D_cubic_variable_c1.py and IV p.9 and V p.2
g30 = 6*np.matmul(Cxy, np.array([complex(1, 0), complex(0, -1), complex(-1, 0), complex(0,  1)]))/8
g21 = 2*np.matmul(Cxy, np.array([complex(3, 0), complex(0, -1), complex( 1, 0), complex(0, -3)]))/8
g12 = 2*np.matmul(Cxy, np.array([complex(3, 0), complex(0,  1), complex( 1, 0), complex(0,  3)]))/8
g03 = 6*np.matmul(Cxy, np.array([complex(1, 0), complex(0,  1), complex(-1, 0), complex(0, -1)]))/8
print ("\ng30 g21 g12 g03    , ", g30, ",", g21, ",", g12, ",", g03)

Cxy = np.array([complex(float(hdr.split(',')[18]), float(hdr.split(',')[18 + hdr_incr])), complex(float(hdr.split(',')[19]), float(hdr.split(',')[19 + hdr_incr])), \
                complex(float(hdr.split(',')[20]), float(hdr.split(',')[20 + hdr_incr])), complex(float(hdr.split(',')[21]), float(hdr.split(',')[21 + hdr_incr])), \
                complex(float(hdr.split(',')[22]), float(hdr.split(',')[22 + hdr_incr]))])
#print ("Cxy quartic,", Cxy) # see Book 'Averaging' p. 57
g40 = 24*np.matmul(Cxy, np.array([complex(1, 0), complex(0,-1), complex(-1, 0), complex(0, 1), complex( 1, 0)]))/16
g31 =  6*np.matmul(Cxy, np.array([complex(4, 0), complex(0,-2), complex( 0, 0), complex(0,-2), complex(-4, 0)]))/16
g22 =  4*np.matmul(Cxy, np.array([complex(6, 0), complex(0, 0), complex( 2, 0), complex(0, 0), complex( 6, 0)]))/16
g13 =  6*np.matmul(Cxy, np.array([complex(4, 0), complex(0, 2), complex( 0, 0), complex(0, 2), complex(-4, 0)]))/16
g04 = 24*np.matmul(Cxy, np.array([complex(1, 0), complex(0, 1), complex(-1, 0), complex(0,-1), complex( 1, 0)]))/16
print ("\ng40 g31 g22 g13 g04, ", g40, ",", g31, ",", g22, ",", g13, ",", g04)

if hdr_incr == 21:                                      # quintic model
    Cxy = np.array([complex(float(hdr.split(',')[23]), float(hdr.split(',')[23 + hdr_incr])), complex(float(hdr.split(',')[24]), float(hdr.split(',')[24 + hdr_incr])), \
                    complex(float(hdr.split(',')[25]), float(hdr.split(',')[25 + hdr_incr])), complex(float(hdr.split(',')[26]), float(hdr.split(',')[26 + hdr_incr])), \
                    complex(float(hdr.split(',')[27]), float(hdr.split(',')[27 + hdr_incr])), complex(float(hdr.split(',')[28]), float(hdr.split(',')[28 + hdr_incr]))])
    print ("Cxy quintic,", Cxy)
    g50 = 120*np.matmul(Cxy, np.array([  1, -1j, -1,  1j,  1,  -1j]))/32
    g41 =  24*np.matmul(Cxy, np.array([  5, -3j, -1, -1j, -3,   5j]))/32
    g32 =  12*np.matmul(Cxy, np.array([ 10, -2j,  2, -2j,  2, -10j]))/32
    g23 =  12*np.matmul(Cxy, np.array([ 10,  2j,  2,  2j,  2,  10j]))/32
    g14 =  24*np.matmul(Cxy, np.array([  5,  3j, -1,  1j, -3,  -5j]))/32
    g05 = 120*np.matmul(Cxy, np.array([  1,  1j, -1, -1j,  1,   1j]))/32
    print ("\ng50 g41 g32 g23 g14 g05, ", g50, ",", g41, ",", g32, ",", g23, ",", g14, ",", g05)

print ("...........................................................")

# see 'Normal Form of Quartic Map', p. 3
# generate inverse of a quadratic z transform : cubic terms (N-S bifurc, Dec 21, 2021, p. 11)

A30 = 0.5*h20*h20 + 0.5*h11*h02.conjugate()             # cubic terms of w(z)
A21 = 1.5*h20*h11 + 0.5*h02*h02.conjugate() + h11*h11.conjugate()
A12 = 0.5*h20*h02 + h11*h11 + 0.5*h11*h20.conjugate() + h02*h11.conjugate()
A03 = 0.5*h11*h02 + 0.5*h02*h20.conjugate()
print ("\nA30 A21 A12 A03    , ", A30, ",", A21, ",", A12, ",", A03)

# generate inverse of a quadratic z transform : quartic terms

A40 = - h20*A30 - h11*A03.conjugate() - h20*h20*h20/8 - h20*h11*h02.conjugate()/4 - h02*h02.conjugate()*h02.conjugate()/8
A31 = - h11*A30 - h20*A21 - h11*A12.conjugate() - h02*A03.conjugate() - h20*h20*h11/2 - h20*h11*h11.conjugate()/2 - h11*h11*h02.conjugate()/2 - h02*h02.conjugate()*h11.conjugate()/2
A22 = - h20*A12 - h02*A12.conjugate() - h11*A21 - h11*A21.conjugate() - h20*h20*h02/4 - h20*h11*h11/2 - h20*h11*h20.conjugate()/4 - h11*h11*h11.conjugate() - h02*h11.conjugate()*h11.conjugate()/2  - h02*h02.conjugate()*h20.conjugate()/4 - h11*h02*h02.conjugate()/4
A13 = - h20*A03 - h11*A30.conjugate() - h11*A12 - h02*A21.conjugate() - h20*h11*h02/2 - h11*h11*h20.conjugate()/2 - h02*h11.conjugate()*h20.conjugate()/2 - h02*h11*h11.conjugate()/2
A04 = - h02*A30.conjugate() - h11*A03 - h20*h02*h02/8 - h02*h20.conjugate()*h20.conjugate()/8 - h11*h02*h20.conjugate()/4
print ("\nA40 A31 A22 A13 A04, ", A40, ",", A31, ",", A22, ",", A13, ",", A04)

# generate inverse of a quadratic z transform : quintic terms

A50 = - h20*A40                       - h11*A04.conjugate()
A41 = - h20*A31 - h02*A04.conjugate() - h11*A13.conjugate() - h11*A40
A32 = - h20*A22 - h02*A13.conjugate() - h11*A22.conjugate() - h11*A31
A23 = - h20*A13 - h02*A22.conjugate() - h11*A31.conjugate() - h11*A22
A14 = - h20*A04 - h02*A31.conjugate() - h11*A40.conjugate() - h11*A13
A05 =           - h02*A40.conjugate()                       - h11*A04

A50 += h20*h20*A30/2
A41 += h20*h20*A21/2 + h20*h11*A30
A32 += h20*h20*A12/2 + h20*h11*A21 + h20*h02*A30/2
A23 += h20*h20*A03/2 + h20*h11*A12 + h20*h02*A21/2
A14 +=                 h20*h11*A03 + h20*h02*A12/2
A05 +=                               h20*h02*A03/2

A50 += h02*h02.conjugate()*A03.conjugate()/2
A41 += h02*h02.conjugate()*A12.conjugate()/2 + h02*h11.conjugate()*A03.conjugate()
A32 += h02*h02.conjugate()*A21.conjugate()/2 + h02*h11.conjugate()*A12.conjugate() + h02*h20.conjugate()*A03.conjugate()/2
A23 += h02*h02.conjugate()*A30.conjugate()/2 + h02*h11.conjugate()*A21.conjugate() + h02*h20.conjugate()*A12.conjugate()/2
A14 +=                                         h02*h11.conjugate()*A30.conjugate() + h02*h20.conjugate()*A21.conjugate()/2
A05 +=                                                                               h02*h20.conjugate()*A30.conjugate()/2

A50 += h11*h02.conjugate()*A30/2
A41 += h11*h02.conjugate()*A21/2 + h11*h11.conjugate()*A30
A32 += h11*h02.conjugate()*A12/2 + h11*h11.conjugate()*A21 + h11*h20.conjugate()*A30/2
A23 += h11*h02.conjugate()*A03/2 + h11*h11.conjugate()*A12 + h11*h20.conjugate()*A21/2
A14 +=                             h11*h11.conjugate()*A03 + h11*h20.conjugate()*A12/2
A05 +=                                                       h11*h20.conjugate()*A03/2

A50 += h11*h20*A03.conjugate()/2
A41 += h11*h20*A12.conjugate()/2 + h11*h11*A03.conjugate()
A32 += h11*h20*A21.conjugate()/2 + h11*h11*A12.conjugate() + h11*h02*A03.conjugate()/2
A23 += h11*h20*A30.conjugate()/2 + h11*h11*A21.conjugate() + h11*h02*A12.conjugate()/2
A14 +=                             h11*h11*A30.conjugate() + h11*h02*A21.conjugate()/2
A05 +=                                                       h11*h02*A30.conjugate()/2
print ("\nA50_41_32_23_14_05 , ", A50, ",", A41, ",", A32, ",", A23, ",", A14, ",", A05)
print ('--------------------------------------------------------------------')

#   generate the product of a quintic z-map and a quadratic transform from z to w
#   to produce a quintic function in w
#   see attachment 'Normal Form of Quartic Map', p. 6-7
#   see attachment 'Final Quintic Map' (Oct 30, 2024), p.10

X0 = np.zeros((6,6), complex)    # final accumulated coeff of w^i*wbar^k
                                 # generated by z-map and z(w)
#print (g10, ',', g10*h20/2, ',', g10*h11, ',', g10*h02/2)
X0 += g10*prod_z_z (0, 1, 0, 0, 0, 0, 1,     1,     0,     0,     0,      0)
X0 += g10*prod_z_z (0, 1, 0, 0, 0, 0, 2, h20/2,   h11, h02/2,     0,      0)

X0 += 0.50*g20*prod_z_z (1,     1,   0,     0, 0, 0, 1,     1,   0,     0,      0,      0)
X0 += 0.50*g20*prod_z_z (1,     2,   0,     0, 0, 0, 2, h20/2, h11, h02/2,      0,      0)
X0 += 0.50*g20*prod_z_z (2, h20/2, h11, h02/2, 0, 0, 2, h20/2, h11, h02/2,      0,      0)

X0 +=      g11*prod_z_zbar (1,     1,   0,     0, 0, 0, 1,     1,   0,     0,   0,      0)
X0 +=      g11*prod_z_zbar (1,     1,   0,     0, 0, 0, 2, h20/2, h11, h02/2,   0,      0)
X0 +=      g11*prod_z_zbar (2, h20/2, h11, h02/2, 0, 0, 1,     1,   0,     0,   0,      0)
X0 +=      g11*prod_z_zbar (2, h20/2, h11, h02/2, 0, 0, 2, h20/2, h11, h02/2,   0,      0)

X0 += 0.50*g02*prod_zbar_zbar (1,     1,   0,     0, 0, 0, 1,     1,   0,     0,  0,    0)
X0 += 0.50*g02*prod_zbar_zbar (1,     2,   0,     0, 0, 0, 2, h20/2, h11, h02/2,  0,    0)
X0 += 0.50*g02*prod_zbar_zbar (2, h20/2, h11, h02/2, 0, 0, 2, h20/2, h11, h02/2,  0,    0)

X0 += 1.0/6.0*g30*prod_z_z (1,       1,     0,       0, 0, 0, 2,     1,   0,     0,  0,  0)
X0 += 1.0/6.0*g30*prod_z_z (2,       3,     0,       0, 0, 0, 2, h20/2, h11, h02/2,  0,  0)
X0 += 1.0/6.0*g30*prod_z_z (3, 3*h20/2, 3*h11, 3*h02/2, 0, 0, 2, h20/2, h11, h02/2,  0,  0)

X0 += 0.5*g21*prod_z_zbar (2,   1,     0,    0,     0, 0, 1,     1,   0,     0,  0,  0)
X0 += 0.5*g21*prod_z_zbar (3, h20, 2*h11,  h02,     0, 0, 1,     1,   0,     0,  0,  0)
X0 += 0.5*g21*prod_z_z    (3,   0, h20/2,  h11, h02/2, 0, 2, h20/2, h11, h02/2,  0,  0)
X0 += 0.5*g21*prod_z_zbar (2,   1,     0,    0,     0, 0, 2, h20/2, h11, h02/2,  0,  0)
X0 += 0.5*g21*prod_z_zbar (3, h20, 2*h11,  h02,     0, 0, 2, h20/2, h11, h02/2,  0,  0)

X0 += 0.5*g12*prod_z_zbar    (1,     1,     0,     0,     0, 0, 2,     1,     0,     0,  0,  0)
X0 += 0.5*g12*prod_z_zbar    (1,     1,     0,     0,     0, 0, 3,   h20, 2*h11,   h02,  0,  0)
X0 += 0.5*g12*prod_zbar_zbar (3,     0, h20/2,   h11, h02/2, 0, 2, h20/2,   h11, h02/2,  0,  0)
X0 += 0.5*g12*prod_z_zbar    (2, h20/2,   h11, h02/2,     0, 0, 2,     1,     0,     0,  0,  0)
X0 += 0.5*g12*prod_z_zbar    (2, h20/2,   h11, h02/2,     0, 0, 3,   h20, 2*h11,   h02,  0,  0)

X0 += 1.0/6.0*g03*prod_zbar_zbar (1,       1,     0,       0, 0, 0, 2,     1,   0,     0,  0,  0)
X0 += 1.0/6.0*g03*prod_zbar_zbar (2,       3,     0,       0, 0, 0, 2, h20/2, h11, h02/2,  0,  0)
X0 += 1.0/6.0*g03*prod_zbar_zbar (3, 3*h20/2, 3*h11, 3*h02/2, 0, 0, 2, h20/2, h11, h02/2,  0,  0)

X0 += 1.0/24.0*g40*prod_z_z   (2,   1,     0,   0, 0, 0, 2,     1,   0,     0,  0,  0)
X0 += 1.0/24.0*g40*prod_z_z   (3,   4,     0,   0, 0, 0, 2, h20/2, h11, h02/2,  0,  0)

X0 += 1.0/6.0*g31*prod_z_zbar (3,     1,     0,     0, 0, 0, 1,     1,   0,     0,  0,  0)
X0 += 1.0/6.0*g31*prod_z_zbar (3,     1,     0,     0, 0, 0, 2, h20/2, h11, h02/2,  0,  0)
X0 += 1.0/6.0*g31*prod_z_zbar (4, h20/2,   h11, h02/2, 0, 0, 1,     3,   0,     0,  0,  0)

X0 += 1.0/4.0*g22*prod_z_zbar (2,     1,     0,     0, 0, 0, 2,     1,   0,     0,  0,  0)
X0 += 1.0/4.0*g22*prod_z_zbar (2,     2,     0,     0, 0, 0, 3, h20/2, h11, h02/2,  0,  0)
X0 += 1.0/4.0*g22*prod_z_zbar (3, h20/2,   h11, h02/2, 0, 0, 2,     2,   0,     0,  0,  0)

X0 += 1.0/6.0*g13*prod_z_zbar (1,     1,     0,     0, 0, 0, 3,     1,   0,     0,  0,  0)
X0 += 1.0/6.0*g13*prod_z_zbar (2, h20/2,   h11, h02/2, 0, 0, 3,     1,   0,     0,  0,  0)
X0 += 1.0/6.0*g13*prod_z_zbar (1,     3,     0,     0, 0, 0, 4, h20/2, h11, h02/2,  0,  0)

X0 += 1.0/24.0*g04*prod_zbar_zbar (2, 1,     0,     0, 0, 0, 2,     1,   0,     0,  0,  0)
X0 += 1.0/24.0*g04*prod_zbar_zbar (3, 4,     0,     0, 0, 0, 2, h20/2, h11, h02/2,  0,  0)

if hdr_incr == 21:                                      # quintic model
    print ("add X0 quintic")
    X0 += 1.0/120.0*g50*prod_z_z       (3,   1,    0,   0, 0, 0, 2,     1,   0,     0,  0,  0)
    X0 +=  1.0/24.0*g41*prod_z_zbar    (4,   1,    0,   0, 0, 0, 1,     1,   0,     0,  0,  0)
    X0 +=  1.0/12.0*g32*prod_z_zbar    (3,   1,    0,   0, 0, 0, 2,     1,   0,     0,  0,  0)
    X0 +=  1.0/12.0*g23*prod_z_zbar    (2,   1,    0,   0, 0, 0, 3,     1,   0,     0,  0,  0)
    X0 +=  1.0/24.0*g14*prod_z_zbar    (1,   1,    0,   0, 0, 0, 4,     1,   0,     0,  0,  0)
    X0 += 1.0/120.0*g05*prod_zbar_zbar (3,   1,    0,   0, 0, 0, 2,     1,   0,     0,  0,  0)

print ('X0 =\n', X0)
print ('arr X0')
print_arr (X0)

# invert the z(w) transform

W = np.zeros((6,6), complex)    # final w-map, based on X0: coeff of w^l*wbar^m
W = W + concat(1, 0) - h20*concat(2, 0)/2 - h11*concat(1, 1) - h02*concat(0, 2)/2
W = W + A30*concat(3, 0) + A21*concat(2, 1) + A12*concat(1, 2) + A03*concat(0, 3)
W = W + A40*concat(4, 0) + A31*concat(3, 1) + A22*concat(2, 2) + A13*concat(1, 3) + A04*concat(0, 4)
W = W + A50*concat(5, 0) + A41*concat(4, 1) + A32*concat(3, 2) + A23*concat(2, 3) + A14*concat(1, 4) + A05*concat(0, 5)

print ('W =\n', W)
print ('arr W')
print_arr (W)

# extract complex gij (non-factorial) notation of 'uniform_quartic.py'

Q3 = (W[3][0], W[2][1], W[1][2], W[0][3])
Q4 = (W[4][0], W[3][1], W[2][2], W[1][3], W[0][4])
Q5 = (W[5][0], W[4][1], W[3][2], W[2][3], W[1][4], W[0][5])
print ("\nQ3, ", Q3)
print ("\nQ4, ", Q4)
print ("\nQ5, ", Q5)

# convert from gij (non-factorial) to complex (transformed) Dxy (non-factorial)
# for conversion matrices, use transpose of 'M' from 'expand_z.py'

D30 = np.matmul(Q3, np.array([  1,   1,   1,   1])) 
D21 = np.matmul(Q3, np.array([ 3j,  1j, -1j, -3j]))
D12 = np.matmul(Q3, np.array([ -3,   1,   1,  -3]))
D03 = np.matmul(Q3, np.array([-1j,  1j, -1j,  1j]))
print ("\nD30 D21 D12 D03    , ", D30, ",", D21, ",", D12, ",", D03)

D40 = np.matmul(Q4, np.array([  1,   1,   1,   1,   1])) 
D31 = np.matmul(Q4, np.array([ 4j,  2j,  0j, -2j, -4j]))
D22 = np.matmul(Q4, np.array([ -6,   0,   2,   0,  -6]))
D13 = np.matmul(Q4, np.array([-4j,  2j,  0j, -2j,  4j]))
D04 = np.matmul(Q4, np.array([  1,  -1,   1,  -1,   1]))
print ("\nD40 D31 D22 D13 D04, ", D40, ",", D31, ",", D22, ",", D13, ",", D04)

D50 = np.matmul(Q5, np.array([   1,   1,   1,   1,   1,   1])) 
D41 = np.matmul(Q5, np.array([  5j,  3j,  1j, -1j, -3j, -5j]))
D32 = np.matmul(Q5, np.array([ -10,  -2,   2,   2,  -2, -10]))
D23 = np.matmul(Q5, np.array([-10j,  2j,  2j, -2j, -2j, 10j]))
D14 = np.matmul(Q5, np.array([   5,  -3,   1,   1,  -3,   5]))
D05 = np.matmul(Q5, np.array([  1j, -1j,  1j, -1j,  1j, -1j]))
print ("\nD50 D41 D32 D23 D14 D05, ", D50, ",", D41, ",", D32, ",", D23, ",", D14, ",", D05)
print ()

#   print transformed Dxy hdr

print ('    private static String hdr = "', end='')
temp = hdr.split(',')
for i in range(8):
    print (temp[i], ",", end='')
print (" 0,", g10.real, ",",-g10.imag, ", 0, 0, 0,", D30.real, ",", D21.real, ",", D12.real, ",", D03.real, ",", D40.real, ",", D31.real, ",", D22.real, ",", D13.real, ",", D04.real, ",", end='')
print (D50.real, ",", D41.real, ",", D32.real, ",", D23.real, ",", D14.real, ",", D05.real, ",", end='')
print (" 0,", g10.imag, ",", g10.real, ", 0, 0, 0,", D30.imag, ",", D21.imag, ",", D12.imag, ",", D03.imag, ",", D40.imag, ",", D31.imag, ",", D22.imag, ",", D13.imag, ",", D04.imag, ",", end='')
print (D50.imag, ",", D41.imag, ",", D32.imag, ",", D23.imag, ",", D14.imag, ",", D05.imag, end='')
print ('";')
