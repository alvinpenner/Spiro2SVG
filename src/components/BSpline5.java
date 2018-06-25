
package components;

import java.io.*;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 5-point B-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// constrain only the slopes at the endpoints and keep P2 arbitrary
// linearize the equations wrt (d1, d2, x2, y2) and solve a 4x4 system of equations
// decompose the B-Spline into 2 Beziers, range (0,1) and (1,2).
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book2, Dec 2016, page 50

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BSpline5.java

public class BSpline5
{
    public static double t1_start = 0;           // Math.PI/3;
    public static final double t1_end = Math.PI/4; // Math.PI/4;
    public static final int N = 100;
    public static double[] Splinex, Spliney;    // 5 point spline
    public static double[][] Bezx;              // 2 Beziers, 4 points each, x component
    public static double[][] Bezy;              // 2 Beziers, 4 points each, y component
    //private static CircleFxn fitted;
    //private static CycloidFxn fitted;       // = new CycloidFxn(.5);           // set c value
    private static epiTrochoidFxn fitted; // = new epiTrochoidFxn(-2);       // set c value
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[4][N+1];        // partial wrt (d1, d2, x2, y2)
    private static double Jacdet = Double.NaN;
    private static double[] eig = new double[] {0,0,0,0};
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 80;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //read_data(80, 0);
        //read_data(1, 16);
        //fitted = new epiTrochoidFxn(10);            // keep this, 25 iterations to converge at c = 10
        //iterate_at_P2(19, 26, 175.2, 62);           // keep this, 25 iterations to converge at c = 10
        fitted = new epiTrochoidFxn(10);
        iterate_at_P2(19.983314292966476, 26.4276333695844, 175.47633731103548, 59.0566819528455);
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //System.out.println("BSpline phi c t1_start t1_end  = ," + phi + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end);
        //System.out.println("solve_at_P2 = " + solve_at_P2(19, 26, 175.2, 62, true));
        //iterate_at_P2(31.80597009761532, 11.910842854679, 162.18087919481738, 86.7856396874945);
        //solve_at_P2(23.84923550198231, 23.84923550197984, 170.525250238704, 70.63387137593989, true);
        //solve_at_P2(9.075207733717743, 48.893417162666445, 191.71824092253422, 34.261684630473695, true);
        //grid_search_at_P2(16.48986402964012, 29.256356995461363, 174.1358694806713, 59.53986312008041);

/*        // exercise matrix functions, to be deleted
        double[] v1 = new double[] {1.1, 3.2, 5.7, 6.8, 9.3};
        double[] v2 = new double[] {11.1, 3.7, 7.2, 9.1, 10.0};
        double[][] m = new double[][] {{1.7, 7.2, 9.7, 12.8, -9.3},
                                       {2.1, 4.2, 5.7,  5.1, 21.4},
                                       {4.1, 5.3, 7.0,  6.8, -3.2}};
        //System.out.println(multvv(v1, v2));
        double[] ret = multmv(m, v2);
        //System.out.println(ret[0] + ", " + ret[1] + ", " + ret[2]);
        m = new double[][] {{1.1, 2.3, 3.1, 5.7, -10.2},
                            {2.1, 4.2, 5.7, 5.1,  17.3},
                            {4.1, 5.3, 7.0, 6.8,   6.1},
                            {8.1, 7.3, 9.0, 2.8,   4.3},
                            {-2.1, 17, -9.2, 2.2, -1.7}};
        //System.out.println("cofactor = " + cofactor(m, 2, 1));
        //m = new double[][] {{1.1, 2.3, 3.1, 5.7},
        //                    {2.1, 4.2, 5.7, 5.1},
        //                    {4.1, 5.3, 7.0, 6.8},
        //                    {8.1, 7.3, 9.0, 2.8}};
        System.out.println("detm     = " + detm(m));
        double[][] minv = invertm2(m);
        System.out.println("invert = " + minv.toString());
        for (int i = 0; i < minv.length; i++)
        {
            for (int j = 0; j < minv.length; j++)
                System.out.print((float)minv[i][j] + ", ");
            System.out.println();
        }
*/
    }

    private static void iterate_at_P2(double d1, double d2, double x2, double y2)
    {
        // calculate a new estimate of (d1, d2, x2, y2) by setting dF = 0
        // include only first-order responses
        // see Spiro2SVG Book 3, page 54 (applied to 5-point cubic B-Spline)
        // setup 4-variable Newton-Raphson iteration

        final int MAXLOOP = 2000;
        double[] f_gx = new double[N+1];
        double[] f_gy = new double[N+1];
        double[] dfxdu = new double[N+1];
        double[] d2fxdudu = new double[N+1];
        double[] dfydu = new double[N+1];
        double[] d2fydudu = new double[N+1];
        double[][] dfxdd = new double[4][N+1];
        double[][] dfydd = new double[4][N+1];

        double[][] Jac = new double[4][4];
        double[] dFdd = new double[4];
        double[] trap_in = new double[N+1];
        double[] deld;                                              // (-Δd1, -Δd2, -Δx2, -Δy2)
        int i, j, k, loop = 0;
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(d1, d2, x2, y2, false)))   // initiallize at (x2, y2)
            {
                System.out.println("fail at " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
                return;
            }

            for (i = 0; i <= N; i++)
            {
                t1 = t1_start + i*(t1_end - t1_start)/N;
                f_gx[i] = multvv(Splinex, N43(t2[i])) - fitted.getx(t1);
                f_gy[i] = multvv(Spliney, N43(t2[i])) - fitted.gety(t1);
                dfxdu[i] = multvv(Splinex, dN43(t2[i]));
                d2fxdudu[i] = multvv(Splinex, d2N43(t2[i]));
                dfydu[i] = multvv(Spliney, dN43(t2[i]));
                d2fydudu[i] = multvv(Spliney, d2N43(t2[i]));
                dfxdd[0][i] = Math.cos(theta_start)*N43(t2[i])[1];
                dfydd[0][i] = Math.sin(theta_start)*N43(t2[i])[1];
                dfxdd[1][i] = -Math.cos(theta_end)*N43(t2[i])[3];
                dfydd[1][i] = -Math.sin(theta_end)*N43(t2[i])[3];
                dfxdd[2][i] = N43(t2[i])[2];
                dfydd[2][i] = 0;
                dfxdd[3][i] = 0;
                dfydd[3][i] = N43(t2[i])[2];
            }

            // calc dFdd[j] at current (d1, d2, x2, y2)

            for (i = 0; i < 4; i++)
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]); // original code
                    //trap_in[k] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];         // new code
                dFdd[i] = t2_vs_t1.integrate(trap_in);
            }

            // calc d2Fdd[i]dd[j] (Jacobean matrix)

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N; k++)
                    {
                        //    System.out.println(k + ", " + ", " + t2[k] + ", " + trap_in[k] + ", " + d2udddd[i][j][k] + ", " + f_gx[k] + ", " + dfxdu[k] + ", " + f_gy[k] + ", " + dfydu[k]);
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k]                // new code
                                   + dfydd[i][k]*dfydd[j][k]
                                   - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k];
                        //trap_in[k] = dfxdd[i][k]*dfxdd[j][k] + (dfxdd[j][k]*dfxdu[k] + f_gx[k]*d2fxdudd[j][k])*t2dd[i][k]   // old code
                        //           + dfydd[i][k]*dfydd[j][k] + (dfydd[j][k]*dfydu[k] + f_gy[k]*d2fydudd[j][k])*t2dd[i][k];
                        //System.out.println(k + ", " + trap_in[k]);
                        //System.out.println(k + ", " + (dfxdd[i][k]*dfxdd[j][k] + dfydd[i][k]*dfydd[j][k]));
                        //System.out.println(k + ", " + ((dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k]));
                    }
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }

            deld = multmv(invertm(Jac), dFdd);  // this is actually the negative of Δd
            d1 -= deld[0];
            d2 -= deld[1];
            x2 -= deld[2];
            y2 -= deld[3];
            Jacdet = detm(Jac);
            // calculate four eigenvalues
            double qua = -Jac[0][0] - Jac[1][1] - Jac[2][2] - Jac[3][3];
            double qub =  Jac[0][0]*Jac[1][1] + Jac[0][0]*Jac[2][2] + Jac[0][0]*Jac[3][3] + Jac[1][1]*Jac[2][2] + Jac[1][1]*Jac[3][3] + Jac[2][2]*Jac[3][3]
                       -  Jac[0][1]*Jac[0][1] - Jac[0][2]*Jac[0][2] - Jac[0][3]*Jac[0][3] - Jac[1][2]*Jac[1][2] - Jac[1][3]*Jac[1][3] - Jac[2][3]*Jac[2][3];
            double quc = -Jac[0][0]*Jac[1][1]*Jac[2][2] - Jac[0][0]*Jac[1][1]*Jac[3][3] - Jac[0][0]*Jac[2][2]*Jac[3][3] - Jac[1][1]*Jac[2][2]*Jac[3][3]
                       +  Jac[0][1]*Jac[0][1]*(Jac[2][2] + Jac[3][3]) + Jac[0][2]*Jac[0][2]*(Jac[1][1] + Jac[3][3]) + Jac[0][3]*Jac[0][3]*(Jac[1][1] + Jac[2][2])
                       +  Jac[1][2]*Jac[1][2]*(Jac[0][0] + Jac[3][3]) + Jac[1][3]*Jac[1][3]*(Jac[0][0] + Jac[2][2]) + Jac[2][3]*Jac[2][3]*(Jac[0][0] + Jac[1][1])
                       - 2*Jac[0][1]*Jac[1][2]*Jac[0][2] - 2*Jac[0][2]*Jac[2][3]*Jac[0][3] - 2*Jac[0][1]*Jac[1][3]*Jac[0][3] - 2*Jac[1][2]*Jac[2][3]*Jac[1][3];
            double qud = Jac[0][0]*Jac[1][1]*Jac[2][2]*Jac[3][3]
                       + Jac[0][1]*Jac[0][1]*Jac[2][3]*Jac[2][3] + Jac[0][2]*Jac[0][2]*Jac[1][3]*Jac[1][3] + Jac[0][3]*Jac[0][3]*Jac[1][2]*Jac[1][2]
                       - Jac[0][1]*Jac[0][1]*Jac[2][2]*Jac[3][3] - Jac[0][2]*Jac[0][2]*Jac[1][1]*Jac[3][3] - Jac[0][3]*Jac[0][3]*Jac[1][1]*Jac[2][2]
                       - Jac[1][2]*Jac[1][2]*Jac[0][0]*Jac[3][3] - Jac[1][3]*Jac[1][3]*Jac[0][0]*Jac[2][2] - Jac[2][3]*Jac[2][3]*Jac[0][0]*Jac[1][1]
                       + 2*Jac[0][0]*Jac[1][2]*Jac[2][3]*Jac[1][3] + 2*Jac[1][1]*Jac[0][2]*Jac[2][3]*Jac[0][3] + 2*Jac[2][2]*Jac[0][1]*Jac[1][3]*Jac[0][3] + 2*Jac[3][3]*Jac[0][1]*Jac[1][2]*Jac[0][2]
                       - 2*Jac[0][1]*Jac[1][2]*Jac[2][3]*Jac[0][3] - 2*Jac[0][1]*Jac[0][2]*Jac[2][3]*Jac[1][3] - 2*Jac[0][2]*Jac[0][3]*Jac[1][3]*Jac[1][2];
            eig = fitymoment.solve_quartic_all(1, qua, qub, quc, qud);
            //System.out.println("eigenvalue = " + eig[0] + ", " + eig[1] + ", " + eig[2] + ", " + eig[3]);
            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + dFdd[2] + ", " + dFdd[3] + ", " + Jacdet);
            System.out.println("deld = " + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3]);
            dump_Jac(Jac);

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                t2[i] -= t2dd[0][i]*deld[0] + t2dd[1][i]*deld[1] + t2dd[2][i]*deld[2] + t2dd[3][i]*deld[3];   // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL) && (Math.abs(deld[2]) < TOL) && (Math.abs(deld[3]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 x2 y2 = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
            solve_at_P2(d1, d2, x2, y2, true);                          // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ")");
    }

    public static void dump_Jac(double[][] J)
    {
        // see 'python_eigenvalues.txt' to calculate eigenvalues using numpy
        System.out.print("a = np.array([");
        for (int i = 0; i < J.length; i++)
        {
            if (i > 0) System.out.print(", ");
            System.out.print("[");
            for (int j = 0; j < J.length; j++)
            {
                if (j > 0) System.out.print(", ");
                System.out.print(J[i][j]);
            }
            System.out.print("]");
        }
        System.out.println("])");
    }
/*
    private static void iterate_at_P2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        // calculate a new estimate of (x2, y2) by setting dF/dx2 = dF/dy2 = 0
        // fit curvature at endpoints to determine d1, d2 as fxns of P2
        // include only first-order responses
        // see Spiro2SVG Book 3, page 52
        // setup two elliptical equations for x2, y2

        final int MAXLOOP = 200;
        double A0, B0, C0, D0, E0, F0;
        double A1, B1, C1, D1, E1, F1;
        double A2, C2, D2, E2, F2;
        double t1, old_x2, old_y2;
        double[] h_gx = new double[N+1];
        double[] h_gy = new double[N+1];
        double[] dhx = new double[N+1];
        double[] dhy = new double[N+1];
        double[] ux = new double[N+1];
        double[] uy = new double[N+1];
        double[] vx = new double[N+1];
        double[] vy = new double[N+1];
        double[] dux = new double[N+1];
        double[] duy = new double[N+1];
        double[] dvx = new double[N+1];
        double[] dvy = new double[N+1];
        double[] trap_in = new double[N+1];
        int i, loop = 0;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(x2, y2, false)))               // initiallize at (x2, y2)
            {
                System.out.println("fail at " + x2 + ", " + y2 + ", " + calc_d1(x2, y2) + ", " + calc_d2(x2, y2));
                return;
            }

            for (i = 0; i <= N; i++)
            {
                t1 = t1_start + i*(t1_end - t1_start)/N;
                h_gx[i] = multvv(Splinex, N43(t2[i])) - fitted.getx(t1);
                h_gy[i] = multvv(Spliney, N43(t2[i])) - fitted.gety(t1);
                dhx[i] = multvv(Splinex, dN43(t2[i]));
                dhy[i] = multvv(Spliney, dN43(t2[i]));
                ux[i] = N43(t2[i])[2] + Math.cos(theta_start)*N43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.cos(theta_end)*N43(t2[i])[3]*calc_d2dx2(x2, y2);
                vx[i] =                 Math.cos(theta_start)*N43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.cos(theta_end)*N43(t2[i])[3]*calc_d2dy2(x2, y2);
                uy[i] =                 Math.sin(theta_start)*N43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.sin(theta_end)*N43(t2[i])[3]*calc_d2dx2(x2, y2);
                vy[i] = N43(t2[i])[2] + Math.sin(theta_start)*N43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.sin(theta_end)*N43(t2[i])[3]*calc_d2dy2(x2, y2);
                dux[i] = dN43(t2[i])[2] + Math.cos(theta_start)*dN43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.cos(theta_end)*dN43(t2[i])[3]*calc_d2dx2(x2, y2);
                dvx[i] =                  Math.cos(theta_start)*dN43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.cos(theta_end)*dN43(t2[i])[3]*calc_d2dy2(x2, y2);
                duy[i] =                  Math.sin(theta_start)*dN43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.sin(theta_end)*dN43(t2[i])[3]*calc_d2dx2(x2, y2);
                dvy[i] = dN43(t2[i])[2] + Math.sin(theta_start)*dN43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.sin(theta_end)*dN43(t2[i])[3]*calc_d2dy2(x2, y2);
                h_gx[i] -= ux[i]*x2 + vx[i]*y2;
                h_gy[i] -= uy[i]*x2 + vy[i]*y2;
                dhx[i] -= dux[i]*x2 + dvx[i]*y2;
                dhy[i] -= duy[i]*x2 + dvy[i]*y2;
            }

            // set dFdx2 = 0

            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dux[i] + uy[i]*duy[i])*t2dx2[i];
            A0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (vx[i]*dvx[i] + vy[i]*dvy[i])*t2dx2[i];
            B0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dvx[i] + dux[i]*vx[i] + uy[i]*dvy[i] + duy[i]*vy[i])*t2dx2[i];
            C0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = ux[i]*ux[i] + (ux[i]*dhx[i] + dux[i]*h_gx[i])*t2dx2[i]
                           + uy[i]*uy[i] + (uy[i]*dhy[i] + duy[i]*h_gy[i])*t2dx2[i];
            D0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = vx[i]*ux[i] + (vx[i]*dhx[i] + dvx[i]*h_gx[i])*t2dx2[i]
                           + vy[i]*uy[i] + (vy[i]*dhy[i] + dvy[i]*h_gy[i])*t2dx2[i];
            E0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = h_gx[i]*(ux[i] + dhx[i]*t2dx2[i])
                           + h_gy[i]*(uy[i] + dhy[i]*t2dx2[i]);
            F0 = integrate(trap_in);

            // set dFdy2 = 0

            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dux[i] + uy[i]*duy[i])*t2dy2[i];
            A1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (vx[i]*dvx[i] + vy[i]*dvy[i])*t2dy2[i];
            B1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dvx[i] + dux[i]*vx[i] + uy[i]*dvy[i] + duy[i]*vy[i])*t2dy2[i];
            C1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = ux[i]*vx[i] + (ux[i]*dhx[i] + dux[i]*h_gx[i])*t2dy2[i]
                           + uy[i]*vy[i] + (uy[i]*dhy[i] + duy[i]*h_gy[i])*t2dy2[i];
            D1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = vx[i]*vx[i] + (vx[i]*dhx[i] + dvx[i]*h_gx[i])*t2dy2[i]
                           + vy[i]*vy[i] + (vy[i]*dhy[i] + dvy[i]*h_gy[i])*t2dy2[i];
            E1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = h_gx[i]*(vx[i] + dhx[i]*t2dy2[i])
                           + h_gy[i]*(vy[i] + dhy[i]*t2dy2[i]);
            F1 = integrate(trap_in);

            // see: https://math.stackexchange.com/questions/1767225/algorithm-intersection-of-two-conics

            A2 = A0*B1 - A1*B0;
            C2 = C0*B1 - C1*B0;
            D2 = D0*B1 - D1*B0;
            E2 = E0*B1 - E1*B0;
            F2 = F0*B1 - F1*B0;
            old_x2 = x2;
            old_y2 = y2;
            x2 = fitymoment.solve_quartic(A0*C2*C2 + B0*A2*A2 - C0*C2*A2,
                                          2*A0*C2*E2 + D0*C2*C2 + 2*B0*A2*D2 - C0*C2*D2 - (C0*E2 + C2*E0)*A2,
                                          A0*E2*E2 + 2*D0*C2*E2 + F0*C2*C2 + 2*B0*A2*F2 + B0*D2*D2
                                          - C0*C2*F2 - (C0*E2 + C2*E0)*D2 - E0*E2*A2,
                                          D0*E2*E2 + 2*F0*C2*E2 + 2*B0*D2*F2 - (C0*E2 + C2*E0)*F2 - E0*E2*D2,
                                          F0*E2*E2 + B0*F2*F2 - E0*E2*F2,
                                          true);
            y2 = -(A2*x2*x2 + D2*x2 + F2)/(C2*x2 + E2);
            //System.out.println("new           d1 d2    = ," + d1 + ", " + d2
            //                 + ", " + (A0*d1*d1 + B0*d2*d2 + C0*d1*d2 + D0*d1 + E0*d2 + F0)
            //                 + ", " + (A1*d1*d1 + B1*d2*d2 + C1*d1*d2 + D1*d1 + E1*d2 + F1));

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                t2[i] += t2dx2[i]*(x2 - old_x2) + t2dy2[i]*(y2 - old_y2);   // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(old_x2 - x2) < TOL) && (Math.abs(old_y2 - y2) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new x2 y2 = , , , , , , " + x2 + ", " + y2 + ", " + Math.abs(old_x2 - x2) + ", " + Math.abs(old_y2 - y2));
            solve_at_P2(x2, y2, true);                  // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops!");
    }
*/
    private static void read_data(int theta, double tempc)
    {
        // read cubic Bezier (d1, d2) data from a file
        // match field 1 = theta for a Cycloid / root number for an epiTrochoid
        // match field 2 = c for an epiTrochoid
        // convert to P2 for a 5 point B-Spline, and initiallize 'solve_at_P2()'

        String str, firstline = "";
        double d1 = 0;                  // Bezier arm length
        double d2 = 0;

        try
        {
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig5_hypocofmd1d2.csv"));
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig7_hypoODFd1d2.csv"));
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig2_Cycloidcofmd1d2.csv"));
            BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig2_CycloidODFd1d2rms.csv"));
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig2_CircleODFd1d2rms.csv"));
            try
            {
                if (istr.ready())
                    firstline = istr.readLine();
                while (istr.ready())
                {
                    str = istr.readLine();
                    //System.out.println(str);
                    if (!str.isEmpty())
                        if ((firstline.startsWith("circle") && (Integer.parseInt(str.split(",")[0].trim()) == theta))
                        ||  (firstline.startsWith("theta") && (Integer.parseInt(str.split(",")[0].trim()) == theta))
                        ||  (firstline.startsWith("root") && (Integer.parseInt(str.split(",")[0].trim()) == theta) && (Double.parseDouble(str.split(",")[1]) == tempc)))
                        {
                            d1 = Double.parseDouble(str.split(",")[3]);
                            d2 = Double.parseDouble(str.split(",")[4]);
                            break;
                        }
                }
            }
            catch (IOException e)
            {
                System.out.println("read error : " + e.getMessage());
                return;
            }
        }
        catch (FileNotFoundException e)
        {
            System.out.println("file not found : " + e.getMessage());
            return;
        }
        if (d1 == 0 && d2 == 0)
        {
            System.out.println("file data match not found, abort");
            return;
        }
        //if (firstline.startsWith("circle"))                     // circle, radius c
        //{
        //    tempc = 4.5;                                        // override for circle
        //    fitted = new CircleFxn(tempc);                      // set c value (radius)
        //    t1_start = Math.PI - theta*Math.PI/180;             // set arc angle
        //}
        if (firstline.startsWith("theta"))                    // cycloid, calculate c
        {
            // calculate the point of maximum curvature of a cycloid from a tangent angle
            tempc = Math.sqrt(1 - .75*Math.cos(theta*Math.PI/180)*Math.cos(theta*Math.PI/180)); // override for cycloid
            //fitted = new CycloidFxn(tempc);                             // set c value
            t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        }
        //if (firstline.startsWith("root"))
        //    fitted = new epiTrochoidFxn(tempc);                           // set c value
        //d1 = 2*9.5;        // fix fix test code
        //d2 = 2*52.9;
        System.out.println("file data at theta c t d1 d2   = ," + theta + ", , " + tempc + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + fitted.getkappa(t1_start) + ", " + fitted.getkappa(t1_end));

        // calculate P2 for a B-Spline

        theta_start = fitted.gettheta(t1_start);                    // possibly redundant
        theta_end = fitted.gettheta(t1_end);
        double[] Bez = new double[] {fitted.getx(t1_start),
                                     fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                     fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                     fitted.getx(t1_end)};
        double P2x = (Bez[1] + Bez[2])/2;       // convert from Bezier to B-Spline
        Bez = new double[] {fitted.gety(t1_start),
                            fitted.gety(t1_start) + d1*Math.sin(theta_start),
                            fitted.gety(t1_end) - d2*Math.sin(theta_end),
                            fitted.gety(t1_end)};
        double P2y = (Bez[1] + Bez[2])/2;
        //P2x = 0.67; // 1.2545; // 2.408125;
        //P2y = 1.39; // 1.5819; // 1.540;
        // convert from cubic Bezier (d1, d2) to B-Spline (d1, d2)
        System.out.println("return code solve_at_P2 = " + solve_at_P2(d1/2, d2/2, P2x, P2y, true));
        //grid_search_at_P2(P2x, P2y);
        //iterate_at_P2(d1/2, d2/2, P2x, P2y);     // default
        //iterate_at_P2(9.5, 52.9, 192, 33.5);      // over-ride
        //iterate_at_P2(16.6, 26, 180, 50);      // over-ride
        //iterate_at_P2(10, 44, 190, 36);      // over-ride
    }

    private static void grid_search_at_P2(double d1, double d2, double x2, double y2)
    {
        double del = 0.0001;

        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
                for (int k = -1; k < 2; k++)
                    for (int l = -1; l < 2; l++)
                    {
                        System.out.println("grid =, " + i + ", " + j + ", " + k + ", " + l);
                        solve_at_P2(d1 + i*del, d2 + j*del, x2 + k*del, y2 + l*del, false);
                    }
    }

/*
    private static void grid_search_at_P2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        final int MAXCOUNT = 1000;
        double del = 0.000001;
        double err, errmin = 999999;
        int iold = 0, jold = 0;
        int imin = 0, jmin = 0;
        int count = 0;

        do
        {
            count++;
            errmin = 999999;
            for (int i = -1; i < 2; i++)
                for (int j = -1; j < 2; j++)
                {
                    err = solve_at_P2(x2 + (i + iold)*del, y2 + (j + jold)*del, false);
                    if (err < errmin)
                    {
                        imin = i;
                        jmin = j;
                        errmin = err;
                    }
                }
            iold += imin;
            jold += jmin;
            System.out.println(",,,,,,,,,,, " + (float)(x2 + iold*del) + ", " + (float)(y2 + jold*del) + ", " + errmin);
        } while ((count < MAXCOUNT) && !(imin == 0 && jmin == 0));
        if (count == MAXCOUNT)
            System.out.println("Not converged!");
        else
            System.out.println("Converged in " + count);
    }
*/
    private static double solve_at_P2(double d1, double d2, double x2, double y2, boolean print)
    {
        // 5-point cubic B-Spline curve
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        Splinex = new double[] {fitted.getx(t1_start),
                                fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                x2,
                                fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                fitted.getx(t1_end)};
        Spliney = new double[] {fitted.gety(t1_start),
                                fitted.gety(t1_start) + d1*Math.sin(theta_start),
                                y2,
                                fitted.gety(t1_end) - d2*Math.sin(theta_end),
                                fitted.gety(t1_end)};
        Bezx = new double[][] {{Splinex[0],
                                Splinex[1],
                                (Splinex[1] + Splinex[2])/2,
                                (Splinex[1] + 2*Splinex[2] + Splinex[3])/4},
                               {(Splinex[1] + 2*Splinex[2] + Splinex[3])/4,
                                (Splinex[2] + Splinex[3])/2,
                                Splinex[3],
                                Splinex[4]}};
        Bezy = new double[][] {{Spliney[0],
                                Spliney[1],
                                (Spliney[1] + Spliney[2])/2,
                                (Spliney[1] + 2*Spliney[2] + Spliney[3])/4},
                               {(Spliney[1] + 2*Spliney[2] + Spliney[3])/4,
                                (Spliney[2] + Spliney[3])/2,
                                Spliney[3],
                                Spliney[4]}};
        if (t2[N] == 0)
            System.out.println("__start B-Spline5 at theta c t d1 d2 = , " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + calc_error());
        //fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2, t2dd1, t2dd2, t2dx2, t2dy2");
        int seg = 0;                // Bezier segment, before or after the splice
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i, seg);
            if (seg == 0 && t2[i] > 1)
            {
                seg++;
                solve_quintic_for_t2(i, seg);   // re-calculate
            }
            if (seg == 1 && t2[i] < 1)
                return Double.NaN;
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                scan_quintic_near_t2(i, seg, t2[i]);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dxy(i, t2[i], "d1");
            t2dd[1][i] = calc_t2dxy(i, t2[i], "d2");
            t2dd[2][i] = calc_t2dxy(i, t2[i], "x2");
            t2dd[3][i] = calc_t2dxy(i, t2[i], "y2");
            if (print)
                System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i]);
        }
        double retVal = calc_error();
        System.out.println("__new t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + (float) retVal + ", " + (float) Jacdet + ", " + (float) eig[1] + ", " + (float) eig[3] + ", " + (float) eig[2] + ", " + (float) eig[0]);
        return retVal;
    }

    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

        double a_b = 180;         // scale factor to make rms error dimensionless
        //double a_b = 1;             // Cycloid only
        double t1 = t1_start;
        double[] trap_in = new double[N+1];
        int seg = 0;                // Bezier segment, before or after the splice

        if ((Math.abs(2 - t2[N]) > TOL) || (Math.abs(t2[0]) > TOL))
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i <= N; i++)
        {
            if (seg == 0 && t2[i] > 1)
                seg++;
            trap_in[i] = (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))*(t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))
                       + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1))*(t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1));
            //System.out.println(i + ", " + seg + ", " + (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1)) + ", " + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1)) + ", " + Math.sqrt(trap_in[i]));
            t1 += (t1_end - t1_start)/N;
        }
        return Math.sqrt(t2_vs_t1.integrate(trap_in))/a_b;
    }

    private static void solve_quintic_for_t2(int i, int seg)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the cubic Bezier t-value

        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double f, fprime, f2prime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double t;
        int loop = 0;

        // initial estimate using quadratic approximation

        if (i == 0) t = 0;
        else t = t2[i-1];
        t -= seg;               // compensate for Bezier segment offset
        f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
        fprime = t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d2fn(Bezx[seg], t) + t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.dfn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d2fn(Bezy[seg], t);
        f2prime = 3*t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.d2fn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d3fn(Bezx[seg], t) + 3*t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.d2fn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d3fn(Bezy[seg], t);
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
            del_t = -fprime/f2prime;
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            if (del_t < 0 || del_t > 2)
            {
                del_t = (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
                if (del_t < 0 || del_t > 2)
                {
                    System.out.println("\nBad init  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
                    t2[i] = Double.NaN;
                    return;
                }
            }
        }
        t += del_t;

        //System.out.println("roots = " + (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime + ", " + (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime);
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
        do
        {
            f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
            fprime = t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d2fn(Bezx[seg], t) + t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.dfn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d2fn(Bezy[seg], t);
            if (loop > 100)
            {
                t2[i] = Double.NaN;
                return;
            }
            if (f == 0 && fprime == 0)
                del_t = 0;
            else
                del_t = -f/fprime;
            t += del_t;
            loop++;
            //System.out.println("         t2 =, " + t + ", " + f + ", " + fprime);
        } while (Math.abs(del_t) > TOL);
        t2[i] = t + seg;               // compensate for Bezier segment offset
    }

    private static void scan_quintic_near_t2(int index, int seg, double t2_bad)
    {
        // if solve_quintic_for_t2 fails, scan the area for other roots
        double t1 = t1_start + index*(t1_end - t1_start)/N;
        double f, t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        System.out.println("\nscanning at " + t1 + ", " + t2_bad);
        for (int i = -20; i <= 200; i++)
        {
            t = -i*t2_bad/10;
            t -= seg;               // compensate for Bezier segment offset
            f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
            System.out.println(t + ", " + f);
        }
        System.out.println();
    }

    private static double calc_t2dxy(int i, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double fnx = multvv(Splinex, N43(t2));
        double fny = multvv(Spliney, N43(t2));
        double dfnx = multvv(Splinex, dN43(t2));
        double dfny = multvv(Spliney, dN43(t2));
        double d2fnx = multvv(Splinex, d2N43(t2));
        double d2fny = multvv(Spliney, d2N43(t2));
        double denom = dfnx*dfnx + (fnx - X)*d2fnx + dfny*dfny + (fny - Y)*d2fny;
        double numer = Double.NaN;
        if (type.equals("x2"))
            numer = dfnx*N43(t2)[2] + (fnx - X)*dN43(t2)[2];
        else if(type.equals("y2"))
            numer = dfny*N43(t2)[2] + (fny - Y)*dN43(t2)[2];
        else if(type.equals("d1"))
            numer = (Math.cos(theta_start)*dfnx + Math.sin(theta_start)*dfny)*N43(t2)[1]
                  + (Math.cos(theta_start)*(fnx - X) + Math.sin(theta_start)*(fny - Y))*dN43(t2)[1];
        else if(type.equals("d2"))
            numer = -(Math.cos(theta_end)*dfnx + Math.sin(theta_end)*dfny)*N43(t2)[3]
                  -  (Math.cos(theta_end)*(fnx - X) + Math.sin(theta_end)*(fny - Y))*dN43(t2)[3];
        //System.out.println("calc_t2dx2 = " + t1 + ", " + t2 + ", " + rhs1 + ", " + rhs2 + ", " + fprime + ", " + ((-3*t2*(1 - t2)*(1 - t2)*rhs1 - 3*(1 - t2)*(1 - 3*t2)*rhs2)/fprime));
        return -numer/denom;
    }

    private static double[] N43(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {(1 - u)*(1 - u)*(1 - u),
                                 u*(12 - 18*u + 7*u*u)/4,
                                 u*u*(3 - 2*u)/2,
                                 u*u*u/4,
                                 0};
        else if (u < 2 + TOL)
            return new double[] {0,
                                 (2 - u)*(2 - u)*(2 - u)/4,
                                 (2 - u)*(2 - u)*(2*u - 1)/2,
                                 (2 - u)*(7*u*u - 10*u + 4)/4,
                                 (u - 1)*(u - 1)*(u - 1)};
        else
            return null;
    }

    private static double[] dN43(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {-3*(1 - u)*(1 - u),
                                 (12 - 36*u + 21*u*u)/4,
                                 3*u*(1 - u),
                                 3*u*u/4,
                                 0};
        else if (u < 2 + TOL)
            return new double[] {0,
                                 -3*(2 - u)*(2 - u)/4,
                                 -3*(2 - u)*(u - 1),
                                 -(21*u*u - 48*u + 24)/4,
                                 3*(u - 1)*(u - 1)};
        else
            return null;
    }

    private static double[] d2N43(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {6*(1 - u),
                                 (-36 + 42*u)/4,
                                 3*(1 - 2*u),
                                 3*u/2,
                                 0};
        else if (u < 2 + TOL)
            return new double[] {0,
                                 3*(2 - u)/2,
                                 3*(2*u - 3),
                                 (-42*u + 48)/4,
                                 6*(u - 1)};
        else
            return null;
    }

    public static double multvv(double[] v1, double[] v2)
    {
        if (v1.length != v2.length) return Double.NaN;
        double retVal = 0;

        for (int i = 0; i < v1.length; i++)
            retVal += v1[i]*v2[i];
        return retVal;
    }

    public static double[] multmv(double[][] m, double[] v)
    {
        if (m[0].length != v.length) return null;
        double[] retVal = new double[m.length];

        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < v.length; j++)
                retVal[i] += m[i][j]*v[j];
        return retVal;
    }

    public static double[][] invertm(double[][] m)
    {
        if (m.length != m[0].length) return null;
        double det = detm(m);
        double[][] retVal = new double[m.length][m.length];
        double[][] sub;
        int sgni = -1, sgnj;

        for (int i = 0; i < m.length; i++)
        {
            sgni *= -1;
            sgnj = -1;
            for (int j = 0; j < m.length; j++)
            {
                sgnj *= -1;
                sub = new double[m.length - 1][m.length - 1];
                for (int ii = 0; ii < m.length; ii++)
                    for (int jj = 0; jj < m.length; jj++)
                        if (ii < i && jj < j)
                            sub[ii][jj] = m[ii][jj];
                        else if (ii < i && jj > j)
                            sub[ii][jj - 1] = m[ii][jj];
                        else if (ii > i && jj < j)
                            sub[ii - 1][jj] = m[ii][jj];
                        else if (ii > i && jj > j)
                            sub[ii - 1][jj - 1] = m[ii][jj];
                retVal[j][i] = sgni*sgnj*detm(sub)/det;
            }
        }
        //System.out.println("invertm");
        //for (int i = 0; i < retVal.length; i++)
        //{
        //    for (int j = 0; j < retVal.length; j++)
        //        System.out.print(retVal[i][j] + ", ");
        //    System.out.println();
        //}
        return retVal;
    }

    public static double detm(double[][] m)
    {
        if (m.length != m[0].length) return Double.NaN;
        if (m.length == 1) return m[0][0];

        int sgn = -1;
        double det = 0;
        for (int k = 0; k < m.length; k++)
        {
            double[][] sub = new double[m.length - 1][m.length - 1];
            sgn *= -1;
            for (int i = 1; i < m.length; i++)
                for (int j = 0; j < m.length; j++)
                    if (j < k)
                        sub[i-1][j] = m[i][j];
                    else if (j > k)
                        sub[i-1][j-1] = m[i][j];
            det += sgn*m[0][k]*detm(sub);
        }
        return det;
    }

    public static double[] gaussj(double[][] m, double[] v)
    {
        // solve m*x = v
        // based on the routine gaussj() in "Numerical Recipes in C", p.39, W.H.Press
        // Gauss-Jordan elimination with full pivoting
        if (m.length != m[0].length) return null;
        if (m.length != v.length) return null;
        int i, j, k, l, ll, icol = 0, irow = 0;
        int n = m.length;
        double big, dum, pivinv;
        double[][] a = new double[n][n];
        double[] b = new double[n];
        int[] indxc = new int[n];
        int[] indxr = new int[n];
        int[] ipiv = new int[n];
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                a[i][j] = m[i][j];
        for (i = 0; i < n; i++)
            b[i] = v[i];
        for (i = 0; i < n; i++)
            ipiv[i] = 0;
        for (i = 0; i < n; i++)
        {
            big = 0;
            for (j = 0; j < n; j++)
                if (ipiv[j] != 1)
                    for (k = 0; k < n; k++)
                        if (ipiv[k] == 0)
                        {
                            if (Math.abs(a[j][k]) >= big)
                            {
                                big = Math.abs(a[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                        else if (ipiv[k] > 1)
                        {
                            System.out.println("gaussj error : Singular Matrix 1");
                            return null;
                        }
            ++(ipiv[icol]);
            if (irow != icol)
            {
                for (l = 0; l < n; l++) swapa (a, irow, l, icol, l);
                swapv(b, irow, icol);
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if (a[icol][icol] == 0)
            {
                System.out.println("gaussj error : Singular Matrix 2");
                return null;
            }
            pivinv = 1/a[icol][icol];
            a[icol][icol] = 1;
            for (l = 0; l < n; l++) a[icol][l] *= pivinv;
            b[icol] *= pivinv;
            for (ll = 0; ll < n; ll++)
                if (ll != icol)
                {
                    dum = a[ll][icol];
                    a[ll][icol] = 0;
                    for (l = 0; l < n; l++) a[ll][l] -= a[icol][l]*dum;
                    b[ll] -= b[icol]*dum;
                }
        }
        for (l = n - 1; l >= 0; l--)
            if (indxr[l] != indxc[l])
                for (k = 0; k < n; k++)
                    swapa(a, k, indxr[l], k, indxc[l]);
        //System.out.println("a inverse");
        //for (i = 0; i < n; i++)
        //{
        //    for (j = 0; j < n; j++)
        //        System.out.print(a[i][j] + ", ");
        //    System.out.println();
        //}
        //return a;
        //System.out.println("b");
        //for (i = 0; i < n; i++)
        //    System.out.println(b[i]);
        return b;
    }

    private static void swapa(double[][] a, int i1, int j1, int i2, int j2)
    {
        double temp = a[i1][j1];
        a[i1][j1] = a[i2][j2];
        a[i2][j2] = temp;
    }

    private static void swapv(double[] v, int i1, int i2)
    {
        double temp = v[i1];
        v[i1] = v[i2];
        v[i2] = temp;
    }
}
