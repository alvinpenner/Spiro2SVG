
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 4-point cubic Bezier (P0 - P3) to it, using parameter 0 < t2 < 1.
// constrain only the slopes at the endpoints
// linearize the equations wrt (d1, d2) and solve a 2x2 system of equations
// Bezier[4] = f(x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// this is a complete rewrite of the code in t2_vs_t1.iterate_at_d1_d2

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BezierCubic.java

import java.io.FileWriter;

public class BezierCubic
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI/4; // Math.PI;
    public static final int N = 100;
    public static double[] Bezx;                // cubic Bezier, 4 points, x component
    public static double[] Bezy;                // cubic Bezier, 4 points, y component
    //private static CircleFxn fitted;
    //private static CycloidFxn fitted;
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[2][N+1];    // partial wrt (d1, d2)
    private static double Jacdet = Double.NaN;
    private static double eig0 = Double.NaN, eig1 = Double.NaN;
    private static double eigangle = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 53;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //fitted = new epiTrochoidFxn(2.5);
        //iterate_at_P2(55.6, 34.3);
        fitted = new epiTrochoidFxn(3.597585);
        iterate_at_P2(57.88557853692204, 30.031054374161535);
        //System.out.println("cubic Bezier solve_at_P2 = " + solve_at_P2(60, 40, true) + "\n");
        //calc_array();               // generate Python 2D contour plot of rms
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
    }

    private static void iterate_at_P2(double d1, double d2)
    {
        // calculate a new estimate of (d1, d2) by setting dF = 0
        // include only first-order responses
        // see Spiro2SVG Book 3, page 54 (applied to quartic Bezier)
        // setup 2-variable Newton-Raphson iteration

        final double gain = 1;                              // fudge factor to reduce gain
        final int MAXLOOP = 1000;
        double[] f_gx = new double[N+1];
        double[] f_gy = new double[N+1];
        double[] dfxdu = new double[N+1];
        double[] dfydu = new double[N+1];
        double[] d2fxdudu = new double[N+1];
        double[] d2fydudu = new double[N+1];
        double[][] dfxdd = new double[2][N+1];
        double[][] dfydd = new double[2][N+1];
        double[][] d2fxdudd = new double[2][N+1];
        double[][] d2fydudd = new double[2][N+1];
        double[] df_gxdc = new double[N+1];                 // used only for calc of dx[]/dc
        double[] df_gydc = new double[N+1];
//        double[] d2fxdudc = new double[N+1];                // used only for calc of dx[]/dc
//        double[] d2fydudc = new double[N+1];

        double[] trap_in = new double[N+1];
        double[][] Jac = new double[2][2];
        double[] dFdd = new double[2];
        double[] d2Fdddc = new double[2];                           // augmented matrix
        double dFdc;
        double d2Fdcdc;                                             // augmented matrix
        double[][] Augment = new double[3][3];
        double[] deld;                                              // (-Δd1, -Δd2)
        int i, j, k, loop = 0;
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(d1, d2, false)))           // initiallize at (x2, y2)
            {
                System.out.println("fail at (d1, d2): " + d1 + ", " + d2);
                return;
            }

            for (i = 0; i <= N; i++)
            {
                t1 = t1_start + i*(t1_end - t1_start)/N;
                f_gx[i] = t2_vs_t1.fn(Bezx, t2[i]) - fitted.getx(t1);
                f_gy[i] = t2_vs_t1.fn(Bezy, t2[i]) - fitted.gety(t1);
                dfxdu[i] = t2_vs_t1.dfn(Bezx, t2[i]);
                d2fxdudu[i] = t2_vs_t1.d2fn(Bezx, t2[i]);
                dfydu[i] = t2_vs_t1.dfn(Bezy, t2[i]);
                d2fydudu[i] = t2_vs_t1.d2fn(Bezy, t2[i]);
                df_gxdc[i] = calc_df_gxdc(t1, t2[i]);
                df_gydc[i] = calc_df_gydc(t1, t2[i]);
//                d2fxdudc[i] = calc_d2fxdudc(t2[i]);           // set to zero as a test
//                d2fydudc[i] = calc_d2fydudc(t2[i]);           // set to zero as a test

                dfxdd[0][i] = Math.cos(theta_start)*N33(t2[i])[1];
                dfydd[0][i] = Math.sin(theta_start)*N33(t2[i])[1];
                dfxdd[1][i] = -Math.cos(theta_end)*N33(t2[i])[2];
                dfydd[1][i] = -Math.sin(theta_end)*N33(t2[i])[2];
                d2fxdudd[0][i] = Math.cos(theta_start)*dN33(t2[i])[1];
                d2fydudd[0][i] = Math.sin(theta_start)*dN33(t2[i])[1];
                d2fxdudd[1][i] = -Math.cos(theta_end)*dN33(t2[i])[2];
                d2fydudd[1][i] = -Math.sin(theta_end)*dN33(t2[i])[2];
            }

            // calc dFdd[j] at current (d1, d2)

            for (i = 0; i < 2; i++)
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]); // original code
                    //trap_in[k] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];                   // new code
                dFdd[i] = t2_vs_t1.integrate(trap_in);
                for (k = 0; k <= N; k++)
//                    trap_in[k] = (df_gxdc[k] + dfxdu[k]*calc_t2dxy(k, t2[k], "c"))*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k])    // original code
//                               + f_gx[k]*(d2fxdudc[k] + d2fxdudu[k]*calc_t2dxy(k, t2[k], "c"))*t2dd[i][k]
//                               + f_gx[k]*d2fxdudd[i][k]*calc_t2dxy(k, t2[k], "c")
//                               + (df_gydc[k] + dfydu[k]*calc_t2dxy(k, t2[k], "c"))*(dfydd[i][k] + dfydu[k]*t2dd[i][k])
//                               + f_gy[k]*(d2fydudc[k] + d2fydudu[k]*calc_t2dxy(k, t2[k], "c"))*t2dd[i][k]
//                               + f_gy[k]*d2fydudd[i][k]*calc_t2dxy(k, t2[k], "c");
                    trap_in[k] = (df_gxdc[k] + dfxdu[k]*calc_t2dxy(k, t2[k], "c"))*dfxdd[i][k]  // new code
                               + f_gx[k]*d2fxdudd[i][k]*calc_t2dxy(k, t2[k], "c")
                               + (df_gydc[k] + dfydu[k]*calc_t2dxy(k, t2[k], "c"))*dfydd[i][k]
                               + f_gy[k]*d2fydudd[i][k]*calc_t2dxy(k, t2[k], "c");
                d2Fdddc[i] = t2_vs_t1.integrate(trap_in);               // augmented matrix
            }

            // calc d2Fdd[i]dd[j] (Jacobean matrix)

            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N; k++)
                    {
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k] + dfydd[i][k]*dfydd[j][k]
                                   - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k];
                        //System.out.println(k + ", " + trap_in[k]);
                        //System.out.println(k + ", " + (dfxdd[i][k]*dfxdd[j][k] + dfydd[i][k]*dfydd[j][k]));
                        //System.out.println(k + ", " + ((dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k]));
                    }
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }

            // calculate determinant of augmented matrix

            for (k = 0; k <= N; k++)
                trap_in[k] = f_gx[k]*df_gxdc[k] + f_gy[k]*df_gydc[k];
            dFdc = t2_vs_t1.integrate(trap_in);
            for (k = 0; k <= N; k++)
                trap_in[k] = df_gxdc[k]*df_gxdc[k]
                           + (df_gxdc[k]*dfxdu[k] + f_gx[k]*calc_d2fxdudc(t2[k]))*calc_t2dxy(k, t2[k], "c")
                           + df_gydc[k]*df_gydc[k]
                           + (df_gydc[k]*dfydu[k] + f_gy[k]*calc_d2fydudc(t2[k]))*calc_t2dxy(k, t2[k], "c");
            d2Fdcdc = t2_vs_t1.integrate(trap_in);
            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                    Augment[i][j] = Jac[i][j];
            for (i = 0; i < 2; i++)
            {
                Augment[i][2] = d2Fdddc[i];
                Augment[2][i] = Augment[i][2];
            }
            Augment[2][2] = d2Fdcdc;
            //System.out.println("dFdc = " + fitted.getc() + ", " + d1 + ", " + d2 + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);

            //deld = BSpline5.multmv(BSpline5.invertm(Jac), dFdd);  // this is actually the negative of Δd
            deld = BSpline5.gaussj(Jac, dFdd);                      // this is actually the negative of Δd
            d1 -= deld[0]/gain;
            d2 -= deld[1]/gain;

            Jacdet = BSpline5.detm(Jac);
            eig0 = (Jac[0][0] + Jac[1][1] - Math.sqrt((Jac[0][0] - Jac[1][1])*(Jac[0][0] - Jac[1][1]) + 4*Jac[0][1]*Jac[0][1]))/2;
            eig1 = (Jac[0][0] + Jac[1][1] + Math.sqrt((Jac[0][0] - Jac[1][1])*(Jac[0][0] - Jac[1][1]) + 4*Jac[0][1]*Jac[0][1]))/2;
            eigangle = Math.atan((Jac[0][0] - eig0)/Jac[0][1]);     // angle of eigenvector transform
            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + Jacdet + ", " + eig0 + ", " + eig1 + ", " + eigangle + ", " + BSpline5.detm(Augment));
            System.out.println("deld = " + deld[0] + ", " + deld[1]);
            //BSpline5.dump_Jac(Jac);
            //double[][] invJac = BSpline5.invertm(Jac);
            //System.out.println("invJac = " + invJac[0][0] + ", " + invJac[0][1] + ", " + invJac[1][0] + ", " + invJac[1][1]);
            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)            // disabled due to crashes
            {
                t2[i] -= t2dd[0][i]*deld[0] + t2dd[1][i]*deld[1];   // first-order response
                // System.out.println("incrementing t2 array : " + i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 = , , , , , , " + d1 + ", " + d2 + ", " + fitted.getc() + ", " + Jac[0][0] + ", " + Jac[1][1] + ", " + Jac[0][1]);
            solve_at_P2(d1, d2, true);                          // final run just for good measure
            // calculate dddc[i]
            double[] dt2dc = BSpline5.gaussj(Jac, d2Fdddc);
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (Math.cos(eigangle)*d2Fdddc[0] - Math.sin(eigangle)*d2Fdddc[1]) + ", " + (float) (Math.sin(eigangle)*d2Fdddc[0] + Math.cos(eigangle)*d2Fdddc[1]));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (eigangle*180.0/Math.PI) + ", " + (float) (Math.atan(d2Fdddc[0]/d2Fdddc[1])*180.0/Math.PI));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) -d2Fdddc[0] + ", " + (float) -d2Fdddc[1]);
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (Jac[1][0]/Jac[0][0]) + ", " + (float) (Jac[1][1]/Jac[0][1]) + ", " + (float) (d2Fdddc[1]/d2Fdddc[0]));
            System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) BSpline5.detm(Augment));
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ")");
    }

    private static double solve_at_P2(double d1, double d2, boolean print)
    {
        // 4-point cubic Bezier curve
        // perform a single calculation of a complete t2[] profile
        // and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        //d2 += 0.01;
        Bezx = new double[] {fitted.getx(t1_start),
                             fitted.getx(t1_start) + d1*Math.cos(theta_start),
                             fitted.getx(t1_end) - d2*Math.cos(theta_end),
                             fitted.getx(t1_end)};
        Bezy = new double[] {fitted.gety(t1_start),
                             fitted.gety(t1_start) + d1*Math.sin(theta_start),
                             fitted.gety(t1_end) - d2*Math.sin(theta_end),
                             fitted.gety(t1_end)};

        if (t2[N] == 0)
            System.out.println("__start cubic Bezier at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + calc_error());
        if (d1 < 0 || d2 < 0)
            System.out.println("WARNING: negative arm length = " + d1 + ", " + d2);
        //fitted.gen_Bezier(new double[] {Bezx[0], Bezy[0], Bezx[1], Bezy[1], Bezx[2], Bezy[2], Bezx[3], Bezy[3]});
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);

        if (print) System.out.println("\n           , t1, t2, t2dd1, t2dd2");
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i);
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(1 - t2[i]) > TOL)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //for (int ii = 0; ii <= i; ii++)
                //    System.out.println(ii + ", " + t2[ii]);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dxy(i, t2[i], "d1");
            t2dd[1][i] = calc_t2dxy(i, t2[i], "d2");
            if (print)
            {
                //double t1 = t1_start + i*(t1_end - t1_start)/N;
                //double f = (t2_vs_t1.fn(Bezx, t2[i]) - fitted.getx(t1))*t2_vs_t1.dfn(Bezx, t2[i]) + (t2_vs_t1.fn(Bezy, t2[i]) - fitted.gety(t1))*t2_vs_t1.dfn(Bezy, t2[i]);
                //System.out.println("cubic Bez, " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + f);
                //System.out.println("cubic Bez, " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + (t2_vs_t1.fn(Bezx, t2[i]) - fitted.getx(t1)) + ", " + (t2_vs_t1.fn(Bezy, t2[i]) - fitted.gety(t1)));
                System.out.println("cubic Bez, " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + calc_t2dxy(i, t2[i], "c"));
            }
        }
        double retVal = calc_error();
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + retVal + ", " + (float) Jacdet + ", " + (float) eig0 + ", " + (float) eig1 + ", " + (float) eigangle);
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

        if ((Math.abs(1 - t2[N]) > TOL) || (Math.abs(t2[0]) > TOL))
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i <= N; i++)
        {
            trap_in[i] = (t2_vs_t1.fn(Bezx, t2[i]) - fitted.getx(t1))*(t2_vs_t1.fn(Bezx, t2[i]) - fitted.getx(t1))
                       + (t2_vs_t1.fn(Bezy, t2[i]) - fitted.gety(t1))*(t2_vs_t1.fn(Bezy, t2[i]) - fitted.gety(t1));
            //System.out.println(i + ", " + ", " + (fn(Bezx, t2[i]) - fitted.getx(t1)) + ", " + (fn(Bezy, t2[i]) - fitted.gety(t1)) + ", " + Math.sqrt(trap_in[i]));
            t1 += (t1_end - t1_start)/N;
        }
        return Math.sqrt(t2_vs_t1.integrate(trap_in))/a_b;
    }

    private static void calc_array()
    {
        // generate data suitable for MATPLOTLIB
        // rms error as function of d1, d2 (copied from Beziererror.java)
        // to be used by \APP\MATPLOTLIB\mycontour.py

        double d1, d2;
        double a_b = 180;                               // header info
        double d1start = 57.8235, d2start = 30.1442;
        double d1end = 57.3829, d2end = 30.9519;
        int steps = 64;                                 // number of length segments (even)
        int buttsteps = 32;                              // +/- steps added on to length
        double width = .0025;                           // dimensionless relative to length

        double length = Math.sqrt((d1end - d1start)*(d1end - d1start) + (d2end - d2start)*(d2end - d2start));
        double theta = Math.atan2(d2end - d2start, d1end - d1start);

        try
        {
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\CubicBezier\\scan_error.txt");
            java.text.DateFormat df = java.text.DateFormat.getInstance();
            out.write(df.format(new java.util.Date()) + " : a-b, c = " + a_b + ", " + fitted.getc() + "\r\n");
            out.write("extent = (, " + (-width/2) + ", " + (width/2) + ", " + (-(float)buttsteps/steps) + ", " + (1 + (float)buttsteps/steps) + ",) step size = " + 1./steps + ", width = " + width + "\r\n");
            out.write(String.format("(%f, %f)-(%f, %f) in %d", d1start, d2start, d1end, d2end, steps));
            out.write("\r\n");
            for (int i = -buttsteps; i <= steps + buttsteps; i++)
            {
                for (int j = -steps/2; j <= steps/2; j++)
                {
                    d1 = d1start + i*length*Math.cos(theta)/steps + j*length*width*Math.cos(theta - Math.PI/2)/steps;
                    d2 = d2start + i*length*Math.sin(theta)/steps + j*length*width*Math.sin(theta - Math.PI/2)/steps;
                    out.write(String.format(" %g", solve_at_P2(d1, d2, false)));
                    //out.write(String.format(" %.8f", solve_at_P2(d1, d2, false)));
                }
                out.write("\r\n");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("calc_array() save error = " + e);}
    }

    private static void solve_quintic_for_t2(int i)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the cubic Bezier t-value

        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double f, fprime, f2prime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double t;
        int loop = 0;
        int success = 0;

        // initial estimate using quadratic approximation

        if (i == 0) t = 0;
        else t = t2[i-1];
        f = (t2_vs_t1.fn(Bezx, t) - X)*t2_vs_t1.dfn(Bezx, t) + (t2_vs_t1.fn(Bezy, t) - Y)*t2_vs_t1.dfn(Bezy, t);
        fprime = t2_vs_t1.dfn(Bezx, t)*t2_vs_t1.dfn(Bezx, t) + (t2_vs_t1.fn(Bezx, t) - X)*t2_vs_t1.d2fn(Bezx, t) + t2_vs_t1.dfn(Bezy, t)*t2_vs_t1.dfn(Bezy, t) + (t2_vs_t1.fn(Bezy, t) - Y)*t2_vs_t1.d2fn(Bezy, t);
        f2prime = 3*t2_vs_t1.dfn(Bezx, t)*t2_vs_t1.d2fn(Bezx, t) + (t2_vs_t1.fn(Bezx, t) - X)*t2_vs_t1.d3fn(Bezx, t) + 3*t2_vs_t1.dfn(Bezy, t)*t2_vs_t1.d2fn(Bezy, t) + (t2_vs_t1.fn(Bezy, t) - Y)*t2_vs_t1.d3fn(Bezy, t);
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
            del_t = -fprime/f2prime;
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            if (del_t < 0 || del_t > 3)
            {
                del_t = (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
                if (del_t < 0 || del_t > 3)
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
            f = (t2_vs_t1.fn(Bezx, t) - X)*t2_vs_t1.dfn(Bezx, t) + (t2_vs_t1.fn(Bezy, t) - Y)*t2_vs_t1.dfn(Bezy, t);
            fprime = t2_vs_t1.dfn(Bezx, t)*t2_vs_t1.dfn(Bezx, t) + (t2_vs_t1.fn(Bezx, t) - X)*t2_vs_t1.d2fn(Bezx, t) + t2_vs_t1.dfn(Bezy, t)*t2_vs_t1.dfn(Bezy, t) + (t2_vs_t1.fn(Bezy, t) - Y)*t2_vs_t1.d2fn(Bezy, t);
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
            if (Math.abs(del_t) < TOL) success++;
        } while (success < 2);
        //} while (Math.abs(del_t) > TOL);
        t2[i] = t;
    }

    private static double calc_t2dxy(int i, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double denom = t2_vs_t1.dfn(Bezx, t2)*t2_vs_t1.dfn(Bezx, t2) + (t2_vs_t1.fn(Bezx, t2) - X)*t2_vs_t1.d2fn(Bezx, t2) + t2_vs_t1.dfn(Bezy, t2)*t2_vs_t1.dfn(Bezy, t2) + (t2_vs_t1.fn(Bezy, t2) - Y)*t2_vs_t1.d2fn(Bezy, t2);
        double numer = Double.NaN;

        if (type.equals("d1"))
            numer = (Math.cos(theta_start)*t2_vs_t1.dfn(Bezx, t2) + Math.sin(theta_start)*t2_vs_t1.dfn(Bezy, t2))*N33(t2)[1]
                  + (Math.cos(theta_start)*(t2_vs_t1.fn(Bezx, t2) - X) + Math.sin(theta_start)*(t2_vs_t1.fn(Bezy, t2) - Y))*dN33(t2)[1];
        else if (type.equals("d2"))
            numer = -(Math.cos(theta_end)*t2_vs_t1.dfn(Bezx, t2) + Math.sin(theta_end)*t2_vs_t1.dfn(Bezy, t2))*N33(t2)[2]
                  -  (Math.cos(theta_end)*(t2_vs_t1.fn(Bezx, t2) - X) + Math.sin(theta_end)*(t2_vs_t1.fn(Bezy, t2) - Y))*dN33(t2)[2];
        else if (type.equals("c"))
            numer = calc_df_gxdc(t1, t2)*t2_vs_t1.dfn(Bezx, t2) + (t2_vs_t1.fn(Bezx, t2) - X)*calc_d2fxdudc(t2)
                  + calc_df_gydc(t1, t2)*t2_vs_t1.dfn(Bezy, t2) + (t2_vs_t1.fn(Bezy, t2) - Y)*calc_d2fydudc(t2);
        return -numer/denom;
    }

    private static double calc_df_gxdc(double t1, double t2)
    {
        return fitted.getdxdc(t1_start)*(N33(t2)[0] + N33(t2)[1]) + fitted.getdxdc(t1_end)*(N33(t2)[2] + N33(t2)[3]) - fitted.getdxdc(t1);
    }

    private static double calc_df_gydc(double t1, double t2)
    {
        return fitted.getdydc(t1_start)*(N33(t2)[0] + N33(t2)[1]) + fitted.getdydc(t1_end)*(N33(t2)[2] + N33(t2)[3]) - fitted.getdydc(t1);
    }

    private static double calc_d2fxdudc(double t2)
    {
        //return 0;                                                   // fix fix test code
        return fitted.getdxdc(t1_start)*(dN33(t2)[0] + dN33(t2)[1]) + fitted.getdxdc(t1_end)*(dN33(t2)[2] + dN33(t2)[3]);
    }

    private static double calc_d2fydudc(double t2)
    {
        //return 0;                                                   // fix fix test code
        return fitted.getdydc(t1_start)*(dN33(t2)[0] + dN33(t2)[1]) + fitted.getdydc(t1_end)*(dN33(t2)[2] + dN33(t2)[3]);
    }

    private static double[] N33(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:   N33 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:   N33 too large u value = " + u);
        return new double[] {(1 - u)*(1 - u)*(1 - u),
                             3*u*(1 - u)*(1 - u),
                             3*u*u*(1 - u),
                             u*u*u};
    }

    private static double[] dN33(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:  dN33 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:  dN33 too large u value = " + u);
        return new double[] {-3*(1 - u)*(1 - u),
                              3*(1 - u)*(1 - 3*u),
                              3*u*(2 - 3*u),
                              3*u*u};
    }
}
