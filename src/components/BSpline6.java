
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 6-point B-Spline (P0 - P5) to it, using parameter 0 < t2 < 3.
// constrain only the slopes at the endpoints and keep P2, P3 arbitrary
// linearize the equations wrt (d1, d2, x2, y2, x3, y3) and solve a 6x6 system of equations
// decompose the B-Spline into 3 Beziers, range (0,1), (1,2), and (2,3).
// Bezier[3][4] = f[3](x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book5, Jan 2018, page 1

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BSpline6.java

public class BSpline6
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI/4; // Math.PI/4;
    public static final int N = 100;
    public static double[] Splinex, Spliney;    // 6 point spline
    public static double[][] Bezx;              // 3 Beziers, 4 points each, x component
    public static double[][] Bezy;              // 3 Beziers, 4 points each, y component
    //private static CircleFxn fitted;
    //private static CycloidFxn fitted;
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[6][N+1];    // partial wrt (d1, d2, x2, y2, x3, y3)
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 90;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //fitted = new epiTrochoidFxn(10);
        //System.out.println("BSpline6 iterate_at_P2_P3 = " + iterate_at_P2_P3(15, 20, 183, 38, 158, 79) + "\n");
        fitted = new epiTrochoidFxn(0.);
        System.out.println("BSpline6 iterate_at_P2_P3 = " +
        iterate_at_P2_P3(16.624489486951397, 14.946247668118907, 176.57737789450186, 48.59041870307481, 155.94233690458145, 91.34718144983118) + "\n");
        //fitted = new epiTrochoidFxn(18.44);
        //System.out.println("BSpline6 convert_to_P2_P3 = " + convert_to_P2_P3(9.016620273713244, 75.83622913947018, true) + "\n");
        //System.out.println("BSpline6 iterate_at_P2_P3 = " + iterate_at_P2_P3(10.432741575200822, 13.925540044687363, 191.17890042589667, 30.996508335507055, 160.1186993640138, 68.3685414189886) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //for (int i = 0; i <= 300; i++)
        //    System.out.println(i/100.0 + ", " + dN53(i/100.)[0] + ", " + dN53(i/100.)[1] + ", " + dN53(i/100.)[2] + ", " + dN53(i/100.)[3] + ", " + dN53(i/100.)[4] + ", " + dN53(i/100.)[5]);
    }

    private static double iterate_at_P2_P3(double d1, double d2, double x2, double y2, double x3, double y3)
    {
        // calculate a new estimate of (d1, d2, x2, y2, x3, y3) by setting dF = 0
        // include only first-order responses
        // setup 6-variable Newton-Raphson iteration

        final double gain = 1;                              // fudge factor to reduce gain
        final int MAXLOOP = 1000;
        double[] f_gx = new double[N+1];
        double[] f_gy = new double[N+1];
        double[] dfxdu = new double[N+1];
        double[] d2fxdudu = new double[N+1];
        double[] dfydu = new double[N+1];
        double[] d2fydudu = new double[N+1];
        double[][] dfxdd = new double[6][N+1];
        double[][] dfydd = new double[6][N+1];
        double[][] Jac = new double[6][6];
        double[] dFdd = new double[6];
        double[] trap_in = new double[N+1];
        double[] deld;                                      // (-Δd1, -Δd2, -Δx2, -Δy2, -Δx3, -Δy3)
        int loop = 0;
        int i, j, k, seg;                                  // Bezier segment, before or after the splice
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2_P3(d1, d2, x2, y2, x3, y3, false)))   // initiallize at (P2, P3)
            {
                System.out.println("fail at " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
                return Double.NaN;
            }

            seg = 0;
            for (i = 0; i <= N; i++)
            {
                if (seg == 0 && t2[i] > 1)
                    seg++;
                if (seg == 1 && t2[i] > 2)
                    seg++;
                t1 = t1_start + i*(t1_end - t1_start)/N;
                f_gx[i] = t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1);
                f_gy[i] = t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1);
                dfxdu[i] = t2_vs_t1.dfn(Bezx[seg], t2[i] - seg);
                d2fxdudu[i] = t2_vs_t1.d2fn(Bezx[seg], t2[i] - seg);
                dfydu[i] = t2_vs_t1.dfn(Bezy[seg], t2[i] - seg);
                d2fydudu[i] = t2_vs_t1.d2fn(Bezy[seg], t2[i] - seg);
                dfxdd[0][i] = calc_dfxdd1(t2[i]);
                dfydd[0][i] = calc_dfydd1(t2[i]);
                dfxdd[1][i] = calc_dfxdd2(t2[i]);
                dfydd[1][i] = calc_dfydd2(t2[i]);
                dfxdd[2][i] = calc_dfxdx2(t2[i]);
                dfydd[2][i] = 0;
                dfxdd[3][i] = 0;
                dfydd[3][i] = calc_dfydy2(t2[i]);
                dfxdd[4][i] = calc_dfxdx3(t2[i]);
                dfydd[4][i] = 0;
                dfxdd[5][i] = 0;
                dfydd[5][i] = calc_dfydy3(t2[i]);
                //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + t2dd[5][i] + ", " + f_gx[i] + ", " + f_gy[i] + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + dfxdd[0][i] + ", " + dfydd[0][i] + ", " + dfxdd[1][i] + ", " + dfydd[1][i] + ", " + dfxdd[2][i] + ", " + dfydd[2][i] + ", " + dfxdd[3][i] + ", " + dfydd[3][i] + ", " + dfxdd[4][i] + ", " + dfydd[4][i] + ", " + dfxdd[5][i] + ", " + dfydd[5][i] + ", " + d2fxdudd[0][i] + ", " + d2fydudd[0][i] + ", " + d2fxdudd[1][i] + ", " + d2fydudd[1][i] + ", " + d2fxdudd[2][i] + ", " + d2fydudd[2][i] + ", " + d2fxdudd[3][i] + ", " + d2fydudd[3][i] + ", " + d2fxdudd[4][i] + ", " + d2fydudd[4][i] + ", " + d2fxdudd[5][i] + ", " + d2fydudd[5][i]);
            }

            // calc dFdd[j] at current (d1, d2, x2, y2, x3, y3)

            for (i = 0; i < 6; i++)
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]); // original code
                    //trap_in[k] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];         // new code
                dFdd[i] = t2_vs_t1.integrate(trap_in);
            }

            // calc d2Fdd[i]dd[j] (symmetric Jacobean matrix)

            for (i = 0; i < 6; i++)
                for (j = 0; j < 6; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N; k++)
                    {
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k] + dfydd[i][k]*dfydd[j][k]
                                   - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k];
                        //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }

            deld = BSpline5.multmv(BSpline5.invertm(Jac), dFdd);  // this is actually the negative of Δd
            d1 -= deld[0]/gain;                 // gain is just a fudge factor to 'improve' convergence
            d2 -= deld[1]/gain;
            x2 -= deld[2]/gain;
            y2 -= deld[3]/gain;
            x3 -= deld[4]/gain;
            y3 -= deld[5]/gain;

            //System.out.println("Jac");
            //for (i = 0; i < Jac.length; i++)
            //{
            //    for (j = 0; j < Jac.length; j++)
            //        System.out.print(Jac[i][j] + ", ");
            //    System.out.println();
            //}
            Jacdet = BSpline5.detm(Jac);
            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + dFdd[2] + ", " + dFdd[3] + ", " + dFdd[4] + ", " + dFdd[5] + ", " + Jacdet);
            System.out.println("deld = " + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ", " + deld[4] + ", " + deld[5]);
            BSpline5.dump_Jac(Jac);

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                for (j = 0; j < 6; j++)
                    t2[i] -= t2dd[j][i]*deld[j];                        // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL) && (Math.abs(deld[2]) < TOL) && (Math.abs(deld[3]) < TOL) && (Math.abs(deld[4]) < TOL) && (Math.abs(deld[5]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 x2 y2 x3 y3 = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
            return solve_at_P2_P3(d1, d2, x2, y2, x3, y3, true);                    // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ", " + deld[4] + ", " + deld[5] + ")");
        return Double.NaN;
    }

    private static double convert_to_P2_P3(double d1, double d2, boolean print)
    {
        // convert a 4-point cubic Bezier to a 6-point B-Spline
        System.out.println("B-Spline6 convert_to_P2_P3: Bezier d1 d2 = , , , , , , " + d1 + ", " + d2);
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        double x0 = fitted.getx(t1_start);
        double x1 = x0 + d1*Math.cos(theta_start);
        double x3 = fitted.getx(t1_end);
        double x2 = x3 - d2*Math.cos(theta_end);
        double y0 = fitted.gety(t1_start);
        double y1 = y0 + d1*Math.sin(theta_start);
        double y3 = fitted.gety(t1_end);
        double y2 = y3 - d2*Math.sin(theta_end);
        //fitted.gen_Bezier(new double[] {x0, y0, x1, y1, x2, y2, x3, y3});

        return iterate_at_P2_P3(d1/3, d2/3, (2*x0 + 5*x1 + 2*x2)/9, (2*y0 + 5*y1 + 2*y2)/9, (2*x1 + 5*x2 + 2*x3)/9, (2*y1 + 5*y2 + 2*y3)/9);
        //return solve_at_P2_P3(d1/3, d2/3, (2*x0 + 5*x1 + 2*x2)/9, (2*y0 + 5*y1 + 2*y2)/9, (2*x1 + 5*x2 + 2*x3)/9, (2*y1 + 5*y2 + 2*y3)/9, print);
    }

    private static double solve_at_P2_P3(double d1, double d2, double x2, double y2, double x3, double y3, boolean print)
    {
        // 6-point cubic B-Spline curve
        // perform a single calculation of a complete t2[] profile
        // at a given P2, P3, and calculate the rms error

        //y3 -= 0.001;
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        Splinex = new double[] {fitted.getx(t1_start),
                                fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                x2,
                                x3,
                                fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                fitted.getx(t1_end)};
        Spliney = new double[] {fitted.gety(t1_start),
                                fitted.gety(t1_start) + d1*Math.sin(theta_start),
                                y2,
                                y3,
                                fitted.gety(t1_end) - d2*Math.sin(theta_end),
                                fitted.gety(t1_end)};
        Bezx = new double[][] {{Splinex[0],
                                Splinex[1],
                                (Splinex[1] + Splinex[2])/2,
                                (3*Splinex[1] + 7*Splinex[2] + 2*Splinex[3])/12},
                               {(3*Splinex[1] + 7*Splinex[2] + 2*Splinex[3])/12,
                                (2*Splinex[2] + Splinex[3])/3,
                                (Splinex[2] + 2*Splinex[3])/3,
                                (2*Splinex[2] + 7*Splinex[3] + 3*Splinex[4])/12},
                               {(2*Splinex[2] + 7*Splinex[3] + 3*Splinex[4])/12,
                                (Splinex[3] + Splinex[4])/2,
                                Splinex[4],
                                Splinex[5]}};
        Bezy = new double[][] {{Spliney[0],
                                Spliney[1],
                                (Spliney[1] + Spliney[2])/2,
                                (3*Spliney[1] + 7*Spliney[2] + 2*Spliney[3])/12},
                               {(3*Spliney[1] + 7*Spliney[2] + 2*Spliney[3])/12,
                                (2*Spliney[2] + Spliney[3])/3,
                                (Spliney[2] + 2*Spliney[3])/3,
                                (2*Spliney[2] + 7*Spliney[3] + 3*Spliney[4])/12},
                               {(2*Spliney[2] + 7*Spliney[3] + 3*Spliney[4])/12,
                                (Spliney[3] + Spliney[4])/2,
                                Spliney[4],
                                Spliney[5]}};
        if (t2[N] == 0)
            System.out.println("__start B-Spline6 theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3 + ", " + calc_error());
        //fitted.gen_Bezier3(Bezx, Bezy);
        //for (int i = 0; i < 3; i++)
        //    for (int j = 0; j < 4; j++)
        //        System.out.println(Bezx[i][j] + "\t " + Bezy[i][j]);

        if (print) System.out.println("\nseg, t1, t2, t2dd1, t2dd2, t2dx2, t2dy2, t2dx3, t2dy3");
        int seg = 0;                // Bezier segment, 0 - 2
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i, seg);
            if (seg == 0 && t2[i] > 1)
            {
                seg++;
                solve_quintic_for_t2(i, seg);   // re-calculate
            }
            if (seg == 1 && t2[i] < 1)
            {
                System.out.println("seg reverted from 1 to 0 at i = " + i + ", " + t2[i]);
                return Double.NaN;
            }
            if (seg == 1 && t2[i] > 2)
            {
                seg++;
                solve_quintic_for_t2(i, seg);   // re-calculate
            }
            if (seg == 2 && t2[i] < 2)
            {
                System.out.println("seg reverted from 2 to 1 at i = " + i + ", " + t2[i]);
                return Double.NaN;
            }
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(3 - t2[i]) > TOL)
            ||  (i == N && seg != 2)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //scan_quintic_near_t2(i, seg, t2[i]);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dd(i, seg, t2[i], "d1");
            t2dd[1][i] = calc_t2dd(i, seg, t2[i], "d2");
            t2dd[2][i] = calc_t2dd(i, seg, t2[i], "x2");
            t2dd[3][i] = calc_t2dd(i, seg, t2[i], "y2");
            t2dd[4][i] = calc_t2dd(i, seg, t2[i], "x3");
            t2dd[5][i] = calc_t2dd(i, seg, t2[i], "y3");
            if (print)
                System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + t2dd[5][i]);
        }
        double retVal = calc_error();
        System.out.println("__new t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3 + ", " + (float) retVal + ", " + (float) Jacdet);
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

        if ((Math.abs(3 - t2[N]) > TOL) || (Math.abs(t2[0]) > TOL))
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i <= N; i++)
        {
            if (seg == 0 && t2[i] > 1)
                seg++;
            if (seg == 1 && t2[i] > 2)
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

    private static double calc_t2dd(int i, int seg, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        double denom = t2_vs_t1.dfn(Bezx[seg], t2 - seg)*t2_vs_t1.dfn(Bezx[seg], t2 - seg) + (t2_vs_t1.fn(Bezx[seg], t2 - seg) - X)*t2_vs_t1.d2fn(Bezx[seg], t2 - seg) + t2_vs_t1.dfn(Bezy[seg], t2 - seg)*t2_vs_t1.dfn(Bezy[seg], t2 - seg) + (t2_vs_t1.fn(Bezy[seg], t2 - seg) - Y)*t2_vs_t1.d2fn(Bezy[seg], t2 - seg);
        double numer = Double.NaN;
        if (type.equals("d1"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2 - seg)*calc_dfxdd1(t2) + t2_vs_t1.dfn(Bezy[seg], t2 - seg)*calc_dfydd1(t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2 - seg) - X)*calc_d2fxdudd1(t2) + (t2_vs_t1.fn(Bezy[seg], t2 - seg) - Y)*calc_d2fydudd1(t2);
        else if (type.equals("d2"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2 - seg)*calc_dfxdd2(t2) + t2_vs_t1.dfn(Bezy[seg], t2 - seg)*calc_dfydd2(t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2 - seg) - X)*calc_d2fxdudd2(t2) + (t2_vs_t1.fn(Bezy[seg], t2 - seg) - Y)*calc_d2fydudd2(t2);
        else if(type.equals("x2"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2 - seg)*calc_dfxdx2(t2) + (t2_vs_t1.fn(Bezx[seg], t2 - seg) - X)*calc_d2fxdudx2(t2);
        else if (type.equals("y2"))
            numer = t2_vs_t1.dfn(Bezy[seg], t2 - seg)*calc_dfydy2(t2) + (t2_vs_t1.fn(Bezy[seg], t2 - seg) - Y)*calc_d2fydudy2(t2);
        else if(type.equals("x3"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2 - seg)*calc_dfxdx3(t2) + (t2_vs_t1.fn(Bezx[seg], t2 - seg) - X)*calc_d2fxdudx3(t2);
        else if (type.equals("y3"))
            numer = t2_vs_t1.dfn(Bezy[seg], t2 - seg)*calc_dfydy3(t2) + (t2_vs_t1.fn(Bezy[seg], t2 - seg) - Y)*calc_d2fydudy3(t2);
        return -numer/denom;
    }

    private static double calc_dfxdd1(double t2)
    {
        return N53(t2)[1]*Math.cos(theta_start);
    }

    private static double calc_dfydd1(double t2)
    {
        return N53(t2)[1]*Math.sin(theta_start);
    }

    private static double calc_d2fxdudd1(double t2)
    {
        return dN53(t2)[1]*Math.cos(theta_start);
    }

    private static double calc_d2fydudd1(double t2)
    {
        return dN53(t2)[1]*Math.sin(theta_start);
    }

    private static double calc_dfxdd2(double t2)
    {
        return -N53(t2)[4]*Math.cos(theta_end);
    }

    private static double calc_dfydd2(double t2)
    {
        return -N53(t2)[4]*Math.sin(theta_end);
    }

    private static double calc_d2fxdudd2(double t2)
    {
        return -dN53(t2)[4]*Math.cos(theta_end);
    }

    private static double calc_d2fydudd2(double t2)
    {
        return -dN53(t2)[4]*Math.sin(theta_end);
    }

    private static double calc_dfxdx2(double t2)
    {
        return N53(t2)[2];
    }

    private static double calc_d2fxdudx2(double t2)
    {
        return dN53(t2)[2];
    }

    private static double calc_dfydy2(double t2)
    {
        return N53(t2)[2];
    }

    private static double calc_d2fydudy2(double t2)
    {
        return dN53(t2)[2];
    }

    private static double calc_dfxdx3(double t2)
    {
        return N53(t2)[3];
    }

    private static double calc_d2fxdudx3(double t2)
    {
        return dN53(t2)[3];
    }

    private static double calc_dfydy3(double t2)
    {
        return N53(t2)[3];
    }

    private static double calc_d2fydudy3(double t2)
    {
        return dN53(t2)[3];
    }

    private static double[] N53(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {(1 - u)*(1 - u)*(1 - u),
                                 u*(12 - 18*u + 7*u*u)/4,
                                 u*u*(18 - 11*u)/12,
                                 u*u*u/6,
                                 0,
                                 0};
        else if (u < 2)
            return new double[] {0,
                                 (2 - u)*(2 - u)*(2 - u)/4,
                                 (7*u*u*u - 36*u*u + 54*u - 18)/12,
                                 (-7*u*u*u + 27*u*u - 27*u + 9)/12,
                                 (u - 1)*(u - 1)*(u - 1)/4,
                                 0};
        else if (u < 3 + TOL)
            return new double[] {0,
                                 0,
                                 (3 - u)*(3 - u)*(3 - u)/6,
                                 (3 - u)*(3 - u)*(11*u - 15)/12,
                                 (3 - u)*(7*u*u - 24*u + 21)/4,
                                 (u - 2)*(u - 2)*(u - 2)};
        else
            return null;
    }

    private static double[] dN53(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {-3*(1 - u)*(1 - u),
                                 (12 - 36*u + 21*u*u)/4,
                                 u*(12 - 11*u)/4,
                                 u*u/2,
                                 0,
                                 0};
        else if (u < 2)
            return new double[] {0,
                                 -3*(2 - u)*(2 - u)/4,
                                 (7*u*u - 24*u + 18)/4,
                                 (-7*u*u + 18*u - 9)/4,
                                 3*(u - 1)*(u - 1)/4,
                                 0};
        else if (u < 3 + TOL)
            return new double[] {0,
                                 0,
                                 -(3 - u)*(3 - u)/2,
                                 -(3 - u)*(11*u - 21)/4,
                                 -3*(7*u*u - 30*u + 31)/4,
                                 3*(u - 2)*(u - 2)};
        else
            return null;
    }
}
