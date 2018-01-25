
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 5-point B-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// choose d1 and d2 to fit the curvature at the endpoints and keep P2 arbitrary
// decompose the B-Spline into 2 Beziers, range (0,1) and (0,2).
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book2, Dec 2016, page 50, revisited in Jan 2018
// this is the second attempt to write this constrained code (Jan 22, 2018)
// the first attempt was in E:components31\BSpline5.java (Oct 29, 2017)

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BSpline5constrain.java

public class BSpline5constrain
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI/4; // Math.PI/4;
    public static final int N = 100;
    public static double[] Splinex, Spliney;    // 5 point spline
    public static double[][] Bezx;              // 2 Beziers, 4 points each, x component
    public static double[][] Bezy;              // 2 Beziers, 4 points each, y component
    //private static CycloidFxn fitted;
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[2][N+1];    // partial wrt (x2, y2)
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 15;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        fitted = new epiTrochoidFxn(6);
        //System.out.println("BSpline5constrain convert_to_P2 = " + convert_to_P2(52.63858334027786, 40.19155607809456, true) + "\n");
        System.out.println("BSpline5constrain iterate_at_P2 = " + iterate_at_P2(183.7343334945213, 44.7816370208017) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
    }

    private static double iterate_at_P2(double x2, double y2)
    {
        // calculate a new estimate of (x2, y2) by setting dF = 0
        // include only first-order responses
        // setup 2-variable Newton-Raphson iteration

        final double gain = 1;                                  // fudge factor to reduce gain
        final int MAXLOOP = 2000;
        double[] f_gx = new double[N+1];
        double[] f_gy = new double[N+1];
        double[] dfxdu = new double[N+1];
        double[] dfydu = new double[N+1];
        double[][] dfxdd = new double[2][N+1];
        double[][] dfydd = new double[2][N+1];
        double[][] d2fxdudd = new double[2][N+1];
        double[][] d2fydudd = new double[2][N+1];
        double[][][] d2fxdddd = new double[2][2][N+1];
        double[][][] d2fydddd = new double[2][2][N+1];
        double[][] Jac = new double[2][2];
        double[] dFdd = new double[2];
        double[] trap_in = new double[N+1];
        double[] deld;                                          // (-Δx2, -Δy2)
        int loop = 0;
        int i, j, k, seg;                                       // Bezier segment, before or after the splice
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(x2, y2, false)))       // initiallize at P2
            {
                System.out.println("fail at " + x2 + ", " + y2);
                return Double.NaN;
            }

            seg = 0;
            for (i = 0; i <= N; i++)
            {
                if (seg == 0 && t2[i] > 1)
                    seg++;
                t1 = t1_start + i*(t1_end - t1_start)/N;
                f_gx[i] = t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1);
                f_gy[i] = t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1);
                dfxdu[i] = t2_vs_t1.dfn(Bezx[seg], t2[i] - seg);
                dfydu[i] = t2_vs_t1.dfn(Bezy[seg], t2[i] - seg);
                dfxdd[0][i] = calc_dfxdx2(t2[i]) + calc_dfxdd1(t2[i])*calc_dd1dx2(x2, y2) + calc_dfxdd2(t2[i])*calc_dd2dx2(x2, y2);
                dfydd[0][i] =                      calc_dfydd1(t2[i])*calc_dd1dx2(x2, y2) + calc_dfydd2(t2[i])*calc_dd2dx2(x2, y2);
                dfxdd[1][i] =                      calc_dfxdd1(t2[i])*calc_dd1dy2(x2, y2) + calc_dfxdd2(t2[i])*calc_dd2dy2(x2, y2);
                dfydd[1][i] = calc_dfydy2(t2[i]) + calc_dfydd1(t2[i])*calc_dd1dy2(x2, y2) + calc_dfydd2(t2[i])*calc_dd2dy2(x2, y2);
                d2fxdudd[0][i] = calc_d2fxdudx2(t2[i]) + calc_d2fxdudd1(t2[i])*calc_dd1dx2(x2, y2) + calc_d2fxdudd2(t2[i])*calc_dd2dx2(x2, y2);
                d2fydudd[0][i] =                         calc_d2fydudd1(t2[i])*calc_dd1dx2(x2, y2) + calc_d2fydudd2(t2[i])*calc_dd2dx2(x2, y2);
                d2fxdudd[1][i] =                         calc_d2fxdudd1(t2[i])*calc_dd1dy2(x2, y2) + calc_d2fxdudd2(t2[i])*calc_dd2dy2(x2, y2);
                d2fydudd[1][i] = calc_d2fydudy2(t2[i]) + calc_d2fydudd1(t2[i])*calc_dd1dy2(x2, y2) + calc_d2fydudd2(t2[i])*calc_dd2dy2(x2, y2);
                d2fxdddd[0][0][i] = calc_dfxdd1(t2[i])*calc_d2d1dx2dx2(x2, y2) + calc_dfxdd2(t2[i])*calc_d2d2dx2dx2(x2, y2);
                d2fxdddd[0][1][i] = calc_dfxdd1(t2[i])*calc_d2d1dx2dy2(x2, y2) + calc_dfxdd2(t2[i])*calc_d2d2dx2dy2(x2, y2);
                d2fxdddd[1][0][i] = d2fxdddd[0][1][i];
                d2fxdddd[1][1][i] = calc_dfxdd1(t2[i])*calc_d2d1dy2dy2(x2, y2) + calc_dfxdd2(t2[i])*calc_d2d2dy2dy2(x2, y2);
                d2fydddd[0][0][i] = calc_dfydd1(t2[i])*calc_d2d1dx2dx2(x2, y2) + calc_dfydd2(t2[i])*calc_d2d2dx2dx2(x2, y2);
                d2fydddd[0][1][i] = calc_dfydd1(t2[i])*calc_d2d1dx2dy2(x2, y2) + calc_dfydd2(t2[i])*calc_d2d2dx2dy2(x2, y2);
                d2fydddd[1][0][i] = d2fydddd[0][1][i];
                d2fydddd[1][1][i] = calc_dfydd1(t2[i])*calc_d2d1dy2dy2(x2, y2) + calc_dfydd2(t2[i])*calc_d2d2dy2dy2(x2, y2);
                //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + f_gx[i] + ", " + f_gy[i] + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + dfxdd[0][i] + ", " + dfydd[0][i] + ", " + dfxdd[1][i] + ", " + dfydd[1][i] + ", " + d2fxdudd[0][i] + ", " + d2fydudd[0][i] + ", " + d2fxdudd[1][i] + ", " + d2fydudd[1][i]);
            }

            // calc dFdd[j] at current (x2, y2)

            for (i = 0; i < 2; i++)
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]);
                dFdd[i] = t2_vs_t1.integrate(trap_in);
            }

            // calc d2Fdd[i]dd[j] (symmetric Jacobean matrix)

            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    for (k = 0; k <= N; k++)
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k] + (dfxdd[j][k]*dfxdu[k] + f_gx[k]*d2fxdudd[j][k])*t2dd[i][k] + f_gx[k]*d2fxdddd[i][j][k]
                                   + dfydd[i][k]*dfydd[j][k] + (dfydd[j][k]*dfydu[k] + f_gy[k]*d2fydudd[j][k])*t2dd[i][k] + f_gy[k]*d2fydddd[i][j][k];
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }

            deld = BSpline5.multmv(BSpline5.invertm(Jac), dFdd);  // this is actually the negative of Δd
            x2 -= deld[0]/gain;                 // gain is just a fudge factor to 'improve' convergence
            y2 -= deld[1]/gain;

            //System.out.println("Jac");
            //for (i = 0; i < Jac.length; i++)
            //{
            //    for (j = 0; j < Jac.length; j++)
            //        System.out.print(Jac[i][j] + ", ");
            //    System.out.println();
            //}

            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + BSpline5.detm(Jac));
            System.out.println("deld = " + deld[0] + ", " + deld[1]);
            //BSpline5.dump_Jac(Jac);

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                for (j = 0; j < 2; j++)
                    t2[i] -= t2dd[j][i]*deld[j];                        // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new x2 y2 = , , , , , , " + x2 + ", " + y2);
            return solve_at_P2(x2, y2, true);                           // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ")");
        return Double.NaN;
    }

    private static double convert_to_P2(double d1, double d2, boolean print)
    {
        // convert a 4-point curvature-fit cubic Bezier to a 5-point constrained B-Spline
        // (the Bezier must be 'curvature-fit' because the endpoint curvature will be constrained at all times)
        System.out.println("B-Spline5constrained convert_to_P2: Bezier d1 d2 = , , , , , , " + d1 + ", " + d2);
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

        return iterate_at_P2((x1 + x2)/2, (y1 + y2)/2);
        //return solve_at_P2((x1 + x2)/2, (y1 + y2)/2, print);
    }

    private static double solve_at_P2(double x2, double y2, boolean print)
    {
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error
        // (calculate d1, d2 to satisfy the curvature at the endpoints)

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);
        double d1 = calc_d1(x2, y2);                // constrain the endpoint curvature
        double d2 = calc_d2(x2, y2);
        if (Double.isNaN(d1) || Double.isNaN(d2))
        {
            System.out.println("Bad d1 d2 = " + d1 + ", " + d2);
            return Double.NaN;
        }

        //y2 += 0.002;
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
            System.out.println("__start at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + calc_error());
        //fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println("d1dx2 = " + calc_dd1dx2(x2, y2) + ", " + calc_dd2dx2(x2, y2) + ", " + calc_dd1dy2(x2, y2) + ", " + calc_dd2dy2(x2, y2));
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2, t2dx2, t2dy2");
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
            {
                System.out.println("seg reverted from 1 to 0 at i = " + i + ", " + t2[i]);
                return Double.NaN;
            }
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
//                scan_quintic_near_t2(t1, t2[i]);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dd(i, seg, t2[i], "x2") + calc_t2dd(i, seg, t2[i], "d1")*calc_dd1dx2(x2, y2) + calc_t2dd(i, seg, t2[i], "d2")*calc_dd2dx2(x2, y2);
            t2dd[1][i] = calc_t2dd(i, seg, t2[i], "y2") + calc_t2dd(i, seg, t2[i], "d1")*calc_dd1dy2(x2, y2) + calc_t2dd(i, seg, t2[i], "d2")*calc_dd2dy2(x2, y2);
            if (print) System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i]);
        }
        double retVal = calc_error();
        System.out.println("__new t2[] at theta c t d1 d2 rms = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + retVal);
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
        return -numer/denom;
    }

    private static double calc_d1(double x2, double y2)
    {
        return Math.sqrt(((y2 - fitted.gety(t1_start))*Math.cos(theta_start) - (x2 - fitted.getx(t1_start))*Math.sin(theta_start))/3/fitted.getkappa(t1_start));
    }

    private static double calc_dd1dx2(double x2, double y2)
    {
        return -Math.sin(theta_start)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start);
    }

    private static double calc_dd1dy2(double x2, double y2)
    {
        return Math.cos(theta_start)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start);
    }

    private static double calc_d2d1dx2dx2(double x2, double y2)
    {
        return Math.sin(theta_start)/calc_d1(x2, y2)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start)*calc_dd1dx2(x2, y2);
    }

    private static double calc_d2d1dy2dy2(double x2, double y2)
    {
        return -Math.cos(theta_start)/calc_d1(x2, y2)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start)*calc_dd1dy2(x2, y2);
    }

    private static double calc_d2d1dx2dy2(double x2, double y2)
    {
        return Math.sin(theta_start)/calc_d1(x2, y2)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start)*calc_dd1dy2(x2, y2);
    }

    private static double calc_d2(double x2, double y2)
    {
        return Math.sqrt(((y2 - fitted.gety(t1_end))*Math.cos(theta_end) - (x2 - fitted.getx(t1_end))*Math.sin(theta_end))/3/fitted.getkappa(t1_end));
    }

    private static double calc_dd2dx2(double x2, double y2)
    {
        return -Math.sin(theta_end)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end);
    }

    private static double calc_dd2dy2(double x2, double y2)
    {
        return Math.cos(theta_end)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end);
    }

    private static double calc_d2d2dx2dx2(double x2, double y2)
    {
        return Math.sin(theta_end)/calc_d2(x2, y2)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end)*calc_dd2dx2(x2, y2);
    }

    private static double calc_d2d2dy2dy2(double x2, double y2)
    {
        return -Math.cos(theta_end)/calc_d2(x2, y2)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end)*calc_dd2dy2(x2, y2);
    }

    private static double calc_d2d2dx2dy2(double x2, double y2)
    {
        return Math.sin(theta_end)/calc_d2(x2, y2)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end)*calc_dd2dy2(x2, y2);
    }

    private static double calc_dfxdd1(double t2)
    {
        return N43(t2)[1]*Math.cos(theta_start);
    }

    private static double calc_dfydd1(double t2)
    {
        return N43(t2)[1]*Math.sin(theta_start);
    }

    private static double calc_d2fxdudd1(double t2)
    {
        return dN43(t2)[1]*Math.cos(theta_start);
    }

    private static double calc_d2fydudd1(double t2)
    {
        return dN43(t2)[1]*Math.sin(theta_start);
    }

    private static double calc_dfxdd2(double t2)
    {
        return -N43(t2)[3]*Math.cos(theta_end);
    }

    private static double calc_dfydd2(double t2)
    {
        return -N43(t2)[3]*Math.sin(theta_end);
    }

    private static double calc_d2fxdudd2(double t2)
    {
        return -dN43(t2)[3]*Math.cos(theta_end);
    }

    private static double calc_d2fydudd2(double t2)
    {
        return -dN43(t2)[3]*Math.sin(theta_end);
    }

    private static double calc_dfxdx2(double t2)
    {
        return N43(t2)[2];
    }

    private static double calc_d2fxdudx2(double t2)
    {
        return dN43(t2)[2];
    }

    private static double calc_dfydy2(double t2)
    {
        return N43(t2)[2];
    }

    private static double calc_d2fydudy2(double t2)
    {
        return dN43(t2)[2];
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
}
