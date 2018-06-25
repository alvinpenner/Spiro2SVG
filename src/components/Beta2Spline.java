
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1.
// fit a 5-point Beta2-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// constrain only the slopes at the endpoints and keep P2 (Q3/Q4) arbitrary.
// the Beta2 spline has an adjustable arm length, d, at the Bezier splice at Q3/Q4.
// psi, ψ, is the angle at the Bezier splice (Q3-Q2/Q5-Q4) = angle of (Q6-Q1).
// linearize the equations wrt (d1, d2, x2, y2, d) and solve a 5x5 system of equations.
// decompose the Beta2-Spline into 2 Beziers, range (0,1) and (1,2).
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2).
// t2 must be chosen to minimize the distance to g(t1).
// see Spiro2SVG Book4, Dec 2017, page 13.

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\Beta2Spline.java

public class Beta2Spline
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI/4; // Math.PI/4;
    public static final int N = 100;
    public static double[][] Bezx;              // 2 Beziers, 4 points each, x component
    public static double[][] Bezy;              // 2 Beziers, 4 points each, y component
    //private static CircleFxn fitted;
    //private static CycloidFxn fitted;
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[5][N+1];        // partial wrt (d1, d2, x2, y2, d)
    private static double psi;
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 10;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //fitted = new epiTrochoidFxn(2.0);
        fitted = new epiTrochoidFxn(10.0);
        //System.out.println("Beta2-Spline convert_at_P2 = " + convert_at_P2(19.983314292966483, 26.42763336958588, 175.47633731103565, 59.05668195284478, true) + "\n");
        //System.out.println("Beta2-Spline solve_at_P2 = " + solve_at_P2(22.742672063451316, 21.951281542489866, 166.71089691010013, 67.4862883523671, 24.830428922781234, true) + "\n");
        //System.out.println("Beta2-Spline iterate_at_P2 = "
        //                  + iterate_at_P2(18.767825828964003, 18.76782582892743, 166.29676657566654, 68.88237609446028, 28.074383011128273) + "\n");
        //System.out.println("Beta2-Spline iterate_at_P2 = "
        //                  + iterate_at_P2(33.0, 11.5, 159.0, 82.8, 20.6) + "\n");
        System.out.println("Beta2-Spline solve_at_P2 = "
                          + solve_at_P2(25.0, 30.0, 166.0, 68.0, 25, true) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
    }

    private static double iterate_at_P2(double d1, double d2, double x2, double y2, double d)
    {
        // calculate a new estimate of (d1, d2, x2, y2, d) by setting dF = 0
        // include second-order responses: d2f/ddi/ddj
        // see Spiro2SVG Book 3, page 54 (applied to 5-point cubic Beta2-Spline)
        // setup 5-variable Newton-Raphson iteration

        final double gain = 1;                          // fudge factor to reduce gain
        final int MAXLOOP = 4000;
        double[] f_gx = new double[N+1];
        double[] f_gy = new double[N+1];
        double[] dfxdu = new double[N+1];
        double[] d2fxdudu = new double[N+1];
        double[] dfydu = new double[N+1];
        double[] d2fydudu = new double[N+1];
        double[][] dfxdd = new double[5][N+1];
        double[][] dfydd = new double[5][N+1];
        double[][] d2fxdudd = new double[5][N+1];
        double[][] d2fydudd = new double[5][N+1];
        double[][][] d2fxdddd = new double[5][5][N+1];
        double[][][] d2fydddd = new double[5][5][N+1];

        double[][] Jac = new double[5][5];
        double[] dFdd = new double[5];
        double[] trap_in = new double[N+1];
        double[] deld;                                  // (-Δd1, -Δd2, -Δx2, -Δy2, -Δd)
        int loop = 0;
        int i, j, k, seg;                               // Bezier segment, before or after the splice
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(d1, d2, x2, y2, d, false)))   // initiallize at (x2, y2)
            {
                System.out.println("fail at " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d);
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
                d2fxdudu[i] = t2_vs_t1.d2fn(Bezx[seg], t2[i] - seg);
                dfydu[i] = t2_vs_t1.dfn(Bezy[seg], t2[i] - seg);
                d2fydudu[i] = t2_vs_t1.d2fn(Bezy[seg], t2[i] - seg);
                dfxdd[0][i] = calc_dfxdd1(seg, t2[i] - seg);
                dfydd[0][i] = calc_dfydd1(seg, t2[i] - seg);
                dfxdd[1][i] = calc_dfxdd2(seg, t2[i] - seg);
                dfydd[1][i] = calc_dfydd2(seg, t2[i] - seg);
                dfxdd[2][i] = calc_dfxdx2(seg, t2[i] - seg);
                dfydd[2][i] = 0;
                dfxdd[3][i] = 0;
                dfydd[3][i] = calc_dfydy2(seg, t2[i] - seg);
                dfxdd[4][i] = calc_dfxdd(seg, t2[i] - seg);
                dfydd[4][i] = calc_dfydd(seg, t2[i] - seg);

                d2fxdudd[0][i] = calc_d2fxdudd1(seg, t2[i] - seg);
                d2fydudd[0][i] = calc_d2fydudd1(seg, t2[i] - seg);
                d2fxdudd[1][i] = calc_d2fxdudd2(seg, t2[i] - seg);
                d2fydudd[1][i] = calc_d2fydudd2(seg, t2[i] - seg);
                d2fxdudd[2][i] = calc_d2fxdudx2(seg, t2[i] - seg);
                d2fydudd[2][i] = 0;
                d2fxdudd[3][i] = 0;
                d2fydudd[3][i] = calc_d2fydudy2(seg, t2[i] - seg);
                d2fxdudd[4][i] = calc_d2fxdudd(seg, t2[i] - seg);
                d2fydudd[4][i] = calc_d2fydudd(seg, t2[i] - seg);

                d2fxdddd[0][0][i] = calc_d2fxdd1dd1(seg, t2[i] - seg);  // non-linear effects
                d2fxdddd[0][1][i] = calc_d2fxdd1dd2(seg, t2[i] - seg);
                d2fxdddd[1][1][i] = calc_d2fxdd2dd2(seg, t2[i] - seg);
                d2fxdddd[0][4][i] = calc_d2fxdd1dd(seg, t2[i] - seg);
                d2fxdddd[1][4][i] = calc_d2fxdd2dd(seg, t2[i] - seg);
                d2fxdddd[1][0][i] = d2fxdddd[0][1][i];
                d2fxdddd[4][0][i] = d2fxdddd[0][4][i];
                d2fxdddd[4][1][i] = d2fxdddd[1][4][i];
                d2fydddd[0][0][i] = calc_d2fydd1dd1(seg, t2[i] - seg);
                d2fydddd[0][1][i] = calc_d2fydd1dd2(seg, t2[i] - seg);
                d2fydddd[1][1][i] = calc_d2fydd2dd2(seg, t2[i] - seg);
                d2fydddd[0][4][i] = calc_d2fydd1dd(seg, t2[i] - seg);
                d2fydddd[1][4][i] = calc_d2fydd2dd(seg, t2[i] - seg);
                d2fydddd[1][0][i] = d2fydddd[0][1][i];
                d2fydddd[4][0][i] = d2fydddd[0][4][i];
                d2fydddd[4][1][i] = d2fydddd[1][4][i];
                //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i] + ", " + t2dx2[i] + ", " + t2dy2[i] + ", " + t2dd[i] + ", " + f_gx[i] + ", " + f_gy[i] + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + dfxdd1[i] + ", " + dfydd1[i] + ", " + dfxdd2[i] + ", " + dfydd2[i] + ", " + dfxdx2[i] + ", " + dfydx2[i] + ", " + dfxdy2[i] + ", " + dfydy2[i] + ", " + dfxdd[i] + ", " + dfydd[i] + ", " + d2fxdudd1[i] + ", " + d2fydudd1[i] + ", " + d2fxdudd2[i] + ", " + d2fydudd2[i] + ", " + d2fxdudx2[i] + ", " + d2fydudx2[i] + ", " + d2fxdudy2[i] + ", " + d2fydudy2[i] + ", " + d2fxdudd[i] + ", " + d2fydudd[i]);
            }

            // calc dFdd[j] at current (d1, d2, x2, y2, d)

            for (i = 0; i < 5; i++)
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]); // original code
                    //trap_in[k] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];         // new code
                dFdd[i] = t2_vs_t1.integrate(trap_in);
            }

            // calc d2Fdd[i]dd[j] (Jacobean matrix)

            for (i = 0; i < 5; i++)
                for (j = 0; j < 5; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N; k++)
                    {
                        //    System.out.println(k + ", " + ", " + t2[k] + ", " + trap_in[k] + ", " + d2udddd[i][j][k] + ", " + f_gx[k] + ", " + dfxdu[k] + ", " + f_gy[k] + ", " + dfydu[k]);
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k]                // new code
                                   + f_gx[k]*d2fxdddd[i][j][k]
                                   + dfydd[i][k]*dfydd[j][k]
                                   + f_gy[k]*d2fydddd[i][j][k]
                                   - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k];
                        //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }

            deld = BSpline5.multmv(BSpline5.invertm(Jac), dFdd);    // this is actually the negative of Δd
            d1 -= deld[0]/gain;                                     // fix fix blatant fudge factor in gain
            d2 -= deld[1]/gain;
            x2 -= deld[2]/gain;
            y2 -= deld[3]/gain;
            d  -= deld[4]/gain;
            //System.out.println("Jac");
            //for (i = 0; i < Jac.length; i++)
            //{
            //    for (int j = 0; j < Jac.length; j++)
            //        System.out.print(Jac[i][j] + ", ");
            //    System.out.println();
            //}

            Jacdet = BSpline5.detm(Jac);
            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + dFdd[2] + ", " + dFdd[3] + ", " + dFdd[4] + ", " + Jacdet);
            System.out.println("deld = " + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ", " + deld[4]);
            BSpline5.dump_Jac(Jac);

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            for (i = 0; i <= N; i++)
                for (j = 0; j < 5; j++)
                    t2[i] -= t2dd[j][i]*deld[j];                        // first-order response
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL) && (Math.abs(deld[2]) < TOL) && (Math.abs(deld[3]) < TOL) && (Math.abs(deld[4]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 x2 y2 d = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d);
            double x1 = fitted.getx(t1_start) + d1*Math.cos(theta_start);   // calculate beta2
            double y1 = fitted.gety(t1_start) + d1*Math.sin(theta_start);
            double x3 = fitted.getx(t1_end) - d2*Math.cos(theta_end);
            double y3 = fitted.gety(t1_end) - d2*Math.sin(theta_end);
            double d0 = 0.25*Math.sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1));
            System.out.println("c, beta2 =, " + fitted.getc() + ", " + 8*(d0 - d)/d);
            return solve_at_P2(d1, d2, x2, y2, d, true);                    // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ", " + deld[4] + ")");
        return Double.NaN;
    }

    private static double convert_at_P2(double d1, double d2, double x2, double y2, boolean print)
    {
        // convert a 5-point B-Spline to a pair of spliced cubic Beziers
        System.out.println("convert_at_P2: B-Spline theta c t d1 d2 = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        double x1 = fitted.getx(t1_start) + d1*Math.cos(theta_start);
        double y1 = fitted.gety(t1_start) + d1*Math.sin(theta_start);
        double x3 = fitted.getx(t1_end) - d2*Math.cos(theta_end);
        double y3 = fitted.gety(t1_end) - d2*Math.sin(theta_end);

        return iterate_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, 0.25*Math.sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1)));
        //return solve_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, 0.25*Math.sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1)), print);
    }

    private static double solve_at_P2(double d1, double d2, double x2, double y2, double d, boolean print)
    {
        // 5-point cubic Beta2-Spline curve, decomposed as two Beziers
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        //d2 += 0.005;
        //y2 += 0.001;                          // temporary code for calculating response functions
        //d += 0.005;

        Bezx = new double[][] {{fitted.getx(t1_start),
                                fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                x2,
                                x2},
                               {x2,
                                x2,
                                fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                fitted.getx(t1_end)}};
        Bezy = new double[][] {{fitted.gety(t1_start),
                                fitted.gety(t1_start) + d1*Math.sin(theta_start),
                                y2,
                                y2},
                               {y2,
                                y2,
                                fitted.gety(t1_end) - d2*Math.sin(theta_end),
                                fitted.gety(t1_end)}};


        psi = Math.atan2(Bezy[1][2] - Bezy[0][1], Bezx[1][2] - Bezx[0][1]);
        //System.out.println("del d = " + (Bezx[1][2] - Bezx[0][1]) + ", " + (Bezy[1][2] - Bezy[0][1]));
        //System.out.println("psi = " + psi);
        Bezx[0][2] -= d*Math.cos(psi);
        Bezx[1][1] += d*Math.cos(psi);
        Bezy[0][2] -= d*Math.sin(psi);
        Bezy[1][1] += d*Math.sin(psi);

        if (t2[N] == 0)
            System.out.println("__start Beta2-Spline theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d);
        else
            System.out.println("__solve new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d + ", " + calc_error());
        //fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2, t2dd1, t2dd2, t2dx2, t2dy2, t2dd");
        int seg = 0;                // Bezier segment, before or after the splice
        for (int i = 0; i <= N; i++)
        {
            if (i == 0)
                t2[i] = solve_quintic_for_t2(i, 0, seg);
            else
                t2[i] = solve_quintic_for_t2(i, t2[i-1], seg);
            if (seg == 0 && t2[i] > 1)
            {
                seg++;
                t2[i] = solve_quintic_for_t2(i, 1, seg);   // re-calculate (initiallize at the splice)
            }
            if (seg == 1 && t2[i] < 1)
            {
                System.out.println("seg reverted back to zero at i = " + i + ", " + t2[i] + ", " + d);
                return Double.NaN;
            }
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //scan_quintic_near_t2(i, seg, t2[i]);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dxy(i, seg, t2[i] - seg, "d1");
            t2dd[1][i] = calc_t2dxy(i, seg, t2[i] - seg, "d2");
            t2dd[2][i] = calc_t2dxy(i, seg, t2[i] - seg, "x2");
            t2dd[3][i] = calc_t2dxy(i, seg, t2[i] - seg, "y2");
            t2dd[4][i] = calc_t2dxy(i, seg, t2[i] - seg, "d");
            if (print)
            {
                double t1 = t1_start + i*(t1_end - t1_start)/N;
                double f = (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))*t2_vs_t1.dfn(Bezx[seg], t2[i] - seg) + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1))*t2_vs_t1.dfn(Bezy[seg], t2[i] - seg);
                System.out.println(seg + ", " + t1 + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + f);
            }
        }
        double retVal = calc_error();
        System.out.println("__new t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d + ", " + (float) retVal + ", " + (float) Jacdet);
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

    private static double solve_quintic_for_t2(double d, double t, int seg)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the cubic Bezier t-value

        double t1 = t1_start + d*(t1_end - t1_start)/N;
        double f, fprime, f2prime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        int loop = 0;
        int success = 0;

        // initial estimate using quadratic approximation

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
                    return Double.NaN;
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
                System.out.println("\ntoo many loops  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + del_t);
                return Double.NaN;
            }
            if (f == 0 && fprime == 0)
                del_t = 0;
            else
                del_t = -f/fprime;
            t += del_t;
            loop++;
            if (Math.abs(del_t) < TOL) success++;
            //System.out.println("         t2 =, " + t + ", " + f + ", " + fprime);
        } while (success < 2);
//        } while (Math.abs(del_t) > TOL);
        return t + seg;                     // compensate for Bezier segment offset
    }

    private static double calc_t2dxy(int i, int seg, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        double denom = t2_vs_t1.dfn(Bezx[seg], t2)*t2_vs_t1.dfn(Bezx[seg], t2) + (t2_vs_t1.fn(Bezx[seg], t2) - X)*t2_vs_t1.d2fn(Bezx[seg], t2) + t2_vs_t1.dfn(Bezy[seg], t2)*t2_vs_t1.dfn(Bezy[seg], t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*t2_vs_t1.d2fn(Bezy[seg], t2);
        double numer = Double.NaN;
        if (type.equals("x2"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdx2(seg, t2) + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudx2(seg, t2);
        else if (type.equals("y2"))
            numer = t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydy2(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudy2(seg, t2);
        else if (type.equals("d1"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdd1(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydd1(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudd1(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudd1(seg, t2);
        else if (type.equals("d2"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdd2(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydd2(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudd2(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudd2(seg, t2);
        else if (type.equals("d"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdd(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydd(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudd(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudd(seg, t2);
        return -numer/denom;
    }

    private static double calc_dfxdx2(int seg, double t2)
    {
        return N33(t2)[2 - 2*seg] + N33(t2)[3 - 2*seg];
    }

    private static double calc_d2fxdudx2(int seg, double t2)
    {
        return dN33(t2)[2 - 2*seg] + dN33(t2)[3 - 2*seg];
    }

    private static double calc_dfydy2(int seg, double t2)
    {
        return N33(t2)[2 - 2*seg] + N33(t2)[3 - 2*seg];
    }

    private static double calc_d2fydudy2(int seg, double t2)
    {
        return dN33(t2)[2 - 2*seg] + dN33(t2)[3 - 2*seg];
    }

    private static double calc_dfxdd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.cos(theta_start)*N33(t2)[1];
        return retVal + sgn*calc_dpsidd("d1")*d*Math.sin(psi)*N33(t2)[2 - seg];
    }

    private static double calc_dfydd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.sin(theta_start)*N33(t2)[1];
        return retVal - sgn*calc_dpsidd("d1")*d*Math.cos(psi)*N33(t2)[2 - seg];
    }

    private static double calc_d2fxdudd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.cos(theta_start)*dN33(t2)[1];
        return retVal + sgn*calc_dpsidd("d1")*d*Math.sin(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fydudd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.sin(theta_start)*dN33(t2)[1];
        return retVal - sgn*calc_dpsidd("d1")*d*Math.cos(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_dfxdd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.cos(theta_end)*N33(t2)[2];
        return retVal + sgn*calc_dpsidd("d2")*d*Math.sin(psi)*N33(t2)[2 - seg];
    }

    private static double calc_dfydd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.sin(theta_end)*N33(t2)[2];
        return retVal - sgn*calc_dpsidd("d2")*d*Math.cos(psi)*N33(t2)[2 - seg];
    }

    private static double calc_d2fxdudd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.cos(theta_end)*dN33(t2)[2];
        return retVal + sgn*calc_dpsidd("d2")*d*Math.sin(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fydudd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.sin(theta_end)*dN33(t2)[2];
        return retVal - sgn*calc_dpsidd("d2")*d*Math.cos(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_dfxdd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.cos(psi)*N33(t2)[2 - seg];
    }

    private static double calc_dfydd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.sin(psi)*N33(t2)[2 - seg];
    }

    private static double calc_d2fxdudd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.cos(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fydudd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.sin(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fxdd1dd1(int seg, double t2)
    {
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*d*(Math.cos(psi)*calc_dpsidd("d1")*calc_dpsidd("d1") + Math.sin(psi)*calc_dpsidd("d1d1"));
    }

    private static double calc_d2fxdd1dd2(int seg, double t2)
    {
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*d*(Math.cos(psi)*calc_dpsidd("d1")*calc_dpsidd("d2") + Math.sin(psi)*calc_dpsidd("d1d2"));
    }

    private static double calc_d2fxdd2dd2(int seg, double t2)
    {
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*d*(Math.cos(psi)*calc_dpsidd("d2")*calc_dpsidd("d2") + Math.sin(psi)*calc_dpsidd("d2d2"));
    }

    private static double calc_d2fxdd1dd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*Math.sin(psi)*calc_dpsidd("d1");
    }

    private static double calc_d2fxdd2dd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*Math.sin(psi)*calc_dpsidd("d2");
    }

    private static double calc_d2fydd1dd1(int seg, double t2)
    {
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*d*(Math.sin(psi)*calc_dpsidd("d1")*calc_dpsidd("d1") - Math.cos(psi)*calc_dpsidd("d1d1"));
    }

    private static double calc_d2fydd1dd2(int seg, double t2)
    {
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*d*(Math.sin(psi)*calc_dpsidd("d1")*calc_dpsidd("d2") - Math.cos(psi)*calc_dpsidd("d1d2"));
    }

    private static double calc_d2fydd2dd2(int seg, double t2)
    {
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return sgn*N33(t2)[2 - seg]*d*(Math.sin(psi)*calc_dpsidd("d2")*calc_dpsidd("d2") - Math.cos(psi)*calc_dpsidd("d2d2"));
    }

    private static double calc_d2fydd1dd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*N33(t2)[2 - seg]*Math.cos(psi)*calc_dpsidd("d1");
    }

    private static double calc_d2fydd2dd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*N33(t2)[2 - seg]*Math.cos(psi)*calc_dpsidd("d2");
    }

    private static double calc_dpsidd(String type)
    {
        double num = Bezy[1][2] - Bezy[0][1];
        double den = Bezx[1][2] - Bezx[0][1];
        double sqr = num*num + den*den;

        if (type.equals("d1"))
            return (-den*Math.sin(theta_start) + num*Math.cos(theta_start))/sqr;
        else if (type.equals("d2"))
            return (-den*Math.sin(theta_end) + num*Math.cos(theta_end))/sqr;
        else if (type.equals("d1d1"))
            return ((num*num - den*den)*Math.sin(2*theta_start) + 2*num*den*Math.cos(2*theta_start))/sqr/sqr;
        else if (type.equals("d2d2"))
            return ((num*num - den*den)*Math.sin(2*theta_end) + 2*num*den*Math.cos(2*theta_end))/sqr/sqr;
        else if (type.equals("d1d2"))
            return ((num*num - den*den)*Math.sin(theta_start + theta_end) + 2*num*den*Math.cos(theta_start + theta_end))/sqr/sqr;
        return Double.NaN;
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
