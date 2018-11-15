
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 4-point cubic Bezier (P0 - P3) to it, using parameter 0 < t2 < 1.
// constrain only the slopes at the endpoints
// constrain the area to be correct
// optimize wrt to d1 subject to this constraint
// Bezier[4] = f(x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// this is based on BezierCubic.java which is from t2_vs_t1.iterate_at_d1_d2

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BezierCubicOneDim.java

public class BezierCubicOneDim
{
    private static final double a_b = 180;          // scale factor to make rms error dimensionless
    private static final double t1_start = 0;
    private static final double t1_end = Math.PI/4;
    private static final int N = 100;
    private static double[] Bezx;                   // cubic Bezier, 4 points, x component
    private static double[] Bezy;                   // cubic Bezier, 4 points, y component
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[2][N+1];    // partial wrt (d1, d2)
    private static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        //fitted = new epiTrochoidFxn(3.61 + 0*0.00001);
        //iterate_at_P2(57.31291807448238 + 0*0.01, 0);
        //fitted = new epiTrochoidFxn(3.6145);
        //iterate_at_P2(57.3, 0);
        fitted = new epiTrochoidFxn(20);
        iterate_at_P2(27, 0);
        //System.out.println("oneDim cubic Bezier solve_at_P2 = " + solve_at_P2(17.889026722854, 92.22096516520845, true) + "\n");
        //System.out.println("oneDim cubic Bezier solve_at_P2 = " + solve_at_P2(56.811, 32.006, true) + "\n");
        //scan_at_P2();
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //for (int i = 0; i <= 10; i++)
        //{
        //    double tempd1 = 60 + 0.01*i;
        //    System.out.println(i + ", " + spiro_area() + ", " + tempd1 + ", " + calc_d2(tempd1) + ", " + calc_dd2dd1(tempd1) + ", " + calc_d2d2dd1dd1(tempd1));
        //}
        //double tempd1 = 18;
        //System.out.println("test , " + fitted.getc() + ", " + tempd1 + ", " + calc_d2(tempd1) + ", " + calc_dd2dd1(tempd1) + ", " + calc_dd2dc(tempd1) + ", " + calc_d2d2dd1dc(tempd1));
        //System.out.println("test , " + fitted.getc() + ", " + tempd1 + ", " + calc_dd2dc_at_d1(tempd1) + ", " + calc_dd1dc_at_d2(tempd1) + ", " + calc_d2d1dcdd2_at_d2(tempd1) + ", " + calc_d2d2dcdd1_at_d1(tempd1));
        //System.out.println("test , " + fitted.getc() + ", " + tempd1 + ", " + calc_d2(tempd1) + ", " + calc_dd2dd1(tempd1) + ", " + calc_d2d2dd1dd1(tempd1));
        //calc_d1_d2_at_h();
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

        double[] trap_in = new double[N+1];
        double[][] Jac = new double[1][1];
        double[] dFdd = new double[1];
        double[] d2Fdddc = new double[1];                           // augmented matrix
        double dFdc;
        double d2Fdh2, d2Fdhdc;                                     // response wrt anti-symmetric variable h
        double d2Fdd2dc;                                            // for testing only, alternative to d2Fdddc[0]
        //double d2Fdcdc;                                           // augmented matrix
        double[][] Augment = new double[2][2];
        double deld;                                                // (-Δd1)
        int i, j, k, loop = 0;
        double t1;

        d2 = calc_d2(d1);                                           // initiallize
        //fitted = new epiTrochoidFxn(fitted.getc() + 0.00001);        // just for numerical d2Fdddc
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
                dfydu[i] = t2_vs_t1.dfn(Bezy, t2[i]);
                d2fxdudu[i] = t2_vs_t1.d2fn(Bezx, t2[i]);
                d2fydudu[i] = t2_vs_t1.d2fn(Bezy, t2[i]);
                df_gxdc[i] = calc_df_gxdc(t1, t2[i]);
                df_gydc[i] = calc_df_gydc(t1, t2[i]);

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

            double dd1dd1, dd1dd2, dd2dd2, dd1dc, dd2dc, dd1, dd2;
            for (i = 0; i < 1; i++)                 // one loop only
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] =  f_gx[k]*(dfxdd[0][k] + dfxdu[k]*t2dd[0][k]) + f_gy[k]*(dfydd[0][k] + dfydu[k]*t2dd[0][k]) // original code
                               + (f_gx[k]*(dfxdd[1][k] + dfxdu[k]*t2dd[1][k]) + f_gy[k]*(dfydd[1][k] + dfydu[k]*t2dd[1][k]))*calc_dd2dd1(d1);
                    //trap_in[k] = f_gx[k]*dfxdd[0][k] + f_gy[k]*dfydd[0][k]
                    //           + (f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k])*calc_dd2dd1(d1);   // new code
                dFdd[i] = t2_vs_t1.integrate(trap_in);
                for (k = 0; k <= N; k++)
                {
                    dd1dd2 = dfxdd[0][k]*dfxdd[1][k] + dfydd[0][k]*dfydd[1][k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[1][k];
                    dd2dd2 = dfxdd[1][k]*dfxdd[1][k] + dfydd[1][k]*dfydd[1][k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*t2dd[1][k];
                    trap_in[k] = dfxdd[0][k]*df_gxdc[k] + dfydd[0][k]*df_gydc[k]
                               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*calc_t2dxy(k, t2[k], "c")
                               + dd1dd2*calc_dd2dc(d1);
                    trap_in[k] += ((dfxdd[1][k]*df_gxdc[k] + dfydd[1][k]*df_gydc[k]
                               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*calc_t2dxy(k, t2[k], "c"))
                               + dd2dd2*calc_dd2dc(d1))*calc_dd2dd1(d1);
                    trap_in[k] += (f_gx[k]*(dfxdd[1][k] + dfxdu[k]*t2dd[1][k]) + f_gy[k]*(dfydd[1][k] + dfydu[k]*t2dd[1][k]))*calc_d2d2dd1dc(d1);
                    //trap_in[k] += (f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k])*calc_d2d2dd1dc(d1);
                    // System.out.println(i + ", " + k + ", " + trap_in[k]);
                }
                d2Fdddc[i] = t2_vs_t1.integrate(trap_in);               // augmented matrix
            }

            // test the hypothesis : dd1dd2 + calc_dd2dd1(d1)*dd2dd2 = 0

            for (k = 0; k <= N; k++)
                trap_in[k] = dfxdd[0][k]*dfxdd[0][k] + dfydd[0][k]*dfydd[0][k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[0][k];
            dd1dd1 = t2_vs_t1.integrate(trap_in);
            for (k = 0; k <= N; k++)
                trap_in[k] = dfxdd[0][k]*dfxdd[1][k] + dfydd[0][k]*dfydd[1][k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[1][k];
            dd1dd2 = t2_vs_t1.integrate(trap_in);
            for (k = 0; k <= N; k++)
                trap_in[k] = dfxdd[1][k]*dfxdd[1][k] + dfydd[1][k]*dfydd[1][k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*t2dd[1][k];
            dd2dd2 = t2_vs_t1.integrate(trap_in);
            System.out.println("test dd1dd2 + calc_dd2dd1(d1)*dd2dd2 = 0 : " + dd1dd2 + ", " + calc_dd2dd1(d1) + ", " + dd2dd2 + ", " + calc_dd2dc_at_d1(d1)*(dd1dd2 + calc_dd2dd1(d1)*dd2dd2));
            System.out.println("test dd1dd1 + calc_dd2dd1(d1)*dd1dd2 = 0 : " + dd1dd1 + ", " + calc_dd2dd1(d1) + ", " + dd1dd2 + ", " + calc_dd1dc_at_d2(d1)*(dd1dd1 + calc_dd2dd1(d1)*dd1dd2));

            // calc d2Fdd[i]dd[j] (Jacobean matrix)

            for (i = 0; i < 1; i++)                 // one loop only
                for (j = 0; j < 1; j++)             // one loop only
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N; k++)
                    {
                        dd1dd1 = dfxdd[0][k]*dfxdd[0][k] + dfydd[0][k]*dfydd[0][k]
                               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[0][k];
                        dd1dd2 = dfxdd[0][k]*dfxdd[1][k] + dfydd[0][k]*dfydd[1][k]
                               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[1][k];
                        dd2dd2 = dfxdd[1][k]*dfxdd[1][k] + dfydd[1][k]*dfydd[1][k]
                               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*t2dd[1][k];
                        dd2 = f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k];
                        trap_in[k] = dd1dd1 + 2*calc_dd2dd1(d1)*dd1dd2 + calc_dd2dd1(d1)*calc_dd2dd1(d1)*dd2dd2
                                   + calc_d2d2dd1dd1(d1)*dd2;
                    }
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }
            //System.out.println("test Jac(d2) : " + Jac[0][0] + ", " + Jac[0][0]/calc_dd2dd1(d1)/calc_dd2dd1(d1));

            // d2Fdd2dc : for testing only, alternative to d2Fdddc[0]

            for (k = 0; k <= N; k++)
            {
                dd1dd1 = dfxdd[0][k]*dfxdd[0][k] + dfydd[0][k]*dfydd[0][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[0][k];
                dd1dd2 = dfxdd[0][k]*dfxdd[1][k] + dfydd[0][k]*dfydd[1][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[1][k];
                trap_in[k] = dfxdd[1][k]*df_gxdc[k] + dfydd[1][k]*df_gydc[k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*calc_t2dxy(k, t2[k], "c")
                           + dd1dd2*calc_dd1dc_at_d2(d1);
                trap_in[k] += ((dfxdd[0][k]*df_gxdc[k] + dfydd[0][k]*df_gydc[k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*calc_t2dxy(k, t2[k], "c"))
                           + dd1dd1*calc_dd1dc_at_d2(d1))/calc_dd2dd1(d1);
                trap_in[k] += (f_gx[k]*(dfxdd[0][k] + dfxdu[k]*t2dd[0][k]) + f_gy[k]*(dfydd[0][k] + dfydu[k]*t2dd[0][k]))*calc_d2d1dcdd2_at_d2(d1);
                //trap_in[k] += (f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k])*calc_d2d2dd1dc(d1);
                // System.out.println(i + ", " + k + ", " + trap_in[k]);
            }
            d2Fdd2dc = t2_vs_t1.integrate(trap_in);               // augmented matrix

            // d2Fdh2: the second-order response wrt the anti-symmetric variable h = ln(d1/d2)

            for (k = 0; k <= N; k++)
            {
                dd1 = f_gx[k]*dfxdd[0][k] + f_gy[k]*dfydd[0][k];
                dd2 = f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k];
                dd1dd1 = dfxdd[0][k]*dfxdd[0][k] + dfydd[0][k]*dfydd[0][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[0][k];
                dd1dd2 = dfxdd[0][k]*dfxdd[1][k] + dfydd[0][k]*dfydd[1][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[1][k];
                dd2dd2 = dfxdd[1][k]*dfxdd[1][k] + dfydd[1][k]*dfydd[1][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*t2dd[1][k];
                trap_in[k] = calc_dd1dh(d1)*calc_dd1dh(d1)*dd1dd1 + 2*calc_dd1dh(d1)*calc_dd2dh(d1)*dd1dd2 + calc_dd2dh(d1)*calc_dd2dh(d1)*dd2dd2
                           + calc_d2d1dh2(d1)*dd1 + calc_d2d2dh2(d1)*dd2;
            }
            d2Fdh2 = t2_vs_t1.integrate(trap_in);

            // d2Fdhdc: the second-order response wrt the anti-symmetric variable h = ln(d1/d2)

            for (k = 0; k <= N; k++)
            {
                dd1 = f_gx[k]*dfxdd[0][k] + f_gy[k]*dfydd[0][k];
                dd2 = f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k];
                dd1dc = dfxdd[0][k]*df_gxdc[k] + dfydd[0][k]*df_gydc[k]
                      - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*calc_t2dxy(k, t2[k], "c");
                dd2dc = dfxdd[1][k]*df_gxdc[k] + dfydd[1][k]*df_gydc[k]
                      - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*calc_t2dxy(k, t2[k], "c");
                dd1dd1 = dfxdd[0][k]*dfxdd[0][k] + dfydd[0][k]*dfydd[0][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[0][k];
                dd1dd2 = dfxdd[0][k]*dfxdd[1][k] + dfydd[0][k]*dfydd[1][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[0][k]*t2dd[1][k];
                dd2dd2 = dfxdd[1][k]*dfxdd[1][k] + dfydd[1][k]*dfydd[1][k]
                       - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[1][k]*t2dd[1][k];
                trap_in[k] = calc_dd1dh(d1)*calc_dd1dc_at_h(d1)*dd1dd1 + (calc_dd1dh(d1)*calc_dd2dc_at_h(d1) + calc_dd2dh(d1)*calc_dd1dc_at_h(d1))*dd1dd2 + calc_dd2dh(d1)*calc_dd2dc_at_h(d1)*dd2dd2
                           + calc_d2d1dcdh(d1)*dd1 + calc_d2d2dcdh(d1)*dd2
                           + calc_dd1dh(d1)*dd1dc + calc_dd2dh(d1)*dd2dc;
            }
            d2Fdhdc = t2_vs_t1.integrate(trap_in);

            // calculate determinant of augmented matrix

            for (k = 0; k <= N; k++)
                trap_in[k] =  f_gx[k]*df_gxdc[k] + f_gy[k]*df_gydc[k]
                           + (f_gx[k]*dfxdd[1][k] + f_gy[k]*dfydd[1][k])*calc_dd2dc(d1);
            dFdc = t2_vs_t1.integrate(trap_in);
            //for (k = 0; k <= N; k++)
            //    trap_in[k] = df_gxdc[k]*df_gxdc[k] + df_gydc[k]*df_gydc[k]
            //               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*calc_t2dxy(k, t2[k], "c")*calc_t2dxy(k, t2[k], "c");
            //d2Fdcdc = t2_vs_t1.integrate(trap_in);        // not needed
            for (i = 0; i < 1; i++)                         // one loop only
                for (j = 0; j < 1; j++)
                    Augment[i][j] = Jac[i][j];
            for (i = 0; i < 1; i++)
            {
                //Augment[i][1] = d2Fdddc[i];
                //Augment[1][i] = Augment[i][1];
            }
            //Augment[1][1] = d2Fdcdc;
            //System.out.println("dFdc = " + fitted.getc() + ", " + d1 + ", " + d2 + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);

            deld = dFdd[0]/Jac[0][0];                               // this is actually the negative of Δd1
            d1 -= deld/gain;
            d2 = calc_d2(d1);                                       // try to prevent error propagation

            System.out.println("dFdd = , " + dFdd[0] + ", " + Jac[0][0] + ", " + BSpline5.detm(Augment));
            System.out.println("deld = , " + deld + ", " + dFdc + ", " + d2Fdddc[0]);
            //System.out.println("\ndFdc = , " + fitted.getc() + ", " + d1 + ", " + d2 + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);
            //BSpline5.dump_Jac(Jac);
            //BSpline5.dump_Jac(Augment);
            //double[][] invJac = BSpline5.invertm(Jac);
            //System.out.println("invJac = " + invJac[0][0] + ", " + invJac[0][1] + ", " + invJac[1][0] + ", " + invJac[1][1]);
            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            //for (i = 0; i <= N; i++)            // disabled due to crashes
            //{
            //    t2[i] -= t2dd[0][i]*deld[0] + t2dd[1][i]*deld[1];   // first-order response
                // System.out.println("incrementing t2 array : " + i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            //}
        } while ((loop < MAXLOOP) && !(Math.abs(deld) < TOL));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 = , , , , , , " + d1 + ", " + d2 + ", " + fitted.getc() + ", " + Jac[0][0]);
            double rms = solve_at_P2(d1, d2, true);                          // final run just for good measure
            // calculate dddc[i]
            //double[] dt2dc = BSpline5.gaussj(Jac, d2Fdddc);
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (Math.cos(eigangle)*d2Fdddc[0] - Math.sin(eigangle)*d2Fdddc[1]) + ", " + (float) (Math.sin(eigangle)*d2Fdddc[0] + Math.cos(eigangle)*d2Fdddc[1]));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (eigangle*180.0/Math.PI) + ", " + (float) (Math.atan(d2Fdddc[0]/d2Fdddc[1])*180.0/Math.PI));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) -d2Fdddc[0] + ", " + (float) -d2Fdddc[1]);
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (Jac[1][0]/Jac[0][0]) + ", " + (float) (Jac[1][1]/Jac[0][1]) + ", " + (float) (d2Fdddc[1]/d2Fdddc[0]));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) BSpline5.detm(Augment));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + eig0 + ", " + rms*rms*180*180);
            //System.out.println("\nfinal oneDim CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + rms + ", " + Jac[0][0] + ", " + d2Fdddc[0] + ", " + d2Fdd2dc);
            System.out.println("d2Fdh2 / d2Fdhdc = , , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + rms + ", " + d2Fdh2 + ", " + d2Fdhdc);
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld + ")");
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
            System.out.println("__start oneDim cubic Bezier at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + calc_error());
        if (d1 < 0 || d2 < 0)
            System.out.println("WARNING: negative arm length = " + d1 + ", " + d2);
        //fitted.gen_Bezier(new double[] {Bezx[0], Bezy[0], Bezx[1], Bezy[1], Bezx[2], Bezy[2], Bezx[3], Bezy[3]});
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);

        if (print) System.out.println("\n                , t1, t2, t2dd1, t2dd2, t2dc");
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
                System.out.println("oneDim cubic Bez, " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + calc_t2dxy(i, t2[i], "c"));
            }
        }
        double retVal = calc_error();
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + retVal);
        //System.out.println("calc_dd2dd1(d1) = " + calc_dd2dd1(d1) + ", " + calc_d2(d1));
        //System.out.println("calc_dd1dh      = " + calc_dd1dh(d1) + ", " + calc_dd2dh(d1) + ", " + calc_d2d1dh2(d1) + ", " + calc_d2d2dh2(d1));
        //System.out.println("calc_dd1dc_at_h = " + calc_dd1dc_at_h(d1) + ", " + calc_dd2dc_at_h(d1) + ", " + calc_d2d1dcdh(d1) + ", " + calc_d2d2dcdh(d1));
        //System.out.println("1 consistency (d/dh)    : " + (calc_dd1dh(d1)/d1 - calc_dd2dh(d1)/d2 - 1));
        //System.out.println("2 consistency (d2/dh2)  : " + (-calc_dd1dh(d1)*calc_dd1dh(d1)/d1/d1 + calc_dd2dh(d1)*calc_dd2dh(d1)/d2/d2 + calc_d2d1dh2(d1)/d1 - calc_d2d2dh2(d1)/d2));
        //System.out.println("3 consistency (d2/dcdh) : " + (-calc_dd1dc_at_h(d1)*calc_dd1dh(d1)/d1/d1 + calc_dd2dc_at_h(d1)*calc_dd2dh(d1)/d2/d2 + calc_d2d1dcdh_alt(d1)/d1 - calc_d2d2dcdh_alt(d1)/d2));
        //System.out.println("4 consistency (d/dc)    : " + (2*calc_dd1dc_at_h(d1)*(l1 - l2/Math.sqrt(2)) + 2*calc_dd2dc_at_h(d1)*(l2 - l1/Math.sqrt(2)) - (d1*calc_dd2dc_at_h(d1) + d2*calc_dd1dc_at_h(d1))/Math.sqrt(2) + 2.*(d1 - d2)*(1. + 1./Math.sqrt(2)) - 20.*dspiro_areadc()/3.));
        //System.out.println("5 consistency (d/dh)    : " + (2*calc_dd1dh(d1)*(l1 - l2/Math.sqrt(2)) + 2*calc_dd2dh(d1)*(l2 - l1/Math.sqrt(2)) - (d1*calc_dd2dh(d1) + d2*calc_dd1dh(d1))/Math.sqrt(2)));
        //System.out.println("6 consistency (d/dc)    : " + (d2*calc_dd1dc_at_h(d1) - d1*calc_dd2dc_at_h(d1)));
        //System.out.println("7 consistency (d2/dh2)  : " + (calc_d2d1dh2(d1)*(2*(l1 - l2/Math.sqrt(2)) - d2/Math.sqrt(2)) + calc_d2d2dh2(d1)*(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2)) - 2*calc_dd1dh(d1)*calc_dd2dh(d1)/Math.sqrt(2)));
        //System.out.println("8 consistency (d2/dcdh) : " + (calc_d2d1dcdh(d1)*(2*(l1 - l2/Math.sqrt(2)) - d2/Math.sqrt(2)) + calc_d2d2dcdh(d1)*(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2))
        //                                                +  calc_dd1dh(d1)*(2*(1 + 1/Math.sqrt(2)) - calc_dd2dc_at_h(d1)/Math.sqrt(2)) + calc_dd2dh(d1)*(-2*(1 + 1/Math.sqrt(2)) - calc_dd1dc_at_h(d1)/Math.sqrt(2))));
        //System.out.println("cancel test = " + calc_dd1dh(d1)*calc_dd2dc_at_h(d1) + ", " + calc_dd2dh(d1)*calc_dd1dc_at_h(d1) + ", " + (calc_dd1dh(d1)*calc_dd2dc_at_h(d1) + calc_dd2dh(d1)*calc_dd1dc_at_h(d1))
        //                   + ", " + (20*dspiro_areadc()/3 + 2*(d2 - d1)*(1 + 1/Math.sqrt(2)))*2*d1*d2*(d2*(l2 - l1/Math.sqrt(2)) - d1*(l1 - l2/Math.sqrt(2)))/2/(d1*(l1 - l2/Math.sqrt(2)) + d2*(l2 - l1/Math.sqrt(2)) - d1*d2/Math.sqrt(2))/2/(d1*(l1 - l2/Math.sqrt(2)) + d2*(l2 - l1/Math.sqrt(2)) - d1*d2/Math.sqrt(2)));
        //System.out.println("verify cancel = " + calc_d2(d1) + ", " + calc_dd1dh_dd2dc(d1));
        //System.out.println("F = , " + (float) fitted.getc() + ", " + d1 + ", " + d2 + ", " + a_b*a_b*retVal*retVal/2);
        return retVal;
    }

    private static void scan_at_P2()
    {
        // calculate F over the whole range of (d1, d2)
        // holding area constant
        // to run this, temporarily comment out lines 248 and 284 to reduce the output
        // keep only line 285 : 'F = , ...'

        double d1;
        //double d_average = 57.1;
        double d1_start = 5; //d_average - 1;
        double d1_end = 50; //d_average + 1;
        int Nd1 = 4;

        fitted = new epiTrochoidFxn(12);
        System.out.println("scan F @ , c, d1, d2, " + fitted.getc());
        for (int i = 0; i <= Nd1; i++)
        {
            t2[N] = 0;                                  // just to control the output
            d1 = d1_start + i*(d1_end - d1_start)/Nd1;
            if (calc_d2(d1) < 0) break;
            solve_at_P2(d1, calc_d2(d1), false);
        }
    }

    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

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
        return fitted.getdxdc(t1_start)*(dN33(t2)[0] + dN33(t2)[1]) + fitted.getdxdc(t1_end)*(dN33(t2)[2] + dN33(t2)[3]);
    }

    private static double calc_d2fydudc(double t2)
    {
        return fitted.getdydc(t1_start)*(dN33(t2)[0] + dN33(t2)[1]) + fitted.getdydc(t1_end)*(dN33(t2)[2] + dN33(t2)[3]);
    }

    private static double spiro_area()
    {
        return a_b*a_b*(Math.PI/2 - Math.sqrt(2))/4                     // <1>
             - fitted.getc()*fitted.getc()*(3*Math.PI/2 - Math.sqrt(2))/4;
    }

    private static double dspiro_areadc()
    {
        return -fitted.getc()*(3*Math.PI/2 - Math.sqrt(2))/2;
    }

    private static void calc_d1_d2_at_h()
    {
        // satisfy area constraint at given c
        // satisfy h = ln(d1/d2)
        double h = -1.64;
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double tempT = l1 - l2/Math.sqrt(2) + Math.exp(-h)*(l2 - l1/Math.sqrt(2));
        double d1 = Math.sqrt(2)*Math.exp(h)*(tempT - Math.sqrt(tempT*tempT - Math.exp(-h)*20*spiro_area()/3/Math.sqrt(2)));
        System.out.println("calc_d1_d2_at_h = " + d1 + ", " + calc_d2(d1));
    }

    private static double calc_d2(double d1)
    {
        // satisfy area constraint
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20*spiro_area()/3 - 2*d1*(l1 - l2/Math.sqrt(2)))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_dd2dd1(double d1)
    {
        // satisfy area constraint (first derivative)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return -(2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_dd2dc(double d1)
    {
        // satisfy area constraint (first derivative)
        // holding d1 fixed
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return ((2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2))*(20./3.*dspiro_areadc() - 2*d1*(1 + 1/Math.sqrt(2)))
               + (20*spiro_area()/3 - 2*d1*(l1 - l2/Math.sqrt(2)))*2*(1 + 1/Math.sqrt(2)))
               /(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_dd2dc_at_d1(double d1)
    {
        // satisfy area constraint (first derivative)
        // holding d1 fixed
        // proposed improvement to calc_dd2dc (tested - works!)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20./3.*dspiro_areadc() + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2)))
              /(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_dd1dc_at_d2(double d1)
    {
        // satisfy area constraint (first derivative)
        // holding d2 fixed (tested - works!)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20./3.*dspiro_areadc() + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2)))
              /(2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2));
    }

    private static double calc_d2d1dcdd2_at_d2(double d1)
    {
        // satisfy area constraint (first derivative)
        // holding d2 fixed (tested - works!)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double Num = 20./3.*dspiro_areadc() + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2));
        double Den = 2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2);
        return (2*Den*(1. - 1./calc_dd2dd1(d1))*(1. + 1./Math.sqrt(2)) + Num/Math.sqrt(2))/Den/Den;
    }

    private static double calc_d2d2dcdd1_at_d1(double d1)
    {
        // satisfy area constraint (first derivative)
        // proposed improvement to calc_dd2dc (tested - works!)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double Num = 20./3.*dspiro_areadc() + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2));
        double Den = 2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2);
        return (2*Den*(calc_dd2dd1(d1) - 1.)*(1. + 1./Math.sqrt(2)) + Num/Math.sqrt(2))/Den/Den;
    }

    private static double calc_d2d2dd1dd1(double d1)
    {
        // satisfy area constraint (second derivative)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return Math.sqrt(2)*calc_dd2dd1(d1)/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_d2d2dd1dc(double d1)
    {
        // satisfy area constraint (second derivative)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return -((2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2))*(2*(1 + 1/Math.sqrt(2)) - calc_dd2dc(d1)/Math.sqrt(2))
               + (2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2))*2*(1 + 1/Math.sqrt(2)))
                /(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_dd1dh(double d1)
    {
        // dimensionless variable h = ln(d1/d2)
        // satisfy area constraint
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return d1*calc_d2(d1)*(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2))/2
              /(calc_d2(d1)*(l2 - l1/Math.sqrt(2)) + d1*(l1 - l2/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2));
    }

    private static double calc_dd2dh(double d1)
    {
        // dimensionless variable h = ln(d1/d2)
        // satisfy area constraint
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return -d1*calc_d2(d1)*(2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2))/2
               /(d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2));
//        return -d1*calc_d2(d1)*(2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2))
//               /(d1*(2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2)) + calc_d2(d1)*(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2)));
    }

    private static double calc_dd1dc_at_h(double d1)
    {
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20*dspiro_areadc()/3 + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2)))*d1/2
               /(d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2));
    }

    private static double calc_dd2dc_at_h(double d1)
    {
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20*dspiro_areadc()/3 + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2)))*calc_d2(d1)/2
               /(d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2));
    }

    private static double calc_dd1dh_dd2dc(double d1)
    {
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        // calculate the cross product : dd1dh*dd2dc + dd2dh*dd1dc
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20*dspiro_areadc()/3 + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2)))*2*d1*calc_d2(d1)*(calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*(l1 - l2/Math.sqrt(2)))/2/(d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2))/2/(d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2));
    }

    private static double calc_d2d1dh2(double d1)
    {
        // dimensionless variable h = ln(d1/d2)
        // satisfy area constraint
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double num1 = 1/d1*(l2 - l1/Math.sqrt(2)) + 1/calc_d2(d1)*(l1 - l2/Math.sqrt(2)) - 1/Math.sqrt(2);
        double num2 = 2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2);
        //return (-num1*calc_dd1dh(d1)/Math.sqrt(2) + num2*(1/d1/d1*calc_dd1dh(d1)*(l2 - l1/Math.sqrt(2)) + 1/calc_d2(d1)/calc_d2(d1)*calc_dd2dh(d1)*(l1 - l2/Math.sqrt(2))))
        //       /2/num1/num1;
        return (2./calc_d2(d1)*calc_dd1dh(d1)*calc_dd2dh(d1)/Math.sqrt(2) + num2*(calc_dd1dh(d1)/d1 + calc_dd2dh(d1)/calc_d2(d1)))/2/num1;
    }

    private static double calc_d2d2dh2(double d1)
    {
        // dimensionless variable h = ln(d1/d2)
        // satisfy area constraint
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double num1 = 1/calc_d2(d1)*(l1 - l2/Math.sqrt(2)) + 1/d1*(l2 - l1/Math.sqrt(2)) - 1/Math.sqrt(2);
        double num2 = 2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2);
        //return -(-num1*calc_dd2dh(d1)/Math.sqrt(2) + num2*(1/calc_d2(d1)/calc_d2(d1)*calc_dd2dh(d1)*(l1 - l2/Math.sqrt(2)) + 1/d1/d1*calc_dd1dh(d1)*(l2 - l1/Math.sqrt(2))))
        //        /2/num1/num1;
        return (2./d1*calc_dd1dh(d1)*calc_dd2dh(d1)/Math.sqrt(2) - num2*(calc_dd1dh(d1)/d1 + calc_dd2dh(d1)/calc_d2(d1)))/2/num1;
    }

    private static double calc_d2d1dcdh(double d1)
    {
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        // calculate d/dh at fixed c
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        //double Num = 20*dspiro_areadc()/3 + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2));
        //double Den = d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2);
        //return Num*calc_dd1dh(d1)/2/Den + d1/2/Den/Den
        //       *(Den*2*(calc_dd2dh(d1) - calc_dd1dh(d1))*(1 + 1/Math.sqrt(2))
        //       - Num*(calc_dd1dh(d1)*(l1 - l2/Math.sqrt(2)) + calc_dd2dh(d1)*(l2 - l1/Math.sqrt(2)) - (d1*calc_dd2dh(d1) + calc_d2(d1)*calc_dd1dh(d1))/Math.sqrt(2)));
        double num1 = 1/d1*(l2 - l1/Math.sqrt(2)) + 1/calc_d2(d1)*(l1 - l2/Math.sqrt(2)) - 1/Math.sqrt(2);
        double num2 = 2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2);
        double num3 = calc_dd1dh(d1)/d1 + calc_dd2dh(d1)/calc_d2(d1);
        return (2./calc_d2(d1)*(1 + 1/Math.sqrt(2))*(calc_dd2dh(d1) - calc_dd1dh(d1)) + calc_dd1dc_at_h(d1)*num3/Math.sqrt(2) + 1/calc_d2(d1)*num2*calc_dd2dc_at_h(d1))/2/num1;
    }

    private static double calc_d2d2dcdh(double d1)
    {
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        // calculate d/dh at fixed c
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        //double Num = 20*dspiro_areadc()/3 + 2*(calc_d2(d1) - d1)*(1 + 1/Math.sqrt(2));
        //double Den = d1*(l1 - l2/Math.sqrt(2)) + calc_d2(d1)*(l2 - l1/Math.sqrt(2)) - d1*calc_d2(d1)/Math.sqrt(2);
        //return Num*calc_dd2dh(d1)/2/Den + calc_d2(d1)/2/Den/Den
        //       *(Den*2*(calc_dd2dh(d1) - calc_dd1dh(d1))*(1 + 1/Math.sqrt(2))
        //       - Num*(calc_dd1dh(d1)*(l1 - l2/Math.sqrt(2)) + calc_dd2dh(d1)*(l2 - l1/Math.sqrt(2)) - (d1*calc_dd2dh(d1) + calc_d2(d1)*calc_dd1dh(d1))/Math.sqrt(2)));
        double num1 = 1/d1*(l2 - l1/Math.sqrt(2)) + 1/calc_d2(d1)*(l1 - l2/Math.sqrt(2)) - 1/Math.sqrt(2);
        double num2 = 2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2);
        double num3 = calc_dd1dh(d1)/d1 + calc_dd2dh(d1)/calc_d2(d1);
        return (2./d1*(1 + 1/Math.sqrt(2))*(calc_dd2dh(d1) - calc_dd1dh(d1)) + calc_dd2dc_at_h(d1)*num3/Math.sqrt(2) - 1/d1*num2*calc_dd1dc_at_h(d1))/2/num1;
    }

    private static double calc_d2d1dcdh_alt(double d1)
    {
        // alternate calc of d2d1dcdh (hopefully not needed)
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        // calculate d/dh at fixed c
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double Num = 2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2);
        double Den = 2*(1/d1*(l2 - l1/Math.sqrt(2)) + 1/calc_d2(d1)*(l1 - l2/Math.sqrt(2)) - 1/Math.sqrt(2));
        return (Den*(-2*(1 + 1/Math.sqrt(2)) - calc_dd1dc_at_h(d1)/Math.sqrt(2))
             - 2*Num*(-calc_dd1dc_at_h(d1)/d1/d1*(l2 - l1/Math.sqrt(2)) - calc_dd2dc_at_h(d1)/calc_d2(d1)/calc_d2(d1)*(l1 - l2/Math.sqrt(2))
             - 1/d1*(1 + 1/Math.sqrt(2)) + 1/calc_d2(d1)*(1 + 1/Math.sqrt(2))))/Den/Den;
    }

    private static double calc_d2d2dcdh_alt(double d1)
    {
        // alternate calc of d2d2dcdh (hopefully not needed)
        // satisfy area constraint
        // calculate d/dc at fixed h = ln(b1/b2)
        // calculate d/dh at fixed c
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        double Num = 2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2);
        double Den = 2*(1/d1*(l2 - l1/Math.sqrt(2)) + 1/calc_d2(d1)*(l1 - l2/Math.sqrt(2)) - 1/Math.sqrt(2));
        return -(Den*(2*(1 + 1/Math.sqrt(2)) - calc_dd2dc_at_h(d1)/Math.sqrt(2))
             - 2*Num*(-calc_dd1dc_at_h(d1)/d1/d1*(l2 - l1/Math.sqrt(2)) - calc_dd2dc_at_h(d1)/calc_d2(d1)/calc_d2(d1)*(l1 - l2/Math.sqrt(2))
             - 1/d1*(1 + 1/Math.sqrt(2)) + 1/calc_d2(d1)*(1 + 1/Math.sqrt(2))))/Den/Den;
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
