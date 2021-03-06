
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1.
// fit a 5-point Beta1-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// constrain only the slopes at the endpoints and keep P2 (Q3/Q4) arbitrary.
// the Beta1 spline has an asymmetric arm length, d- and d+, at the Bezier splice at Q3/Q4.
// linearize the equations wrt (d1, d2, x2, y2, beta1) and solve a 5x5 system of equations.
// decompose the Beta1-Spline into 2 Beziers, range (0,1) and (1,2).
// note that t2dd[i][j] will normally be discontinuous at t2 = 1 (if beta1 != 1)
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2).
// t2 must be chosen to minimize the distance to g(t1).
// see Spiro2SVG Book 4, Dec 2017, page 24.

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\Beta1Spline.java

public class Beta1Spline
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
    private static double[][] t2dd = new double[5][N+3];    // partial wrt (d1, d2, x2, y2, beta1); [N+1][N+2] is response before and after splice
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static double beta1;
    private static int splicei;                             // t1 index i before splice
    private static double spliced;                          // fractional portion of index i at splice
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 85;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //fitted = new epiTrochoidFxn(0.9045);    // keep for finite difference 2nd partial difference of F
        //fitted = new epiTrochoidFxn(4.175);
        //fitted = new epiTrochoidFxn(19.8 - .0025);
        fitted = new epiTrochoidFxn(0.0);
        //System.out.println("Beta1-Spline convert_at_P2 = " + convert_at_P2(11.244692865076667, 50.78732161191724, 190.0813110769242, 36.1022521075289, true) + "\n");
        //System.out.println("Beta1-Spline solve_at_P2 = " + solve_at_P2(25, 30, 166, 68, 1.5, true) + "\n");
        //System.out.println("Beta1-Spline solve_at_P2 = " + solve_at_P2(20.23, 24.41, 164.22, 72.57, 0.6227, true) + "\n");
        //System.out.println("Beta1-Spline iterate_at_P2 = "
        //                  + iterate_at_P2(19.894534447557504, 15.42257264455448, 172.73626099162183, 54.37109910205891, 2.085275) + "\n");
        System.out.println("Beta1-Spline solve_at_P2 = " 
                          + iterate_at_P2(23.849235501951565, 23.849235502011112, 166.29841931017188, 68.88306067938679, 1.0) + "\n");
        //System.out.println("Beta1-Spline iterate_at_P2 = " + iterate_at_P2(22.56, 19.97, 156.42, 79.28, 0.40) + "\n");
        //System.out.println("Beta1-Spline iterate_at_P2 = "
        //                  + iterate_at_P2(23.686962633324175, 18.326697796941833, 139.88954856635974, 86.58082601420564, 0.1675789980868262) + "\n");
        //System.out.println("Beta1-Spline iterate_at_P2 = "
        //                  + iterate_at_P2(20.238531117463033, 24.415493401132725, 164.22656338664095, 72.5702320294945, 0.6227750517813624) + "\n");
        //System.out.println("Beta1-Spline solve_at_P2 = "
        //                  + solve_at_P2(11.095840455163636, 35.76777395541423, 177.6310736380369, 32.57345106824122, 3.421498759870758, true) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
    }

    private static double iterate_at_P2(double d1, double d2, double x2, double y2, double m_beta1)
    {
        // calculate a new estimate of (d1, d2, x2, y2, m_beta1) by setting dF = 0
        // include second-order responses: d2f/ddi/ddj
        // see Spiro2SVG Book 4, page 24 (applied to 5-point cubic Beta1-Spline)
        // setup 5-variable Newton-Raphson iteration

        final double gain = 1;                              // fudge factor to reduce gain
        final int MAXLOOP = 2000;
        double[] f_gx = new double[N+3];                    // [N+1][N+2] is fxn before and after splice
        double[] f_gy = new double[N+3];
        double[] dfxdu = new double[N+3];
        double[] dfydu = new double[N+3];
        double[] d2fxdudu = new double[N+3];
        double[] d2fydudu = new double[N+3];
        double[][] dfxdd = new double[5][N+3];
        double[][] dfydd = new double[5][N+3];
        double[][] d2fxdudd = new double[5][N+3];
        double[][] d2fydudd = new double[5][N+3];
        double[][][] d2fxdddd = new double[5][5][N+3];
        double[][][] d2fydddd = new double[5][5][N+3];
        double[] df_gxdc = new double[N+1];                 // used only for calc of dx[]/dc
        double[] df_gydc = new double[N+1];
        double[][] d2fxdcdd = new double[5][N+1];
        double[][] d2fydcdd = new double[5][N+1];
        double[] t2ddc = new double[N+1];

        double[] trap_in = new double[N+3];
        double[][] Jac = new double[5][5];
        double[] dFdd = new double[5];
        double[] d2Fdddc = new double[5];                   // augmented matrix
        double dFdc;
        double d2Fdcdc;                                     // augmented matrix
        double[][] Augment = new double[6][6];
        double[] deld;                                      // (-Δd1, -Δd2, -Δx2, -Δy2, -Δβ1)
        int loop = 0;
        int i, j, k, seg;                                   // Bezier segment, before or after the splice
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(d1, d2, x2, y2, m_beta1, false)))   // initiallize at (x2, y2)
            {
                System.out.println("fail at " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + m_beta1);
                return Double.NaN;
            }

            seg = 0;
            //System.out.println("i , seg , t2[i], f_gx[i], f_gy[i], df_gxdc[i] , df_gydc[i] , dfxdd[4][i] , dfydd[4][i] , d2fxdcdd[4][i] , d2fydcdd[4][i]");
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
                df_gxdc[i] = calc_df_gxdc(seg, t1, t2[i] - seg);
                df_gydc[i] = calc_df_gydc(seg, t1, t2[i] - seg);
                d2fxdcdd[4][i] = calc_d2fxdcdbeta1(seg, t2[i] - seg);
                d2fydcdd[4][i] = calc_d2fydcdbeta1(seg, t2[i] - seg);
                dfxdd[0][i] = calc_dfxdd1(seg, t2[i] - seg);
                dfydd[0][i] = calc_dfydd1(seg, t2[i] - seg);
                dfxdd[1][i] = calc_dfxdd2(seg, t2[i] - seg);
                dfydd[1][i] = calc_dfydd2(seg, t2[i] - seg);
                dfxdd[2][i] = calc_dfxdx2(seg, t2[i] - seg);
                dfydd[2][i] = 0;
                dfxdd[3][i] = 0;
                dfydd[3][i] = calc_dfydy2(seg, t2[i] - seg);
                dfxdd[4][i] = calc_dfxdbeta1(seg, t2[i] - seg);
                dfydd[4][i] = calc_dfydbeta1(seg, t2[i] - seg);
                d2fxdudd[0][i] = calc_d2fxdudd1(seg, t2[i] - seg);
                d2fydudd[0][i] = calc_d2fydudd1(seg, t2[i] - seg);
                d2fxdudd[1][i] = calc_d2fxdudd2(seg, t2[i] - seg);
                d2fydudd[1][i] = calc_d2fydudd2(seg, t2[i] - seg);
                d2fxdudd[2][i] = calc_d2fxdudx2(seg, t2[i] - seg);
                d2fydudd[2][i] = 0;
                d2fxdudd[3][i] = 0;
                d2fydudd[3][i] = calc_d2fydudy2(seg, t2[i] - seg);
                d2fxdudd[4][i] = calc_d2fxdudbeta1(seg, t2[i] - seg);
                d2fydudd[4][i] = calc_d2fydudbeta1(seg, t2[i] - seg);

                d2fxdddd[0][4][i] = calc_d2fxdd1dbeta1(seg, t2[i] - seg);
                d2fxdddd[1][4][i] = calc_d2fxdd2dbeta1(seg, t2[i] - seg);
                d2fxdddd[2][4][i] = calc_d2fxdx2dbeta1(seg, t2[i] - seg);
                d2fxdddd[3][4][i] = 0;
                d2fxdddd[4][4][i] = calc_d2fxdbeta1dbeta1(seg, t2[i] - seg);
                d2fxdddd[4][0][i] = d2fxdddd[0][4][i];
                d2fxdddd[4][1][i] = d2fxdddd[1][4][i];
                d2fxdddd[4][2][i] = d2fxdddd[2][4][i];
                d2fxdddd[4][3][i] = d2fxdddd[3][4][i];
                d2fydddd[0][4][i] = calc_d2fydd1dbeta1(seg, t2[i] - seg);
                d2fydddd[1][4][i] = calc_d2fydd2dbeta1(seg, t2[i] - seg);
                d2fydddd[2][4][i] = 0;
                d2fydddd[3][4][i] = calc_d2fydy2dbeta1(seg, t2[i] - seg);
                d2fydddd[4][4][i] = calc_d2fydbeta1dbeta1(seg, t2[i] - seg);
                d2fydddd[4][0][i] = d2fydddd[0][4][i];
                d2fydddd[4][1][i] = d2fydddd[1][4][i];
                d2fydddd[4][2][i] = d2fydddd[2][4][i];
                d2fydddd[4][3][i] = d2fydddd[3][4][i];
                //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + f_gx[i] + ", " + f_gy[i] + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + dfxdd[0][i] + ", " + dfydd[0][i] + ", " + dfxdd[1][i] + ", " + dfydd[1][i] + ", " + dfxdd[2][i] + ", " + dfydd[2][i] + ", " + dfxdd[3][i] + ", " + dfydd[3][i] + ", " + dfxdd[4][i] + ", " + dfydd[4][i] + ", " + d2fxdudd[0][i] + ", " + d2fydudd[0][i] + ", " + d2fxdudd[1][i] + ", " + d2fydudd[1][i] + ", " + d2fxdudd[2][i] + ", " + d2fydudd[2][i] + ", " + d2fxdudd[3][i] + ", " + d2fydudd[3][i] + ", " + d2fxdudd[4][i] + ", " + d2fydudd[4][i]);
                //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + f_gx[i] + ", " + f_gy[i] + ", " + df_gxdc[i] + ", " + df_gydc[i] + ", " + dfxdd[4][i] + ", " + dfydd[4][i] + ", " + d2fxdcdd[4][i] + ", " + d2fydcdd[4][i]);

                t2ddc[i] = calc_t2dd(i, seg, t2[i] - seg, "c");
            }
            for (i = N + 1; i <= N + 2; i++)    // [N+1][N+2] is data before and after splice
            {
                seg = i - N - 1;
                t1 = t1_start + (splicei + spliced)*(t1_end - t1_start)/N;
                f_gx[i] = t2_vs_t1.fn(Bezx[seg], 1 - seg) - fitted.getx(t1);    // use this method preferentially because it retains the symmetry of Jac
                f_gy[i] = t2_vs_t1.fn(Bezy[seg], 1 - seg) - fitted.gety(t1);    // see hypoBeta1Splined1d2rms_org.csv
                dfxdu[i] = t2_vs_t1.dfn(Bezx[seg], 1 - seg);
                d2fxdudu[i] = t2_vs_t1.d2fn(Bezx[seg], 1 - seg);
                dfydu[i] = t2_vs_t1.dfn(Bezy[seg], 1 - seg);
                d2fydudu[i] = t2_vs_t1.d2fn(Bezy[seg], 1 - seg);
                dfxdd[0][i] = calc_dfxdd1(seg, 1 - seg);
                dfydd[0][i] = calc_dfydd1(seg, 1 - seg);
                dfxdd[1][i] = calc_dfxdd2(seg, 1 - seg);
                dfydd[1][i] = calc_dfydd2(seg, 1 - seg);
                dfxdd[2][i] = calc_dfxdx2(seg, 1 - seg);
                dfydd[2][i] = 0;
                dfxdd[3][i] = 0;
                dfydd[3][i] = calc_dfydy2(seg, 1 - seg);
                dfxdd[4][i] = calc_dfxdbeta1(seg, 1 - seg);
                dfydd[4][i] = calc_dfydbeta1(seg, 1 - seg);
                d2fxdudd[0][i] = calc_d2fxdudd1(seg, 1 - seg);
                d2fydudd[0][i] = calc_d2fydudd1(seg, 1 - seg);
                d2fxdudd[1][i] = calc_d2fxdudd2(seg, 1 - seg);
                d2fydudd[1][i] = calc_d2fydudd2(seg, 1 - seg);
                d2fxdudd[2][i] = calc_d2fxdudx2(seg, 1 - seg);
                d2fydudd[2][i] = 0;
                d2fxdudd[3][i] = 0;
                d2fydudd[3][i] = calc_d2fydudy2(seg, 1 - seg);
                d2fxdudd[4][i] = calc_d2fxdudbeta1(seg, 1 - seg);
                d2fydudd[4][i] = calc_d2fydudbeta1(seg, 1 - seg);

                d2fxdddd[0][4][i] = calc_d2fxdd1dbeta1(seg, 1 - seg);
                d2fxdddd[1][4][i] = calc_d2fxdd2dbeta1(seg, 1 - seg);
                d2fxdddd[2][4][i] = calc_d2fxdx2dbeta1(seg, 1 - seg);
                d2fxdddd[3][4][i] = 0;
                d2fxdddd[4][4][i] = calc_d2fxdbeta1dbeta1(seg, 1 - seg);
                d2fxdddd[4][0][i] = d2fxdddd[0][4][i];
                d2fxdddd[4][1][i] = d2fxdddd[1][4][i];
                d2fxdddd[4][2][i] = d2fxdddd[2][4][i];
                d2fxdddd[4][3][i] = d2fxdddd[3][4][i];
                d2fydddd[0][4][i] = calc_d2fydd1dbeta1(seg, 1 - seg);
                d2fydddd[1][4][i] = calc_d2fydd2dbeta1(seg, 1 - seg);
                d2fydddd[2][4][i] = 0;
                d2fydddd[3][4][i] = calc_d2fydy2dbeta1(seg, 1 - seg);
                d2fydddd[4][4][i] = calc_d2fydbeta1dbeta1(seg, 1 - seg);
                d2fydddd[4][0][i] = d2fydddd[0][4][i];
                d2fydddd[4][1][i] = d2fydddd[1][4][i];
                d2fydddd[4][2][i] = d2fydddd[2][4][i];
                d2fydddd[4][3][i] = d2fydddd[3][4][i];
                //System.out.println(i + ", " + seg + ", " + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + f_gx[i] + ", " + f_gy[i] + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + dfxdd[0][i] + ", " + dfydd[0][i] + ", " + dfxdd[1][i] + ", " + dfydd[1][i] + ", " + dfxdd[2][i] + ", " + dfydd[2][i] + ", " + dfxdd[3][i] + ", " + dfydd[3][i] + ", " + dfxdd[4][i] + ", " + dfydd[4][i] + ", " + d2fxdudd[0][i] + ", " + d2fydudd[0][i] + ", " + d2fxdudd[1][i] + ", " + d2fydudd[1][i] + ", " + d2fxdudd[2][i] + ", " + d2fydudd[2][i] + ", " + d2fxdudd[3][i] + ", " + d2fydudd[3][i] + ", " + d2fxdudd[4][i] + ", " + d2fydudd[4][i]);
            }
            //for (i = N + 1; i <= N + 2; i++)    // [N+1][N+2] is data before and after splice
            //    System.out.println("before/after = " + t2dd[0][i] + ", " + d2fxdddd[2][4][i]);
            //System.out.println("before/after = " + t2dd[0][N+2]/t2dd[0][N+1] + ", " + d2fydddd[0][4][N+2]/d2fydddd[0][4][N+1]);

            // calc dFdd[j] at current (d1, d2, x2, y2, beta1)

            for (i = 0; i < 5; i++)
            {
                //System.out.println("dFdd " + i);
                for (k = 0; k <= N + 2; k++)
                {
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]); // original code
                    //trap_in[k] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];         // new code
                    //System.out.println(i + ", " + k + ", " + trap_in[k]);
                    //System.out.println(i + ", " + k + ", " + trap_in[k] + ", " + (f_gx[k]*dfxdu[k] + f_gy[k]*dfydu[k]) + ", " + t2dd[i][k]);
                }
                dFdd[i] = integrate(trap_in);
                for (k = 0; k <= N; k++)
                    trap_in[k] = dfxdd[i][k]*df_gxdc[k] + dfydd[i][k]*df_gydc[k]
                               + f_gx[k]*d2fxdcdd[i][k] + f_gy[k]*d2fydcdd[i][k]
                               - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2ddc[k];
                d2Fdddc[i] = integrate(trap_in);                        // augmented matrix
            }

            // calc d2Fdd[i]dd[j] (symmetric Jacobean matrix)

            for (i = 0; i < 5; i++)
                for (j = 0; j < 5; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N + 2; k++)
                    {
                        //trap_in[k] = (dfxdd[i][k] + dfxdu[k]*t2dd[i][k])*(dfxdd[j][k] + dfxdu[k]*t2dd[j][k])    // original code
                        //           + f_gx[k]*(d2fxdddd[i][j][k] + dfxdu[k]*d2udddd[i][j][k]
                        //           + d2fxdudd[i][k]*t2dd[j][k] + d2fxdudd[j][k]*t2dd[i][k] + d2fxdudu[k]*t2dd[i][k]*t2dd[j][k])
                        //           + (dfydd[i][k] + dfydu[k]*t2dd[i][k])*(dfydd[j][k] + dfydu[k]*t2dd[j][k])
                        //           + f_gy[k]*(d2fydddd[i][j][k] + dfydu[k]*d2udddd[i][j][k]
                        //           + d2fydudd[i][k]*t2dd[j][k] + d2fydudd[j][k]*t2dd[i][k] + d2fydudu[k]*t2dd[i][k]*t2dd[j][k]);
                        //if (i == 2 && j == 2 && k <= N && true)           // fix fix use only for response calculations of F
                        //    System.out.println(k + ", " + ", " + t2[k] + ", " + trap_in[k] + ", " + d2udddd[i][j][k] + ", " + f_gx[k] + ", " + dfxdu[k] + ", " + f_gy[k] + ", " + dfydu[k]);
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k]                // new code
                                   + f_gx[k]*d2fxdddd[i][j][k]
                                   + dfydd[i][k]*dfydd[j][k]
                                   + f_gy[k]*d2fydddd[i][j][k]
                                   - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k];
                        //if (j == 4)
                            //System.out.println(k + ", " + (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k]);
                            //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = integrate(trap_in);
                }

            // calculate determinant of augmented matrix

            for (k = 0; k <= N; k++)
                trap_in[k] = f_gx[k]*df_gxdc[k] + f_gy[k]*df_gydc[k];
            dFdc = t2_vs_t1.integrate(trap_in);
            for (k = 0; k <= N; k++)
                trap_in[k] = df_gxdc[k]*df_gxdc[k] + df_gydc[k]*df_gydc[k]
                           - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2ddc[k]*t2ddc[k];
            d2Fdcdc = t2_vs_t1.integrate(trap_in);
            for (i = 0; i < 5; i++)
                for (j = 0; j < 5; j++)
                    Augment[i][j] = Jac[i][j];
            for (i = 0; i < 5; i++)
            {
                Augment[i][5] = d2Fdddc[i];
                Augment[5][i] = Augment[i][5];
            }
            Augment[5][5] = d2Fdcdc;
            //System.out.println("dFdc = " + fitted.getc() + ", " + d1 + ", " + d2 + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);

            //deld = BSpline5.multmv(BSpline5.invertm(Jac), dFdd);  // this is actually the negative of Δd
            //deld = BSpline5.multmv(BSpline5.gaussj(Jac), dFdd);  // this is actually the negative of Δd
            deld = BSpline5.gaussj(Jac, dFdd);  // this is actually the negative of Δd
            d1 -= deld[0]/gain;                 // gain is just a fudge factor to 'improve' convergence
            d2 -= deld[1]/gain;
            x2 -= deld[2]/gain;
            y2 -= deld[3]/gain;
            m_beta1 -= deld[4]/gain;

            //System.out.println("Jac");
            //for (i = 0; i < Jac.length; i++)
            //{
            //    for (j = 0; j < Jac.length; j++)
            //        System.out.print(Jac[i][j] + ", ");
            //    System.out.println();
            //}
            //splicei = 35;                        // debug integrate() routine
            //spliced = 1.0;
            //for (k = 0; k <= splicei; k++)
            //    trap_in[k] = 5 + 2*k;
            //for (k = splicei + 1; k <= N; k++)
            //    trap_in[k] = 2 + 3*k;
            //trap_in[N + 1] = 5 + 2*(splicei + spliced);
            //trap_in[N + 2] = 2 + 3*(splicei + spliced);
            //System.out.println("test integrate , " + splicei + ", " + spliced + ", " + integrate(trap_in));

            Jacdet = BSpline5.detm(Jac);
            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + dFdd[2] + ", " + dFdd[3] + ", " + dFdd[4] + ", " + Jacdet);
            System.out.println("deld = " + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ", " + deld[4]);
            System.out.println("\ndFdc = " + fitted.getc() + ", " + d1 + ", " + d2 + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);
            BSpline5.dump_Jac(Jac);
            BSpline5.dump_Jac(Augment);
            //BSpline5.gaussj(new double[][] {{1,2,1,4,3}, {2,3,1,5,2}, {1,1,2,4,3}, {4,5,4,5,5}, {3,2,3,5,4}}, new double[] {3.1, 3.2, 3.9, 7.2, 2.3});       // fix fix test code

            // temporary code to evaluate behavior at splice
/*            for (int jspl = 0; jspl < 5; jspl++)
            {
                System.out.println("jspl = " + jspl);
                System.out.println("x before ," + f_gx[N + 1] + ", " + dfxdd[jspl][N + 1] + ", " + dfxdu[N + 1] + ", " + t2dd[jspl][N + 1]);
                System.out.println("x after  ," + f_gx[N + 2] + ", " + dfxdd[jspl][N + 2] + ", " + dfxdu[N + 2] + ", " + t2dd[jspl][N + 2]);
                System.out.println("y before ," + f_gy[N + 1] + ", " + dfydd[jspl][N + 1] + ", " + dfydu[N + 1] + ", " + t2dd[jspl][N + 1]);
                System.out.println("y after  ," + f_gy[N + 2] + ", " + dfydd[jspl][N + 2] + ", " + dfydu[N + 2] + ", " + t2dd[jspl][N + 2]);
            }*/
            // end of temporary code to evaluate behavior at splice

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                for (j = 0; j < 5; j++)
                    t2[i] -= t2dd[j][i]*deld[j];                        // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL) && (Math.abs(deld[2]) < TOL) && (Math.abs(deld[3]) < TOL) && (Math.abs(deld[4]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 x2 y2 m_beta1 = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + m_beta1);
            // calculate dddc[i]
            double[] dt2dc = BSpline5.gaussj(Jac, d2Fdddc);
            System.out.println("\nfinal Beta1Spline, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) -dt2dc[2] + ", " + (float) -dt2dc[3] + ", " + (float) -dt2dc[4] + ", " + (float) BSpline5.detm(Augment));
            return solve_at_P2(d1, d2, x2, y2, m_beta1, true);                    // final run just for good measure
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

        return iterate_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, 1);
        //return solve_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, 1, print);
    }

    private static double solve_at_P2(double d1, double d2, double x2, double y2, double m_beta1, boolean print)
    {
        // 5-point cubic Beta1-Spline curve, decomposed as two Beziers
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error

        beta1 = m_beta1;
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        //d1 += 0.002;
        //d2 += 0.002;
        //x2 += 0.001;
        //y2 += 0.001;                          // temporary code for calculating response functions
        //beta1 += 0.00002;
        //System.out.println("Beta1 Spline u(t)");
        //System.out.println("d1, " + d1);
        //System.out.println("d2, " + d2);
        //System.out.println("x2, " + x2);
        //System.out.println("y2, " + y2);
        //System.out.println("beta1, " + beta1);

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

        double delx = (beta1*beta1*(Bezx[0][3] - Bezx[0][1]) + Bezx[1][2] - Bezx[1][0])/2/(beta1 + 1);
        double dely = (beta1*beta1*(Bezy[0][3] - Bezy[0][1]) + Bezy[1][2] - Bezy[1][0])/2/(beta1 + 1);
        Bezx[0][2] -= delx/beta1;
        Bezx[1][1] += delx;
        Bezy[0][2] -= dely/beta1;
        Bezy[1][1] += dely;
        //double dmins = Math.sqrt((Bezx[0][3] - Bezx[0][2])*(Bezx[0][3] - Bezx[0][2]) + (Bezy[0][3] - Bezy[0][2])*(Bezy[0][3] - Bezy[0][2]));
        //double dplus = Math.sqrt((Bezx[1][1] - Bezx[1][0])*(Bezx[1][1] - Bezx[1][0]) + (Bezy[1][1] - Bezy[1][0])*(Bezy[1][1] - Bezy[1][0]));
        //System.out.println("beta1 = ," + beta1 + ", " + dmins + ", " + dplus);
        //System.out.println("test coordinates = ," + beta1 + ", " + Bezx[0][1] + ", " + Bezy[0][1] + ", " + Bezx[0][3] + ", " + Bezy[0][3] + ", " + Bezx[1][2] + ", " + Bezy[1][2]);

        if (t2[N] == 0)
            System.out.println("__start Beta1-Spline at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1 + ", " + calc_error());
        //print_data_at_splice();
        //fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2, t2dd1, t2dd2, t2dx2, t2dy2, t2dbeta1, t2dc");
        int seg = 0;                // Bezier segment, before or after the splice
        for (int i = 0; i <= N; i++)
        {
            if (i == 0)
                t2[i] = solve_quintic_for_t2(i, 0, seg);
            else
                t2[i] = solve_quintic_for_t2(i, t2[i-1], seg);
            if (seg == 0 && ((t2[i] > 1) || Double.isNaN(t2[i])))
            {
                if (Double.isNaN(t2[i]))
                {
                    //scan_quintic_simple(i - 1);
                    //scan_quintic_simple(i);
                    System.out.println("restart due to Double.isNaN");
                }
                seg++;
                //t2[i] = solve_quintic_for_t2(i, t2[i-1], seg);   // re-calculate
                t2[i] = solve_quintic_for_t2(i, 1, seg);   // re-calculate (initiallize at the splice)
                splicei = i - 1;
                //spliced = (1 - t2[i - 1])/(t2[i] - t2[i - 1]);
                spliced = (1 - t2[i - 1])/(beta1*(t2[i] - 1) + 1 - t2[i - 1]);    // test code improved estimate
                if (splicei < 2)
                {
                    System.out.println("Abort: Bad splicei = " + seg + ", " + i + ", " + splicei);
                    //scan_quintic_simple(i);
                    return Double.NaN;
                }
                if (splicei > N - 3)
                {
                    //System.out.println("Abort: Bad splicei = " + splicei);
                    //return Double.NaN;
                }
            }
            if (seg == 1 && t2[i] < 1)
                return Double.NaN;
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //scan_quintic_simple(i - 1);
                //scan_quintic_simple(i);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dd(i, seg, t2[i] - seg, "d1");
            t2dd[1][i] = calc_t2dd(i, seg, t2[i] - seg, "d2");
            t2dd[2][i] = calc_t2dd(i, seg, t2[i] - seg, "x2");
            t2dd[3][i] = calc_t2dd(i, seg, t2[i] - seg, "y2");
            t2dd[4][i] = calc_t2dd(i, seg, t2[i] - seg, "beta1");
            if (print)
            {
                //double t1 = t1_start + i*(t1_end - t1_start)/N;
                //double f = (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))*t2_vs_t1.dfn(Bezx[seg], t2[i] - seg) + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1))*t2_vs_t1.dfn(Bezy[seg], t2[i] - seg);
                //System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + f);
                System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + calc_t2dd(i, seg, t2[i] - seg, "c"));
            }
        }
        System.out.println("start splice at c = " + fitted.getc() + ", " + splicei + ", " + spliced + ", " + solve_quintic_for_t2(splicei + spliced, t2[splicei], 0)); // + ", " + solve_quintic_for_t2(splicei + spliced, 1, 1)); // test code
/*  temporarily disable refine splice because it is probably not needed
        if (t2[splicei - 2] - 2*t2[splicei - 1] + t2[splicei] < 0)    // ensure negative curvature below
            spliced = refine_splice_t1_LHS(splicei);
        else if (t2[splicei + 1] - 2*t2[splicei + 2] + t2[splicei + 3] > 0) // ensure positive curvature above
            spliced = refine_splice_t1_RHS(splicei);
        else
        {
            System.out.println("opposite curvatures at the splice = " + (t2[splicei - 2] - 2*t2[splicei - 1] + t2[splicei])
                                                               + ", " + (t2[splicei + 1] - 2*t2[splicei + 2] + t2[splicei + 3]));
            spliced = refine_splice_t1_LHS(splicei);
        }
*/
        //System.out.println("end   splice at " + fitted.getc() + ", " + splicei + ", " + spliced); // test code
        //print_fx_fy(d1, d2, x2, y2);          // calculate fx and fy response functions
        //print_f_g_squared(d1, d2, x2, y2);    // calculate response of F error function
        for (int i = N + 1; i <= N + 2; i++)    // [N+1][N+2] is data before and after splice
        {
            seg = i - N - 1;
            t2dd[0][i] = calc_t2dd(i, seg, 1 - seg, "d1");
            t2dd[1][i] = calc_t2dd(i, seg, 1 - seg, "d2");
            t2dd[2][i] = calc_t2dd(i, seg, 1 - seg, "x2");
            t2dd[3][i] = calc_t2dd(i, seg, 1 - seg, "y2");
            t2dd[4][i] = calc_t2dd(i, seg, 1 - seg, "beta1");
            if (print)
            {
                //double t1 = t1_start + i*(t1_end - t1_start)/N;
                //double f = (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))*t2_vs_t1.dfn(Bezx[seg], t2[i] - seg) + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1))*t2_vs_t1.dfn(Bezy[seg], t2[i] - seg);
                //System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i] + ", " + f);
                System.out.println(i + ", " + (t1_start + (splicei + spliced)*(t1_end - t1_start)/N) + ", " + "'1'" + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i] + ", " + t2dd[4][i]);
            }
        }
        double retVal = calc_error();
        //System.out.println("__new t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1 + ", " + (float) retVal + ", " + (float) Jacdet + ", " + (float) (splicei + spliced) + ", " + (float) (Math.sqrt(delx*delx + dely*dely)/beta1));
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1 + ", " + retVal + ", " + (float) Jacdet + ", " + (float) (splicei + spliced));
        return retVal;
    }

    private static void print_f_g_squared(double d1, double d2, double x2, double y2)
    {
        // numerically calculate the function (f - g)**2
        // this will be the full response function, including variations in t2[i]
        // this was done for the purpose of independently confirming the Jac matrix
        int seg = 0;
        double rx, ry, t1;

        System.out.println("Beta1 Spline (f - g)**2 : c = " + fitted.getc());
        System.out.println("d1, " + d1);
        System.out.println("d2, " + d2);
        System.out.println("x2, " + x2);
        System.out.println("y2, " + y2);
        System.out.println("beta1, " + beta1);

        for (int i = 0; i <= N; i++)
        {
            if (seg == 0 && t2[i] > 1)
                seg++;
            t1 = t1_start + i*(t1_end - t1_start)/N;
            rx = t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1);
            ry = t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1);
            System.out.println(i + ", " + seg + ", " + t2[i] + ", " + 0.5*(rx*rx + ry*ry));
        }
        System.out.println();
    }

    private static void print_fx_fy(double d1, double d2, double x2, double y2)
    {
        // numerically calculate the response functions dfxddj, dfyddj and second derivatives
        // holding u(t) fixed (i.e. t2[i] fixed)
        // this was done for the purpose of independently confirming the functions calc_d2fxddjdbeta1 and calc_d2fxdbeta1dbeta1
        int seg = 0;
        //fitted = new epiTrochoidFxn(3.0 + 0.001);
        //y2 += 0.005;
        //beta1 += 0.00005;

        System.out.println("Beta1 Spline fx fy : c =, " + fitted.getc());
        System.out.println("d1, " + d1);
        System.out.println("d2, " + d2);
        System.out.println("x2, " + x2);
        System.out.println("y2, " + y2);
        System.out.println("beta1, " + beta1);

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

        double delx = (beta1*beta1*(Bezx[0][3] - Bezx[0][1]) + Bezx[1][2] - Bezx[1][0])/2/(beta1 + 1);
        double dely = (beta1*beta1*(Bezy[0][3] - Bezy[0][1]) + Bezy[1][2] - Bezy[1][0])/2/(beta1 + 1);
        Bezx[0][2] -= delx/beta1;
        Bezx[1][1] += delx;
        Bezy[0][2] -= dely/beta1;
        Bezy[1][1] += dely;

        for (int i = 0; i <= N; i++)
        {
            if (seg == 0 && t2[i] > 1)
                seg++;
            System.out.println(i + ", " + seg + ", " + t2[i] + ", " + t2_vs_t1.fn(Bezx[seg], t2[i] - seg) + ", " + t2_vs_t1.fn(Bezy[seg], t2[i] - seg));
            //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + t2_vs_t1.dfn(Bezx[seg], t2[i] - seg) + ", " + t2_vs_t1.dfn(Bezy[seg], t2[i] - seg) + ", " + t2_vs_t1.d2fn(Bezx[seg], t2[i] - seg) + ", " + t2_vs_t1.d2fn(Bezy[seg], t2[i] - seg));
            //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + calc_d2fxdbeta1dbeta1(seg, t2[i] - seg) + ", " + calc_d2fydbeta1dbeta1(seg, t2[i] - seg));
            //System.out.println(i + ", " + seg + ", " + t2[i] + ", " + calc_d2fxdcdbeta1(seg, t2[i] - seg) + ", " + calc_d2fydcdbeta1(seg, t2[i] - seg));
        }
        System.out.println();
    }

    private static void print_data_at_splice()
    {
        System.out.println("printing data before and after splice");
        System.out.println("fx[i]     = " + t2_vs_t1.fn(Bezx[0], 1) + ", " + t2_vs_t1.fn(Bezx[1], 0));
        System.out.println("fy[i]     = " + t2_vs_t1.fn(Bezy[0], 1) + ", " + t2_vs_t1.fn(Bezy[1], 0));
        System.out.println("dfxdu[i]  = " + t2_vs_t1.dfn(Bezx[0], 1) + ", " + t2_vs_t1.dfn(Bezx[1], 0));
        System.out.println("dfydu[i]  = " + t2_vs_t1.dfn(Bezy[0], 1) + ", " + t2_vs_t1.dfn(Bezy[1], 0));
        System.out.println("dfxdd1[i] = " + calc_dfxdd1(0, 1) + ", " + calc_dfxdd1(1, 0));
        System.out.println("dfydd1[i] = " + calc_dfydd1(0, 1) + ", " + calc_dfydd1(1, 0));
        System.out.println("dfxdd2[i] = " + calc_dfxdd2(0, 1) + ", " + calc_dfxdd2(1, 0));
        System.out.println("dfydd2[i] = " + calc_dfydd2(0, 1) + ", " + calc_dfydd2(1, 0));
        System.out.println("dfxdx2[i] = " + calc_dfxdx2(0, 1) + ", " + calc_dfxdx2(1, 0));
        System.out.println("dfydx2[i] = 0");
        System.out.println("dfxdy2[i] = 0");
        System.out.println("dfydy2[i] = " + calc_dfydy2(0, 1) + ", " + calc_dfydy2(1, 0));
        System.out.println("dfxdd[i]  = " + calc_dfxdbeta1(0, 1) + ", " + calc_dfxdbeta1(1, 0));
        System.out.println("dfydd[i]  = " + calc_dfydbeta1(0, 1) + ", " + calc_dfydbeta1(1, 0));
        System.out.println("d2fxdudd1[i] = " + calc_d2fxdudd1(0, 1) + ", " + calc_d2fxdudd1(1, 0));
        System.out.println("d2fydudd1[i] = " + calc_d2fydudd1(0, 1) + ", " + calc_d2fydudd1(1, 0));
        System.out.println("d2fxdudd2[i] = " + calc_d2fxdudd2(0, 1) + ", " + calc_d2fxdudd2(1, 0));
        System.out.println("d2fydudd2[i] = " + calc_d2fydudd2(0, 1) + ", " + calc_d2fydudd2(1, 0));
        System.out.println("d2fxdudx2[i] = " + calc_d2fxdudx2(0, 1) + ", " + calc_d2fxdudx2(1, 0));
        System.out.println("d2fydudx2[i] = 0");
        System.out.println("d2fxdudy2[i] = 0");
        System.out.println("d2fydudy2[i] = " + calc_d2fydudy2(0, 1) + ", " + calc_d2fydudy2(1, 0));
        System.out.println("d2fxdudd[i]  = " + calc_d2fxdudbeta1(0, 1) + ", " + calc_d2fxdudbeta1(1, 0));
        System.out.println("d2fydudd[i]  = " + calc_d2fydudbeta1(0, 1) + ", " + calc_d2fydudbeta1(1, 0));
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

    private static double integrate(double[] trap)
    {
        // trapezoidal rule integration of a fxn of t1 (N+1 points)
        // compensate for discontinuity at i = splicei + spliced

        double ret = (trap[0] + trap[N])/2;
        for (int i = 1; i < N; i++)
            ret += trap[i];
        //ret -= (trap[splicei] + trap[splicei + 1])/2;
        //ret += spliced*(trap[splicei] + trap[N + 1])/2 + (1 - spliced)*(trap[N + 2] + trap[splicei + 1])/2;
        return ret/N;
    }

    private static double solve_quintic_for_t2(double d, double t, int seg)
    {
        // d = integer index of t1 value (0-N)
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
        if (d == N && seg == 1)
            t = 1;              // special initiallization for endpoint
        f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
        fprime = t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d2fn(Bezx[seg], t) + t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.dfn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d2fn(Bezy[seg], t);
        f2prime = 3*t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.d2fn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d3fn(Bezx[seg], t) + 3*t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.d2fn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d3fn(Bezy[seg], t);
        //System.out.println("\ninit  t1 t2 =, " + (int) d + ", " + seg + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
        {
            del_t = -fprime/f2prime;
            //System.out.println("no roots del_t = " + del_t + ", " + (t + del_t));
        }
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            //System.out.println("first del_t = " + del_t + ", " + (t + del_t));
            if (del_t < 0 || del_t > 4) // used to be 2)
            {
                del_t = (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
                //System.out.println("second del_t = " + del_t + ", " + (t + del_t));
                if (del_t < 0 || del_t > 2)
                {
                    System.out.println("\nBad init  t1 t2 =, " + d + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
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
            if (loop > 200)
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
            //System.out.println("         t2 =, " + (int) d + ", " + seg + ", " + t + ", " + f + ", " + fprime);
        //} while (Math.abs(del_t) > TOL);
        } while (success < 2);
        //} while ((Math.abs(del_t) > TOL) || (Math.abs(f) > TOL/1000));
        return t + seg;                                 // compensate for Bezier segment offset
    }

    private static void scan_quintic_simple(int index)
    {
        // if solve_quintic_for_t2 fails, scan the area for other roots
        double t1 = t1_start + index*(t1_end - t1_start)/N;
        double f, t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        int seg = 0;

        System.out.println("\nscanning at " + index + ", " + seg + ", " + t1);
        for (int i = 0; i <= 200; i++)
        {
            if (i == 100) seg++;
            t = i/100.0;
            t -= seg;               // compensate for Bezier segment offset
            f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
            System.out.println(seg + ", " + (t + seg) + ", " + f);
        }
        System.out.println();
    }

    private static double refine_splice_t1_LHS(int i)
    {
        // improve the initial estimate of variable spliced using regula falsi
        // i = splicei, the index of the original estimate of the t1 value
        double xold = i;
        double xnew = i + 1;
        double fold = t2[i];
        double fnew = t2[i + 1];
        int seg;
        int loop = 0;

        System.out.println("LHS xold = ," + xold + ", " + fold + ", " + fnew);
        do
        {
            loop++;
            seg = 0;
            if (fnew > 1)
                xnew = xold + (xnew - xold)/(1 + beta1*(fnew - 1)/(1 - fold));  // if we bracket the splice
            else
                xnew = xold + (xnew - xold)*(1 - fold)/(fnew - fold);           // both points below the splice
            fnew = solve_quintic_for_t2(xnew, t2[i], seg);
            if (fnew > 1)
            {
                seg++;
                System.out.println("recheck = ," + xnew + ", " + fold + ", " + fnew + ", " + seg);
                fnew = solve_quintic_for_t2(xnew, t2[i], seg);
                if (fnew < 1)
                {
                    System.out.println("error: reversal of xnew = ," + xnew + ", " + fold + ", " + fnew);
                    return Double.NaN;
                }
            }
            System.out.println("xnew    = ," + xnew + ", " + fold + ", " + fnew + ", " + seg);
        } while ((loop < 10) && (Math.abs(fnew - 1) > 0.00000001));
        if (loop == 10)
        {
            System.out.println("too many loops in refine_splice_t1");
            return Double.NaN;
        }
        if (xnew < i) return Double.NaN;
        if (xnew > i + 1) return Double.NaN;
        return xnew - (int) xnew;
    }

    private static double refine_splice_t1_RHS(int i)
    {
        // improve the initial estimate of variable spliced using regula falsi
        // i = splicei, the index of the original estimate of the t1 value
        double xold = i;
        double xnew = i + 1;
        double fold = t2[i];
        double fnew = t2[i + 1];
        int seg;
        int loop = 0;

        System.out.println("RHS xold = ," + xold + ", " + fold + ", " + fnew);
        do
        {
            loop++;
            seg = 1;
            if (fold < 1)
                xold = xold + (xnew - xold)/(1 + beta1*(fnew - 1)/(1 - fold));  // if we bracket the splice
            else
                xold = xold + (xnew - xold)*(1 - fold)/(fnew - fold);           // both points above the splice
            fold = solve_quintic_for_t2(xold, t2[i], seg);
            if (fold < 1)
            {
                seg--;
                System.out.println("recheck = ," + xold + ", " + fold + ", " + fnew + ", " + seg);
                fold = solve_quintic_for_t2(xold, t2[i], seg);
                if (fold > 1)
                {
                    System.out.println("error: reversal of xold = ," + xold + ", " + fold + ", " + fnew);
                    return Double.NaN;
                }
            }
            System.out.println("xold    = ," + xold + ", " + fold + ", " + fnew + ", " + seg);
        } while ((loop < 10) && (Math.abs(fold - 1) > 0.00000001));
        if (loop == 10)
        {
            System.out.println("too many loops in refine_splice_t1");
            return Double.NaN;
        }
        if (xold < i) return Double.NaN;
        if (xold > i + 1) return Double.NaN;
        return xold - (int) xold;
    }

    private static double calc_t2dd(int i, int seg, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        if (i > N)                                                      // calculate exact position of splice
            t1 = t1_start + (splicei + spliced)*(t1_end - t1_start)/N;  // interpolate t1, not g
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double denom = t2_vs_t1.dfn(Bezx[seg], t2)*t2_vs_t1.dfn(Bezx[seg], t2) + (t2_vs_t1.fn(Bezx[seg], t2) - X)*t2_vs_t1.d2fn(Bezx[seg], t2) + t2_vs_t1.dfn(Bezy[seg], t2)*t2_vs_t1.dfn(Bezy[seg], t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*t2_vs_t1.d2fn(Bezy[seg], t2);
        //if (type.equals("x2")) System.out.println(" , , , , , , , , " + denom);
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
        else if (type.equals("beta1"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdbeta1(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydbeta1(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudbeta1(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudbeta1(seg, t2);
        else if (type.equals("c"))
            numer = calc_df_gxdc(seg, t1, t2)*t2_vs_t1.dfn(Bezx[seg], t2) + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudc(seg, t2)
                  + calc_df_gydc(seg, t1, t2)*t2_vs_t1.dfn(Bezy[seg], t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudc(seg, t2);
        return -numer/denom;
    }

    private static double calc_dfxdx2(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[3] + N33(t2)[2]*(1 + beta1)/2/beta1;
        return N33(t2)[0] + N33(t2)[1]*(1 + beta1)/2;
    }

    private static double calc_d2fxdudx2(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[3] + dN33(t2)[2]*(1 + beta1)/2/beta1;
        return dN33(t2)[0] + dN33(t2)[1]*(1 + beta1)/2;
    }

    private static double calc_dfydy2(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[3] + N33(t2)[2]*(1 + beta1)/2/beta1;
        return N33(t2)[0] + N33(t2)[1]*(1 + beta1)/2;
    }

    private static double calc_d2fydudy2(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[3] + dN33(t2)[2]*(1 + beta1)/2/beta1;
        return dN33(t2)[0] + dN33(t2)[1]*(1 + beta1)/2;
    }

    private static double calc_dfxdd1(int seg, double t2)
    {
        if (seg == 0) return (N33(t2)[1] + N33(t2)[2]*beta1/2/(1 + beta1))*Math.cos(theta_start);
        return -N33(t2)[1]*beta1*beta1/2/(1 + beta1)*Math.cos(theta_start);
    }

    private static double calc_dfydd1(int seg, double t2)
    {
        if (seg == 0) return (N33(t2)[1] + N33(t2)[2]*beta1/2/(1 + beta1))*Math.sin(theta_start);
        return -N33(t2)[1]*beta1*beta1/2/(1 + beta1)*Math.sin(theta_start);
    }

    private static double calc_d2fxdudd1(int seg, double t2)
    {
        if (seg == 0) return (dN33(t2)[1] + dN33(t2)[2]*beta1/2/(1 + beta1))*Math.cos(theta_start);
        return -dN33(t2)[1]*beta1*beta1/2/(1 + beta1)*Math.cos(theta_start);
    }

    private static double calc_d2fydudd1(int seg, double t2)
    {
        if (seg == 0) return (dN33(t2)[1] + dN33(t2)[2]*beta1/2/(1 + beta1))*Math.sin(theta_start);
        return -dN33(t2)[1]*beta1*beta1/2/(1 + beta1)*Math.sin(theta_start);
    }

    private static double calc_dfxdd2(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[2]/2/beta1/(1 + beta1)*Math.cos(theta_end);
        return -(N33(t2)[2] + N33(t2)[1]/2/(1 + beta1))*Math.cos(theta_end);
    }

    private static double calc_dfydd2(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[2]/2/beta1/(1 + beta1)*Math.sin(theta_end);
        return -(N33(t2)[2] + N33(t2)[1]/2/(1 + beta1))*Math.sin(theta_end);
    }

    private static double calc_d2fxdudd2(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[2]/2/beta1/(1 + beta1)*Math.cos(theta_end);
        return -(dN33(t2)[2] + dN33(t2)[1]/2/(1 + beta1))*Math.cos(theta_end);
    }

    private static double calc_d2fydudd2(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[2]/2/beta1/(1 + beta1)*Math.sin(theta_end);
        return -(dN33(t2)[2] + dN33(t2)[1]/2/(1 + beta1))*Math.sin(theta_end);
    }

    private static double calc_dfxdbeta1(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[2]*(-beta1*(Bezx[0][3] - Bezx[0][1]) + (1 + 2*beta1)*(Bezx[0][3] - Bezx[0][2]))/beta1/(1 + beta1);
        return N33(t2)[1]*(beta1*(Bezx[0][3] - Bezx[0][1]) - (Bezx[1][1] - Bezx[1][0]))/(1 + beta1);
    }

    private static double calc_dfydbeta1(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[2]*(-beta1*(Bezy[0][3] - Bezy[0][1]) + (1 + 2*beta1)*(Bezy[0][3] - Bezy[0][2]))/beta1/(1 + beta1);
        return N33(t2)[1]*(beta1*(Bezy[0][3] - Bezy[0][1]) - (Bezy[1][1] - Bezy[1][0]))/(1 + beta1);
    }

    private static double calc_d2fxdudbeta1(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[2]*(-beta1*(Bezx[0][3] - Bezx[0][1]) + (1 + 2*beta1)*(Bezx[0][3] - Bezx[0][2]))/beta1/(1 + beta1);
        return dN33(t2)[1]*(beta1*(Bezx[0][3] - Bezx[0][1]) - (Bezx[1][1] - Bezx[1][0]))/(1 + beta1);
    }

    private static double calc_d2fydudbeta1(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[2]*(-beta1*(Bezy[0][3] - Bezy[0][1]) + (1 + 2*beta1)*(Bezy[0][3] - Bezy[0][2]))/beta1/(1 + beta1);
        return dN33(t2)[1]*(beta1*(Bezy[0][3] - Bezy[0][1]) - (Bezy[1][1] - Bezy[1][0]))/(1 + beta1);
    }

//  calculate d2f/ddi/ddj (at least one index must be beta1)

    private static double calc_d2fxdd1dbeta1(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[2]/2/(1 + beta1)/(1 + beta1)*Math.cos(theta_start);
        return -N33(t2)[1]*beta1*(2 + beta1)/2/(1 + beta1)/(1 + beta1)*Math.cos(theta_start);
    }

    private static double calc_d2fydd1dbeta1(int seg, double t2)
    {
        if (seg == 0) return N33(t2)[2]/2/(1 + beta1)/(1 + beta1)*Math.sin(theta_start);
        return -N33(t2)[1]*beta1*(2 + beta1)/2/(1 + beta1)/(1 + beta1)*Math.sin(theta_start);
    }

    private static double calc_d2fxdd2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -N33(t2)[2]*(1 + 2*beta1)/2/beta1/beta1/(1 + beta1)/(1 + beta1)*Math.cos(theta_end);
        return N33(t2)[1]/2/(1 + beta1)/(1 + beta1)*Math.cos(theta_end);
    }

    private static double calc_d2fydd2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -N33(t2)[2]*(1 + 2*beta1)/2/beta1/beta1/(1 + beta1)/(1 + beta1)*Math.sin(theta_end);
        return N33(t2)[1]/2/(1 + beta1)/(1 + beta1)*Math.sin(theta_end);
    }

    private static double calc_d2fxdx2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -N33(t2)[2]/2/beta1/beta1;
        return N33(t2)[1]/2;
    }

    private static double calc_d2fydy2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -N33(t2)[2]/2/beta1/beta1;
        return N33(t2)[1]/2;
    }

    private static double calc_d2fxdbeta1dbeta1(int seg, double t2)
    {
        if (seg == 0)
        {
            double dQ2dbeta1 = (-beta1*(Bezx[0][3] - Bezx[0][1]) + (1 + 2*beta1)*(Bezx[0][3] - Bezx[0][2]))/beta1/(1 + beta1);
            return N33(t2)[2]*(Bezx[0][1] - 2*Bezx[0][2] + Bezx[0][3]
                            - 2*(1 + 2*beta1)*dQ2dbeta1)/beta1/(1 + beta1);
        }
        double dQ5dbeta1 = (beta1*(Bezx[0][3] - Bezx[0][1]) - (Bezx[1][1] - Bezx[1][0]))/(1 + beta1);
        return N33(t2)[1]*(Bezx[0][3] - Bezx[0][1] - 2*dQ5dbeta1)/(1 + beta1);
    }

    private static double calc_d2fydbeta1dbeta1(int seg, double t2)
    {
        if (seg == 0)
        {
            double dQ2dbeta1 = (-beta1*(Bezy[0][3] - Bezy[0][1]) + (1 + 2*beta1)*(Bezy[0][3] - Bezy[0][2]))/beta1/(1 + beta1);
            return N33(t2)[2]*(Bezy[0][1] - 2*Bezy[0][2] + Bezy[0][3]
                            - 2*(1 + 2*beta1)*dQ2dbeta1)/beta1/(1 + beta1);
        }
        double dQ5dbeta1 = (beta1*(Bezy[0][3] - Bezy[0][1]) - (Bezy[1][1] - Bezy[1][0]))/(1 + beta1);
        return N33(t2)[1]*(Bezy[0][3] - Bezy[0][1] - 2*dQ5dbeta1)/(1 + beta1);
    }

//  calculate Augmented matrix

    private static double calc_df_gxdc(int seg, double t1, double t2)
    {
        if (seg == 0)
            return fitted.getdxdc(t1_start)*(N33(t2)[0] + N33(t2)[1] + N33(t2)[2]*beta1/2/(1 + beta1)) - fitted.getdxdc(t1_end)*N33(t2)[2]/2/beta1/(1 + beta1) - fitted.getdxdc(t1);
        return fitted.getdxdc(t1_end)*(N33(t2)[3] + N33(t2)[2] + N33(t2)[1]/2/(1 + beta1)) - fitted.getdxdc(t1_start)*N33(t2)[1]*beta1*beta1/2/(1 + beta1) - fitted.getdxdc(t1);
    }

    private static double calc_df_gydc(int seg, double t1, double t2)
    {
        if (seg == 0)
            return fitted.getdydc(t1_start)*(N33(t2)[0] + N33(t2)[1] + N33(t2)[2]*beta1/2/(1 + beta1)) - fitted.getdydc(t1_end)*N33(t2)[2]/2/beta1/(1 + beta1) - fitted.getdydc(t1);
        return fitted.getdydc(t1_end)*(N33(t2)[3] + N33(t2)[2] + N33(t2)[1]/2/(1 + beta1)) - fitted.getdydc(t1_start)*N33(t2)[1]*beta1*beta1/2/(1 + beta1) - fitted.getdydc(t1);
    }

    private static double calc_d2fxdudc(int seg, double t2)
    {
        if (seg == 0)
            return fitted.getdxdc(t1_start)*(dN33(t2)[0] + dN33(t2)[1] + dN33(t2)[2]*beta1/2/(1 + beta1)) - fitted.getdxdc(t1_end)*dN33(t2)[2]/2/beta1/(1 + beta1);
        return fitted.getdxdc(t1_end)*(dN33(t2)[3] + dN33(t2)[2] + dN33(t2)[1]/2/(1 + beta1)) - fitted.getdxdc(t1_start)*dN33(t2)[1]*beta1*beta1/2/(1 + beta1);
    }

    private static double calc_d2fydudc(int seg, double t2)
    {
        if (seg == 0)
            return fitted.getdydc(t1_start)*(dN33(t2)[0] + dN33(t2)[1] + dN33(t2)[2]*beta1/2/(1 + beta1)) - fitted.getdydc(t1_end)*dN33(t2)[2]/2/beta1/(1 + beta1);
        return fitted.getdydc(t1_end)*(dN33(t2)[3] + dN33(t2)[2] + dN33(t2)[1]/2/(1 + beta1)) - fitted.getdydc(t1_start)*dN33(t2)[1]*beta1*beta1/2/(1 + beta1);
    }

    private static double calc_d2fxdcdbeta1(int seg, double t2)
    {
        if (seg == 0)
            return fitted.getdxdc(t1_start)*N33(t2)[2]/2/(1 + beta1)/(1 + beta1) + fitted.getdxdc(t1_end)*N33(t2)[2]*(1 + 2*beta1)/2/beta1/(1 + beta1)/beta1/(1 + beta1);
        return -fitted.getdxdc(t1_end)*N33(t2)[1]/2/(1 + beta1)/(1 + beta1) - fitted.getdxdc(t1_start)*N33(t2)[1]*beta1*(2 + beta1)/2/(1 + beta1)/(1 + beta1);
    }

    private static double calc_d2fydcdbeta1(int seg, double t2)
    {
        if (seg == 0)
            return fitted.getdydc(t1_start)*N33(t2)[2]/2/(1 + beta1)/(1 + beta1) + fitted.getdydc(t1_end)*N33(t2)[2]*(1 + 2*beta1)/2/beta1/(1 + beta1)/beta1/(1 + beta1);
        return -fitted.getdydc(t1_end)*N33(t2)[1]/2/(1 + beta1)/(1 + beta1) - fitted.getdydc(t1_start)*N33(t2)[1]*beta1*(2 + beta1)/2/(1 + beta1)/(1 + beta1);
    }

//  calculate d2u/ddi/ddj

    private static double calc_d3fxdududd1(int seg, double t2)
    {
        if (seg == 0) return (d2N33(t2)[1] + d2N33(t2)[2]*beta1/2/(1 + beta1))*Math.cos(theta_start);
        return -d2N33(t2)[1]*beta1*beta1/2/(1 + beta1)*Math.cos(theta_start);
    }

    private static double calc_d3fydududd1(int seg, double t2)
    {
        if (seg == 0) return (d2N33(t2)[1] + d2N33(t2)[2]*beta1/2/(1 + beta1))*Math.sin(theta_start);
        return -d2N33(t2)[1]*beta1*beta1/2/(1 + beta1)*Math.sin(theta_start);
    }

    private static double calc_d3fxdududd2(int seg, double t2)
    {
        if (seg == 0) return d2N33(t2)[2]/2/beta1/(1 + beta1)*Math.cos(theta_end);
        return -(d2N33(t2)[2] + d2N33(t2)[1]/2/(1 + beta1))*Math.cos(theta_end);
    }

    private static double calc_d3fydududd2(int seg, double t2)
    {
        if (seg == 0) return d2N33(t2)[2]/2/beta1/(1 + beta1)*Math.sin(theta_end);
        return -(d2N33(t2)[2] + d2N33(t2)[1]/2/(1 + beta1))*Math.sin(theta_end);
    }

    private static double calc_d3fxdududx2(int seg, double t2)
    {
        if (seg == 0) return d2N33(t2)[3] + d2N33(t2)[2]*(1 + beta1)/2/beta1;
        return d2N33(t2)[0] + d2N33(t2)[1]*(1 + beta1)/2;
    }

    private static double calc_d3fydududy2(int seg, double t2)
    {
        if (seg == 0) return d2N33(t2)[3] + d2N33(t2)[2]*(1 + beta1)/2/beta1;
        return d2N33(t2)[0] + d2N33(t2)[1]*(1 + beta1)/2;
    }

    private static double calc_d3fxdududbeta1(int seg, double t2)
    {
        if (seg == 0) return d2N33(t2)[2]*(-beta1*(Bezx[0][3] - Bezx[0][1]) + (1 + 2*beta1)*(Bezx[0][3] - Bezx[0][2]))/beta1/(1 + beta1);
        return d2N33(t2)[1]*(beta1*(Bezx[0][3] - Bezx[0][1]) - (Bezx[1][1] - Bezx[1][0]))/(1 + beta1);
    }

    private static double calc_d3fydududbeta1(int seg, double t2)
    {
        if (seg == 0) return d2N33(t2)[2]*(-beta1*(Bezy[0][3] - Bezy[0][1]) + (1 + 2*beta1)*(Bezy[0][3] - Bezy[0][2]))/beta1/(1 + beta1);
        return d2N33(t2)[1]*(beta1*(Bezy[0][3] - Bezy[0][1]) - (Bezy[1][1] - Bezy[1][0]))/(1 + beta1);
    }

    private static double calc_d3fxdudd1dbeta1(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[2]/2/(1 + beta1)/(1 + beta1)*Math.cos(theta_start);
        return -dN33(t2)[1]*beta1*(2 + beta1)/2/(1 + beta1)/(1 + beta1)*Math.cos(theta_start);
    }

    private static double calc_d3fydudd1dbeta1(int seg, double t2)
    {
        if (seg == 0) return dN33(t2)[2]/2/(1 + beta1)/(1 + beta1)*Math.sin(theta_start);
        return -dN33(t2)[1]*beta1*(2 + beta1)/2/(1 + beta1)/(1 + beta1)*Math.sin(theta_start);
    }

    private static double calc_d3fxdudd2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -dN33(t2)[2]*(1 + 2*beta1)/2/beta1/beta1/(1 + beta1)/(1 + beta1)*Math.cos(theta_end);
        return dN33(t2)[1]/2/(1 + beta1)/(1 + beta1)*Math.cos(theta_end);
    }

    private static double calc_d3fydudd2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -dN33(t2)[2]*(1 + 2*beta1)/2/beta1/beta1/(1 + beta1)/(1 + beta1)*Math.sin(theta_end);
        return dN33(t2)[1]/2/(1 + beta1)/(1 + beta1)*Math.sin(theta_end);
    }

    private static double calc_d3fxdudx2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -dN33(t2)[2]/2/beta1/beta1;
        return dN33(t2)[1]/2;
    }

    private static double calc_d3fydudy2dbeta1(int seg, double t2)
    {
        if (seg == 0) return -dN33(t2)[2]/2/beta1/beta1;
        return dN33(t2)[1]/2;
    }

    private static double calc_d3fxdudbeta1dbeta1(int seg, double t2)
    {
        if (seg == 0)
        {
            double dQ2dbeta1 = (-beta1*(Bezx[0][3] - Bezx[0][1]) + (1 + 2*beta1)*(Bezx[0][3] - Bezx[0][2]))/beta1/(1 + beta1);
            return dN33(t2)[2]*(Bezx[0][1] - 2*Bezx[0][2] + Bezx[0][3]
                            - 2*(1 + 2*beta1)*dQ2dbeta1)/beta1/(1 + beta1);
        }
        double dQ5dbeta1 = (beta1*(Bezx[0][3] - Bezx[0][1]) - (Bezx[1][1] - Bezx[1][0]))/(1 + beta1);
        return dN33(t2)[1]*(Bezx[0][3] - Bezx[0][1] - 2*dQ5dbeta1)/(1 + beta1);
    }

    private static double calc_d3fydudbeta1dbeta1(int seg, double t2)
    {
        if (seg == 0)
        {
            double dQ2dbeta1 = (-beta1*(Bezy[0][3] - Bezy[0][1]) + (1 + 2*beta1)*(Bezy[0][3] - Bezy[0][2]))/beta1/(1 + beta1);
            return dN33(t2)[2]*(Bezy[0][1] - 2*Bezy[0][2] + Bezy[0][3]
                            - 2*(1 + 2*beta1)*dQ2dbeta1)/beta1/(1 + beta1);
        }
        double dQ5dbeta1 = (beta1*(Bezy[0][3] - Bezy[0][1]) - (Bezy[1][1] - Bezy[1][0]))/(1 + beta1);
        return dN33(t2)[1]*(Bezy[0][3] - Bezy[0][1] - 2*dQ5dbeta1)/(1 + beta1);
    }

//  calculate d2f/du/du

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

    private static double[] d2N33(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:  dN33 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:  dN33 too large u value = " + u);
        return new double[] { 6*(1 - u),
                             -6*(2 - 3*u),
                              6*(1 - 3*u),
                              6*u};
    }
 }
