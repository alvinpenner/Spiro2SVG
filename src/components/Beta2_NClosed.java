
package components;

import java.io.*;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a closed N-point Beta2-Spline (P0 -> PN-1) to it, using parameter 0 < t2 < N, where N = 8.
// express the beta2 spline as N Beziers with no constraints on the N Bezier endpoints
// and no constraints on the N beta2 values.
// linearize the equations wrt Pi and betai and solve a 3Nx3N system of equations.
// t2 must be chosen to minimize the distance to g(t1)
// define a = (x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7, beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7)
// see Spiro2SVG Book9, Sept 2019, page 8

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\Beta2_NClosed.java

public class Beta2_NClosed
{
    public static final double t1_start = 0;
    public static final double t1_end = 2*Math.PI;
    public static final int N = 800;                // normally 100*number of segments
    private static double[] ax;                     // x coord of Bezier endpoints
    private static double[] ay;                     // y coord of Bezier endpoints
    private static double[] abeta;                  // beta2 values at Bezier endpoints
    public static double[] Bezx = new double[24];   // cubic Bezier, 8 splices x 3 points, x component
    public static double[] Bezy = new double[24];   // cubic Bezier, 8 splices x 3 points, y component
    //private static epiTrochoidFxn fitted;
    //private static HippopedeFxn fitted;
    private static SuperEllipse fitted;
    private static double[] t2 = new double[N+1];
    private static double[] denom = new double[N+1];                    // = E(u)
    private static double[][] t2ddx = new double[8][N+1];               // partial of u wrt {xj}
    private static double[][] t2ddy = new double[8][N+1];               // partial of u wrt {yj}
    private static double[][] t2ddbeta = new double[8][N+1];            // partial of u wrt {betaj}
    private static double[][] dddx = new double[8][8];                  // dd[i]/dx[j] (see Book 9, p. 54)
    private static double[][][] d2ddxdbeta = new double[8][8][8];       // d2d[i]/dx[j]/dbeta[k]
    private static double[][] ddxdbeta = new double[8][8];              // ddx[i]dbeta[j]
    private static double[][] ddydbeta = new double[8][8];              // ddy[i]dbeta[j]
    private static double[][][] d2dxdbetaibetaj = new double[8][8][8];  // d2dx[i]/dbeta[j]/dbeta[k]
    private static double[][][] d2dydbetaibetaj = new double[8][8][8];  // d2dy[i]/dbeta[j]/dbeta[k]
    private static double[][] Jac;                                      // defined here only so we can share it with 'scan_steepest'
    private static double[] dFda;
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.00000001;      // loosened by 00

    public static void main (String[] args)
    {
        //fitted = new epiTrochoidFxn(-5.5);
        //fitted = new HippopedeFxn(0.5);
        //fitted = new SuperEllipse(0.75);
        //scan_fourfold_symmetry();
        //scan_steepest();
        get_Bezier_endpoints_from_Beta_Spline();
        System.out.println("Beta2_Spline8 iterate_at_x_y = " + iterate_at_x_y(null) + "\n");       // normally include this line
        //System.out.println("Beta2_Spline8 solve_at_x_y = " + solve_at_x_y(true) + "\n");
        //write_Beta2_SplineN_data();                              // generate eigenvalue data
        //for (int i = 0; i < 5; i++)
        //    read_one_line(-0.0002);                                       // initiallize using previous data
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //gen_spline_points();
        //for (int i = 0; i <= 100; i++)
        //    System.out.println(i + ", " + Bi3(i/100.0)[3] + ", " + dBi3(i/100.0)[3] + ", " + d2Bi3(i/100.0)[3] + ", " + d3Bi3(i/100.0)[3]);
    }

    private static double iterate_at_x_y(PrintWriter out)
    {
        // calculate a new estimate of (xi, yi, betai) by setting dF = 0
        // setup 3N-variable Newton-Raphson iteration

        final double gain = 1;                                  // factor to reduce gain
        final int MAXLOOP = 500;
        double[] f_gx = new double[N];
        double[] f_gy = new double[N];
        double[] dfxdu = new double[N];
        double[] dfydu = new double[N];

        double[] trap_in = new double[N];
        Jac = new double[3*ax.length][3*ax.length];
        dFda = new double[3*ax.length];
        double[] d2Fdadc = new double[3*ax.length];             // augmented matrix
        //double dFdc;
        double d2Fdcdc;                                         // augmented matrix
        double[][] Augment = new double[3*ax.length + 1][3*ax.length + 1];
        double[] dela;                                          // (-Δxi, -Δyi)
        int loop = 0;
        int i, j, k;
        boolean outside;

        if (out != null && false)   // temporary code to dump d0 - d3 to a file (disabled)
        {
            double retVal = solve_at_x_y(false);
            out.print("d1d2 = ," + fitted.getc());
            for (i = 0; i < 2; i++)
                out.print(", " + Math.sqrt((Bezx[3*i + 1] - Bezx[3*i])*(Bezx[3*i + 1] - Bezx[3*i]) + (Bezy[3*i + 1] - Bezy[3*i])*(Bezy[3*i + 1] - Bezy[3*i])));
            out.println(", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2);
            return Double.NaN;
        }
        do
        {
            loop++;
            if (Double.isNaN(solve_at_x_y(false)))              // initiallize at (xi, yi)
            {
                System.out.println("fail at " + print_coord(ax) + print_coord(ay) + print_coord(abeta));
                return Double.NaN;
            }

            for (i = 0; i < N; i++)
            {
                f_gx[i] = multvv(Bezx, (int) t2[i], Bi3(t2[i] % 1)) - fitted.getx(t1_start + i*(t1_end - t1_start)/N);
                f_gy[i] = multvv(Bezy, (int) t2[i], Bi3(t2[i] % 1)) - fitted.gety(t1_start + i*(t1_end - t1_start)/N);
                dfxdu[i] = multvv(Bezx, (int) t2[i], dBi3(t2[i] % 1));
                dfydu[i] = multvv(Bezy, (int) t2[i], dBi3(t2[i] % 1));
                //System.out.println(i + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + denom[i]);
                //double[] testf = new double[3*ax.length];               // test for oval branch
                //double testtheta;
                //for (j = 0; j < ax.length; j++)                         // use xi, yi, betai as independent variables
                //{                                                       // This is the wrong way to do it.
                //    testf[j] = calc_dfdai(j, t2[i])*f_gx[i];
                //    testf[j + ax.length] = calc_dfdai(j, t2[i])*f_gy[i];
                //    testf[j + 2*ax.length] = calc_dfdbetai(ddxdbeta, j, t2[i])*f_gx[i]
                //                           + calc_dfdbetai(ddydbeta, j, t2[i])*f_gy[i];
                //}
/*                                                                     // temporary code for testing for submanifold
                for (j = 0; j < ax.length; j++)                         // use xi, yi, di as independent variables
                {                                                       // This is the right way to do it
                    testtheta = Math.PI*(1 + j)/4;                      // assuming angular symmetry (-45, 0)
                    //testf[j] = calc_dfdai_at_di(j, t2[i])*f_gx[i];        // estimate betai
                    //testf[j + ax.length] = calc_dfdai_at_di(j, t2[i])*f_gy[i];
                    //testf[j + 2*ax.length] = calc_dfddi_at_betai(j, t2[i])*f_gx[i]*Math.cos(testtheta)
                    //                       + calc_dfddi_at_betai(j, t2[i])*f_gy[i]*Math.sin(testtheta);
                    testf[j] = calc_d2fdudai_at_di(j, t2[i])*f_gx[i];       // estimete gammai
                    testf[j + ax.length] = calc_d2fdudai_at_di(j, t2[i])*f_gy[i];
                    testf[j + 2*ax.length] = calc_d2fduddi_at_betai(j, t2[i])*f_gx[i]*Math.cos(testtheta)
                                           + calc_d2fduddi_at_betai(j, t2[i])*f_gy[i]*Math.sin(testtheta);
                    //testf[j] = calc_dfdai_at_di(j, t2[i]);      // temporary code only
                    //testf[j + ax.length] = calc_dfdai_at_di(j, t2[i]);
                    //testf[j + 2*ax.length] = calc_dfddi_at_betai(j, t2[i]);
                }
                System.out.println("testf = ," + i + print_coord(testf));
*/
                //System.out.println("testf = ," + (t2[i]%8) + print_coord(testf));
                //System.out.print(t2[i] % 8);
                //for (j = 0; j < ax.length; j++)
                //    System.out.print(", " + calc_dfdai(j, t2[i]));
                //for (j = 0; j < ax.length; j++)
                //    System.out.print(", " + calc_dfdbetai(ddxdbeta, j, t2[i]));
                //for (j = 0; j < ax.length; j++)
                //    System.out.print(", " + calc_dfdbetai(ddydbeta, j, t2[i]));
                //System.out.println();
            }

            // calc dFda[j] at current (xi, yi, betai)

            for (i = 0; i < 3*ax.length; i++)
            {
                for (k = 0; k < N; k++)
                    if (i < ax.length)                  // ax variable
                        trap_in[k] = f_gx[k]*calc_dfdai(i, t2[k])                       // original code
                                   + (f_gx[k]*dfxdu[k] + f_gy[k]*dfydu[k])*t2ddx[i][k];
                        //trap_in[k] = f_gx[k]*calc_dfdai(i, t2[k]);                    // new code
                    else if (i < 2*ax.length)           // ay variable
                        trap_in[k] = f_gy[k]*calc_dfdai(i - ax.length, t2[k])           // original code
                                   + (f_gx[k]*dfxdu[k] + f_gy[k]*dfydu[k])*t2ddy[i - ax.length][k];
                        //trap_in[k] = f_gy[k]*calc_dfdai(i - ax.length, t2[k]);        // new code
                    else                                // abeta variable
                        trap_in[k] = f_gx[k]*calc_dfdbetai(ddxdbeta, i - 2*ax.length, t2[k])    // original code
                                   + f_gy[k]*calc_dfdbetai(ddydbeta, i - 2*ax.length, t2[k])
                                   + (f_gx[k]*dfxdu[k] + f_gy[k]*dfydu[k])*t2ddbeta[i - 2*ax.length][k];
                dFda[i] = integrate(trap_in);
                //System.out.println("dFda ," + i + ", " + dFda[i]);
                for (k = 0; k < N; k++)                                     // fix fix incomplete Hessian ??????
                    if (i < ax.length)                  // ax variable
                        trap_in[k] = -calc_dfdai(i, t2[k])*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)
                                   - denom[k]*t2ddx[i][k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
                    else if (i < 2*ax.length)           // ay variable
                        trap_in[k] = -calc_dfdai(i - ax.length, t2[k])*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)
                                   - denom[k]*t2ddy[i - ax.length][k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
                    else                                // abeta variable
                        trap_in[k] = -calc_dfdbetai(ddxdbeta, i - 2*ax.length, t2[k])*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)
                                   -  calc_dfdbetai(ddydbeta, i - 2*ax.length, t2[k])*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)
                                   - denom[k]*t2ddbeta[i - 2*ax.length][k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
                d2Fdadc[i] = integrate(trap_in);                       // augmented matrix
            }

            // calc d2Fda[i]da[j] (symmetric Jacobean matrix)

            for (i = 0; i < 3*ax.length; i++)
                for (j = 0; j < 3*ax.length; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k < N; k++)
                    {
                        if (i < ax.length)
                            if (j < ax.length)
                                trap_in[k] = calc_dfdai(i, t2[k])*calc_dfdai(j, t2[k])
                                           - denom[k]*t2ddx[i][k]*t2ddx[j][k];
                            else if (j < 2*ax.length)
                                trap_in[k] = -denom[k]*t2ddx[i][k]*t2ddy[j - ax.length][k];
                            else
                                trap_in[k] = calc_dfdai(i, t2[k])*calc_dfdbetai(ddxdbeta, j - 2*ax.length, t2[k])
                                           - denom[k]*t2ddx[i][k]*t2ddbeta[j - 2*ax.length][k]
                                           + f_gx[k]*calc_d2fdaidbetaj(i, j - 2*ax.length, t2[k]);
                        else if (i < 2*ax.length)
                            if (j < ax.length)
                                trap_in[k] = -denom[k]*t2ddy[i - ax.length][k]*t2ddx[j][k];
                            else if (j < 2*ax.length)
                                trap_in[k] = calc_dfdai(i - ax.length, t2[k])*calc_dfdai(j - ax.length, t2[k])
                                           - denom[k]*t2ddy[i - ax.length][k]*t2ddy[j - ax.length][k];
                            else
                                trap_in[k] = calc_dfdai(i - ax.length, t2[k])*calc_dfdbetai(ddydbeta, j - 2*ax.length, t2[k])
                                           - denom[k]*t2ddy[i - ax.length][k]*t2ddbeta[j - 2*ax.length][k]
                                           + f_gy[k]*calc_d2fdaidbetaj(i - ax.length, j - 2*ax.length, t2[k]);
                        else
                            if (j < ax.length)
                                trap_in[k] = calc_dfdbetai(ddxdbeta, i - 2*ax.length, t2[k])*calc_dfdai(j, t2[k])
                                           - denom[k]*t2ddbeta[i - 2*ax.length][k]*t2ddx[j][k]
                                           + f_gx[k]*calc_d2fdaidbetaj(j, i - 2*ax.length, t2[k]);
                            else if (j < 2*ax.length)
                                trap_in[k] = calc_dfdbetai(ddydbeta, i - 2*ax.length, t2[k])*calc_dfdai(j - ax.length, t2[k])
                                           - denom[k]*t2ddbeta[i - 2*ax.length][k]*t2ddy[j - ax.length][k]
                                           + f_gy[k]*calc_d2fdaidbetaj(j - ax.length, i - 2*ax.length, t2[k]);
                            else
                                trap_in[k] = calc_dfdbetai(ddxdbeta, i - 2*ax.length, t2[k])*calc_dfdbetai(ddxdbeta, j - 2*ax.length, t2[k])
                                           + calc_dfdbetai(ddydbeta, i - 2*ax.length, t2[k])*calc_dfdbetai(ddydbeta, j - 2*ax.length, t2[k])
                                           - denom[k]*t2ddbeta[i - 2*ax.length][k]*t2ddbeta[j - 2*ax.length][k]
                                           + f_gx[k]*calc_d2fdbetaibetaj(d2dxdbetaibetaj, i - 2*ax.length, j - 2*ax.length, t2[k])
                                           + f_gy[k]*calc_d2fdbetaibetaj(d2dydbetaibetaj, i - 2*ax.length, j - 2*ax.length, t2[k]);
                        //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = integrate(trap_in);
                }
            //if (true) return Double.NaN;                        // use only for scan_steepest()

            // calculate determinant of augmented matrix

            //for (i = 0; i < 3*ax.length; i++)
            //{
            //    for (j = 0; j < 3*ax.length; j++)
            //        System.out.print(" ," + Jac[i][j]);
            //    System.out.println();
            //}
            for (k = 0; k < N; k++)
                trap_in[k] = -f_gx[k]*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N) - f_gy[k]*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N);
            //dFdc = integrate(trap_in);
            //int checkj = 21;
            //System.out.println("check Aug ," + fitted.getc() + ", " + checkj + ", " + dFdc + ", " + dFda[checkj] + ", " + d2Fdadc[checkj]);
            for (k = 0; k < N; k++)
                trap_in[k] = fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)
                           + fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)
                           - f_gx[k]*fitted.getd2xdc2(t1_start + k*(t1_end - t1_start)/N)
                           - f_gy[k]*fitted.getd2ydc2(t1_start + k*(t1_end - t1_start)/N)
                           - denom[k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c")*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
            d2Fdcdc = integrate(trap_in);
            //System.out.println("check d2Fdc2 , " + fitted.getc() + ", " + dFdc + ", " + d2Fdcdc);
            for (i = 0; i < 3*ax.length; i++)
                for (j = 0; j < 3*ax.length; j++)
                    Augment[i][j] = Jac[i][j];
            for (i = 0; i < 3*ax.length; i++)
            {
                Augment[i][3*ax.length] = d2Fdadc[i];
                Augment[3*ax.length][i] = Augment[i][3*ax.length];
            }
            Augment[3*ax.length][3*ax.length] = d2Fdcdc;
            //System.out.println("dFdc = " + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);

            if (loop == 1)  // perform one iteration of steepest descent as a special initialization
            {
                //BSpline5.dump_Jac(Jac);
                System.out.println("steepest descent at loop = " + loop);
                dela = new double[dFda.length];
                double steepnum = BSpline5.multvv(dFda, dFda);
                double steepdenom = BSpline5.multvv(dFda, BSpline5.multmv(Jac, dFda));
                for (j = 0; j < dFda.length; j++)
                    dela[j] = gain*dFda[j]*steepnum/steepdenom;
            }
            else
                dela = BSpline5.gaussj(Jac, dFda);          // this is actually the negative of Δa
            for (i = 0; i < ax.length; i++)
            {
                ax[i] -= dela[i]/gain;             // gain is just a fudge factor to 'improve' convergence
                ay[i] -= dela[i + ax.length]/gain;
                abeta[i] -= dela[i + 2*ax.length]/gain;
            }

            //Jacdet = BSpline5.detm(Jac);
            System.out.println("dFda = " + print_coord(dFda) + ", " + Jacdet);
            System.out.println("dela = " + print_coord(dela));
            //System.out.println("d2Fdadc    = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(d2Fdadc));
            outside = false;
            for (k = 0; k < dela.length; k++)
                if (Math.abs(dela[k]) > TOL)
                    outside = true;
            //outside = false;                              // force convergence
        } while ((loop < MAXLOOP) && outside);
        if (loop < MAXLOOP)
        {
            System.out.println("\nt2 = [" + t2[0] + ", " + t2[N] + "]");
            //System.out.println("oval det =  ," + fitted.getc() + ", " + calc_oval_det());
            if (Math.abs(t2[N] - t2[0] - ax.length) > TOL
            ||  t2[0] < -TOL)
            {
                System.out.println("bad converged t2[i]: ABORT");
                return Double.NaN;
            }
            System.out.println("\n__converged in " + loop + " at new xi yi betai = , , , , , , " + print_coord(ax) + print_coord(ay) + print_coord(abeta));
            //Jacdet = BSpline5.detm(Jac);
            //System.out.println("dFdc = " + fitted.getc() + print_radial() + print_coord(abeta));     // this is just a header for dump_Jac file data
            //BSpline5.dump_Jac(Jac);
            //BSpline5.dump_Jac(Augment);
            if (out != null && false)                                             // dump oval_det to a file
                out.println("oval det =  ," + fitted.getc() + print_radial() + print_coord(abeta) + ", " + calc_oval_det(null));
            if (out != null && false)                                            // dump d0 - d3 to a file
            {
                out.print("d1d2 = ," + fitted.getc());
                for (i = 0; i < 4; i++)
                    out.print(", " + Math.sqrt((Bezx[3*i + 1] - Bezx[3*i])*(Bezx[3*i + 1] - Bezx[3*i]) + (Bezy[3*i + 1] - Bezy[3*i])*(Bezy[3*i + 1] - Bezy[3*i])));
                out.println();
            }
            if (out != null && true)                                            // dump Jac to a file
            {
                out.println("dFdc = " + fitted.getc() + print_radial() + print_coord(abeta));     // this is just a header for dump_Jac file data
                out.print("a = np.array([");
                for (i = 0; i < Augment.length; i++)
                {
                    if (i > 0) out.print(", ");
                    out.print("[");
                    for (j = 0; j < Augment.length; j++)
                    {
                        if (j > 0) out.print(", ");
                        out.print(Augment[i][j]);
                    }
                    out.print("]");
                }
                out.println("])");
                out.println();
            }
            double origin_y = 250;                                                  // used only for N-point Bezier
            System.out.print("M " + Bezx[0] + ", " + (origin_y + Bezy[0]) + " C");  // dump N-point Bezier
            for (i = 1; i < 3*ax.length; i++)
                System.out.print(" " + Bezx[i] + ", " + (origin_y + Bezy[i]));
            System.out.println(" " + Bezx[0] + ", " + (origin_y + Bezy[0]));
            //write_BSplineN_data(null);                      // use only for dumping to file from dump_Jac()
            //write_BSplineN_data(Jac);
            //write_BSplineN_data(Augment);
            //for (i = 0; i < N; i++)                         // test sum of t2ddx and t2ddy
            //    System.out.println(i + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + denom[i]);
            //double[][] dot_product = new double [N][8];     // test code just for oval test
            //j = 4;                                          // segment index (0 - Splinex.length)
            //int l = 0;
            //for (k = 0; k < N; k++)
            //    if ((t2[k] > j) && (t2[k] < j + 1))
            //    {
            //        for (i = 0; i < 4; i++)                 // only need 4 control pts
            //        {
            //            dot_product[l][i] = f_gx[k]*calc_dfdai((j + i) % ax.length, t2[k]);
            //            dot_product[l][i + 4] = f_gy[k]*calc_dfdai((j + i) % ax.length, t2[k]);
            //        }
            //        l++;
            //    }
            //System.out.println("dot_product:");
            //BSpline5.dump_Jac(dot_product);
            //correlate_f_dfda(dot_product, j);                                  // end of test code
            double retVal = solve_at_x_y(false);                // final run just for good measure
            //System.out.println("d2Fdadc    = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(d2Fdadc));
            //double[] dadc = BSpline5.gaussj(Jac, d2Fdadc);     // -response to change in c
            //System.out.println("dadc       = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(dadc));
            //if (out != null)                                    // use only for replicating file 'betarun8.csv'
            //    out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + print_coord(abeta) + ", " + (float) retVal + ", " + (float) Jacdet);
            //System.out.println("oval det = ," + calc_oval_det(d2Fdadc));
            write_one_line(retVal);                                   // if success, write output to file
            return retVal;
        }
        else
        {
            if (out != null)                                    // use only for replicating file 'betarun8.csv'
                out.println("convergence failed at c = " + fitted.getc());
            System.out.println("\nNOT converged after " + loop + " loops! (" + print_coord(dela) + ")");
        }
        return Double.NaN;
    }

    private static double solve_at_x_y(boolean print)
    {
        // 8-point, uniform, closed, cubic Beta2-Spline
        // perform a single calculation of a complete t2[] profile
        // at a given set {x, y, beta2}, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);
        //System.out.println("theta_start = " + theta_start + ", " + theta_end);

        convert_to_Bezier();                                    // refresh the Bezier representation
        if (t2[N] == 0)
            System.out.println("__start Beta2-Spline8 theta c t x y = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + print_coord(ax) + print_coord(ay) + print_coord(abeta));
        else
        {
            System.out.println("__solve at new d1 d2 rms = , , , , , " + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + calc_error());
            System.out.println("__solve at new d1 d2 rms = , , , , , " + print_radial() + print_coord(abeta) + ", " + calc_error());
        }

        if (print) System.out.println("\n , t1, t2, t2ddx0, t2ddx1, t2ddx2, t2ddx3, t2ddx4, t2ddx5, t2ddx6, t2ddx7, t2ddy0, t2ddy1, t2ddy2, t2ddy3, t2ddy4, t2ddy5, t2ddy6, t2ddy7, t2ddbeta0, t2ddbeta1, t2ddbeta2, t2ddbeta3, t2ddbeta4, t2ddbeta5, t2ddbeta6, t2ddbeta7");
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i);
            if ((i == N && Math.abs(t2[N] - t2[0] - ax.length) > TOL)
            ||  (t2[i] < -5000*TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                return Double.NaN;
            }
            for (int j = 0; j < ax.length; j++)
                t2ddx[j][i] = calc_t2dd(j, i, t2[i], "x");
            for (int j = 0; j < ay.length; j++)
                t2ddy[j][i] = calc_t2dd(j, i, t2[i], "y");
            for (int j = 0; j < abeta.length; j++)
                t2ddbeta[j][i] = calc_t2dd(j, i, t2[i], "beta");
            if (print)
                System.out.println(i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + print_coord(t2ddx, i) + print_coord(t2ddy, i) + print_coord(t2ddbeta, i));
                //System.out.println(i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + calc_t2dd(Integer.MAX_VALUE, i, t2[i], "c"));
        }
        double retVal = calc_error();
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + (float) retVal + ", " + (float) Jacdet);
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2);
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + retVal*retVal*180.0*180.0/2.0);
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + print_coord(abeta) + ", " + (float) retVal + ", " + (float) Jacdet);
        return retVal;
    }

    private static void scan_steepest()
    {
        // calculate a streamline using steepest-descent
        // starting at a midpoint between a saddle point and a minimum
        // use method of 'Chong and Zak', 'An Introduction to Optimization', p.121
        // move in increments of 'steepest' in the direction of the gradient -dFda[]
        // to run this, temporarily comment out lines 380, 408 and 410 to reduce the output
        // also       , temporarily enable line 222 to abort 'iterate_at_x_y()'
        // see Book 9, page 42

        final int Nstr = 30;
        double retVal;
        double steepnum, steepdenom;

        System.out.println("scan steepest, " + Nstr);
        iterate_at_x_y(null);                                       // initiallize dFdd
        retVal = calc_error();
        System.out.println("scan = , , , " + (float) fitted.getc() + ", , 0" + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2);

        for (int i = 0; i < Nstr; i++)
        {
            t2[N] = 0;                                          // just to control the output
            steepnum = BSpline5.multvv(dFda, dFda);
            steepdenom = BSpline5.multvv(dFda, BSpline5.multmv(Jac, dFda));
            for (int j = 0; j < ax.length; j++)
            {
                ax[j] -= dFda[j]*steepnum/steepdenom;
                ay[j] -= dFda[j + ax.length]*steepnum/steepdenom;
                abeta[j] -= dFda[j + 2*ax.length]*steepnum/steepdenom;
            }
            iterate_at_x_y(null);                                   // re-calculate dFdd
            retVal = calc_error();
            System.out.println("scan = , , , " + (float) fitted.getc() + ", , " + (i + 1) + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2 + ", " + steepnum/steepdenom);
        }
    }

    private static void get_Bezier_endpoints_from_Beta_Spline()
    {
        //fitted = new epiTrochoidFxn(5.4);              // fix fix temporary location
        fitted = new SuperEllipse(0.90955);
        //String str = "gauss t2[] @ , 90.0, 90.0, 17.1, 800, , 181.11768, 191.3406, 168.46349, 191.3981, 181.11768, 191.3406, 168.46349, 191.3981, -73.178055, -8.765734, 29.248032, 81.26961, 106.821945, 171.23427, -150.75197, -98.73039, 4.3080115, 1.7127049, 0.0033888216, 1.1912025, 4.3080115, 1.7127049, 0.0033888216, 1.1912025, 6.869E-4, NaN";
        // convert a 8-point Beta-Spline to 8 spliced cubic Beziers (endpoints only)
        double[] beta_in = new double[] {
//168.08003, 197.48953, 168.08003, 197.48953, 168.08003, 197.48953, 168.08003, 197.48953, -63.464863, -4.8974257, 26.535137, 85.10258, 116.53513, 175.10257, -153.46486, -94.89742, 4.023081, 2.4261758, 4.023081, 2.4261758, 4.023081, 2.4261758, 4.023081, 2.4261758
//168.08003, 197.48953, 168.08003, 197.48953, 168.08003, 197.48953, 168.08003, 197.48953, -63.464863, -4.8974257, 26.535137, 85.10258, 116.53513, 175.10257, -153.46486, -94.89742, 4.023081, 2.4261758, 4.023081, 2.4261758, 4.023081, 2.4261758, 4.023081, 2.4261758
//172.97676, 187.03522, 172.97676, 187.03522, 172.97676, 187.03522, 172.97676, 187.03522, -46.46174, -0.96474224, 43.53826, 89.035255, 133.53827, 179.03526, -136.46173, -90.964745, 0.68979716, -0.1626248, 0.68979716, -0.1626248, 0.68979716, -0.1626248, 0.68979716, -0.1626248
//170.5169, 180.04228, 170.5169, 150.57526, 170.5169, 180.04228, 170.5169, 150.57526, -35.763977, 2.3606814E-13, 35.763977, 90.0, 144.23602, 180.0, -144.23602, -90.0, -0.69682413, 2.4031181, -0.69682413, -0.46950108, -0.69682413, 2.4031181, -0.69682413, -0.46950108
185.6804, 180.11008, 185.6804, 180.11008, 185.6804, 180.11008, 185.6804, 180.11008, -45.0, 4.0410976E-14, 45.0, 90.0, 135.0, 180.0, -135.0, -90.0, -0.32908234, -0.79496855, -0.32908234, -0.79496855, -0.32908234, -0.79496855, -0.32908234, -0.79496855
//179.80605, 181.36041, 180.14049, 182.09087, 179.5185, 184.15178, 170.6738, 188.99484, -69.78527, -17.851955, 19.662558, 73.21765, 110.6087, 166.3585, -138.03001, -93.310196, -0.03717412, 0.085394986, 0.21248694, 0.22758241, 0.4741566, 0.3797251, -2.5248592, 7.7442045
//179.85231, 180.1501, 179.85231, 180.1501, 179.85231, 180.1501, 179.85231, 180.1501, -62.850533, -17.501976, 27.149467, 72.498024, 117.14947, 162.49802, -152.85052, -107.501976, -0.13609967, -0.18932837, -0.13609967, -0.18932837, -0.13609967, -0.18932837, -0.13609967, -0.18932837
//159.98505, 199.99086, 159.98505, 199.99086, 159.98505, 199.99086, 159.98505, 199.99086, -45.0, -1.5553355E-14, 45.0, 90.0, 135.0, 180.0, -135.0, -90.0, -2.1210384, -0.12768523, -2.1210384, -0.12768523, -2.1210384, -0.12768523, -2.1210384, -0.12768523
        };
        //if (fitted.getc() < 0)                                          // use only for negative c
        //    reflect_at_45_degrees(beta_in);
        if (beta_in.length != 24)
        {
            System.out.println("convert_beta_to_Bezier: incorrect length = " + beta_in.length);
            return;
        }
        //double mix = 0;                                               // use only for interpolation
        //double[] rad_in1 = new double[] {180.1914, 180.1914, 180.1914, 180.1914, 180.1914, 180.1914, 180.1914, 180.1914, -69.35823, -20.641771, 20.641771, 69.35823, 110.64177, 159.35823, -159.35823, -110.64177, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702};
        //double[] rad_in1 = new double[] {175.49046, 184.52057, 175.49046, 184.52057, 175.49046, 184.52057, 175.49046, 184.52057, -45.0, -2.737488E-13, 45.0, 90.0, 135.0, -180.0, -135.0, -90.0, 3.0205624, -0.6741331, 3.0205624, -0.6741331, 3.0205624, -0.6741331, 3.0205624, -0.6741331};
        //double[] rad_in2 = new double[] {175.4421, 184.44278, 175.4421, 184.44278, 175.4421, 184.44278, 175.4421, 184.44278, -45.0, -2.2157824E-13, 45.0, 90.0, 135.0, -180.0, -135.0, -90.0, -2.3263218, 22.012825, -2.3263218, 22.012825, -2.3263218, 22.012825, -2.3263218, 22.012825};
        //for (int i = 0; i < beta_in.length; i++)
        //    beta_in[i] = mix*rad_in1[i] + (1 - mix)*rad_in2[i];
        System.out.println("convert_Beta_Spline_to_Bezier radial at c = ," + fitted.getc() + print_coord(beta_in));
        abeta = new double[8];
        ax = new double[8];
        ay = new double[8];
        for (int i = 0; i < ax.length; i++)
        {
            ax[i] = beta_in[i]*Math.cos(beta_in[i + ax.length]*Math.PI/180);
            ay[i] = beta_in[i]*Math.sin(beta_in[i + ax.length]*Math.PI/180);
            abeta[i] = beta_in[i + 2*ax.length];
        }
        double incr = 0; //-4.68; //20; // magic factor = 72;
        //double theta = 30;
        //double[] eig_vec_angular0 = new double[] {0.27658725,0.21717687,-0.21717687,-0.27658725,-0.27658725,-0.21717687,0.21717687,0.27658725,0.21717687,0.27658725,0.27658725,0.21717687,-0.21717687,-0.27658725,-0.27658725,-0.21717687,0.03651984,-0.03651984,0.03651984,-0.03651984,0.03651984,-0.03651984,0.03651984,-0.03651984};
        //double[] eig_vec_angular1 = new double[] {0.19880945,0.23317641,-0.26448369,0.2548129,0.19880945,0.23317641,-0.26448369,0.2548129,0.15533342,0.29835982,0.33826469,-0.19927081,0.15533342,0.29835982,0.33826469,-0.19927081,0.04314422,-0.03735079,-0.01050166,-0.02401306,-0.04314422,0.03735079,0.01050166,0.02401306};
        //double[] eig_vec_angular2 = new double[] {0.33826469,-0.19927081,0.15533342,0.29835982,0.33826469,-0.19927081,0.15533342,0.29835982,0.26448369,-0.2548129,-0.19880945,-0.23317641,0.26448369,-0.2548129,-0.19880945,-0.23317641,-0.01050166,-0.02401306,-0.04314422,0.03735079,0.01050166,0.02401306,0.04314422,-0.03735079};
        //double[] eig_vec_angular1 = new double[] {0.2548468501217222 , -0.0005358546507131523 , -0.2501227493954531 , -0.3496846604956383 , -0.25484685019326037 , 0.0005358546505479047 , 0.2501227492800065 , 0.34968466021419686 , 0.25484685012201974 , 0.3496846604329374 , 0.25012274939535156 , 0.000535854650730883 , -0.25484685019313963 , -0.349684660360145 , -0.2501227492801641 , -0.0005358546502497118 , -1.3171685964152857e-12 , 0.014561744013351984 , -3.8062955565187906e-13 , -0.014561743954782047 , 8.22205476445903e-13 , 0.014561743856429827 , 8.22830356368015e-13 , -0.01456174388611318};
        //double[] eig_vec_angular2 = new double[] {-0.25012274943127843 , -0.0005358546505758871 , 0.25484684997240187 , 0.34968465988587116 , 0.25012274924104877 , 0.0005358546508734316 , -0.25484685034491217 , -0.34968466072136284 , -0.2501227494308299 , -0.3496846603134654 , -0.2548468499729285 , 0.00053585464990715 , 0.2501227492419466 , 0.349684660581498 , 0.25484685034472954 , -0.0005358546510560035 , -2.4388269181940814e-12 , 0.014561744160439887 , -2.5197240127727838e-12 , -0.014561743978799007 , 3.997952197163507e-12 , 0.014561743744442983 , 8.377652683990572e-13 , -0.014561743860631209};
        double[] eig_vec_angular  = new double[] {0.29909275880246133 , 0.0005081598453183095 , 0.2990927589829023 , -6.431784553735664e-11 , -0.29909275904993343 , -0.00050815984526378 , -0.2990927588367384 , -5.4796016682456816e-11 , 0.20961899355446123 , -1.3871301329793084e-10 , -0.2096189936807442 , 0.002453100897947564 , -0.20961899372751047 , -1.6704270521601845e-10 , 0.20961899357819866 , -0.002453100897073815 , -0.0032926804965980744 , 0.011781790317967104 , -0.003292680497931598 , -0.48274219642777527 , -0.0032926804974798855 , 0.011781790317261915 , -0.003292680495699311 , -0.48274219627279313};
        for (int i = 0; i < ax.length; i++)                     // increment along the lowest eigenvector
        {
            //ax[i] += incr*(Math.cos(theta*Math.PI/180)*eig_vec_angular2[i] + Math.sin(theta*Math.PI/180)*eig_vec_angular1[i]);
            //ay[i] += incr*(Math.cos(theta*Math.PI/180)*eig_vec_angular2[i + ax.length] + Math.sin(theta*Math.PI/180)*eig_vec_angular1[i + ax.length]);
            //abeta[i] += incr*(Math.cos(theta*Math.PI/180)*eig_vec_angular2[i + 2*ax.length] + Math.sin(theta*Math.PI/180)*eig_vec_angular1[i + 2*ax.length]);
            ax[i] += incr*eig_vec_angular[i];
            ay[i] += incr*eig_vec_angular[i + ax.length];
            abeta[i] += incr*eig_vec_angular[i + 2*ax.length];
        }
        //double[] cart_in1 = new double[] {35.550457, 165.3004, 165.3004, 35.550457, -35.550457, -165.3004, -165.3004, -35.550457, -161.02158, -96.76103, 96.76103, 161.02158, 161.02158, 96.76103, -96.76103, -161.02158, 4.4096084, -1.0530666, -1.0530666, 4.4096084, 4.4096084, -1.0530666, -1.0530666, 4.4096084};
        //double[] cart_in2 = new double[] {143.71088, 198.06859, 143.71088, 1.13499325E-13, -143.71088, -198.06859, -143.71088, -2.0809191E-13, -119.76606, 5.3121654E-14, 119.76606, 161.96565, 119.76606, 6.35461E-13, -119.76606, -161.96565, -1.3941962, 4.584869, -1.3941962, 1.3657236, -1.3941962, 4.584869, -1.3941962, 1.3657236};
        //double[] cart_in = new double[] {42.725063, 191.89647, 152.34802, 27.31546, -42.725063, -191.89647, -152.34802, -27.31546, -180.40105, -26.47261, 73.33422, 191.52939, 180.40105, 26.47261, -73.33422, -191.52939, 9.370194, 3.2129233, 1.2039822, 2.505883, 9.370194, 3.2129233, 1.2039822, 2.505883};
        //double[] cart_in = new double[] {127.330395,179.85931,127.54462,0.30967167,-127.329524,-179.85948,-127.114946,0.31372904,-127.330395,-0.31372904,127.114946,179.85948,127.329524,-0.30967167,-127.54462,-179.85931,-1.629431,28.549601,-1.632839,28.682253,-1.6361683,28.682253,-1.632839,28.549601};
        //double[] cart_in = new double[] {128.91849,180.02619,128.91849,-1.01977E-12,-128.91849,-180.02619,-128.91849,7.1936E-13,-128.91849,5.0848E-13,128.91849,180.02619,128.91849,4.06193E-13,-128.91849,-180.02619,2.6112098,-1.4800188,2.6112098,-1.4800188,2.6112098,-1.4800188,2.6112098,-1.4800188};
        //double mix = 0.999;                                               // use only for interpolation
        for (int i = 0; i < ax.length; i++)                             // interpolate between 2 initial estimates
        {
            //ax[i] = mix*cart_in1[i] + (1 - mix)*cart_in2[i];
            //ay[i] = mix*cart_in1[i + ax.length] + (1 - mix)*cart_in2[i + ax.length];
            //abeta[i] = mix*cart_in1[i + 2*ax.length] + (1 - mix)*cart_in2[i + 2*ax.length];
            //ax[i] = cart_in[i];
            //ay[i] = cart_in[i + ax.length];
            //abeta[i] = cart_in[i + 2*ax.length];
        }
        //ax[3] += 0.01;                                                  // increment for dFdai calculation
        //ay[4] += 0.01;                                                  // increment for dFdai calculation
        //abeta[5] += 0.001;                                               // increment for dFdai calculation
        //abeta[6] += 0.01;                                               // increment for dFdai calculation
        System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta));
        System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_radial() + print_coord(abeta));
    }

    private static void get_Bezier_endpoints_from_B_Spline()
    {
        // convert a 8-point B-Spline to 8 spliced cubic Beziers (endpoints only)
        int j;
        double[] rin = new double[] {181.62039, 217.36208};
        double[] thin = new double[] {-47.694225, -1.0211561};
        //rin = new double[] {202, 202};
        //thin = new double[] {-84.34531, -5.6546926};
        double[] rad = new double[] {rin[0], rin[1], rin[0], rin[1], rin[0], rin[1], rin[0], rin[1],
                                     thin[0], thin[1], thin[0] + 90, thin[1] + 90, thin[0] + 180, thin[1] + 180, thin[0] + 270, thin[1] + 270};
        //rad = new double[] {181.62039, 217.36208, 181.62039, 217.36208, 181.62039, 217.36208, 181.62039, 217.36208, -47.694225, -1.0211561, 42.305775, 88.97884, 132.30577, 178.97885, -137.69423, -91.02116};
        //rad = new double[] {135.0592463,219.317122,192.3538406,259.4224354,254.6998233,72.80109889,212.6029163,121.6552506,44.70002484,24.22774532,8.972626615,62.44718842,46.90915243,164.0546041,48.81407483,9.462322208};
        if (rad.length != 16)
        {
            System.out.println("convert_to_Bezier: incorrect length = " + rad.length);
            return;
        }
        System.out.println("convert_B_Spline_to_Bezier radial at c = ," + fitted.getc() + print_coord(rad));
        //abeta = new double[] {-1.8, 3, -1.8, 3, -1.8, 3, -1.8, 3};  // beta2 values of Bezier endpoints
        abeta = new double[] {-1, 1, -1, 1, -1, 1, -1, 1};
        ax = new double[8];
        ay = new double[8];
        for (int i = 0; i < ax.length; i++)
        {
            ax[i] = 4*rad[i]*Math.cos(rad[i + 8]*Math.PI/180)/6;
            ay[i] = 4*rad[i]*Math.sin(rad[i + 8]*Math.PI/180)/6;
            j = (i + 1) % 8;
            ax[i] += rad[j]*Math.cos(rad[j + 8]*Math.PI/180)/6;
            ay[i] += rad[j]*Math.sin(rad[j + 8]*Math.PI/180)/6;
            j = (i + 7) % 8;
            ax[i] += rad[j]*Math.cos(rad[j + 8]*Math.PI/180)/6;
            ay[i] += rad[j]*Math.sin(rad[j + 8]*Math.PI/180)/6;
        }

        //ax = new double[] {126.57209, 180.99977, 126.57209, -3.6856538E-12, -126.57209, -180.99977, -126.57209, 5.4690805E-12};                      // overwrite
        //ay = new double[] {-126.57209, 4.6605267E-12, 126.57209, 180.99977, 126.57209, -3.0797032E-12, -126.57209, -180.99977};                      // overwrite
        //abeta = new double[] {-2, 0, -2, 0, -2, 0, -2, 0};                      // overwrite
        System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta));
        System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_radial() + print_coord(abeta));
    }

    private static void convert_to_Bezier()
    {
        int i, j;
        // see Spiro2SVG Book 9, p. 53
        double[][] convert = new double[8][8];
        double[] vec = new double[8];
        double[] lhs;
        double[] d;
        //int k = 7;                                  // index of coord to change - test only
        //ay[k] += 0.01;                              // increment Bezier coord   - test only
        //abeta[k] += 0.01;
        //System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta));
        //System.out.println(" , " + k + ": " + ax[k] + " : "+ abeta[k]);
        //int d2ddj = 1;
        //int d2ddk = 1;
        //abeta[d2ddj] += 0*0.01;                                      // calc numer d2ddbetajdbetak
        //abeta[d2ddk] += 0*0.01;                                      // calc numer d2ddxdbeta
        for (i = 0; i < convert.length; i++)
            for (j = 0; j < convert.length; j++)
                convert[i][j] = 0;
        for (i = 0; i < convert.length; i++)
        {
            convert[i][i] = 1;
            convert[i][(i + 1) % convert.length] = 4 + abeta[(i + 1) % convert.length]/2;
            convert[i][(i + 2) % convert.length] = 1;
        }
        //System.out.println("convert[i][j]:");
        //for (i = 0; i < convert.length; i++)
        //{
        //    for (j = 0; j < convert.length; j++)
        //        System.out.print(convert[i][j] + " ,");
        //    System.out.println();
        //}
        for (i = 0; i < vec.length; i++)                        // calculate dx
            vec[i] = ax[(i + 2) % vec.length] - ax[i];
        d = BSpline5.gaussj(convert, vec);
//        System.out.println("dx:");
//        for (i = 0; i < convert.length; i++)
//            System.out.println(i + ", " + d[i]);
        for (i = 0; i < convert.length; i++)
        {
            Bezx[3*i] = ax[i];
            Bezx[3*i + 1] = ax[i] + d[i];
            Bezx[3*i + 2] = ax[(i + 1) % convert.length] - d[(i + 1) % convert.length];
        }
        for (j = 0; j < vec.length; j++)                        // calculate ddxdbeta
        {
            for (i = 0; i < vec.length; i++)
                vec[i] = 0;
            vec[(j + vec.length - 1) % vec.length] = -d[j]/2;
            //for (i = 0; i < vec.length; i++)
            //    System.out.println("vec " + j + ", " + i + ", " + vec[i]);
            lhs = BSpline5.gaussj(convert, vec);
            for (i = 0; i < vec.length; i++)
                ddxdbeta[i][j] = lhs[i];
        }
        for (j = 0; j < vec.length; j++)                        // calculate d2dx[i]/dbeta[j]/dbeta[k]
            for (int k = 0; k < vec.length; k++)
            {
                for (i = 0; i < vec.length; i++)
                    vec[i] = 0;
                vec[(k + vec.length - 1) % vec.length] += -ddxdbeta[k][j]/2;
                vec[(j + vec.length - 1) % vec.length] += -ddxdbeta[j][k]/2;
                lhs = BSpline5.gaussj(convert, vec);
                for (i = 0; i < vec.length; i++)
                    d2dxdbetaibetaj[i][j][k] = lhs[i];
            }
        //ay[d2ddj] += 0.01;                                      // calc numer d2ddxdbeta
        for (i = 0; i < vec.length; i++)                        // calculate dy
            vec[i] = ay[(i + 2) % vec.length] - ay[i];
        d = BSpline5.gaussj(convert, vec);
//        System.out.println("dy:");
//        for (i = 0; i < convert.length; i++)
//            System.out.println(i + ", " + d[i]);
        for (i = 0; i < convert.length; i++)
        {
            Bezy[3*i] = ay[i];
            Bezy[3*i + 1] = ay[i] + d[i];
            Bezy[3*i + 2] = ay[(i + 1) % convert.length] - d[(i + 1) % convert.length];
        }
        for (j = 0; j < vec.length; j++)                        // calculate ddydbeta
        {
            for (i = 0; i < vec.length; i++)
                vec[i] = 0;
            vec[(j + vec.length - 1) % vec.length] = -d[j]/2;
            lhs = BSpline5.gaussj(convert, vec);
            for (i = 0; i < vec.length; i++)
                ddydbeta[i][j] = lhs[i];
        }
        for (j = 0; j < vec.length; j++)                        // calculate d2dy[i]/dbeta[j]/dbeta[k]
            for (int k = 0; k < vec.length; k++)
            {
                for (i = 0; i < vec.length; i++)
                    vec[i] = 0;
                vec[(k + vec.length - 1) % vec.length] += -ddydbeta[k][j]/2;
                vec[(j + vec.length - 1) % vec.length] += -ddydbeta[j][k]/2;
                lhs = BSpline5.gaussj(convert, vec);
                for (i = 0; i < vec.length; i++)
                    d2dydbetaibetaj[i][j][k] = lhs[i];
            }
        for (i = 0; i < vec.length; i++)                        // calculate dddx & dddy
        {
            for (j = 0; j < vec.length; j++)
                vec[j] = 0;
            vec[(i + vec.length - 2) % vec.length] = 1;
            vec[i] = -1;
            lhs = BSpline5.gaussj(convert, vec);
            for (j = 0; j < vec.length; j++)
                dddx[j][i] = lhs[j];
        }
        for (j = 0; j < vec.length; j++)                        // calculate d2d[i]/dx[j]/dbeta[k] & d2d[i]/dy[j]/dbeta[k]
            for (int k = 0; k < vec.length; k++)
            {
                for (i = 0; i < vec.length; i++)
                    vec[i] = 0;
                vec[(k + vec.length - 1) % vec.length] = -dddx[k][j]/2;
                lhs = BSpline5.gaussj(convert, vec);
                for (i = 0; i < vec.length; i++)
                    d2ddxdbeta[i][j][k] = lhs[i];
            }
/*
        System.out.println("ddxdbeta[i][j]:");
        for (i = 0; i < convert.length; i++)
        {
            for (j = 0; j < convert.length; j++)
                System.out.print(ddxdbeta[i][j] + " ,");
            System.out.println();
        }
        System.out.println("ddydbeta[i][j]:");
        for (i = 0; i < convert.length; i++)
        {
            for (j = 0; j < convert.length; j++)
                System.out.print(ddydbeta[i][j] + " ,");
            System.out.println();
        }
*/
        //System.out.println("dddx[i][j]:");
        //for (i = 0; i < convert.length; i++)
        //{
        //    for (j = 0; j < convert.length; j++)
        //        System.out.print(dddx[i][j] + " ,");
        //    System.out.println();
        //}
        //System.out.println(", y[j], beta[k], " + ay[d2ddj] + ", " + abeta[d2ddk]);
        //System.out.println(", beta[j], beta[k], " + abeta[d2ddj] + ", " + abeta[d2ddk]);
        //System.out.println("i, d, ddydbeta[" + d2ddj + "], ddydbeta[" + d2ddk + "], d2ddxdbeta[" + d2ddj + "][" + d2ddk + "]");
        //System.out.println("i, d, dddx[" + d2ddj + "], ddydbeta[" + d2ddk + "], d2ddxdbeta[" + d2ddj + "][" + d2ddk + "]");
        //for (i = 0; i < vec.length; i++)
        //    System.out.println(i + ", " + d[i] + ", " + ddydbeta[i][d2ddj] + ", " + ddydbeta[i][d2ddk] + ", " + d2dydbetaibetaj[i][d2ddj][d2ddk]);
            //System.out.println(i + ", " + d[i] + ", " + dddx[i][d2ddj] + ", " + ddydbeta[i][d2ddk] + ", " + d2ddxdbeta[i][d2ddj][d2ddk]);
        //for (i = 0; i < N; i++)
        //{
        //    double di = i/100.0;
            //System.out.println(i + ", " + multvv(Bezx, (int) di, Bi3(di % 1)) + ", " + calc_dfdai(k, di));
        //    System.out.println(i + ", " + multvv(Bezy, (int) di, Bi3(di % 1)) + ", " + calc_dfdbetai(ddxdbeta, k, di) + ", " + calc_dfdbetai(ddydbeta, k, di));
        //}
    }

    private static double calc_oval_det(double[] d2Fdadc0)
    {
        // see Spiro2SVG Book 10, p.11
        double d1 = Math.sqrt((Bezx[1] - Bezx[0])*(Bezx[1] - Bezx[0]) + (Bezy[1] - Bezy[0])*(Bezy[1] - Bezy[0]));
        double d2 = Math.sqrt((Bezx[3] - Bezx[2])*(Bezx[3] - Bezx[2]) + (Bezy[3] - Bezy[2])*(Bezy[3] - Bezy[2]));
        double theta_s = Math.PI/4;
        double theta_e = Math.PI/2;
        double[][] oval = new double[][] {{2*(Bezy[3] - Bezy[0])*Math.cos(theta_s) - 2*(Bezx[3] - Bezx[0])*Math.sin(theta_s) - 2*d2*Math.sin(theta_e - theta_s), d1*Math.sin(theta_e - theta_s)},
                                          {d2*Math.sin(theta_e - theta_s), -2*(Bezy[3] - Bezy[0])*Math.cos(theta_e) + 2*(Bezx[3] - Bezx[0])*Math.sin(theta_e) - 2*d1*Math.sin(theta_e - theta_s)}};
        //System.out.println("oval det =  ," + fitted.getc() + ", " + d1 + ", " + d2 + ", " + BSpline5.detm(oval));
        //System.out.println("calc_oval_det = " + oval[0][0] + ", " + oval[0][1] + ", " + oval[1][0] + ", " + oval[1][1]);
        //System.out.println("beta , " + -oval[0][1]/Math.sqrt(oval[0][0]*oval[0][0] + oval[0][1]*oval[0][1]) + ", " + oval[0][0]/Math.sqrt(oval[0][0]*oval[0][0] + oval[0][1]*oval[0][1]));
/*
        System.out.println("d2Fdadc," + d2Fdadc0[16] + ", " + d2Fdadc0[17] + ", " + d2Fdadc0[18] + ", " + d2Fdadc0[19] + ", " + d2Fdadc0[20] + ", " + d2Fdadc0[21] + ", " + d2Fdadc0[22] + ", " + d2Fdadc0[23]);
        double[] dx = new double[8];
        double[] dy = new double[8];
        for (int i = 0; i < dx.length; i++)
        {
            dx[i] = Bezx[3*i + 1] - Bezx[3*i];
            dy[i] = Bezy[3*i + 1] - Bezy[3*i];
        }
        System.out.println("beta   " + print_coord(abeta));
        System.out.println("dx     " + print_coord(dx));
        System.out.println("dy     " + print_coord(dy));
*/
        return BSpline5.detm(oval);
    }

    private static void reflect_at_45_degrees(double[] a)
    {
        // transform radial input parameters a[] by reflecting about 45 degree axis
        // use when c is negative
        int i, j, k;
        double[] b = new double[24];
        int offset = 1;
        for (k = 7; k < 11; k++)
        {
            i = k % 8;
            j = 13 - k;
            //System.out.println(i + ", " + j);
            b[i] = a[j];
            b[j] = a[i];
            b[i + 8] = 45 - a[j + 8];
            b[j + 8] = 45 - a[i + 8];
            b[i + 16] = a[j + 16];
            b[j + 16] = a[i + 16];
        }
        for (i = 0; i < 8; i++)
        {
            j = (i + offset) % 8;
            a[i] = b[j];
            a[i + 8] = b[j + 8];
            a[i + 16] = b[j + 16];
        }
    }

    private static void read_one_line(double incr)
    {
        // read one line of Beta2_SplineN data (c, ri, thetai, betai) from a file
        //System.out.println("line_in = " + line_in[0] + ", " + line_in[23]);
        try
        {
            //BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Hippopede\\one_line.txt"));
            BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_SuperEllipse\\one_line.txt"));
            try
            {
                double[] line_in = new double[24];
                String str_in = "";
                while (instr.ready())
                    str_in = instr.readLine();
                if (!str_in.isEmpty())
                {
                    String[] str_split = str_in.split(",");
                    for (int i = 0; i < line_in.length; i++)
                        line_in[i] = Double.parseDouble(str_split[i + 6]);
                    //fitted = new epiTrochoidFxn(Double.parseDouble(str_split[3]) + incr);   // increment c
                    fitted = new SuperEllipse(Double.parseDouble(str_split[3]) + incr);   // increment c
                    abeta = new double[8];
                    ax = new double[8];
                    ay = new double[8];
                    for (int i = 0; i < ax.length; i++)
                    {
                        ax[i] = line_in[i]*Math.cos(line_in[i + ax.length]*Math.PI/180);
                        ay[i] = line_in[i]*Math.sin(line_in[i + ax.length]*Math.PI/180);
                        abeta[i] = line_in[i + 2*ax.length];
                    }
                    System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta));
                    System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_radial() + print_coord(abeta));
                    System.out.println("Beta2_Spline8 iterate_at_x_y = " + iterate_at_x_y(null) + "\n");
                }
                else
                    System.out.println("input file is empty!");
                instr.close();

            }
            catch (IOException e)
            {
                System.out.println("read input error : " + e.getMessage());
                return;
            }
        }
        catch (FileNotFoundException e)
        {
            System.out.println("input file not found : " + e.getMessage());
            return;
        }
    }

    private static void write_one_line(double rms)
    {
        // write one line of Beta2_SplineN data (c, ri, thetai, betai)
        // append to previous runs
        try
        {
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Epi_C2v\\one_line.txt", true);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Hippopede\\one_line.txt", true);
            FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_SuperEllipse\\one_line.txt", true);
            PrintWriter out = new PrintWriter(fw);
            out.print("one_line   @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", , ");
            for (int i = 0; i < ax.length; i++)
                out.print((float) Math.sqrt(ax[i]*ax[i] + ay[i]*ay[i]) + ", ");
            for (int i = 0; i < ax.length; i++)
                out.print((float) (Math.atan2(ay[i], ax[i])*180/Math.PI) + ", ");
            for (double c: abeta)
                out.print((float) c + ", ");
            out.println((float) rms + ", " + (float) Jacdet);
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("write_one_line() save error = " + e);}
    }

    private static void write_Beta2_SplineN_data()
    {
        // generate a file consisting of Jacobian matrices for eigenvalue calculation
        try
        {
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\eig_augment_raw_8.txt", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\eig_augment_temp.txt", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\flip_to_minus_c_output.txt", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\eig_d1d2_temp.txt", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\betarun8.csv", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\steepest\\junk.txt", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Limacon8\\eig_Jac_temp.txt", false);
            FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Hippopede\\eig_Jac_temp.txt", false);
            PrintWriter out = new PrintWriter(fw);
            read_Beta2_SplineN_data(out);
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("read_Beta2_SplineN_data() save error = " + e);}
    }

    private static void read_Beta2_SplineN_data(PrintWriter out)
    {
        // read Beta2_SplineN data (c, ri, thetai, betai) from a file
        // generate Jacobian matrix for egenvalue calculation
        String str_in = "";
        String[] str_split;
        double[] beta_in = new double[24];
        int line_in = 0;
        abeta = new double[8];
        ax = new double[8];
        ay = new double[8];

        try
        {
            //BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\betarun8.csv"));
            //BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\flip_to_minus_c.csv"));
            //BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\steepest\\5_points_final_edited.csv"));
            //BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Limacon8\\Limacon8.csv"));
            //BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Epi_C2v\\Epi_C2v.csv"));
            BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Hippopede\\Hippo_run.csv"));
            try
            {
                while (instr.ready()) // && line_in < 5)
                {
                    line_in++;
                    str_in = instr.readLine();
                    if (line_in >= 0 && line_in <= 3000)
                    {
                        if (str_in.isEmpty())
                            out.println("....................");
                        else if (str_in.startsWith("gauss t2[] @ ,") || str_in.startsWith("one_line   @ ,") || str_in.startsWith("scan = , , , "))
                        {
                            str_split = str_in.split(",");
                            //out.println("line_in = " + line_in + ", " + str_split[3] + ", " + str_in);
                            //fitted = new epiTrochoidFxn(Double.parseDouble(str_split[3]));  // read c value
                            fitted = new SuperEllipse(Double.parseDouble(str_split[3]));
                            for (int i = 0; i < beta_in.length; i++)                        // read input parameters
                                beta_in[i] = Double.parseDouble(str_split[i + 6]);
                            for (int i = 0; i < ax.length; i++)
                            {
                                ax[i] = beta_in[i]*Math.cos(beta_in[i + ax.length]*Math.PI/180);  // read in radial coords
                                ay[i] = beta_in[i]*Math.sin(beta_in[i + ax.length]*Math.PI/180);
                                //ax[i] = beta_in[i];                                 // read in Cartesian coords
                                //ay[i] = beta_in[i + ax.length];
                                abeta[i] = beta_in[i + 2*ax.length];
                            }
                            out.print(line_in + " - ");
                            //out.println(ax[0] + ", " + ay[0] + ", " + abeta[0]);
                            System.out.println("Beta2_Spline8 iterate_at_x_y = " + iterate_at_x_y(out) + "\n");
                        }
                    }
                }
            }
            catch (IOException e)
            {
                System.out.println("read input error : " + e.getMessage());
                return;
            }
        }
        catch (FileNotFoundException e)
        {
            System.out.println("input file not found : " + e.getMessage());
            return;
        }
    }

    private static void scan_fourfold_symmetry()
    {
        //fitted = new epiTrochoidFxn(4.43);                       // fix fix temporary location
        fitted = new SuperEllipse(0.1);             // fix fix bug bug
        double[] r_in_org = new double[] {175.51227, 184.37256};
        double[] th_in_org = new double[] {-45, 0};
        double[] beta_in_org = new double[] {-2.3176925, 22.393124};
        double del_r = 0;
        double del_th = 0;
        double del_beta = .5;
        double[] r_in = new double[2];
        double[] th_in = new double[2];
        double[] beta_in = new double[2];

        abeta = new double[8];
        ax = new double[8];
        ay = new double[8];
        for (int i = -1; i < 2; i += 2)
        for (int j = -1; j < 2; j += 2)
        for (int k = -1; k < 2; k += 2)
        for (int l = -1; l < 2; l += 2)
        for (int m = -1; m < 2; m += 2)
        for (int n = -1; n < 2; n += 2)
        {
            r_in[0] = r_in_org[0] + del_r*i;
            r_in[1] = r_in_org[1] + del_r*j;
            th_in[0] = th_in_org[0] + del_th*k;
            th_in[1] = th_in_org[1] + del_th*l;
            beta_in[0] = beta_in_org[0] + del_beta*m;
            beta_in[1] = beta_in_org[1] + del_beta*n;
            for (int p = 0; p < ax.length; p++)
            {
                ax[p] = r_in[p % 2]*Math.cos(((p/2)*90 + th_in[p % 2])*Math.PI/180);
                ay[p] = r_in[p % 2]*Math.sin(((p/2)*90 + th_in[p % 2])*Math.PI/180);
                abeta[p] = beta_in[p % 2];
            }
            System.out.println("final = " + i + ", " + j + ", " + k + ", " + l + ", " + m + ", " + n + ", " + iterate_at_x_y(null) + ", " + print_radial() + "\n");
        }
    }

    private static void gen_spline_points()
    {
        double Npt = 800;
        double d;
        double origin_x = 0;
        double origin_y = 0;
        double scale = 1;

        System.out.printf("M");
        for (int i = 0; i < Npt; i++)
        {
            d = i/100.0;
            System.out.printf(" %f, %f", origin_x + scale*multvv(Bezx, (int) d, Bi3(d % 1)),
                                         origin_y + scale*multvv(Bezy, (int) d, Bi3(d % 1)));
        }
        System.out.println("\n");
    }

    private static String print_coord(double[] coord)
    {
        String str = "";
        for (double c: coord)
            str += ", " + (float) c; // (float)
        return str;
    }

    private static String print_coord(double[][] coord, int index)
    {
        String str = "";
        for (int i = 0; i < coord.length; i++)
            str += ", " + (float) coord[i][index];
        return str;
    }

    private static String print_radial()
    {
        String str = "";
        for (int i = 0; i < ax.length; i++)
            str += ", " + (float) Math.sqrt(ax[i]*ax[i] + ay[i]*ay[i]);
        for (int i = 0; i < ax.length; i++)
            str += ", " + (float) (Math.atan2(ay[i], ax[i])*180/Math.PI);
        return str;
    }

    private static double integrate(double[] trap)
    {
        // trapezoidal rule integration of a cyclic fxn of t1 (N points)
        double ret = 0;
        for (double t: trap)
            ret += t;
        return ret/trap.length;
    }

    private static double calc_error()
    {
        double t1 = t1_start;
        double[] trap_in = new double[N];

        if (Math.abs(t2[N] - t2[0] - ax.length) > TOL)
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i < N; i++)
        {
            trap_in[i] = (multvv(Bezx, (int) t2[i], Bi3(t2[i] % 1)) - fitted.getx(t1))*(multvv(Bezx, (int) t2[i], Bi3(t2[i] % 1)) - fitted.getx(t1))
                       + (multvv(Bezy, (int) t2[i], Bi3(t2[i] % 1)) - fitted.gety(t1))*(multvv(Bezy, (int) t2[i], Bi3(t2[i] % 1)) - fitted.gety(t1));
            t1 += (t1_end - t1_start)/N;
        }
        //return Math.sqrt(integrate(trap_in))/(fitted.a + fitted.b);
        return Math.sqrt(integrate(trap_in))/180;       // fudge to scale the HippopedeFxn only
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

        // initial estimate using quadratic approximation

        //if (i == 0) t = 0;                                // original clamped Bezier default
        if (i == 0) t = 1;                                  // Oct 18, 2019, use this for the -45 degree sol'n for beta2-spline
        else t = t2[i-1];
        f = (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, dBi3(t % 1))
          + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, dBi3(t % 1));
        fprime = multvv(Bezx, (int) t, dBi3(t % 1))*multvv(Bezx, (int) t, dBi3(t % 1)) + (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, d2Bi3(t % 1))
               + multvv(Bezy, (int) t, dBi3(t % 1))*multvv(Bezy, (int) t, dBi3(t % 1)) + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, d2Bi3(t % 1));
        f2prime = 3*multvv(Bezx, (int) t, dBi3(t % 1))*multvv(Bezx, (int) t, d2Bi3(t % 1)) + (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, d3Bi3(t % 1))
                + 3*multvv(Bezy, (int) t, dBi3(t % 1))*multvv(Bezy, (int) t, d2Bi3(t % 1)) + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, d3Bi3(t % 1));
        //System.out.println("\ninit 1 t1 t2 =, " + i + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
            del_t = -fprime/f2prime;
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            if (del_t < -500000000*TOL || del_t > 3)
            {
                System.out.println("\nBad init1 at i t1 t2 =, " + i + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t + ", " + -5000*TOL);
                del_t = (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
                if (del_t < -5000*TOL || del_t > 3)
                {
                    System.out.println("Bad init2 at i t1 t2 =, " + i + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
                    t2[i] = Double.NaN;
                    return;
                }
            }
        }
        t += del_t;

        //System.out.println("roots = " + (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime + ", " + (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime);
        //System.out.println("init 2 t1 t2 =, " + i + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
        do
        {
            f = (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, dBi3(t % 1))
              + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, dBi3(t % 1));
            fprime = multvv(Bezx, (int) t, dBi3(t % 1))*multvv(Bezx, (int) t, dBi3(t % 1)) + (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, d2Bi3(t % 1))
                   + multvv(Bezy, (int) t, dBi3(t % 1))*multvv(Bezy, (int) t, dBi3(t % 1)) + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, d2Bi3(t % 1));
            if (loop > 100)
            {
                System.out.println("too many loops = " + loop);
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
        if (i == 0 && t < -5000*TOL)
        {
            System.out.println("t2 is negative at start : " + i + ", " + t);
            t2[i] = Double.NaN;
        }
        else if (i > 0 && t < t2[i - 1] - TOL)
        {
            System.out.println("t2 is decreasing : " + i + ", " + t + ", " + t2[i - 1]);
            t2[i] = Double.NaN;
        }
        else
            t2[i] = t;
    }

    private static double calc_t2dd(int j, int i, double u, String type)
    {
        // i is the time index for t1
        // j is the coefficient index for aj, betaj
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        denom[i] = multvv(Bezx, (int) u, dBi3(u % 1))*multvv(Bezx, (int) u, dBi3(u % 1)) + (multvv(Bezx, (int) u, Bi3(u % 1)) - X)*multvv(Bezx, (int) u, d2Bi3(u % 1))
                 + multvv(Bezy, (int) u, dBi3(u % 1))*multvv(Bezy, (int) u, dBi3(u % 1)) + (multvv(Bezy, (int) u, Bi3(u % 1)) - Y)*multvv(Bezy, (int) u, d2Bi3(u % 1));
        //System.out.println("calc_t2dd = " + j + ", " + i + ", " + u + ", " + type + ", " + denom[i]);
        double numer = Double.NaN;
        if (type.equals("x"))
            numer = multvv(Bezx, (int) u, dBi3(u % 1))*calc_dfdai(j, u)
                  + (multvv(Bezx, (int) u, Bi3(u % 1)) - X)*calc_d2fdudai(j, u);
        else if (type.equals("y"))
            numer = multvv(Bezy, (int) u, dBi3(u % 1))*calc_dfdai(j, u)
                  + (multvv(Bezy, (int) u, Bi3(u % 1)) - Y)*calc_d2fdudai(j, u);
        else if (type.equals("beta"))
            numer = multvv(Bezx, (int) u, dBi3(u % 1))*calc_dfdbetai(ddxdbeta, j, u)
                  + (multvv(Bezx, (int) u, Bi3(u % 1)) - X)*calc_d2fdudbetai(ddxdbeta, j, u)
                  + multvv(Bezy, (int) u, dBi3(u % 1))*calc_dfdbetai(ddydbeta, j, u)
                  + (multvv(Bezy, (int) u, Bi3(u % 1)) - Y)*calc_d2fdudbetai(ddydbeta, j, u);
        else if (type.equals("c"))
            numer = - multvv(Bezx, (int) u, dBi3(u % 1))*fitted.getdxdc(t1)
                    - multvv(Bezy, (int) u, dBi3(u % 1))*fitted.getdydc(t1);
        return -numer/denom[i];
    }

    private static double calc_dfdai(int i, double u)
    {
        // this works for both calc_dfdaxi and calc_dfdayi
        //System.out.println("calc_dfdai " + t2[0] + ", " + u);
        double ret = 0;
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;
        if (ui == i)
            ret = Bi3(u % 1)[0] + Bi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = Bi3(u % 1)[2] + Bi3(u % 1)[3];
        ret += dddx[ui][i]*Bi3(u % 1)[1];
        ret -= dddx[(ui + 1) % 8][i]*Bi3(u % 1)[2];
        //System.out.println(i + ", " + u + ", " + ret + ", " + dddx[ui][i] + ", " + dddx[(ui + 1) % 8][i]);
        return ret;
    }

    private static double calc_d2fdudai(int i, double u)
    {
        // this works for both calc_d2fdudaxi and calc_d2fdudayi
        double ret = 0;
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;
        if (ui == i)
            ret = dBi3(u % 1)[0] + dBi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = dBi3(u % 1)[2] + dBi3(u % 1)[3];
        ret += dddx[ui][i]*dBi3(u % 1)[1];
        ret -= dddx[(ui + 1) % 8][i]*dBi3(u % 1)[2];
        return ret;
    }

    private static double calc_dfdbetai(double[][] dddbeta, int i, double u)
    {
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;

        return dddbeta[ui][i]*Bi3(u % 1)[1]
             - dddbeta[(ui + 1) % 8][i]*Bi3(u % 1)[2];
    }

    private static double calc_d2fdudbetai(double[][] dddbeta, int i, double u)
    {
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;

        return dddbeta[ui][i]*dBi3(u % 1)[1]
             - dddbeta[(ui + 1) % 8][i]*dBi3(u % 1)[2];
    }

    private static double calc_d2fdaidbetaj(int i, int j, double u)
    {
        // this works for both calc_d2fdaxidbetaj and calc_d2fdayidbetaj
        //System.out.println("calc_d2fdaidbetaj " + t2[0] + ", " + u);
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;

        return d2ddxdbeta[ui][i][j]*Bi3(u % 1)[1]
             - d2ddxdbeta[(ui + 1) % 8][i][j]*Bi3(u % 1)[2];
    }

    private static double calc_d2fdbetaibetaj(double[][][] d2ddbetaibetaj, int i, int j, double u)
    {
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;

        return d2ddbetaibetaj[ui][i][j]*Bi3(u % 1)[1]
             - d2ddbetaibetaj[(ui + 1) % 8][i][j]*Bi3(u % 1)[2];
    }

    private static double[] Bi3(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:   Bi3 negative u value = " + u);
        if (u > 1 + TOL)
            System.out.println("WARNING:   Bi3 too large u value = " + u);
        return new double[] {(1 - u)*(1 - u)*(1 - u),
                             3*u*(1 - u)*(1 - u),
                             3*u*u*(1 - u),
                             u*u*u};
    }

    private static double[] dBi3(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:  dBi3 negative u value = " + u);
        if (u > 1 + TOL)
            System.out.println("WARNING:  dBi3 too large u value = " + u);
        return new double[] {-3*(1 - u)*(1 - u),
                              3*(1 - u)*(1 - 3*u),
                              3*u*(2 - 3*u),
                              3*u*u};
    }

    private static double[] d2Bi3(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING: d2Bi3 negative u value = " + u);
        if (u > 1 + TOL)
            System.out.println("WARNING: d2Bi3 too large u value = " + u);
        return new double[] { 6*(1 - u),
                              3*(6*u - 4),
                              3*(2 - 6*u),
                              6*u};
    }

    private static double[] d3Bi3(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING: d3Bi3 negative u value = " + u);
        if (u > 1 + TOL)
            System.out.println("WARNING: d3Bi3 too large u value = " + u);
        return new double[] { -6,
                              18,
                             -18,
                               6};
    }

    public static double multvv(double[] v1, int v1Pos, double[] v2)
    {
        if (v1.length != Bezx.length) return Double.NaN;
        if (v2.length != 4) return Double.NaN;
        double retVal = 0;

        if (v1Pos < 0) System.out.println("error in multvv: v1Pos = " + v1Pos);
        for (int i = 0; i < v2.length; i++)
            retVal += v1[(3*v1Pos + i) % v1.length]*v2[i];
        return retVal;
    }

    private static double calc_dfdai_at_di(int i, double u)
    {
        // this works for both calc_dfdaxi and calc_dfdayi
        // this is a simplified calc of dfxdaxi, or dfydayi, holding di constant
        // normally we would hold betai fixed, which would introduce coupling due to changing dxi or dyi

        double ret = 0;
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;
        if (ui == i)
            ret = Bi3(u % 1)[0] + Bi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = Bi3(u % 1)[2] + Bi3(u % 1)[3];
        //ret += dddx[ui][i]*Bi3(u % 1)[1];
        //ret -= dddx[(ui + 1) % 8][i]*Bi3(u % 1)[2];
        //System.out.println(i + ", " + u + ", " + ret + ", " + dddx[ui][i] + ", " + dddx[(ui + 1) % 8][i]);
        return ret;
    }

    private static double calc_dfddi_at_betai(int i, double u)
    {
        // this is a simplified calc of dfxddi, holding betai constant
        // it returns the response to the scalar magnitude of di, must subsequently decompose into x, y, components
        double ret = 0;
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;
        if (ui == i)
            ret = Bi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = -Bi3(u % 1)[2];
        return ret;
    }

    private static double calc_d2fdudai_at_di(int i, double u)
    {
        // this works for both calc_dfdaxi and calc_dfdayi
        // this is a simplified calc of d2fxdudaxi, or d2fydudayi, holding di constant
        // normally we would hold betai fixed, which would introduce coupling due to changing dxi or dyi

        double ret = 0;
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;
        if (ui == i)
            ret = dBi3(u % 1)[0] + dBi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = dBi3(u % 1)[2] + dBi3(u % 1)[3];
        return ret;
    }

    private static double calc_d2fduddi_at_betai(int i, double u)
    {
        // this is a simplified calc of d2fxduddi, holding betai constant
        // it returns the response to the scalar magnitude of di, must subsequently decompose into x, y, components
        double ret = 0;
        if (u < t2[0]) return Double.NaN;
        if (u > t2[0] + 8) return Double.NaN;
        int ui = (int) u;
        ui = ui % 8;
        if (ui == i)
            ret = dBi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = -dBi3(u % 1)[2];
        return ret;
    }
}
