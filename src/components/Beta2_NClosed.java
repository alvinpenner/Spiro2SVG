
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
    private static epiTrochoidFxn fitted;
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
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.0000000001;

    public static void main (String[] args)
    {
        //fitted = new epiTrochoidFxn(-5.5);
        //scan_fourfold_symmetry();
        //get_Bezier_endpoints_from_Beta_Spline();
        //System.out.println("Beta2_Spline8 iterate_at_x_y = " + iterate_at_x_y(null) + "\n");       // normally include this line
        //System.out.println("Beta2_Spline8 solve_at_x_y = " + solve_at_x_y(true) + "\n");
        write_Beta2_SplineN_data();                              // generate eigenvalue data
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
        final int MAXLOOP = 1000;
        double[] f_gx = new double[N];
        double[] f_gy = new double[N];
        double[] dfxdu = new double[N];
        double[] dfydu = new double[N];

        double[] trap_in = new double[N];
        double[][] Jac = new double[3*ax.length][3*ax.length];
        double[] dFda = new double[3*ax.length];
        double[] d2Fdadc = new double[3*ax.length];             // augmented matrix
        double dFdc;
        double d2Fdcdc;                                         // augmented matrix
        double[][] Augment = new double[3*ax.length + 1][3*ax.length + 1];
        double[] dela;                                          // (-Δxi, -Δyi)
        int loop = 0;
        int i, j, k;
        boolean outside;

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

            // calculate determinant of augmented matrix

            //for (i = 0; i < 3*ax.length; i++)
            //{
            //    for (j = 0; j < 3*ax.length; j++)
            //        System.out.print(" ," + Jac[i][j]);
            //    System.out.println();
            //}
            for (k = 0; k < N; k++)
                trap_in[k] = -f_gx[k]*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N) - f_gy[k]*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N);
            dFdc = integrate(trap_in);
            //System.out.println("check Aug ," + "c" + ", " + checkj + ", " + dFdc + ", " + dFda[checkj] + ", " + d2Fdadc[checkj]);
            for (k = 0; k < N; k++)
                trap_in[k] = fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)
                           + fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)
                           - denom[k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c")*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
            d2Fdcdc = integrate(trap_in);
            //System.out.println("check d2Fdc2 , " + ", " + dFdc + ", " + d2Fdcdc);
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

            if (loop == -1)  // perform one iteration of steepest descent as a special initialization
            {
                BSpline5.dump_Jac(Jac);
                System.out.println("loop == 1");
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
            outside = false;
            for (k = 0; k < dela.length; k++)
                if (Math.abs(dela[k]) > TOL)
                    outside = true;
            //outside = false;                            // force convergence
        } while ((loop < MAXLOOP) && outside);
        if (loop < MAXLOOP)
        {
            System.out.println("\nt2 = [" + t2[0] + ", " + t2[N] + "]");
            if (Math.abs(t2[N] - t2[0] - ax.length) > TOL
            ||  t2[0] < -TOL)
            {
                System.out.println("bad converged t2[i]: ABORT");
                return Double.NaN;
            }
            System.out.println("\n__converged in " + loop + " at new xi yi betai = , , , , , , " + print_coord(ax) + print_coord(ay) + print_coord(abeta));
            //Jacdet = BSpline5.detm(Jac);
            System.out.println("dFdc = " + fitted.getc() + print_radial() + print_coord(abeta));     // this is just a header for dump_Jac file data
            BSpline5.dump_Jac(Jac);
            //BSpline5.dump_Jac(Augment);
            if (out != null)                                            // dump Jac to a file
            {
                out.println("dFdc = " + fitted.getc() + print_radial() + print_coord(abeta));     // this is just a header for dump_Jac file data
                out.print("a = np.array([");
                for (i = 0; i < Jac.length; i++)
                {
                    if (i > 0) out.print(", ");
                    out.print("[");
                    for (j = 0; j < Jac.length; j++)
                    {
                        if (j > 0) out.print(", ");
                        out.print(Jac[i][j]);
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
            //if (out != null)                                    // use only for replicating file 'betarun8.csv'
            //    out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + print_coord(abeta) + ", " + (float) retVal + ", " + (float) Jacdet);
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
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + print_coord(abeta) + ", " + (float) retVal + ", " + (float) Jacdet);
        return retVal;
    }

    private static void get_Bezier_endpoints_from_Beta_Spline()
    {
        fitted = new epiTrochoidFxn(6.27);              // fix fix temporary location
        //String str = "gauss t2[] @ , 90.0, 90.0, 17.1, 800, , 181.11768, 191.3406, 168.46349, 191.3981, 181.11768, 191.3406, 168.46349, 191.3981, -73.178055, -8.765734, 29.248032, 81.26961, 106.821945, 171.23427, -150.75197, -98.73039, 4.3080115, 1.7127049, 0.0033888216, 1.1912025, 4.3080115, 1.7127049, 0.0033888216, 1.1912025, 6.869E-4, NaN";
        // convert a 8-point Beta-Spline to 8 spliced cubic Beziers (endpoints only)
        double[] beta_in = new double[] {
//185.39139, 193.71384, 169.07936, 193.46742, 185.39139, 193.71384, 169.07936, 193.46742, -76.67592, -7.854524, 25.704296, 81.88335, 103.32408, 172.14548, -154.2957, -98.11665, 9.370194, 3.2129233, 1.2039822, 2.505883, 9.370194, 3.2129233, 1.2039822, 2.505883
//175.51227, 184.37256, 175.51227, 184.37256, 175.51227, 184.37256, 175.51227, 184.37256, -45.0, -6.140414E-11, 45.0, 90.0, 135.0, 180.0, -135.0, -90.0, -2.3176925, 22.393124, -2.3176925, 22.393124, -2.3176925, 22.393124, -2.3176925, 22.393124
//176.4047, 183.60616, 176.4047, 183.60616, 176.4047, 183.60616, 176.4047, 183.60616, -45.0, -2.6008996E-13, 45.0, 90.0, 135.0, -180.0, -135.0, -90.0, 3.6070511, -0.8296387, 3.6070511, -0.8296387, 3.6070511, -0.8296387, 3.6070511, -0.8296387
//175.51227, 184.37256, 175.51227, 184.37256, 175.51227, 184.37256, 175.51227, 184.37256, -45.0, -6.140414E-11, 45.0, 90.0, 135.0, 180.0, -135.0, -90.0, -2.3176925, 22.393124, -2.3176925, 22.393124, -2.3176925, 22.393124, -2.3176925, 22.393124
//187.07251, 187.07251, 187.07251, 187.07251, 187.07251, 187.07251, 187.07251, 187.07251, -77.75456, -12.245439, 12.245439, 77.75456, 102.24544, 167.75456, -167.75456, -102.24544, 4.0902343, 4.0902343, 4.0902343, 4.0902343, 4.0902343, 4.0902343, 4.0902343, 4.0902343
//168.08003, 197.48953, 168.08003, 197.48953, 168.08003, 197.48953, 168.08003, 197.48953, -63.464863, -4.8974257, 26.535137, 85.10258, 116.53513, 175.10257, -153.46486, -94.89742, 4.023081, 2.4261758, 4.023081, 2.4261758, 4.023081, 2.4261758, 4.023081, 2.4261758
174.02264, 186.3201, 174.02264, 186.29512, 174.02264, 186.3201, 174.02264, 186.29512, -50.11968, 1.6953382E-13, 50.11968, 90.0, 129.88033, 180.0, -129.88033, -90.0, 2.1633155, -0.38572678, 2.1633155, -0.35647652, 2.1633155, -0.38572678, 2.1633155, -0.35647652
//179.85231, 180.1501, 179.85231, 180.1501, 179.85231, 180.1501, 179.85231, 180.1501, -62.850533, -17.501976, 27.149467, 72.498024, 117.14947, 162.49802, -152.85052, -107.501976, -0.13609967, -0.18932837, -0.13609967, -0.18932837, -0.13609967, -0.18932837, -0.13609967, -0.18932837
//159.98505, 199.99086, 159.98505, 199.99086, 159.98505, 199.99086, 159.98505, 199.99086, -45.0, -1.5553355E-14, 45.0, 90.0, 135.0, 180.0, -135.0, -90.0, -2.1210384, -0.12768523, -2.1210384, -0.12768523, -2.1210384, -0.12768523, -2.1210384, -0.12768523
        };
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
        double incr = 0; //20; // magic factor = 72;
        //double theta = 120;
        //double[] eig_vec_angular0 = new double[] {0.27658725,0.21717687,-0.21717687,-0.27658725,-0.27658725,-0.21717687,0.21717687,0.27658725,0.21717687,0.27658725,0.27658725,0.21717687,-0.21717687,-0.27658725,-0.27658725,-0.21717687,0.03651984,-0.03651984,0.03651984,-0.03651984,0.03651984,-0.03651984,0.03651984,-0.03651984};
        //double[] eig_vec_angular1 = new double[] {0.19880945,0.23317641,-0.26448369,0.2548129,0.19880945,0.23317641,-0.26448369,0.2548129,0.15533342,0.29835982,0.33826469,-0.19927081,0.15533342,0.29835982,0.33826469,-0.19927081,0.04314422,-0.03735079,-0.01050166,-0.02401306,-0.04314422,0.03735079,0.01050166,0.02401306};
        //double[] eig_vec_angular2 = new double[] {0.33826469,-0.19927081,0.15533342,0.29835982,0.33826469,-0.19927081,0.15533342,0.29835982,0.26448369,-0.2548129,-0.19880945,-0.23317641,0.26448369,-0.2548129,-0.19880945,-0.23317641,-0.01050166,-0.02401306,-0.04314422,0.03735079,0.01050166,0.02401306,0.04314422,-0.03735079};
        double[] eig_vec_angular = new double[] {-0.35355294812152965 , 0.0004545543395921281 , -0.35355294800178955 , 1.9828270188954056e-10 , 0.3535529486267195 , -0.0004545543401452256 , 0.3535529488062767 , 2.40356135602221e-10 , -0.3535529481219878 , -4.71939709534297e-11 , 0.35355294800200066 , -0.00045455433982321885 , 0.3535529486263335 , -6.760268430119898e-11 , -0.35355294880602967 , 0.00045455433992618727 , -4.4906411922340794e-11 , -0.0006471098200360459 , -3.340369647553132e-11 , 0.0006471098235461131 , 4.2019007739403647e-11 , -0.0006471098275880648 , 2.6353145360369368e-11 , 0.0006471098249856699};
        //double[] eig_vec_radial  = new double[] {0.31665681,0.15725157,-0.15725157,-0.31665681,-0.31665681,-0.15725157,0.15725157,0.31665681,0.15725157,0.31665681,0.31665681,0.15725157,-0.15725157,-0.31665681,-0.31665681,-0.15725157,0.00063766,-0.00063766,0.00063766,-0.00063766,0.00063766,-0.00063766,0.00063766,-0.00063766};
        for (int i = 0; i < ax.length; i++)                     // increment along the lowest eigenvector
        {
            //ax[i] += incr*(Math.cos(theta*Math.PI/180)*eig_vec_angular2[i] + Math.sin(theta*Math.PI/180)*eig_vec_angular1[i]);
            //ay[i] += incr*(Math.cos(theta*Math.PI/180)*eig_vec_angular2[i + ax.length] + Math.sin(theta*Math.PI/180)*eig_vec_angular1[i + ax.length]);
            //abeta[i] += incr*(Math.cos(theta*Math.PI/180)*eig_vec_angular2[i + 2*ax.length] + Math.sin(theta*Math.PI/180)*eig_vec_angular1[i + 2*ax.length]);
            ax[i] += incr*eig_vec_angular[i];
            ay[i] += incr*eig_vec_angular[i + ax.length];
            abeta[i] += incr*eig_vec_angular[i + 2*ax.length];
        }
        //double[] cart_in1 = new double[] {63.521793, 168.62363, 168.62363, 63.521793, -63.521793, -168.62363, -168.62363, -63.521793, -168.62363, -63.521793, 63.521793, 168.62363, 168.62363, 63.521793, -63.521793, -168.62363, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702, -0.08902702};
        //double[] cart_in2 = new double[] {124.05629, 184.44278, 124.05629, -2.2495281E-13, -124.05629, -184.44278, -124.05629, -1.8626032E-13, -124.05629, -7.1329E-13, 124.05629, 184.44278, 124.05629, -3.3561272E-13, -124.05629, -184.44278, -2.3263218, 22.012825, -2.3263218, 22.012825, -2.3263218, 22.012825, -2.3263218, 22.012825};
        //double[] cart_in = new double[] {42.725063, 191.89647, 152.34802, 27.31546, -42.725063, -191.89647, -152.34802, -27.31546, -180.40105, -26.47261, 73.33422, 191.52939, 180.40105, 26.47261, -73.33422, -191.52939, 9.370194, 3.2129233, 1.2039822, 2.505883, 9.370194, 3.2129233, 1.2039822, 2.505883};
        //double[] cart_in = new double[] {159.231966,186.29883,159.231966,-1.30592E-12,-159.231966,-186.29883,-159.231966,1.05438E-12,-85.74218,3.75189E-12,85.74218,186.37727,85.74218,2.85274E-12,-85.74218,-186.37727,2.2298417,-0.31485569,2.2298417,-0.4207099,2.2298417,-0.31485569,2.2298417,-0.4207099};
        //double[] cart_in = new double[] {26.47261,180.40105,191.52939,73.33422,-26.47261,-180.40105,-191.52939,-73.33422,-191.89647,-42.725063,27.31546,152.34802,191.89647,42.725063,-27.31546,-152.34802,3.2129233,9.370194,2.505883,1.2039822,3.2129233,9.370194,2.505883,1.2039822};
        //double mix = 1.;                                               // use only for interpolation
        for (int i = 0; i < ax.length; i++)                             // interpolate between 2 initial estimates
        {
        //    int j = (i + ax.length) % ax.length;
        //    int k = (i + ax.length) % ax.length;
            //ax[i] = mix*cart_in1[j] + (1 - mix)*cart_in2[k];
            //ay[i] = mix*cart_in1[j + ax.length] + (1 - mix)*cart_in2[k + ax.length];
            //abeta[i] = mix*cart_in1[j + 2*ax.length] + (1 - mix)*cart_in2[k + 2*ax.length];
            //ax[i] = cart_in[i];
            //ay[i] = cart_in[i + ax.length];
            //abeta[i] = cart_in[i + 2*ax.length];
        }
        //ax[6] += 0.01;                                                  // increment for dFdai calculation
        //ay[6] += 0.01;                                                  // increment for dFdai calculation
        //abeta[1] += 0.001;                                               // increment for dFdai calculation
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

    private static void write_Beta2_SplineN_data()
    {
        // generate a file consisting of Jacobian matrices for egenvalue calculation
        try
        {
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\eig_augment_raw_8.txt", false);
            FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\eig_augment_temp.txt", false);
            //FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\betarun8.csv", false);
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
            BufferedReader instr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\Beta2_Spline8\\betarun8.csv"));
            try
            {
                while (instr.ready()) // && line_in < 50)
                {
                    line_in++;
                    str_in = instr.readLine();
                    if (line_in >= 1505 && line_in <= 2000)
                    {
                        if (str_in.isEmpty())
                            out.println("....................");
                        else if (str_in.startsWith("gauss t2[] @ ,"))
                        {
                            str_split = str_in.split(",");
                            //out.println("line_in = " + line_in + ", " + str_split[3] + ", " + str_in);
                            fitted = new epiTrochoidFxn(Double.parseDouble(str_split[3]));  // read c value
                            for (int i = 0; i < beta_in.length; i++)                        // read input parameters
                                beta_in[i] = Double.parseDouble(str_split[i + 6]);
                            for (int i = 0; i < ax.length; i++)
                            {
                                ax[i] = beta_in[i]*Math.cos(beta_in[i + ax.length]*Math.PI/180);
                                ay[i] = beta_in[i]*Math.sin(beta_in[i + ax.length]*Math.PI/180);
                                abeta[i] = beta_in[i + 2*ax.length];
                            }
                            out.print(line_in + " - ");
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
        fitted = new epiTrochoidFxn(4.43);                       // fix fix temporary location
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
            str += ", " + (float) c;
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
        return Math.sqrt(integrate(trap_in))/(fitted.a + fitted.b);
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
}
