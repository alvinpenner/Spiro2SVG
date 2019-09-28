
package components;

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
    private static double[][] t2ddx = new double[8][N+1];       // partial of u wrt {xj}
    private static double[][] t2ddy = new double[8][N+1];       // partial of u wrt {yj}
    private static double[][] t2ddbeta = new double[8][N+1];    // partial of u wrt {betaj}
    private static double[][] dddx = new double[8][8];          // see Book 9, p. 54
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        fitted = new epiTrochoidFxn(7.7);
        get_Bezier_endpoints();
        convert_to_Bezier();
        //gen_spline_points();
        //System.out.println("Beta2_SplineN iterate_at_x_y = " + iterate_at_x_y() + "\n");       // normally include this line
        System.out.println("Beta2_SplineN solve_at_x_y = " + solve_at_x_y(true) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //for (int i = 0; i <= 100; i++)
        //    System.out.println(i + ", " + Bi3(i/100.0)[3] + ", " + dBi3(i/100.0)[3] + ", " + d2Bi3(i/100.0)[3] + ", " + d3Bi3(i/100.0)[3]);
    }
/*
    private static double iterate_at_x_y()
    {
        // calculate a new estimate of (xi, yi) by setting dF = 0
        // include only first-order responses
        // setup 2N-variable Newton-Raphson iteration

        final double gain = 1;                                  // factor to reduce gain
        final int MAXLOOP = 500;
        double[] f_gx = new double[N];
        double[] f_gy = new double[N];
        double[] dfxdu = new double[N];
        double[] dfydu = new double[N];
        double[] denom = new double[N];                         // = E(u)

        double[] trap_in = new double[N];
        double[][] Jac = new double[2*Splinex.length][2*Splinex.length];
        double[] dFda = new double[2*Splinex.length];
        double[] d2Fdadc = new double[2*Splinex.length];        // augmented matrix
        double dFdc;
        double d2Fdcdc;                                         // augmented matrix
        double[][] Augment = new double[2*Splinex.length + 1][2*Splinex.length + 1];
        double[] dela;                                          // (-Δxi, -Δyi)
        int loop = 0;
        int i, j, k;
        boolean outside;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_x_y(false)))              // initiallize at (xi, yi)
            {
                System.out.println("fail at " + print_coord(Splinex) + print_coord(Spliney));
                return Double.NaN;
            }

            for (i = 0; i < N; i++)
            {
                f_gx[i] = BSpline5.multvv(permute(Splinex, t2[i]), Ni3(t2[i] % 1)) - fitted.getx(t1_start + i*(t1_end - t1_start)/N);
                f_gy[i] = BSpline5.multvv(permute(Spliney, t2[i]), Ni3(t2[i] % 1)) - fitted.gety(t1_start + i*(t1_end - t1_start)/N);
                dfxdu[i] = BSpline5.multvv(permute(Splinex, t2[i]), dNi3(t2[i] % 1));
                dfydu[i] = BSpline5.multvv(permute(Spliney, t2[i]), dNi3(t2[i] % 1));
                denom[i] = BSpline5.multvv(permute(Splinex, t2[i]), dNi3(t2[i] % 1))*BSpline5.multvv(permute(Splinex, t2[i]), dNi3(t2[i] % 1))
                         + f_gx[i]*BSpline5.multvv(permute(Splinex, t2[i]), d2Ni3(t2[i] % 1))
                         + BSpline5.multvv(permute(Spliney, t2[i]), dNi3(t2[i] % 1))*BSpline5.multvv(permute(Spliney, t2[i]), dNi3(t2[i] % 1))
                         + f_gy[i]*BSpline5.multvv(permute(Spliney, t2[i]), d2Ni3(t2[i] % 1));
                //System.out.println(i + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + denom[i]);
            }

            // calc dFda[j] at current (xi, yi)

            for (i = 0; i < 2*Splinex.length; i++)
            {
                for (k = 0; k < N; k++)
                    if (i < Splinex.length)             // x variable
                        trap_in[k] = f_gx[k]*calc_dfdai(t2[k], i)                       // original code
                                   + (f_gx[k]*dfxdu[k] + f_gy[k]*dfydu[k])*t2ddx[i][k];
                        //trap_in[k] = f_gx[k]*calc_dfdai(t2[k], i);                    // new code
                    else                                // y variable
                        trap_in[k] = f_gy[k]*calc_dfdai(t2[k], i - Splinex.length)      // original code
                                   + (f_gx[k]*dfxdu[k] + f_gy[k]*dfydu[k])*t2ddy[i - Splinex.length][k];
                        //trap_in[k] = f_gy[k]*calc_dfdai(t2[k], i - Splinex.length);   // new code
                dFda[i] = integrate(trap_in);
                //System.out.println("dFda ," + i + ", " + dFda[i]);
                for (k = 0; k < N; k++)
                    if (i < Splinex.length)             // x variable
                        trap_in[k] = -calc_dfdai(t2[k], i)*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)
                                   - denom[k]*t2ddx[i][k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
                    else
                        trap_in[k] = -calc_dfdai(t2[k], i - Splinex.length)*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)
                                   - denom[k]*t2ddy[i - Splinex.length][k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
                d2Fdadc[i] = integrate(trap_in);                       // augmented matrix
            }

            // calc d2Fda[i]da[j] (symmetric Jacobean matrix)

            for (i = 0; i < 2*Splinex.length; i++)
                for (j = 0; j < 2*Splinex.length; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k < N; k++)
                    {
                        if (i < Splinex.length)
                            if (j < Splinex.length)
                                trap_in[k] = calc_dfdai(t2[k], i)*calc_dfdai(t2[k], j)
                                           - denom[k]*t2ddx[i][k]*t2ddx[j][k];
                            else
                                trap_in[k] = -denom[k]*t2ddx[i][k]*t2ddy[j - Splinex.length][k];
                        else if (j < Splinex.length)
                                trap_in[k] = -denom[k]*t2ddy[i - Splinex.length][k]*t2ddx[j][k];
                            else
                                trap_in[k] = calc_dfdai(t2[k], i - Splinex.length)*calc_dfdai(t2[k], j - Splinex.length)
                                           - denom[k]*t2ddy[i - Splinex.length][k]*t2ddy[j - Splinex.length][k];
                        //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = integrate(trap_in);
                }
            //if (true) return Double.NaN;                        // use only for scan_steepest()
            //int checki = 1;
            //int checkj = 7;
            //System.out.println("check M ," + checki + ", " + checkj + ", " + dFda[checki] + ", " + dFda[checkj] + ", " + Jac[checki][checkj]);

            // calculate determinant of augmented matrix

            for (k = 0; k < N; k++)
                trap_in[k] = -f_gx[k]*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N) - f_gy[k]*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N);
            dFdc = integrate(trap_in);
            //System.out.println("check Aug ," + "c" + ", " + checkj + ", " + dFdc + ", " + dFda[checkj] + ", " + d2Fdadc[checkj]);
            for (k = 0; k < N; k++)
                trap_in[k] = fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)*fitted.getdxdc(t1_start + k*(t1_end - t1_start)/N)
                           + fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)*fitted.getdydc(t1_start + k*(t1_end - t1_start)/N)
                           - denom[k]*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c")*calc_t2dd(Integer.MAX_VALUE, k, t2[k], "c");
            d2Fdcdc = integrate(trap_in);
            System.out.println("check d2Fdc2 , " + ", " + dFdc + ", " + d2Fdcdc);
            for (i = 0; i < 2*Splinex.length; i++)
                for (j = 0; j < 2*Splinex.length; j++)
                    Augment[i][j] = Jac[i][j];
            for (i = 0; i < 2*Splinex.length; i++)
            {
                Augment[i][2*Splinex.length] = d2Fdadc[i];
                Augment[2*Splinex.length][i] = Augment[i][2*Splinex.length];
            }
            Augment[2*Splinex.length][2*Splinex.length] = d2Fdcdc;
            System.out.println("dFdc = " + fitted.getc() + print_coord(Splinex) + print_coord(Spliney) + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);

            dela = BSpline5.gaussj(Jac, dFda);          // this is actually the negative of Δa
            for (i = 0; i < Splinex.length; i++)
            {
                Splinex[i] -= dela[i]/gain;             // gain is just a fudge factor to 'improve' convergence
                Spliney[i] -= dela[i + Splinex.length]/gain;
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
            if (Math.abs(t2[N] - t2[0] - Splinex.length) > TOL
            ||  t2[0] < -100*TOL)
            {
                System.out.println("bad converged t2[i]: ABORT");
                return Double.NaN;
            }
            System.out.println("\n__converged in " + loop + " at new xi yi = , , , , , , " + print_coord(Splinex) + print_coord(Spliney));
            //Jacdet = BSpline5.detm(Jac);
            System.out.println("dFdc = " + fitted.getc() + print_radial());     // this is just a header for dump_Jac file data
            BSpline5.dump_Jac(Jac);
            BSpline5.dump_Jac(Augment);
            //write_BSplineN_data(null);                      // use only for dumping to file from dump_Jac()
            //write_BSplineN_data(Jac);
            //write_BSplineN_data(Augment);
            //for (i = 0; i < N; i++)                         // test sum of t2ddx and t2ddy
            //    System.out.println(i + ", " + dfxdu[i] + ", " + dfydu[i] + ", " + denom[i]);
            double[][] dot_product = new double [N][8];     // test code just for oval test
            j = 4;                                          // segment index (0 - Splinex.length)
            int l = 0;
            for (k = 0; k < N; k++)
                if ((t2[k] > j) && (t2[k] < j + 1))
                {
                    for (i = 0; i < 4; i++)                 // only need 4 control pts
                    {
                        dot_product[l][i] = f_gx[k]*calc_dfdai(t2[k], (j + i) % Splinex.length);
                        dot_product[l][i + 4] = f_gy[k]*calc_dfdai(t2[k], (j + i) % Splinex.length);
                    }
                    l++;
                }
            //System.out.println("dot_product:");
            //BSpline5.dump_Jac(dot_product);
            //correlate_f_dfda(dot_product, j);                                  // end of test code
            return solve_at_x_y(true);                              // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + print_coord(dela) + ")");
        return Double.NaN;
    }
*/
    private static double solve_at_x_y(boolean print)
    {
        // 8-point, uniform, closed, cubic Beta2-Spline
        // perform a single calculation of a complete t2[] profile
        // at a given set {x, y, beta2}, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        if (t2[N] == 0)
            System.out.println("__start Beta2-Spline8 theta c t x y = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + print_coord(ax) + print_coord(ay) + print_coord(abeta));
        else
        {
            System.out.println("__solve at new d1 d2 rms = , , , , , " + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + calc_error());
            System.out.println("__solve at new d1 d2 rms = , , , , , " + print_radial() + print_coord(abeta) + ", " + calc_error());
        }

        if (print) System.out.println("\n , t1, t2, t2ddx0, t2ddx1, t2ddx2, t2ddx3, t2ddx4, t2ddx5, t2ddx6, t2ddx7");
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i);
            if ((i == N && Math.abs(t2[N] - t2[0] - ax.length) > TOL)
            ||  (t2[i] < -5000*TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                return Double.NaN;
            }
//            for (int j = 0; j < ax.length; j++)
//                t2ddx[j][i] = calc_t2dd(j, i, t2[i], "x");
//            for (int j = 0; j < ay.length; j++)
//                t2ddy[j][i] = calc_t2dd(j, i, t2[i], "y");
            if (print)
                System.out.println(i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]); // + print_coord(t2ddx, i) + print_coord(t2ddy, i));
                //System.out.println(i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + calc_t2dd(Integer.MAX_VALUE, i, t2[i], "c"));
        }
        double retVal = calc_error();
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(ax) + print_coord(ay) + print_coord(abeta) + ", " + (float) retVal + ", " + (float) Jacdet);
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(Splinex) + print_coord(Spliney) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2);
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + ", " + (float) retVal + ", " + (float) Jacdet);
        return retVal;
    }

    private static void get_Bezier_endpoints()
    {
        // convert a 8-point B-Spline to 8 spliced cubic Beziers (endpoints only)
        double[] rad = new double[] {181.62039, 217.36208, 181.62039, 217.36208, 181.62039, 217.36208, 181.62039, 217.36208, -47.694225, -1.0211561, 42.305775, 88.97884, 132.30577, 178.97885, -137.69423, -91.02116};
        int j;
        //rad = new double[] {135.0592463,219.317122,192.3538406,259.4224354,254.6998233,72.80109889,212.6029163,121.6552506,44.70002484,24.22774532,8.972626615,62.44718842,46.90915243,164.0546041,48.81407483,9.462322208};
        if (rad.length != 16)
        {
            System.out.println("convert_to_Bezier: incorrect length = " + rad.length);
            return;
        }
        System.out.println("convert_to_Bezier radial at c = ," + fitted.getc() + print_coord(rad));
        abeta = new double[] {0, 0, 0, 0, 0, 0, 0, 0};  // beta2 values of Bezier endpoints
        //abeta = new double[] {-1, 1, -1, 1, -1, 1, -1, 1};
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
        System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta));
    }

    private static void convert_to_Bezier()
    {
        int i, j;
        // see Spiro2SVG Book 9, p. 53
        double[][] convert = new double[8][8];
        double[] vec = new double[8];
        double[] lhs;
        double[] d;
        //int k = 0;                                  // index of coord to change - test only
        //ay[k] += 0.01;                              // increment Bezier coord   - test only
        //System.out.println("Bezier endpoints at c = ," + fitted.getc() + print_coord(ax) + print_coord(ay) + print_coord(abeta));
        //System.out.println(" , " + k + ": " + .00);
        for (i = 0; i < convert.length; i++)
            for (j = 0; j < convert.length; j++)
                convert[i][j] = 0;
        for (i = 0; i < convert.length; i++)
        {
            convert[i][i] = 1;
            convert[i][(i + 1) % convert.length] = 4 + abeta[i]/2;
            convert[i][(i + 2) % convert.length] = 1;
        }
//        System.out.println("convert[i][j]:");
//        for (i = 0; i < convert.length; i++)
//        {
//            for (j = 0; j < convert.length; j++)
//                System.out.print(convert[i][j] + " ,");
//            System.out.println();
//        }
        for (i = 0; i < convert.length; i++)
            vec[i] = ax[(i + 2) % convert.length] - ax[i];
        d = BSpline5.gaussj(convert, vec);
//        for (i = 0; i < convert.length; i++)
//            System.out.println(i + ", " + d[i]);
        for (i = 0; i < convert.length; i++)
        {
            Bezx[3*i] = ax[i];
            Bezx[3*i + 1] = ax[i] + d[i];
            Bezx[3*i + 2] = ax[(i + 1) % convert.length] - d[(i + 1) % convert.length];
        }
        for (i = 0; i < convert.length; i++)
            vec[i] = ay[(i + 2) % convert.length] - ay[i];
        d = BSpline5.gaussj(convert, vec);
        for (i = 0; i < convert.length; i++)
        {
            Bezy[3*i] = ay[i];
            Bezy[3*i + 1] = ay[i] + d[i];
            Bezy[3*i + 2] = ay[(i + 1) % convert.length] - d[(i + 1) % convert.length];
        }
        for (i = 0; i < convert.length; i++)                // calculate dddx
        {
            for (j = 0; j < convert.length; j++)
                vec[j] = 0;
            vec[(i + vec.length - 2) % vec.length] = 1;
            vec[i] = -1;
            lhs = BSpline5.gaussj(convert, vec);
            for (j = 0; j < convert.length; j++)
                dddx[j][i] = lhs[j];
        }
        //System.out.println("dddx[i][j]:");
        //for (i = 0; i < convert.length; i++)
        //{
        //    for (j = 0; j < convert.length; j++)
        //        System.out.print(dddx[i][j] + " ,");
        //    System.out.println();
        //}
        //System.out.print("M " + Bezx[0] + ", " + Bezy[0] + " C");
        //for (i = 1; i < 3*convert.length; i++)
        //    System.out.print(" " + Bezx[i] + ", " + Bezy[i]);
        //System.out.println(" " + Bezx[0] + ", " + Bezy[0]);
        //for (i = 0; i < N; i++)
        //{
        //    double di = i/100.0;
        //    System.out.println(i + ", " + multvv(Bezy, (int) di, Bi3(di % 1)) + ", " + calc_dfdaxi(k, di));
        //}
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

    private static String print_radial()
    {
        String str = "";
        for (int i = 0; i < ax.length; i++)
            str += ", " + (double) Math.sqrt(ax[i]*ax[i] + ay[i]*ay[i]);
        for (int i = 0; i < ax.length; i++)
            str += ", " + (double) (Math.atan2(ay[i], ax[i])*180/Math.PI);
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

        if (i == 0) t = 0;
        else t = t2[i-1];
        f = (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, dBi3(t % 1))
          + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, dBi3(t % 1));
        fprime = multvv(Bezx, (int) t, dBi3(t % 1))*multvv(Bezx, (int) t, dBi3(t % 1)) + (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, d2Bi3(t % 1))
               + multvv(Bezy, (int) t, dBi3(t % 1))*multvv(Bezy, (int) t, dBi3(t % 1)) + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, d2Bi3(t % 1));
        f2prime = 3*multvv(Bezx, (int) t, dBi3(t % 1))*multvv(Bezx, (int) t, d2Bi3(t % 1)) + (multvv(Bezx, (int) t, Bi3(t % 1)) - X)*multvv(Bezx, (int) t, d3Bi3(t % 1))
                + 3*multvv(Bezy, (int) t, dBi3(t % 1))*multvv(Bezy, (int) t, d2Bi3(t % 1)) + (multvv(Bezy, (int) t, Bi3(t % 1)) - Y)*multvv(Bezy, (int) t, d3Bi3(t % 1));
        //System.out.println("\ninit i t1 t2 =, " + i + ", " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
            del_t = -fprime/f2prime;
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            if (del_t < -5000*TOL || del_t > 3)
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
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
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
/*
    private static double calc_t2dd(int j, int i, double u, String type)
    {
        // i is the time index for t1
        // j is the coefficient index for aj
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        //System.out.println("calc_t2dd = " + j + ", " + i + ", " + u + ", " + type);
        double denom = BSpline5.multvv(permute(Splinex, u), dNi3(u % 1))*BSpline5.multvv(permute(Splinex, u), dNi3(u % 1)) + (BSpline5.multvv(permute(Splinex, u), Ni3(u % 1)) - X)*BSpline5.multvv(permute(Splinex, u), d2Ni3(u % 1))
                     + BSpline5.multvv(permute(Spliney, u), dNi3(u % 1))*BSpline5.multvv(permute(Spliney, u), dNi3(u % 1)) + (BSpline5.multvv(permute(Spliney, u), Ni3(u % 1)) - Y)*BSpline5.multvv(permute(Spliney, u), d2Ni3(u % 1));
        double numer = Double.NaN;
        if (type.equals("x"))
            numer = BSpline5.multvv(permute(Splinex, u), dNi3(u % 1))*calc_dfdai(u, j)
                  + (BSpline5.multvv(permute(Splinex, u), Ni3(u % 1)) - X)*calc_d2fdudai(u, j);
        else if (type.equals("y"))
            numer = BSpline5.multvv(permute(Spliney, u), dNi3(u % 1))*calc_dfdai(u, j)
                  + (BSpline5.multvv(permute(Spliney, u), Ni3(u % 1)) - Y)*calc_d2fdudai(u, j);
        else if (type.equals("c"))
            numer = - BSpline5.multvv(permute(Splinex, u), dNi3(u % 1))*fitted.getdxdc(t1)
                    - BSpline5.multvv(permute(Spliney, u), dNi3(u % 1))*fitted.getdydc(t1);
        return -numer/denom;
    }
*/
    private static double calc_dfdaxi(int i, double u)
    {
        // this works for both calc_dfdaxi and calc_dfdayi
        double ret = 0;
        int ui = (int) u;
        if (ui < 0) return Double.NaN;
        if (ui > 7) return Double.NaN;
        if (ui == i)
            ret = Bi3(u % 1)[0] + Bi3(u % 1)[1];
        else if ((ui + 1) % 8 == i)
            ret = Bi3(u % 1)[2] + Bi3(u % 1)[3];
        ret += dddx[ui][i]*Bi3(u % 1)[1];
        ret -= dddx[(ui + 1) % 8][i]*Bi3(u % 1)[2];
        return ret;
    }
/*
    private static double calc_d2fdudai(double u, int i)
    {
        int istart = (int) u;
        if ((i + Splinex.length - istart) % Splinex.length > 3)
            return 0;
        else
            return dNi3(u % 1)[(i + Splinex.length - istart) % Splinex.length];
    }
*/
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

        for (int i = 0; i < v2.length; i++)
            retVal += v1[(3*v1Pos + i) % v1.length]*v2[i];
        return retVal;
    }
}
