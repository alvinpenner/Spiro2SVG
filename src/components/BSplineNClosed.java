
package components;

import java.io.*;
//import org.jblas.*;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a closed N-point B-Spline (P0 -> PN-1) to it, using parameter 0 < t2 < N.
// with no contraints on the control points or endpoints.
// linearize the equations wrt Pi and solve a 2Nx2N system of equations.
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book9, May 2019, page 1

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BSplineNClosed.java

public class BSplineNClosed
{
    public static final double t1_start = 0;
    public static final double t1_end = 2*Math.PI;
    public static final int N = 800;                // normally 100*number of segments
    public static double[] Splinex, Spliney;        // N-point spline
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2ddx;                // partial of u wrt {xj}
    private static double[][] t2ddy;                // partial of u wrt {yj}
    private static double[][] Jac;                  // defined here only so we can share it with 'scan_steepest'
    private static double[] dFda;
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        //fitted = new epiTrochoidFxn(8.45);
        fitted = new epiTrochoidFxn(7.7);
        //double rad0 = 195.5;
        //double rad1 = 195;
        //double rad2 = 195.5;
        //double rad3 = 182;
        //double rad4 = 182;
        //double rad5 = 180;
        //double[] rad = new double[] {rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2};
        //double[] rad = new double[] {rad1, rad2, rad3, rad1, rad2, rad3, rad1, rad2, rad3, rad1, rad2, rad3};
        //double[] rad = new double[] {rad0, rad1, rad2, rad3, rad4, rad0, rad1, rad2, rad3, rad4, rad0, rad1, rad2, rad3, rad4, rad0, rad1, rad2, rad3, rad4};
        //double[] rad = new double[] {rad0, rad1, rad2, rad0, rad1, rad2, rad0, rad1, rad2, rad0, rad1, rad2};
        //double[] rad = new double[] {rad0, rad1, rad2, rad3, rad4, rad5, rad0, rad1, rad2, rad3, rad4, rad5, rad0, rad1, rad2, rad3, rad4, rad5, rad0, rad1, rad2, rad3, rad4, rad5};
        double[] rin = new double[] {181.62039, 217.36208};
        double[] thin = new double[] {-47.694225, -1.0211561};
        //double[] thin = new double[] {-76.79729772, -17.43070974};
        //double[] rin = new double[] {174.84395, 200.71751, 192.32411};
        //double[] thin = new double[] {-33.770927, -5.76914, 15.113813};
        //double[] rin = new double[] {195.9031, 202.15785, 195.9031, 178.2428, 160.45967, 178.2428};
        //double[] thin = new double[] {-7.6173887, 5.1160803E-12, 7.6173887, 19.259514, 45.0, 70.740486};
        double[] rad = new double[] {rin[0], rin[1], rin[0], rin[1], rin[0], rin[1], rin[0], rin[1]};
        //rad = new double[] {rin[0], rin[1], rin[2], rin[0], rin[1], rin[2], rin[0], rin[1], rin[2], rin[0], rin[1], rin[2]};
        //rad = new double[] {rin[0], rin[1], rin[2], rin[3], rin[4], rin[0], rin[1], rin[2], rin[3], rin[4], rin[0], rin[1], rin[2], rin[3], rin[4], rin[0], rin[1], rin[2], rin[3], rin[4]};
        //rad = new double[] {rin[0], rin[1], rin[2], rin[3], rin[4], rin[5], rin[0], rin[1], rin[2], rin[3], rin[4], rin[5], rin[0], rin[1], rin[2], rin[3], rin[4], rin[5], rin[0], rin[1], rin[2], rin[3], rin[4], rin[5]};
        //double th0 =  -5.8;
        //double th1 =    0;
        //double th2 =   5.8;
        //double th3 =   20.5;
        //double th4 =   69.5;
        //double th5 =   60;
        //double[] theta = new double[] {th1, th2, th1 + 90, th2 + 90, th1 + 180, th2 + 180, th1 + 270, th2 + 270};
        //double[] theta = new double[] {th1, th2, th3, th1 + 90, th2 + 90, th3 + 90, th1 + 180, th2 + 180, th3 + 180, th1 + 270, th2 + 270, th3 + 270};
        //double[] theta = new double[] {th0, th1, th2, th3, th4, th0 + 90, th1 + 90, th2 + 90, th3 + 90, th4 + 90, th0 + 180, th1 + 180, th2 + 180, th3 + 180, th4 + 180, th0 + 270, th1 + 270, th2 + 270, th3 + 270, th4 + 270};
        //double[] theta = new double[] {th0, th1, th2, th3, th4, th5, th0 + 90, th1 + 90, th2 + 90, th3 + 90, th4 + 90, th5 + 90, th0 + 180, th1 + 180, th2 + 180, th3 + 180, th4 + 180, th5 + 180, th0 + 270, th1 + 270, th2 + 270, th3 + 270, th4 + 270, th5 + 270};
        double[] theta = new double[] {thin[0], thin[1], thin[0] + 90, thin[1] + 90, thin[0] + 180, thin[1] + 180, thin[0] + 270, thin[1] + 270};
        //theta = new double[] {thin[0], thin[1], thin[2], thin[0] + 90, thin[1] + 90, thin[2] + 90, thin[0] + 180, thin[1] + 180, thin[2] + 180, thin[0] + 270, thin[1] + 270, thin[2] + 270};
        //theta = new double[] {thin[0], thin[1], thin[2], thin[3], thin[4], thin[5], thin[0] + 90, thin[1] + 90, thin[2] + 90, thin[3] + 90, thin[4] + 90, thin[5] + 90, thin[0] + 180, thin[1] + 180, thin[2] + 180, thin[3] + 180, thin[4] + 180, thin[5] + 180, thin[0] + 270, thin[1] + 270, thin[2] + 270, thin[3] + 270, thin[4] + 270, thin[5] + 270};
        Splinex = new double[rad.length];
        Spliney = new double[rad.length];
        dFda = new double[2*Splinex.length];
        Jac = new double[2*Splinex.length][2*Splinex.length];
        for (int i = 0; i < Splinex.length; i++)              // normally include these four lines
            Splinex[i] = rad[i]*Math.cos(Math.PI*theta[i]/180 - 0.0*2*Math.PI/Splinex.length);
        for (int i = 0; i < Splinex.length; i++)
            Spliney[i] = rad[i]*Math.sin(Math.PI*theta[i]/180 - 0.0*2*Math.PI/Splinex.length);
        //for (int i = 0; i < Splinex.length; i++)
        //    Spliney[i] += 5;                                    // fix fix tremporary code
        //double rad = 200; // 194.0*Math.sqrt(2);
        //Splinex = new double[] {131.49011, 213.08975, 131.49011, -2.5748248E-9, -131.49011, -213.08975, -131.49011, -1.5923756E-9};
        //Spliney = new double[] {-131.49011, 3.116379E-10, 131.49011, 213.08975, 131.49011, -6.7329753E-10, -131.49011, -213.08975};
        //Splinex = new double[8];
        //Spliney = new double[8];
        //Splinex = new double[] {96,200,190,120,174,-70,140,120};     // fix fix override previous
        //Spliney = new double[] {95,90,30,230,186,20,160,20};      // fix fix override previous
        t2ddx = new double[Splinex.length][N+1];
        t2ddy = new double[Spliney.length][N+1];
        //for (int i = 0; i < Splinex.length; i++)
        //    Splinex[(i+1)%Splinex.length] = rad*Math.cos(i*Math.PI*2/Splinex.length - 0.0*Math.PI/4);
        //for (int i = 0; i < Spliney.length; i++)
        //    Spliney[(i+1)%Splinex.length] = rad*Math.sin(i*Math.PI*2/Spliney.length - 0.0*Math.PI/4);
        //System.out.println("BSplineN iterate_at_x_y = " + iterate_at_x_y() + "\n");       // normally include this line
        //scan_steepest();
        System.out.println("BSplineN solve_at_x_y = " + solve_at_x_y(true) + "\n");
        //System.out.println("BSplineN circle = ," + Splinex.length + ", " + rad + ", " + solve_at_x_y(true) + "\n");
        //read_BSplineN_data(133);
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //for (int i = 0; i < 100; i++)
        //    System.out.println(i/100.0 + ", " + d2Ni3(i/100.0)[0] + ", " + d2Ni3(i/100.0)[1] + ", " + d2Ni3(i/100.0)[2] + ", " + d2Ni3(i/100.0)[3]);
        //double[] perm = permute(Splinex, 4.5);
        //System.out.println("permute = " + perm[0] + ", " + perm[1] + ", " + perm[2] + ", " + perm[3]);
        //gen_spline_points(100);
        //for (int i = 0; i <= 60; i++)               // demo printout
        //{
        //    double d = i/90.0;
        //    System.out.println(i + ", " + fitted.getx(t1_start + i*(t1_end - t1_start)/N) + ", " + fitted.gety(t1_start + i*(t1_end - t1_start)/N)
        //                         + ", " + BSpline5.multvv(permute(Splinex, d), Ni3(d % 1)) + ", " + BSpline5.multvv(permute(Spliney, d), Ni3(d % 1))
        //                         + ", " + BSpline5.multvv(permute(Splinex, d), dNi3(d % 1)) + ", " + BSpline5.multvv(permute(Spliney, d), dNi3(d % 1)));
        //}
        //Splinex = new double[] {10, 100, 125, 110+.01, 20};
        //for (int i = 0; i <= N; i++)
            //System.out.println(i + ", " + calc_dfdai(t2[i], 1) + ", " + calc_dfdai(t2[i], 3));
            //System.out.println(i + ", " + BSpline5.multvv(permute(Splinex, t2[i]), Ni3(t2[i] % 1)));
        //confirm_LHS();
        //DoubleMatrix d = new DoubleMatrix(2, 2, 5, 1.3, 1.3, 6);
        //Eigen.symmetricEigenvalues(d);
        //System.out.println(d);
        //System.out.println(Eigen.symmetricEigenvalues(d));
        //solve_for_beta();
    }

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
//        double[][] Jac = new double[2*Splinex.length][2*Splinex.length];
//        double[] dFda = new double[2*Splinex.length];
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

    private static void read_BSplineN_data(int ln)
    {
        // read BSplineN data (c, ri, thetai) from a file
        // use only one line (ln)

        String str = "";
        double c = Double.NaN;
        int row = 0;

        try
        {
            //BufferedReader istr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\B_SplineN\\run16.csv"));
            BufferedReader istr = new BufferedReader(new FileReader("\\APP\\Java\\SpiroGraph\\B_SplineN\\run8.csv"));
            try
            {
                while (istr.ready() && row < ln)
                {
                    str = istr.readLine();
                    row++;
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

        if (!str.startsWith("gauss t2[] @"))
        {
            System.out.println("valid data not found: '" + str + "'");
            return;
        }
        c = Double.parseDouble(str.split(",")[3]);
        fitted = new epiTrochoidFxn(c);
        theta_start = fitted.gettheta(t1_start);                        // possibly redundant
        theta_end = fitted.gettheta(t1_end);
        for (int i = 0; i < Splinex.length; i++)
            Splinex[i] = Double.parseDouble(str.split(",")[6 + i])*Math.cos(Math.PI*Double.parseDouble(str.split(",")[Splinex.length + 6 + i])/180);
        for (int i = 0; i < Splinex.length; i++)
            Spliney[i] = Double.parseDouble(str.split(",")[6 + i])*Math.sin(Math.PI*Double.parseDouble(str.split(",")[Splinex.length + 6 + i])/180);
        System.out.println("file data at line " + ln + ", " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial());
        System.out.println("BSplineN iterate_at_x_y = " + iterate_at_x_y() + "\n");
    }

    private static void write_BSplineN_data(double[][] J)
    {
        // write BSplineN data to file as in BSpline5.dump_Jac(J)
        // see 'python_eigenvalues.txt' to calculate eigenvalues using numpy
        try
        {
            FileWriter fw = new FileWriter("\\APP\\Java\\SpiroGraph\\B_SplineN\\dump_augment_raw_8.txt", true);
            PrintWriter out = new PrintWriter(fw);
            if (J == null)
                out.println("\ndFdc = " + fitted.getc() + print_radial());    // write header for data
            else
            {
                out.print("a = np.array([");
                for (int i = 0; i < J.length; i++)
                {
                    if (i > 0) out.print(", ");
                    out.print("[");
                    for (int j = 0; j < J.length; j++)
                    {
                        if (j > 0) out.print(", ");
                        out.print(J[i][j]);
                    }
                    out.print("]");
                }
                out.println("])");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("write_BSplineN_data() save error = " + e);}
    }

    private static double solve_at_x_y(boolean print)
    {
        // N-point uniform, closed, cubic B-Spline
        // perform a single calculation of a complete t2[] profile
        // at a given set {x, y}, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        if (t2[N] == 0)
            System.out.println("__start B-SplineN theta c t x y = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + print_coord(Splinex) + print_coord(Spliney));
        else
        {
            System.out.println("__solve at new d1 d2 rms = , , , , , " + print_coord(Splinex) + print_coord(Spliney) + ", " + calc_error());
            System.out.println("__solve at new d1 d2 rms = , , , , , " + print_radial() + ", " + calc_error());
        }
        //fitted.gen_BezierN(Splinex, Spliney);

        if (print) System.out.println("\n , t1, t2, t2ddx0, t2ddx1, t2ddx2, t2ddx3, t2ddx4, t2ddx5, t2ddx6, t2ddx7");
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i);
            if ((i == N && Math.abs(t2[N] - t2[0] - Splinex.length) > TOL)
            ||  (t2[i] < -5000*TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                return Double.NaN;
            }
            for (int j = 0; j < Splinex.length; j++)
                t2ddx[j][i] = calc_t2dd(j, i, t2[i], "x");
            for (int j = 0; j < Spliney.length; j++)
                t2ddy[j][i] = calc_t2dd(j, i, t2[i], "y");
            if (print)
                System.out.println(i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + print_coord(t2ddx, i) + print_coord(t2ddy, i));
                //System.out.println(i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + calc_t2dd(Integer.MAX_VALUE, i, t2[i], "c"));
        }
        //print_one_point(155);
        double retVal = calc_error();
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(Splinex) + print_coord(Spliney) + ", " + (float) retVal + ", " + (float) Jacdet);
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_coord(Splinex) + print_coord(Spliney) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2);
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + ", " + (float) retVal + ", " + (float) Jacdet);
        return retVal;
    }

    private static void scan_steepest()
    {
        // calculate a streamline using steepest-descent starting at a saddle point
        // use method of 'Chong and Zak', 'An Introduction to Optimization', p.121
        // move in increments of 'steepest' in the direction of the gradient -dFda[]
        // to run this, temporarily comment out lines 401, 430 and 431 to reduce the output
        // to run this, temporarily enable line 223 to abort 'iterate_at_x_y()'
        // see Book 9, page 42

        final int Nstr = 50;
        double steepest = Double.NaN;
        double retVal;
        double num, denom;

        System.out.println("scan steepest, " + Nstr);
        iterate_at_x_y();                                       // initiallize dFdd
        retVal = calc_error();
        System.out.println("scan = , , , " + (float) fitted.getc() + ", " + print_coord(Splinex) + print_coord(Spliney) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2 + ", " + steepest);

        for (int i = 0; i < Nstr; i++)
        {
            t2[N] = 0;                                          // just to control the output
            num = BSpline5.multvv(dFda, dFda);
            denom = BSpline5.multvv(dFda, BSpline5.multmv(Jac, dFda));
            steepest = num/denom;
            for (int j = 0; j < Splinex.length; j++)
            {
                Splinex[j] -= steepest*dFda[j];
                Spliney[j] -= steepest*dFda[j + Splinex.length];
            }
            iterate_at_x_y();                                   // re-calculate dFdd
            retVal = calc_error();
            System.out.println("scan = , , , " + (float) fitted.getc() + ", " + print_coord(Splinex) + print_coord(Spliney) + ", " + (fitted.a + fitted.b)*(fitted.a + fitted.b)*retVal*retVal/2 + ", " + steepest);
        }
    }

    private static void confirm_LHS()
    {
        // calculate LHS of dfydu * dfxdai
        double u;
        int j = 3;                          // index of betax[j]
        System.out.println("confirm_LHS = " + print_coord(Spliney));
        double[] D = new double[6];
        double[] coeff = new double[] {1*(Spliney[2] - Spliney[0]),
                                       4*(Spliney[2] - Spliney[1]),
                                       1*(Spliney[3] - Spliney[1]),
                                       0,
                                       0,
                                       0};
        coeff = new double[] { 4*(Spliney[2] - Spliney[0]),
                              12*(Spliney[2] - Spliney[0]) + 16*(Spliney[2] - Spliney[1]),
                               6*(Spliney[2] - Spliney[0]) + 48*(Spliney[2] - Spliney[1]) +  4*(Spliney[3] - Spliney[1]),
                               1*(Spliney[2] - Spliney[0]) + 24*(Spliney[2] - Spliney[1]) + 12*(Spliney[3] - Spliney[1]),
                               4*(Spliney[2] - Spliney[1]) +  6*(Spliney[3] - Spliney[1]),
                               1*(Spliney[3] - Spliney[1])};
        coeff = new double[] { 1*(Spliney[2] - Spliney[0]),
                               6*(Spliney[2] - Spliney[0]) +  4*(Spliney[2] - Spliney[1]),
                              12*(Spliney[2] - Spliney[0]) + 24*(Spliney[2] - Spliney[1]) + 1*(Spliney[3] - Spliney[1]),
                               4*(Spliney[2] - Spliney[0]) + 48*(Spliney[2] - Spliney[1]) + 6*(Spliney[3] - Spliney[1]),
                              16*(Spliney[2] - Spliney[1]) + 12*(Spliney[3] - Spliney[1]),
                               4*(Spliney[3] - Spliney[1])};
        coeff = new double[] { 0,
                               0,
                               0,
                               1*(Spliney[2] - Spliney[0]),
                               4*(Spliney[2] - Spliney[1]),
                               1*(Spliney[3] - Spliney[1])};
        for (int i = 0; i <= 100; i++)
        {
            u = i/100.0;
            for (int k = 0; k < 6; k++)
                D[k] = Math.pow(u, k)*Math.pow(1 - u, 5 - k);
            System.out.println(i + ", " + BSpline5.multvv(permute(Spliney, u), dNi3(u))*Ni3(u)[j]
                                 + ", " + BSpline5.multvv(D, coeff)/12);
        }
    }

    private static void correlate_f_dfda(double[][] m, int seg)
    {
        System.out.println("correlate_f_dfda: seg = " + seg);
        for (int k = 0; k < N; k++)
        {
            System.out.print("m = " + k);
            for (int i = 0; i < 8; i++)
                System.out.print(", " + m[k][i]);
            System.out.println();
        }
        double[][] corr = new double[8][8];
//        double vec[] = new double[7];
//        double ret[] = new double[7];
//        double retN[] = new double[8];

        for (int i = 0; i < 8; i++)
        {
            //for (int k = 0; k < N; k++)
            //    vec[i - 1] += m[k][0]*m[k][i];
            for (int j = i; j < 8; j++)
            {
                for (int k = 0; k < N; k++)
                    corr[i][j] += m[k][i]*m[k][j];
                corr[j][i] = corr[i][j];
                //System.out.print(" (" + i + ", " + j + ") " + corr[i][j]);
            }
            //System.out.println();
        }
        System.out.println("\ncorr:");
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 8; j++)
                System.out.print((float) corr[i][j] + ", ");
            System.out.println();
        }
        System.out.println("\nseg det corr = , " + seg + ", " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + N + ", " + print_radial() + ", " + (float) BSpline5.detm(corr));
        BSpline5.dump_Jac(corr);
        //ret = BSpline5.gaussj(corr, vec);
        //retN[0] = 1;
        //for (int i = 1; i < 8; i++)
        //    retN[i] = -ret[i - 1];
        //System.out.println("correlate vec = " + print_coord(ret));
        //System.out.println("error     vec = " + print_coord(BSpline5.multmv(m, retN)));
    }

    private static void solve_for_beta()
    {
        int i, j;
        double rad1 = 1.;
        double rad2 = 1.123;
        double[] rad  = new double[] {rad1, rad2, rad1, rad2, rad1, rad2, rad1, rad2};
        double th1 = -40.;
        double th2 = -30.;
        double[] th = new double[] {th1, th2, th1 + 90, th2 + 90, th1 + 180, th2 + 180, th1 + 270, th2 + 270};
        double[] x = new double[2*rad.length];
        double[] y = new double[2*rad.length];

        for (i = 0; i < rad.length; i++)
        {
            x[i] = rad[i]*Math.cos(Math.PI*th[i]/180);
            x[i + rad.length] = x[i];                   // allow for cyclic permutations
            y[i] = rad[i]*Math.sin(Math.PI*th[i]/180);
            y[i + rad.length] = y[i];
        }
        System.out.print("x = ");
        for (i = 0; i < rad.length; i++)
            System.out.print(", " + x[i]);
        System.out.println();
        System.out.print("y = ");
        for (i = 0; i < rad.length; i++)
            System.out.print(", " + y[i]);
        System.out.println();
        double[][] sol = new double[8*rad.length][8*rad.length];

        for (i = 0; i < rad.length; i++)
        {
            j = i + 4;
            sol[i][8*i]     =     y[j + 2] - y[j];
            sol[i][8*i + 1] =  4*(y[j + 2] - y[j]);
            sol[i][8*i + 2] =     y[j + 2] - y[j];
            sol[i][8*i + 4] =   -(x[j + 2] - x[j]);
            sol[i][8*i + 5] = -4*(x[j + 2] - x[j]);
            sol[i][8*i + 6] =   -(x[j + 2] - x[j]);
        }

        for (i = 0; i < rad.length; i++)
        {
            j = i + 4;
            sol[rad.length + i][8*i]     =   4*(y[j + 2] - y[j + 1]);
            sol[rad.length + i][8*i + 1] =  16*(y[j + 2] - y[j + 1]) + 12*(y[j + 2] - y[j]);
            sol[rad.length + i][8*i + 2] =   4*(y[j + 2] - y[j + 1]) +  6*(y[j + 2] - y[j]);
            sol[rad.length + i][8*i + 4] =  -4*(x[j + 2] - x[j + 1]);
            sol[rad.length + i][8*i + 5] = -16*(x[j + 2] - x[j + 1]) - 12*(x[j + 2] - x[j]);
            sol[rad.length + i][8*i + 6] =  -4*(x[j + 2] - x[j + 1]) -  6*(x[j + 2] - x[j]);
        }

        for (i = 0; i < rad.length; i++)
        {
            j = i + 4;
            sol[2*rad.length + i][8*i]     =     y[j + 3] - y[j + 1];
            sol[2*rad.length + i][8*i + 1] =  4*(y[j + 3] - y[j + 1]) + 48*(y[j + 2] - y[j + 1]) +  6*(y[j + 2] - y[j]);
            sol[2*rad.length + i][8*i + 2] =     y[j + 3] - y[j + 1]  + 24*(y[j + 2] - y[j + 1]) + 12*(y[j + 2] - y[j]);
            sol[2*rad.length + i][8*i + 4] =   -(x[j + 3] - x[j + 1]);
            sol[2*rad.length + i][8*i + 5] = -4*(x[j + 3] - x[j + 1]) - 48*(x[j + 2] - x[j + 1]) -  6*(x[j + 2] - x[j]);
            sol[2*rad.length + i][8*i + 6] =   -(x[j + 3] - x[j + 1]) - 24*(x[j + 2] - x[j + 1]) - 12*(x[j + 2] - x[j]);
        }

        for (i = 0; i < rad.length; i++)
        {
            j = i + 4;
            sol[3*rad.length + i][8*i + 1] =  12*(y[j + 3] - y[j + 1]) + 24*(y[j + 2] - y[j + 1]) +    y[j + 2] - y[j];
            sol[3*rad.length + i][8*i + 2] =   6*(y[j + 3] - y[j + 1]) + 48*(y[j + 2] - y[j + 1]) + 4*(y[j + 2] - y[j]);
            sol[3*rad.length + i][8*i + 3] =                                                           y[j + 2] - y[j];
            sol[3*rad.length + i][8*i + 5] = -12*(x[j + 3] - x[j + 1]) - 24*(x[j + 2] - x[j + 1]) -   (x[j + 2] - x[j]);
            sol[3*rad.length + i][8*i + 6] =  -6*(x[j + 3] - x[j + 1]) - 48*(x[j + 2] - x[j + 1]) - 4*(x[j + 2] - x[j]);
            sol[3*rad.length + i][8*i + 7] =                                                         -(x[j + 2] - x[j]);
        }

        for (i = 0; i < rad.length; i++)
        {
            j = i + 4;
            sol[4*rad.length + i][8*i + 1] =   6*(y[j + 3] - y[j + 1]) +  4*(y[j + 2] - y[j + 1]);
            sol[4*rad.length + i][8*i + 2] =  12*(y[j + 3] - y[j + 1]) + 16*(y[j + 2] - y[j + 1]);
            sol[4*rad.length + i][8*i + 3] =                              4*(y[j + 2] - y[j + 1]);
            sol[4*rad.length + i][8*i + 5] =  -6*(x[j + 3] - x[j + 1]) -  4*(x[j + 2] - x[j + 1]);
            sol[4*rad.length + i][8*i + 6] = -12*(x[j + 3] - x[j + 1]) - 16*(x[j + 2] - x[j + 1]);
            sol[4*rad.length + i][8*i + 7] =                             -4*(x[j + 2] - x[j + 1]);
        }

        for (i = 0; i < rad.length; i++)
        {
            j = i + 4;
            sol[5*rad.length + i][8*i + 1] =     y[j + 3] - y[j + 1];
            sol[5*rad.length + i][8*i + 2] =  4*(y[j + 3] - y[j + 1]);
            sol[5*rad.length + i][8*i + 3] =     y[j + 3] - y[j + 1];
            sol[5*rad.length + i][8*i + 5] =   -(x[j + 3] - x[j + 1]);
            sol[5*rad.length + i][8*i + 6] = -4*(x[j + 3] - x[j + 1]);
            sol[5*rad.length + i][8*i + 7] =   -(x[j + 3] - x[j + 1]);
        }

        for (i = 0; i < rad.length; i++)                // splice constraint on betaxi
        {
            sol[6*rad.length + i][8*i + 1] = 1;         // coeff of D50 segment i-1
            sol[6*rad.length + i][8*i + 2] = 4;
            sol[6*rad.length + i][8*i + 3] = 1;
            sol[6*rad.length + i][(8*i + 0 + 8) % (8*rad.length)] = -1;    // coeff of D05 segment i
            sol[6*rad.length + i][(8*i + 1 + 8) % (8*rad.length)] = -4;
            sol[6*rad.length + i][(8*i + 2 + 8) % (8*rad.length)] = -1;
        }

        for (i = 0; i < rad.length; i++)                // normalize the betay0 (one per segment)
            sol[7*rad.length + i][8*i] = 1;

        System.out.println("c r0 r1 th0 th1 = , " + fitted.getc() + ", " + Math.sqrt(x[0]*x[0] + y[0]*y[0]) + ", " + Math.sqrt(x[1]*x[1] + y[1]*y[1])
                                           + ", " + Math.atan(y[0]/x[0])*180/Math.PI + ", " + Math.atan(y[1]/x[1])*180/Math.PI);

        double[] rhs = new double[8*rad.length];
        for (i = 0; i < rad.length; i++)                // initiallize the betay0 = 1
            rhs[7*rad.length + i] = 1;
        double[] lhs = BSpline5.gaussj(sol, rhs);
        for (i = 0; i < rad.length; i++)
            System.out.println((float) lhs[8*i + 0] + ", " + (float) lhs[8*i + 1] + ", " + (float) lhs[8*i + 2] + ", " + (float) lhs[8*i + 3] + ", " + (float) lhs[8*i + 4] + ", " + (float) lhs[8*i + 5] + ", " + (float) lhs[8*i + 6] + ", " + (float) lhs[8*i + 7]);
    }

    private static void print_one_point(int i)
    {
        //t2[i] += 0.001;
        System.out.println("\n" + i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]
                                    + ", " + fitted.getx(t1_start + i*(t1_end - t1_start)/N) + ", " + fitted.gety(t1_start + i*(t1_end - t1_start)/N)
                                    + ", " + BSpline5.multvv(permute(Splinex, t2[i]), Ni3(t2[i] % 1)) + ", " + BSpline5.multvv(permute(Spliney, t2[i]), Ni3(t2[i] % 1)));
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
        for (int i = 0; i < Splinex.length; i++)
            str += ", " + (double) Math.sqrt(Splinex[i]*Splinex[i] + Spliney[i]*Spliney[i]);
        for (int i = 0; i < Splinex.length; i++)
            str += ", " + (double) (Math.atan2(Spliney[i], Splinex[i])*180/Math.PI);
        return str;
    }

    private static void gen_spline_points(double Npt)
    {
        double origin_x = 0; // 50;
        double origin_y = 0; // 600;
        double scale = 1; // 16;

        System.out.printf("M");
        for (int i = 0; i <= Npt*Splinex.length; i++)
        {
            double d = i/Npt;
            System.out.printf(" %f, %f", origin_x + scale*BSpline5.multvv(permute(Splinex, d), Ni3(d % 1)), origin_y - scale*BSpline5.multvv(permute(Spliney, d), Ni3(d % 1)));
        }
        System.out.println("\n");
    }

    private static double[] permute(double[] arr, double u)
    {
        if (arr.length < 4) return null;
        if (u < -5000*TOL) return null;

        double[] ret = new double[4];
        int istart = (int) u;
        if (istart > arr.length)
        {
            System.out.println("u is too large: " + u + " > " + arr.length);
            return null;
        }
        for (int i = 0; i < ret.length; i++)
            ret[i] = arr[(i + istart)%arr.length];
        return ret;
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
        //double a_b = 180;           // scale factor to make rms error dimensionless
        //double a_b = 1;             // Cycloid only
        double t1 = t1_start;
        double[] trap_in = new double[N];

        if (Math.abs(t2[N] - t2[0] - Splinex.length) > TOL)
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i < N; i++)
        {
            trap_in[i] = (BSpline5.multvv(permute(Splinex, t2[i]), Ni3(t2[i] % 1)) - fitted.getx(t1))*(BSpline5.multvv(permute(Splinex, t2[i]), Ni3(t2[i] % 1)) - fitted.getx(t1))
                       + (BSpline5.multvv(permute(Spliney, t2[i]), Ni3(t2[i] % 1)) - fitted.gety(t1))*(BSpline5.multvv(permute(Spliney, t2[i]), Ni3(t2[i] % 1)) - fitted.gety(t1));
            //if (i == 155)
            //    System.out.println(i + ", , " + (BSpline5.multvv(permute(Splinex, t2[i]), Ni3(t2[i] % 1)) - fitted.getx(t1)) + ", " + (BSpline5.multvv(permute(Spliney, t2[i]), Ni3(t2[i] % 1)) - fitted.gety(t1)) + ", " + Math.sqrt(trap_in[i]));
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
        f = (BSpline5.multvv(permute(Splinex, t), Ni3(t % 1)) - X)*BSpline5.multvv(permute(Splinex, t), dNi3(t % 1))
          + (BSpline5.multvv(permute(Spliney, t), Ni3(t % 1)) - Y)*BSpline5.multvv(permute(Spliney, t), dNi3(t % 1));
        fprime = BSpline5.multvv(permute(Splinex, t), dNi3(t % 1))*BSpline5.multvv(permute(Splinex, t), dNi3(t % 1)) + (BSpline5.multvv(permute(Splinex, t), Ni3(t % 1)) - X)*BSpline5.multvv(permute(Splinex, t), d2Ni3(t % 1))
               + BSpline5.multvv(permute(Spliney, t), dNi3(t % 1))*BSpline5.multvv(permute(Spliney, t), dNi3(t % 1)) + (BSpline5.multvv(permute(Spliney, t), Ni3(t % 1)) - Y)*BSpline5.multvv(permute(Spliney, t), d2Ni3(t % 1));
        f2prime = 3*BSpline5.multvv(permute(Splinex, t), dNi3(t % 1))*BSpline5.multvv(permute(Splinex, t), d2Ni3(t % 1)) + (BSpline5.multvv(permute(Splinex, t), Ni3(t % 1)) - X)*BSpline5.multvv(permute(Splinex, t), d3Ni3(t % 1))
                + 3*BSpline5.multvv(permute(Spliney, t), dNi3(t % 1))*BSpline5.multvv(permute(Spliney, t), d2Ni3(t % 1)) + (BSpline5.multvv(permute(Spliney, t), Ni3(t % 1)) - Y)*BSpline5.multvv(permute(Spliney, t), d3Ni3(t % 1));
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
            f = (BSpline5.multvv(permute(Splinex, t), Ni3(t % 1)) - X)*BSpline5.multvv(permute(Splinex, t), dNi3(t % 1))
              + (BSpline5.multvv(permute(Spliney, t), Ni3(t % 1)) - Y)*BSpline5.multvv(permute(Spliney, t), dNi3(t % 1));
            fprime = BSpline5.multvv(permute(Splinex, t), dNi3(t % 1))*BSpline5.multvv(permute(Splinex, t), dNi3(t % 1)) + (BSpline5.multvv(permute(Splinex, t), Ni3(t % 1)) - X)*BSpline5.multvv(permute(Splinex, t), d2Ni3(t % 1))
                   + BSpline5.multvv(permute(Spliney, t), dNi3(t % 1))*BSpline5.multvv(permute(Spliney, t), dNi3(t % 1)) + (BSpline5.multvv(permute(Spliney, t), Ni3(t % 1)) - Y)*BSpline5.multvv(permute(Spliney, t), d2Ni3(t % 1));
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

    private static double calc_dfdai(double u, int i)
    {
        int istart = (int) u;
        if ((i + Splinex.length - istart) % Splinex.length > 3)
            return 0;
        else
            return Ni3(u % 1)[(i + Splinex.length - istart) % Splinex.length];
    }

    private static double calc_d2fdudai(double u, int i)
    {
        int istart = (int) u;
        if ((i + Splinex.length - istart) % Splinex.length > 3)
            return 0;
        else
            return dNi3(u % 1)[(i + Splinex.length - istart) % Splinex.length];
    }

    private static double[] Ni3(double u)
    {
        if (u < -5000*TOL)
            return null;
        else if (u < 1 + TOL)
            return new double[] {(1 - u)*(1 - u)*(1 - u)/6,
                                 ((u + 2)*(1 - u)*(1 - u) + (2 - u)*((u + 1)*(1 - u) + u*(2 - u)))/6,
                                 ((u + 1)*((u + 1)*(1 - u) + u*(2 - u)) + u*u*(3 - u))/6,
                                 u*u*u/6};
        else
            return null;
    }

    private static double[] dNi3(double u)
    {
        if (u < -5000*TOL)
            return null;
        else if (u < 1 + TOL)
            return new double[] {-(1 - u)*(1 - u)/2,
                                 ((1 - u)*(1 - u) - 2*(u + 2)*(1 - u)
                               - (u + 1)*(1 - u) + (2 - u)*(1 - u) - (2 - u)*(u + 1)
                               - 2*(2 - u)*u + (2 - u)*(2 - u))/6,
                                 (2*(u + 1)*(1 - u) - (u + 1)*(u + 1)
                               + (2 - u)*u - (u + 1)*u + (u + 1)*(2 - u)
                               - u*u + 2*(3 - u)*u)/6,
                                 u*u/2};
        else
            return null;
    }

    private static double[] d2Ni3(double u)
    {
        if (u < -5000*TOL)
            return null;
        else if (u < 1 + TOL)
            return new double[] {1 - u,
                                 (-6*(1 - u) - 6*(2 - u)
                               + 2*(u + 2) + 2*(u + 1) + 2*u)/6,
                                 (2*(1 - u) + 2*(2 - u) + 2*(3 - u)
                               - 6*(u + 1) - 6*u)/6,
                                 u};
        else
            return null;
    }

    private static double[] d3Ni3(double u)
    {
        if (u < -5000*TOL)
            return null;
        else if (u < 1 + TOL)
            return new double[] {-1,
                                  3,
                                 -3,
                                  1};
        else
            return null;
    }
}
