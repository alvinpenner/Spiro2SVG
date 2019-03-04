
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
    private static double[] dFdd = new double[2];           // defined here so we can share it with 'scan_streanline_at_P2'
    //private static double del_d, del_angle;                 // size and direction of iteration for streamlines
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
        fitted = new epiTrochoidFxn(3.592);
        //iterate_at_P2(59, 29);     // normal optimization
        //iterate_at_P2(14, 38, 2);                     // area-constrained optimization
        //scan_streamline_at_P2(56, 34.5);
        //scan_streamline_at_P2(56.1, 34.56);
        //scan_streamline_at_P2(56.053 + .005, 34.557);
        //scan_streamline_at_P2(58.94795711624154, 29.62577804214821);
        //scan_streamline_at_P2(56.05, 34.53);
        //scan_streamline_at_P2(56.06272 - .1, 34.54501);
//        scan_streamline_at_P2(56.05399864376891, 34.559340683666015);   // intermediate position
        //scan_streamline_at_P2(53.3, 38.9);
        //iterate_at_P2(36.063996220825075, 63.23990306508096, 36.28674943705403, 62.840719240125544);    // linear-constrained optimization
        //fitted = new epiTrochoidFxn(3.59611);
        //iterate_at_P2(14.959259818763828, 92.24198234164545);
        scan_streamline_at_P2(57.13628153609134, 31.41210991725839);
        //fitted = new epiTrochoidFxn(3.5);
        //iterate_at_P2(57.6415663538542, 30.721902560475314);
        //fitted = new epiTrochoidFxn(3.5965);
        //iterate_at_P2(57.42308520622449 ,  30.879968);
        //solve_at_P2(60, 40, false);
        //System.out.println("cubic Bezier solve_at_P2 = " + solve_at_P2(60, 40, true) + "\n");
        //calc_array();               // generate Python 2D contour plot of rms
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
    }

    //private static void iterate_at_P2(double d1, double d2)
    private static void iterate_at_P2(double... values)
    {
        // calculate a new estimate of (d1, d2) by setting dF = 0
        // include only first-order responses
        // see Spiro2SVG Book 3, page 54 (applied to quartic Bezier)
        // setup 2-variable Newton-Raphson iteration
        // see Core Java Vol I, p. 189, for use of 'variable number of parameters'

        final double gain = 1;                              // fudge factor to reduce gain
        final int MAXLOOP = 1000;
        double d1, d2, theta_perp = Double.NaN;             // initial (d1,d2), angle of scan
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
        double[] df_gxdc = new double[N+1];                         // used only for calc of dx[]/dc
        double[] df_gydc = new double[N+1];
        double[] h_inv = new double[N+1];                           // Book 7, page 40
        double alpha = 1.81418686001726;                             // Book 7, page 40

        double[] trap_in = new double[N+1];
        double[][] Jac = new double[2][2];
        //double[] dFdd = new double[2];
        double[] d2Fdddc = new double[2];                           // augmented matrix
        double dFdc;
        double d2Fdcdc;                                             // augmented matrix
        double[][] Augment = new double[3][3];
        double[] deld;                                              // (-Δd1, -Δd2)
        int i, j, k, loop = 0;
        double t1;

        if (values.length == 2)                 // normal 2D optimization of (d1, d2)
        {
            d1 = values[0];
            d2 = values[1];
        }
        else if (values.length == 3)
        {
            double index = values[2];           // should be (0 - 4) or Double.NaN
            if (Double.isNaN(index))
            {
                d1 = values[0];                 // steepest-descent optimization
                d2 = values[1];
            }
            else                                // area-constrained optimization
            {
                // interpolate d1 linearly between values[0] and values[1]
                // calculate d2 using area constraint
                // define scan direction as perpendicular to dd2/dd1
                d1 = values[0] + index*(values[1] - values[0])/4;
                d2 = calc_d2(d1);
                theta_perp = Math.atan(calc_dd2dd1(d1)) - Math.PI/2;
                System.out.println("area-constrained optimization at , " + values[0] + ", " + values[1] + ", " + values[2] + ", " + d1 + ", " + d2 + ", " + theta_perp*180/Math.PI);
            }
        }
        else if (values.length == 4)            // linear-constrained 1D optimization
        {                                       // moving perpendicular to the midpoint between two solutions
            d1 = (values[0] + values[2])/2;     // calculate midpoint
            d2 = (values[1] + values[3])/2;
            theta_perp = Math.atan2(values[3] - values[1], values[2] - values[0]) - Math.PI/2;
            System.out.println("linear-constrained optimization at , " + values[0] + ", " + values[1] + ", " + values[2] + ", " + values[3] + ", " + d1 + ", " + d2 + ", " + theta_perp*180/Math.PI);
        }
        else
        {
            System.out.println("Wrong number of arguments in 'iterate_at_P2'");
            return;
        }
        //System.out.println("args = " + values.length + ", " + d1 + ", " + d2 + ", " + theta_perp);
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

                //double testf0 = (dfxdd[0][i] + dfxdu[i]*t2dd[0][i])*fitted.getdxdc(t1)
                //              + (dfydd[0][i] + dfydu[i]*t2dd[0][i])*fitted.getdydc(t1);    // temporary test code fix fix
                //double testf1 = (dfxdd[1][i] + dfxdu[i]*t2dd[1][i])*fitted.getdxdc(t1)
                //              + (dfydd[1][i] + dfydu[i]*t2dd[1][i])*fitted.getdydc(t1);    // temporary test code fix fix
                //double testf0 = (dfxdd[0][i] + 0*dfxdu[i]*t2dd[0][i])*df_gxdc[i]
                //              + (dfydd[0][i] + 0*dfydu[i]*t2dd[0][i])*df_gydc[i];    // temporary test code fix fix
                //double testf1 = (dfxdd[1][i] + 0*dfxdu[i]*t2dd[1][i])*df_gxdc[i]
                //              + (dfydd[1][i] + 0*dfydu[i]*t2dd[1][i])*df_gydc[i];    // temporary test code fix fix
                //double testf0 = (dfxdd[0][i] + 0*dfxdu[i]*t2dd[0][i])*f_gx[i]
                //              + (dfydd[0][i] + 0*dfydu[i]*t2dd[0][i])*f_gy[i];    // temporary test code fix fix
                //double testf1 = (dfxdd[1][i] + 0*dfxdu[i]*t2dd[1][i])*f_gx[i]
                //              + (dfydd[1][i] + 0*dfydu[i]*t2dd[1][i])*f_gy[i];    // temporary test code fix fix
                //double testf0 = (dfxdd[0][i] + 0*dfxdu[i]*t2dd[0][i])*fitted.getx(t1)
                //              + (dfydd[0][i] + 0*dfydu[i]*t2dd[0][i])*fitted.gety(t1);    // temporary test code fix fix
                //double testf1 = (dfxdd[1][i] + 0*dfxdu[i]*t2dd[1][i])*fitted.getx(t1)
                //              + (dfydd[1][i] + 0*dfydu[i]*t2dd[1][i])*fitted.gety(t1);    // temporary test code fix fix
                double testf0 = dfxdd[0][i]*f_gx[i] + dfydd[0][i]*f_gy[i];    // temporary test code fix fix
                double testf1 = dfxdd[1][i]*f_gx[i] + dfydd[1][i]*f_gy[i];    // temporary test code fix fix
                //double perp0 = df_gxdc[i]*f_gx[i] + df_gydc[i]*f_gy[i];
                //double perp1 = df_gxdc[1]*f_gx[i] + df_gydc[1]*f_gy[i];
                //double dot0 = f_gx[i]*dfxdu[i] + f_gy[i]*dfydu[i];
                //double testf0 = dfxdd[0][i]*f_gx[i] + dfydd[0][i]*f_gy[i];    // temporary test code fix fix
                //double testf1 = df_gxdc[i]*f_gx[i] + df_gydc[i]*f_gy[i];    // temporary test code fix fix
                //double alpha = 1.7625079209276;
                //double testdfdd1 = (dfydd[0][i] - alpha*dfydd[1][i])/(dfxdd[0][i] - alpha*dfxdd[1][i]);
                //double testd2fdudd1 = (d2fydudd[0][i] - alpha*d2fydudd[1][i])/(d2fxdudd[0][i] - alpha*d2fxdudd[1][i]);
                //double magf_g = Math.sqrt(f_gx[i]*f_gx[i] + f_gy[i]*f_gy[i]);
                //double mag0 = Math.sqrt(dfxdd[0][i]*dfxdd[0][i] + dfydd[0][i]*dfydd[0][i]);
                //double mag1 = Math.sqrt(dfxdd[1][i]*dfxdd[1][i] + dfydd[1][i]*dfydd[1][i]);
                //System.out.println("testf = ," + t1 + ", " + magf_g + ", " + mag0 + ", " + mag1 + ", " + testf0 + ", " + testf1 + ", " + f_gy[i]/f_gx[i]);
                //System.out.println("testf = ," + (float) t1 + ", " + (float) testf0 + ", " + (float) testf1 + ", " + testf0/testf1); // + ", " + (float) perp0 + ", " + perp0/testf1);
                // calculate the normalization function h(u) - for the oval branch
                double num = Math.sqrt(dfxdu[i]*dfxdu[i] + dfydu[i]*dfydu[i]);
                double den = Math.sqrt((dfxdd[0][i] - alpha*dfxdd[1][i])*(dfxdd[0][i] - alpha*dfxdd[1][i]) + (dfydd[0][i] - alpha*dfydd[1][i])*(dfydd[0][i] - alpha*dfydd[1][i]));
                h_inv[i] = -(t2dd[0][i] - alpha*t2dd[1][i]);
                //System.out.println("normalization = , " + i + ", " + num/den + ", " + 1/h_inv[i]);
                //System.out.println("components  = , " + (dfxdd[0][i]*dfxdd[0][i] + dfydd[0][i]*dfydd[0][i])
                //                               + ", " + (dfxdu[i]*dfxdu[i] + dfydu[i]*dfydu[i] + f_gx[i]*d2fxdudu[i] + f_gy[i]*d2fydudu[i])*t2dd[0][i]*t2dd[0][i]
                //                               + ", " + (dfxdd[0][i]*dfxdd[1][i] + dfydd[0][i]*dfydd[1][i])
                //                               + ", " + (dfxdu[i]*dfxdu[i] + dfydu[i]*dfydu[i] + f_gx[i]*d2fxdudu[i] + f_gy[i]*d2fydudu[i])*t2dd[0][i]*t2dd[1][i]);
                //System.out.println("test d2fdudd = " + (d2fydudd[0][i] + d2fydudu[i]*t2dd[0][i])/(d2fxdudd[0][i] + d2fxdudu[i]*t2dd[0][i]) + ", "
                //                                     + (d2fydudd[1][i] + d2fydudu[i]*t2dd[1][i])/(d2fxdudd[1][i] + d2fxdudu[i]*t2dd[1][i]));
                //System.out.println("alpha CubicBezier," + i + " , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + testf0/testf1);
                //System.out.println("testf = ," + t1 + ", " + testf0 + ", " + testf1 + ", " + testf0/testf1 + ", " + d2fydudu[i]/d2fxdudu[i] + ", " + testd2fdudd1);
                //System.out.println("testf = ," + t1 + ", " + testf0 + ", " + testf1 + ", " + testf0/testf1 + ", " + Math.sqrt(dfydu[i]*dfydu[i] + dfxdu[i]*dfxdu[i])
                //                        + ", " + Math.sqrt((dfydd[0][i] - alpha*dfydd[1][i])*(dfydd[0][i] - alpha*dfydd[1][i]) + (dfxdd[0][i] - alpha*dfxdd[1][i])*(dfxdd[0][i] - alpha*dfxdd[1][i])));
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

            Jacdet = BSpline5.detm(Jac);
            eig0 = (Jac[0][0] + Jac[1][1] - Math.sqrt((Jac[0][0] - Jac[1][1])*(Jac[0][0] - Jac[1][1]) + 4*Jac[0][1]*Jac[0][1]))/2;
            eig1 = (Jac[0][0] + Jac[1][1] + Math.sqrt((Jac[0][0] - Jac[1][1])*(Jac[0][0] - Jac[1][1]) + 4*Jac[0][1]*Jac[0][1]))/2;
            eigangle = -Math.atan((Jac[0][0] - eig0)/Jac[0][1]);        // angle of eigenvector transform
            //del_angle = Math.atan2(-dFdd[1], -dFdd[0]);                 // we wish to decrease F (zeroth-order estimate, obsolete)
            //del_d = Math.sqrt(dFdd[0]*dFdd[0] + dFdd[1]*dFdd[1])/(Math.cos(del_angle)*Math.cos(del_angle)*Jac[0][0] + 2*Math.cos(del_angle)*Math.sin(del_angle)*Jac[0][1] + Math.sin(del_angle)*Math.sin(del_angle)*Jac[1][1]);
            //if (del_d > .01) del_d = .01;
            //del_d = del_d/4;
            //del_d = .0004;
            if (values.length == 3 && Double.isNaN(values[2]))          // steepest-descent optimization
                return;

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

            // calculate the "oval RHS residual", Book 7, page 40

            for (k = 0; k <= N; k++)
                trap_in[k] = h_inv[k]*(f_gx[k]*calc_d2fxdudc(t2[k]) + f_gy[k]*calc_d2fydudc(t2[k]));
            double rhsc = -t2_vs_t1.integrate(trap_in);
            System.out.println("oval residual c = ," + fitted.getc() + ", " + d1 + ", " + d2 + ", " + alpha + ", " + (d2Fdddc[0] - alpha*d2Fdddc[1]) + ", " + rhsc);
            for (k = 0; k <= N; k++)
                trap_in[k] = h_inv[k]*(f_gx[k]*d2fxdudd[0][k] + f_gy[k]*d2fydudd[0][k]);
            double rhs0 = -t2_vs_t1.integrate(trap_in);
            System.out.println("oval residual 0 = ," + fitted.getc() + ", " + d1 + ", " + d2 + ", " + alpha + ", " + (Jac[0][0] - alpha*Jac[0][1]) + ", " + rhs0);
            for (k = 0; k <= N; k++)
                trap_in[k] = h_inv[k]*(f_gx[k]*d2fxdudd[1][k] + f_gy[k]*d2fydudd[1][k]);
            double rhs1 = -t2_vs_t1.integrate(trap_in);
            System.out.println("oval residual 1 = ," + fitted.getc() + ", " + d1 + ", " + d2 + ", " + alpha + ", " + (Jac[1][0] - alpha*Jac[1][1]) + ", " + rhs1);
            //for (k = 0; k <= N; k++)
            //    System.out.println(k + ", " + h_inv[k]*(f_gx[k]*calc_d2fxdudc(t2[k]) + f_gy[k]*calc_d2fydudc(t2[k]))
            //                         + ", " + h_inv[k]*(f_gx[k]*d2fxdudd[0][k] + f_gy[k]*d2fydudd[0][k])
            //                         + ", " + h_inv[k]*(f_gx[k]*d2fxdudd[1][k] + f_gy[k]*d2fydudd[1][k]));

            if (values.length == 2)                         // normal 2D optimization
                deld = BSpline5.gaussj(Jac, dFdd);          // this is actually the negative of Δd
            else                                            // constrained 1D optimization
            {
                double fp = dFdd[0]*Math.cos(theta_perp) + dFdd[1]*Math.sin(theta_perp);
                double f2p = Jac[0][0]*Math.cos(theta_perp)*Math.cos(theta_perp) + 2*Jac[0][1]*Math.cos(theta_perp)*Math.sin(theta_perp) + Jac[1][1]*Math.sin(theta_perp)*Math.sin(theta_perp);
                deld = new double[] {fp/f2p*Math.cos(theta_perp), fp/f2p*Math.sin(theta_perp)};
            }
            d1 -= deld[0]/gain;
            d2 -= deld[1]/gain;

            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + Jacdet + ", " + eig0 + ", " + eig1 + ", " + (float) (eigangle*180/Math.PI) + ", " + BSpline5.detm(Augment));
            System.out.println("deld = " + deld[0] + ", " + deld[1]);
            System.out.println("\ndFdc = " + fitted.getc() + ", " + d1 + ", " + d2 + ", , " + (float) dFdc + ", " + (float) d2Fdcdc);
            BSpline5.dump_Jac(Jac);
            BSpline5.dump_Jac(Augment);
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
            double rms = solve_at_P2(d1, d2, true);                          // final run just for good measure
            // calculate dddc[i]
            double[] dt2dc = BSpline5.gaussj(Jac, d2Fdddc);
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (Math.cos(eigangle)*d2Fdddc[0] - Math.sin(eigangle)*d2Fdddc[1]) + ", " + (float) (Math.sin(eigangle)*d2Fdddc[0] + Math.cos(eigangle)*d2Fdddc[1]));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (eigangle*180.0/Math.PI) + ", " + (float) (Math.atan(d2Fdddc[0]/d2Fdddc[1])*180.0/Math.PI));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) -d2Fdddc[0] + ", " + (float) -d2Fdddc[1]);
            //System.out.println("final Det, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + Jac[0][0] + ", " + Jac[0][1] + ", " + Jac[1][1]);
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) (Jac[1][0]/Jac[0][0]) + ", " + (float) (Jac[1][1]/Jac[0][1]) + ", " + (float) (d2Fdddc[1]/d2Fdddc[0]));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + (float) -dt2dc[0] + ", " + (float) -dt2dc[1] + ", " + (float) eig0 + ", " + (float) BSpline5.detm(Augment));
            if (values.length == 2)
                System.out.println("\nfinal CubicBezier , , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + rms*rms*180.0*180.0/2.0 + ", " + eig0 + ", " + eigangle + ", " + eig1);
            else if (values.length == 3)
                System.out.println("\narea-constrained Cubic, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + rms*rms*180.0*180.0/2.0 + ", , " + (theta_perp + Math.PI/2));
            else if (values.length == 4)
                System.out.println("\ninterpolated Cubic, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + rms*rms*180.0*180.0/2.0 + ", , " + (theta_perp + Math.PI/2));
            //System.out.println("\nfinal CubicBezier, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + eig0 + ", " + 180.0*eigangle/Math.PI);
            //System.out.println(  "collinearity test, , , " + fitted.getc() + ", " + d1 + ", " + d2 + ", " + Jac[1][0]/Jac[0][0] + ", " + Jac[1][1]/Jac[0][1] + ", " + d2Fdddc[1]/d2Fdddc[0]);
            //System.out.println("Jac ratio  = ," + fitted.getc() + ", " + d1 + ", " + d2 + ", " + d2Fdddc[0]/d2Fdddc[1] + ", " + Jac[0][0]/Jac[0][1] + ", " + Jac[1][0]/Jac[1][1]);
            //calc_alpha(d1, d2);
            //scan_at_P2(d1, d2);             // scan 5 points to get F response function
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
        double a_b = 180;         // scale factor to make rms error dimensionless

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
            ; //System.out.println("__start cubic Bezier at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2);
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
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + retVal + ", " + (float) Jacdet + ", " + (float) eig0 + ", " + (float) eig1 + ", " + (float) (eigangle*180/Math.PI));
        System.out.println("F = , , , " + (float) fitted.getc() + ", " + d1 + ", " + d2 + ", " + a_b*a_b*retVal*retVal/2 + ", " + eig0 + ", " + eigangle + ", " + dFdd[0] + ", " + dFdd[1]);
        return retVal;
    }

    private static void scan_streamline_at_P2(double d1_org, double d2_org)
    {
        // calculate a streamline of steepest-descent starting at a saddle point
        // perform one normal optimization to determine a saddle point
        // then move in increments of del_d at an angle determined by dFdd[1]/dFdd[0]
        // to run this, temporarily comment out lines 393 and 429 to reduce the output
        // keep only line 430 : "F = , ..."
        // see Book 8, page 33

        int Nstr = 200;
        double del_d = 0.002;
        double old_angle = Double.NaN;
//        double a0p, a1p, C;

        System.out.println("streamline, " + d1_org + ", " + d2_org + ", " + Nstr + ", " + del_d);
        //iterate_at_P2(d1_org, d2_org);  // normal optimization
                                        // this MUST be a previously converged result
        iterate_at_P2(d1_org, d2_org, Double.NaN);  // initiallize dFdd
        old_angle = eigangle;
        //old_angle = Math.atan2(-dFdd[1], -dFdd[0]);

        for (int i = 0; i < Nstr; i++)
        {
            t2[N] = 0;                                  // just to control the output
            //System.out.println("dFdd , " + dFdd[0] + ", " + dFdd[1]);
//            if (Math.abs(dFdd[0]) < TOL/100. && Math.abs(dFdd[1]) < TOL/1000.)
//            {
//                System.out.println("override");
                //angle = eigangle + (eigangle - old_angle)/2; // + Math.PI;
                //old_angle = eigangle;
//            }
//            else
                //angle = Math.atan2(-dFdd[1], -dFdd[0]);  // we wish to decrease F (zeroth-order estimate, obsolete)
            //eigangle = Math.atan2(dFdd[1], dFdd[0]); // test uphill
            //a0p = ( Math.cos(eigangle)*dFdd[0] + Math.sin(eigangle)*dFdd[1])/eig0;   // current position in
            //a1p = (-Math.sin(eigangle)*dFdd[0] + Math.cos(eigangle)*dFdd[1])/eig1;   // normalized coordinates (page 32)
            //C = Math.pow(Math.abs(a1p), 1.0/eig1)/Math.pow(Math.abs(a0p), 1.0/eig0);
            //System.out.println("eigen data: " + i + ", " + Math.pow(Math.abs(a1p), 1.0/eig1) + ", " + Math.pow(Math.abs(a0p), 1.0/eig0));
            //System.out.println("eigen data: " + i + ", " + a0p + ", " + a1p + ", " + C + ", " + eig0 + ", " + eig1 + ", " + eigangle*180/Math.PI);

            d1_org += del_d*Math.cos(eigangle + (eigangle - old_angle)/2);
            d2_org += del_d*Math.sin(eigangle + (eigangle - old_angle)/2);
            old_angle = eigangle;
            iterate_at_P2(d1_org, d2_org, Double.NaN);  // steepest-descent optimization
        }
        //System.out.println("....................");
        //for (int i = 0; i < 360; i+= 5)
        //{
        //    iterate_at_P2(d1_org + del_d*Math.cos(i*Math.PI/180), d2_org + del_d*Math.sin(i*Math.PI/180), Double.NaN);
        //}
    }

    private static void scan_at_P2(double d1_org, double d2_org)
    {
        // calculate F over a local range around (d1_org, d2_org)
        // move in direction of eigenvector corresponding to the lowest eigenvalue
        // to run this, temporarily comment out lines 341 and 377 to reduce the output
        // keep only line 378 : "F = , ..."

        double d_inc = 0.5;
        int N_inc = 7;

        //eigangle = 1.071734874;       // overwrite the angle
        System.out.println("del =," + d_inc + ", eig0 =," + (float) eig0 + ", eigangle = " + (float) (eigangle*180/Math.PI));
        //d_inc /= 2;                                     // interpolate between data points to do a least squares fit
        System.out.println("scan F:, c, d1, d2, F");
        for (int i = -(N_inc - 1)/2; i <= (N_inc - 1)/2; i++)
        {
            t2[N] = 0;                                  // just to control the output
            solve_at_P2(d1_org + i*d_inc*Math.cos(eigangle), d2_org + i*d_inc*Math.sin(eigangle), false);
            //solve_at_P2(d1_org + i*d_inc*Math.cos(-eigangle), d2_org + i*d_inc*Math.sin(-eigangle), false); // test code fix fix
            //solve_at_P2(d1_org + i*d_inc*Math.cos(eigangle + Math.PI/2), d2_org + i*d_inc*Math.sin(eigangle + Math.PI/2), false);
        }
    }

    private static void calc_alpha(double d1, double d2_org)
    {
        // express: A + B*d1 + C*d2 + D*d1*d2 = 0 (oval solution only)
        double A = 4*((Bezy[3] - Bezy[0])*Math.cos(theta_start) - (Bezx[3] - Bezx[0])*Math.sin(theta_start))
                    *((Bezy[3] - Bezy[0])*Math.cos(theta_end) - (Bezx[3] - Bezx[0])*Math.sin(theta_end));
        double B = -4*((Bezy[3] - Bezy[0])*Math.cos(theta_start) - (Bezx[3] - Bezx[0])*Math.sin(theta_start))*Math.sin(theta_start - theta_end);
        double C =  4*((Bezy[3] - Bezy[0])*Math.cos(theta_end) - (Bezx[3] - Bezx[0])*Math.sin(theta_end))*Math.sin(theta_start - theta_end);
        double D = -3*Math.sin(theta_start - theta_end)*Math.sin(theta_start - theta_end);
        double calc_d2 = -(A + B*d1)/(C + D*d1);                // re-calculate d2 as a double-check
        double alpha1 = (B/2/Math.sin(theta_start - theta_end)/Math.sin(theta_start - theta_end) - 2*d2_org)/d1;
        double alpha2 = d2_org/(C/2/Math.sin(theta_start - theta_end)/Math.sin(theta_start - theta_end) - 2*d1);
        System.out.println("calc_alpha = ," + fitted.getc() + ", " + d1 + ", " + d2_org + ", (" + calc_d2 + "), " + alpha1 + ", (" + alpha2 + ")");

        double[][] oval = new double[][] {{2*(Bezy[3] - Bezy[0])*Math.cos(theta_start) - 2*(Bezx[3] - Bezx[0])*Math.sin(theta_start) - 2*d2_org*Math.sin(theta_end - theta_start), -d1*Math.sin(theta_end - theta_start)},
                                          {d2_org*Math.sin(theta_end - theta_start), 2*(Bezy[3] - Bezy[0])*Math.cos(theta_end) - 2*(Bezx[3] - Bezx[0])*Math.sin(theta_end) + 2*d1*Math.sin(theta_end - theta_start)}};
        System.out.println("oval  det =  ," + fitted.getc() + ", " + d1 + ", " + d2_org + ", " + BSpline5.detm(oval));
        double[][] ovalc = new double[][] {{d1*(1 + Math.cos(theta_end - theta_start)), (Bezx[3] - Bezx[0])*Math.sin(theta_start) - (Bezy[3] - Bezy[0])*Math.cos(theta_start) + d2_org*Math.sin(theta_end - theta_start), -d1*Math.sin(theta_end - theta_start)},
                                           {4*(Bezx[3] - Bezx[0])*(Math.cos(theta_start) + Math.cos(theta_end)) + 4*(Bezy[3] - Bezy[0])*(Math.sin(theta_start) + Math.sin(theta_end)), 3*d2_org*Math.sin(theta_end - theta_start), -3*d1*Math.sin(theta_end - theta_start)},
                                           {d2_org*(1 + Math.cos(theta_end - theta_start)), d2_org*Math.sin(theta_end - theta_start), (Bezx[3] - Bezx[0])*Math.sin(theta_end) - (Bezy[3] - Bezy[0])*Math.cos(theta_end) - d1*Math.sin(theta_end - theta_start)}};
        System.out.println("ovalc det =  ," + fitted.getc() + ", " + d1 + ", " + d2_org + ", " + BSpline5.detm(ovalc));
        double[][] suboval = new double[][] {{ovalc[1][1], ovalc[1][2]},
                                             {ovalc[2][1], ovalc[2][2]}};
        double[] subvec = new double[] {-ovalc[1][0], -ovalc[2][0]};
        double[] gamma = BSpline5.gaussj(suboval, subvec);
        System.out.println("gamma vector = ,1 ," + gamma[0] + ", " + gamma[1]);

        //BSpline5.dump_Jac(oval);
        //BSpline5.dump_Jac(ovalc);

        // calculate det(oval) and det(ovalc) analytically

        double xs1 = (Bezx[3] - Bezx[0])*Math.sin(theta_start);
        double xs2 = (Bezx[3] - Bezx[0])*Math.sin(theta_end);
        double xc1 = (Bezx[3] - Bezx[0])*Math.cos(theta_start);
        double xc2 = (Bezx[3] - Bezx[0])*Math.cos(theta_end);
        double ys1 = (Bezy[3] - Bezy[0])*Math.sin(theta_start);
        double ys2 = (Bezy[3] - Bezy[0])*Math.sin(theta_end);
        double yc1 = (Bezy[3] - Bezy[0])*Math.cos(theta_start);
        double yc2 = (Bezy[3] - Bezy[0])*Math.cos(theta_end);
        double sdel = Math.sin(theta_end - theta_start);
        double cdel = Math.cos(theta_end - theta_start);

        A = 4*(yc1 - xs1)*(yc2 - xs2);
        B = 4*sdel*(yc1 - xs1);
        C = -4*sdel*(yc2 - xs2);
        D = -3*sdel*sdel;
        System.out.println("anal oval  det = ," + fitted.getc() + ", " + A + ", " + B + ", " + C + ", " + D + ", " + (A + B*d1 + C*d2_org + D*d1*d2_org));

        A = -4*(xc1 + xc2 + ys1 + ys2)*(xs1 - yc1)*(xs2 - yc2);
        B =  4*(xc1 + xc2 + ys1 + ys2)*sdel*(xs1 - yc1);
        C = -4*(xc1 + xc2 + ys1 + ys2)*sdel*(xs2 - yc2);
        D =  3*(1 + cdel)*sdel*(xs2 - yc2 - xs1 + yc1);
        System.out.println("anal ovalc det = ," + fitted.getc() + ", " + A + ", " + B + ", " + C + ", " + D + ", " + (A + B*d1 + C*d2_org + D*d1*d2_org));
        System.out.println("alpha vs gamma = ," + fitted.getc() + ", " + d1 + ", " + alpha1 + ", " + gamma[0] + ", " + gamma[1]);
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

    private static double spiro_area()
    {
        // copied from BezierCubicOneDim.java
        double a_b = 180;         // scale factor to make rms error dimensionless
        return a_b*a_b*(Math.PI/2 - Math.sqrt(2))/4                     // <1>
             - fitted.getc()*fitted.getc()*(3*Math.PI/2 - Math.sqrt(2))/4;
    }

    private static double calc_d2(double d1)
    {
        // copied from BezierCubicOneDim.java
        double a_b = 180;         // scale factor to make rms error dimensionless
        // satisfy area constraint
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return (20*spiro_area()/3 - 2*d1*(l1 - l2/Math.sqrt(2)))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }

    private static double calc_dd2dd1(double d1)
    {
        // copied from BezierCubicOneDim.java
        double a_b = 180;         // scale factor to make rms error dimensionless
        // satisfy area constraint (first derivative)
        double l1 = a_b + fitted.getc();       // distance to start point (l1, 0)
        double l2 = a_b - fitted.getc();       // distance to end   point (l2/√2, l2/√2)
        return -(2*(l1 - l2/Math.sqrt(2)) - calc_d2(d1)/Math.sqrt(2))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2));
    }
}
