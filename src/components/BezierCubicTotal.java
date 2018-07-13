
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 4-point cubic Bezier (P0 - P3) to it, using parameter 0 < t2 < 1.
// constrain only the slopes at the endpoints
// Bezier[4] = f(x0, x1, x2, x3, t2)

// this is an implementation of the "Total" method discussed by Ahn
// in which all the variables, including all t2(t1) are simultaneously solved for

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BezierCubicTotal.java

//import java.io.FileWriter;

public class BezierCubicTotal
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI/4;
    public static final int N = 10;
    public static double[] Bezx;            // cubic Bezier, 4 x-coordinates
    public static double[] Bezy;            // cubic Bezier, 4 y-coordinates
    private static epiTrochoidFxn fitted;
    private static double[] t2;             // t2[0] = d1
                                            // t2[1] = d2
                                            // t2[i = 2-N] = u(t[i-1])
    private static double Jacdet = Double.NaN;
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        //long startTime = System.currentTimeMillis();
        fitted = new epiTrochoidFxn(9.1);
        //t2 = new double[] {15.9753, 92.2577, 0.1945, 0.3346, 0.4501, 0.5510, 0.6418, 0.7250, 0.8018, 0.8729, 0.9389};
        t2 = new double[] {7.6036, 79.2178, 0.2367, 0.3672, 0.4729, 0.5657, 0.6502, 0.7288, 0.8024, 0.8720, 0.9377};
        //t2 = new double[] {65, 25, 0.2296, 0.3599, 0.4660, 0.5595, 0.6448, 0.7241, 0.7987, 0.8693, 0.9363};
        iterate_at_P2();
        //System.out.println("cubic Bezier (Total) solve_at_P2 = " + solve_at_P2(true) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //System.out.println("Elapsed Time = " + (System.currentTimeMillis() - startTime));
    }

    private static void iterate_at_P2()
    {
        // calculate a new estimate of (d1, d2, t2[]) by setting dF = 0
        // see Spiro2SVG Book 6, page 6 (applied to cubic Bezier)
        // setup N+1 variables in Newton-Raphson iteration

        final double gain = 1;                          // fudge factor to reduce gain
        final int MAXLOOP = 100;
        double[] f_gx = new double[N-1];
        double[] f_gy = new double[N-1];
        double[] dfxdu = new double[N-1];
        double[] dfydu = new double[N-1];
        double[] d2fxdudu = new double[N-1];
        double[] d2fydudu = new double[N-1];
        double[][] dfxdd = new double[2][N-1];
        double[][] dfydd = new double[2][N-1];
        double[][] d2fxdudd = new double[2][N-1];
        double[][] d2fydudd = new double[2][N-1];
        double[] df_gxdc = new double[N-1];             // used only for calc of dx[]/dc
        double[] df_gydc = new double[N-1];
        double[] d2fxdudc = new double[N-1];            // used only for calc of dx[]/dc
        double[] d2fydudc = new double[N-1];

        double[] trap_in = new double[N+1];
        double[][] Jac = new double[N+1][N+1];
        double[] dFdd = new double[N+1];
        double[] deld;                                  // (-Δd1, -Δd2)
        //double[][] dJacdc = new double[N+1][N+1];       // used only for calc of dx[]/dc
        double[] d2Fdddc = new double[N+1];
        int i, j, k, loop = 0;
        double t1, maxdel;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(false)))       // initiallize at (x2, y2)
            {
                System.out.println("fail at (d1, d2): " + t2[0] + ", " + t2[1]);
                return;
            }

            for (i = 0; i < N - 1; i++)
            {
                t1 = t1_start + (i + 1)*(t1_end - t1_start)/N;
                f_gx[i] = t2_vs_t1.fn(Bezx, t2[i + 2]) - fitted.getx(t1);
                f_gy[i] = t2_vs_t1.fn(Bezy, t2[i + 2]) - fitted.gety(t1);
                dfxdu[i] = t2_vs_t1.dfn(Bezx, t2[i + 2]);
                dfydu[i] = t2_vs_t1.dfn(Bezy, t2[i + 2]);
                d2fxdudu[i] = t2_vs_t1.d2fn(Bezx, t2[i + 2]);
                d2fydudu[i] = t2_vs_t1.d2fn(Bezy, t2[i + 2]);
                df_gxdc[i] = fitted.getdxdc(t1_start)*(N33(t2[i + 2])[0] + N33(t2[i + 2])[1]) + fitted.getdxdc(t1_end)*(N33(t2[i + 2])[2] + N33(t2[i + 2])[3]) - fitted.getdxdc(t1);
                df_gydc[i] = fitted.getdydc(t1_start)*(N33(t2[i + 2])[0] + N33(t2[i + 2])[1]) + fitted.getdydc(t1_end)*(N33(t2[i + 2])[2] + N33(t2[i + 2])[3]) - fitted.getdydc(t1);
                d2fxdudc[i] = fitted.getdxdc(t1_start)*(dN33(t2[i + 2])[0] + dN33(t2[i + 2])[1]) + fitted.getdxdc(t1_end)*(dN33(t2[i + 2])[2] + dN33(t2[i + 2])[3]);
                d2fydudc[i] = fitted.getdydc(t1_start)*(dN33(t2[i + 2])[0] + dN33(t2[i + 2])[1]) + fitted.getdydc(t1_end)*(dN33(t2[i + 2])[2] + dN33(t2[i + 2])[3]);

                dfxdd[0][i] = Math.cos(theta_start)*N33(t2[i + 2])[1];
                dfydd[0][i] = Math.sin(theta_start)*N33(t2[i + 2])[1];
                dfxdd[1][i] = -Math.cos(theta_end)*N33(t2[i + 2])[2];
                dfydd[1][i] = -Math.sin(theta_end)*N33(t2[i + 2])[2];
                d2fxdudd[0][i] = Math.cos(theta_start)*dN33(t2[i + 2])[1];
                d2fydudd[0][i] = Math.sin(theta_start)*dN33(t2[i + 2])[1];
                d2fxdudd[1][i] = -Math.cos(theta_end)*dN33(t2[i + 2])[2];
                d2fydudd[1][i] = -Math.sin(theta_end)*dN33(t2[i + 2])[2];
            }

            // calc dFdd[j] at current (d1, d2)

            trap_in[0] = 0;
            trap_in[N] = 0;
            for (i = 0; i < 2; i++)
            {
                for (k = 0; k < N - 1; k++)
                    trap_in[k + 1] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];
                dFdd[i] = sum(trap_in);
                for (k = 0; k < N - 1; k++)
                    trap_in[k + 1] = df_gxdc[k]*dfxdd[i][k] + df_gydc[k]*dfydd[i][k];
                d2Fdddc[i] = sum(trap_in);
            }
            for (i = 0; i < N - 1; i++)
            {
                dFdd[i + 2] = f_gx[i]*dfxdu[i] + f_gy[i]*dfydu[i];
                d2Fdddc[i + 2] = df_gxdc[i]*dfxdu[i] + df_gydc[i]*dfydu[i]
                               + f_gx[i]*d2fxdudc[i] + f_gy[i]*d2fydudc[i];
            }

            // calc d2Fdd[i]dd[j] (Jacobean matrix)

            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k < N - 1; k++)
                    {
                        trap_in[k + 1] = dfxdd[i][k]*dfxdd[j][k] + dfydd[i][k]*dfydd[j][k];
                        //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = sum(trap_in);
                    //dJacdc[i][j] = 0;
                }
            for (i = 0; i < N - 1; i++)
                for (j = 0; j < N - 1; j++)
                    if (i == j)
                    {
                        Jac[i + 2][i + 2] = dfxdu[i]*dfxdu[i] + f_gx[i]*d2fxdudu[i]
                                          + dfydu[i]*dfydu[i] + f_gy[i]*d2fydudu[i];
                        //dJacdc[i + 2][j + 2] = -dgxdc[i]*d2fxdudu[i] - dgydc[i]*d2fydudu[i];
                    }
                    else
                    {
                        Jac[i + 2][j + 2] = 0;
                        //dJacdc[i + 2][j + 2] = 0;
                    }
            for (i = 0; i < 2; i++)
                for (j = 0; j < N - 1; j++)
                {
                    Jac[i][j + 2] = dfxdd[i][j]*dfxdu[j] + f_gx[j]*d2fxdudd[i][j]
                                  + dfydd[i][j]*dfydu[j] + f_gy[j]*d2fydudd[i][j];
                    Jac[j + 2][i] = Jac[i][j + 2];
                    //dJacdc[i][j + 2] = -dgxdc[j]*d2fxdudd[i][j] - dgydc[j]*d2fydudd[i][j];
                    //dJacdc[j + 2][i] = dJacdc[i][j + 2];
                }

            deld = BSpline5.gaussj(Jac, dFdd);      // this is actually the negative of Δd
            for (i = 0; i <= N; i++)
                t2[i] -= deld[i]/gain;

            //Jacdet = BSpline5.detm(Jac);
            System.out.print("dFdd = ");
            for (i = 0; i <= N; i++) System.out.print(", " + dFdd[i]);
            System.out.println();
            System.out.print("deld = ");
            for (i = 0; i <= N; i++) System.out.print(", " + deld[i]);
            System.out.println();
            maxdel = 0;
            for (i = 0; i <= N; i++)
                if (Math.abs(deld[i]) > maxdel)
                    maxdel = Math.abs(deld[i]);
        } while ((loop < MAXLOOP) && (maxdel > TOL));
        if (loop < MAXLOOP)
        {
            BSpline5.dump_Jac(Jac);
            //Jacdet = BSpline5.detm(Jac);
            System.out.println("\n__converged in " + loop + " at new d1 d2 = , , , , , , " + t2[0] + ", " + t2[1] + ", " + fitted.getc());
            for (i = 0; i < N; i++) System.out.printf("%.4f, ", t2[i]);
            System.out.printf("%.4f\n", t2[N]);
            solve_at_P2(true);                          // final run just for good measure
            // calculate -dL/dc - (dM/dc)*x
            //BSpline5.dump_Jac(dJacdc);
            double[] dt2dc = BSpline5.gaussj(Jac, d2Fdddc);
            System.out.println("\ncheck -dL/dc - (dM/dc)*x at c = ," + fitted.getc());
            for (i = 0; i <= N; i++)
                //System.out.println(i + ", " + dFdd[i] + ", " + d2Fdddc[i]);
                System.out.println(i + ", " + t2[i] + ", " + d2Fdddc[i] + ", " + -dt2dc[i]);
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ")");
    }

    private static double solve_at_P2(boolean print)
    {
        // 4-point cubic Bezier curve
        // perform a single calculation of a complete (d1, d2, t2[]) profile
        // and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        Bezx = new double[] {fitted.getx(t1_start),
                             fitted.getx(t1_start) + t2[0]*Math.cos(theta_start),
                             fitted.getx(t1_end) - t2[1]*Math.cos(theta_end),
                             fitted.getx(t1_end)};
        Bezy = new double[] {fitted.gety(t1_start),
                             fitted.gety(t1_start) + t2[0]*Math.sin(theta_start),
                             fitted.gety(t1_end) - t2[1]*Math.sin(theta_end),
                             fitted.gety(t1_end)};

        if (t2[0] < 0 || t2[1] < 0)
            System.out.println("WARNING: negative arm length = " + t2[0] + ", " + t2[1]);
        //fitted.gen_Bezier(new double[] {Bezx[0], Bezy[0], Bezx[1], Bezy[1], Bezx[2], Bezy[2], Bezx[3], Bezy[3]});
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);

        if (print)
        {
            System.out.println("\n           , t1, t2, del x, del y");
            System.out.println("cubic Bez (Total), " + t1_start + ", " + "0" + ", " + (t2_vs_t1.fn(Bezx, 0) - fitted.getx(t1_start)) + ", " + (t2_vs_t1.fn(Bezy, 0) - fitted.gety(t1_start)));
            for (int i = 1; i < N; i++)
            {
                double t1 = t1_start + i*(t1_end - t1_start)/N;
                System.out.println("cubic Bez (Total), " + t1 + ", " + t2[i+1] + ", " + (t2_vs_t1.fn(Bezx, t2[i+1]) - fitted.getx(t1)) + ", " + (t2_vs_t1.fn(Bezy, t2[i+1]) - fitted.gety(t1)));
            }
            System.out.println("cubic Bez (Total), " + t1_end + ", " + "1" + ", " + (t2_vs_t1.fn(Bezx, 1) - fitted.getx(t1_end)) + ", " + (t2_vs_t1.fn(Bezy, 1) - fitted.gety(t1_end)));
        }
        double retVal = calc_error();
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + t2[0] + ", " + t2[1] + ", " + retVal + ", " + (float) Jacdet);
        return retVal;
    }

    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

        double a_b = 180;               // scale factor to make rms error dimensionless
        //double a_b = 1;               // Cycloid only
        double t1 = t1_start;
        double[] trap_in = new double[N+1];

        trap_in[0] = 0;
        for (int i = 1; i < N; i++)
        {
            t1 += (t1_end - t1_start)/N;
            trap_in[i] = (t2_vs_t1.fn(Bezx, t2[i+1]) - fitted.getx(t1))*(t2_vs_t1.fn(Bezx, t2[i+1]) - fitted.getx(t1))
                       + (t2_vs_t1.fn(Bezy, t2[i+1]) - fitted.gety(t1))*(t2_vs_t1.fn(Bezy, t2[i+1]) - fitted.gety(t1));
        }
        trap_in[N] = 0;
        return Math.sqrt(t2_vs_t1.integrate(trap_in))/a_b;
    }

    private static double sum(double[] trap)
    {
        // based on trapezoidal rule integration of a fxn of t1 (N+1 points)

        //System.out.println("trap length = " + trap.length);
        double ret = (trap[0] + trap[trap.length - 1])/2;
        for (int i = 1; i < trap.length - 1; i++)
            ret += trap[i];
        return ret;
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
