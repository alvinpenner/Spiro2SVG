
package components;

// Model a saddle point and a minimum in two dimensions.
// Error functional F is cubic in x, quadratic in y.
// Coefficients are Cxy where x is the power of x, y is the power of y

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\streamline.java

import java.io.FileWriter;

public class streamline
{
    private static final double C00 = 1, C10 = 0, C30 = 1, C20 = -C30;
    private static final double C01 = 0, C02 = 100;
    private static final double C11 = Math.sqrt(10);

    public static void main (String[] args)
    {
        //calc_array("F");
        //calc_array("x");
        //calc_array("y");
        //calc_array("dFdx");
        //calc_array("dFdy");
        //solve_initial(0, -0.01);
        //solve_boundary();
        solve_initial_delx(0.1, -0.0035);
        //solve_runge_kutta(0.1, -0.01);
        //solve_theo_y();

        //double alpha = Math.tan(diagonalize(0, 0));         // test alpha and beta at (0,0)
        //double beta = C20*alpha/(C20 + 1.5*alpha*C11 - 0.5*C02);
        //System.out.println("\ninit slope = " + alpha + ", " + beta);
        //for (int i = 0; i <= 5; i++)
        //{
        //    double x = 0.001*i;
        //    double y = alpha*x + 0.5*beta*x*x;
        //    System.out.println(x + ", " + calc_dFdy(x, y)/calc_dFdx(x, y) + ", " + calc_d2ydx2(x, y));
        //}
    }

    private static void calc_array(String type)
    {
        // generate F(x, y) data suitable for MATPLOTLIB
        // to be used by \APP\MATPLOTLIB\streamlines\mystreamline.py

        double min_x = -1, max_x = 2;
        double min_y = -0.3, max_y = 0.3;
        int N_x = 50, N_y = 50;
        double x, y;

        try
        {
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\streamlines\\scan_" + type + ".txt");
            if (type.equals("F"))
            {
                java.text.DateFormat df = java.text.DateFormat.getInstance();
                out.write(df.format(new java.util.Date()) + " : C02, C11 = " + C02 + ", " + C11 + "\r\n");
                out.write("extent = (, " + min_x + ", " + max_x + ", " + min_y + ", " + max_y + ",) N_x = " + N_x + ", N_y = " + N_y + "\r\n");
                double A = C02*C30;
                double B = C02*C20 - C11*C11;
                double C = C02*C10 - C11*C01;
                double x0 = (-B - Math.sqrt(B*B - 4*A*C))/2/A;
                double y0 = (-C01 - C11*x0)/C02;
                double x1 = (-B + Math.sqrt(B*B - 4*A*C))/2/A;
                double y1 = (-C01 - C11*x1)/C02;
                out.write("(x0 y0) (x1 y1), " + x0 + ", " + y0 + ", " + x1 + ", " + y1 + "\r\n");
                diagonalize(x0, y0, true);
                diagonalize(x1, y1, true);
            }
            for (int i = 0; i <= N_y; i++)
            {
                y = min_y + i*(max_y - min_y)/N_y;
                for (int j = 0; j <= N_x; j++)
                {
                    x = min_x + j*(max_x - min_x)/N_x;
                    if (type.equals("F"))
                        out.write(String.format(" %g", calc_F(x, y)));
                    else if (type.equals("x"))
                        out.write(String.format(" %g", x));
                    else if (type.equals("y"))
                        out.write(String.format(" %g", y));
                    else if (type.equals("dFdx"))
                        out.write(String.format(" %g", calc_dFdx(x, y)));
                    else if (type.equals("dFdy"))
                        out.write(String.format(" %g", calc_dFdy(x, y)));
                }
                out.write("\r\n");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("calc_array() save error = " + e);}
    }

    private static void solve_initial(double x, double y)
    {
        // constant increment of s, arc length
        // first-order extrapolation of angle

        int Nstr = 500;
        double del = 0.005;
        double angle, vel;

        //angle = diagonalize(x, y, true);                      // initial angle
        angle = Math.atan2(-calc_dFdy(x, y), -calc_dFdx(x, y));  // we wish to decrease F
        java.text.DateFormat df = java.text.DateFormat.getInstance();
        System.out.println(" " + df.format(new java.util.Date()));
        System.out.println("C20  C02  C11 = , " + C20 + ", " + C02 + ", " + C11);
        System.out.println("initial value = , " + x + ", " + y + ", " + Nstr + ", " + del);
        for (int i = 0; i < Nstr; i++)
        {
            //vel = Math.sqrt(calc_dFdx(x, y)*calc_dFdx(x, y) + calc_dFdy(x, y)*calc_dFdy(x, y));
            //del = 0.0001*vel;
            //if (del < 0.00001) del = 0.00001;
            if (i > 0)
                angle = Math.atan2(-calc_dFdy(x, y), -calc_dFdx(x, y));  // we wish to decrease F
            //System.out.println(i + ", " + angle);
            x += del*Math.cos(angle);
            y += del*Math.sin(angle);
            System.out.println("                , " + x + ", " + y);
        }
    }

    private static void solve_initial_delx(double x, double y)
    {
        // constant increment of x
        // second-order extrapolation of y
        // see Froberg, p. 266
        // see Spiro2SVG Book 8, page 56

        int Nstr = 400;
        double delx= 0.0025;
        double slope = Math.tan(diagonalize(x, y, true));                   // initial slope
        double curv = C20*slope/(C20 + 1.5*slope*C11 - 0.5*C02);            // initial curvature
        if (x != 0 || y != 0)
        {
            slope = calc_dFdy(x, y)/calc_dFdx(x, y);
            curv = calc_d2ydx2(x, y);
        }
        //double yold, ynew;

        java.text.DateFormat df = java.text.DateFormat.getInstance();
        System.out.println("initial delx - " + df.format(new java.util.Date()));
        System.out.println("C20  C02  C11 Nstr delx = , " + C20 + ", " + C02 + ", " + C11 + ", " + Nstr + ", " + delx);
        //System.out.println("initial value = , " + x + ", " + y + ", " + slope + ", " + curv + ", " + Nstr + ", " + delx);
        System.out.println(" , x, y, slope, tan theta0, curv, F");
        System.out.println("initial value = , " + x + ", " + y + ", " + slope + ", " + Math.tan(diagonalize(x, y, false)) + ", " + curv + ", " + calc_F(x, y));
        //yold = y;
        //y = yold + slope*delx;
        y += slope*delx + 0.5*curv*delx*delx;
        x += delx;
        for (int i = 1; i < Nstr; i++)
        {
            slope = calc_dFdy(x, y)/calc_dFdx(x, y);
            curv = calc_d2ydx2(x, y);
            System.out.println("                , " + x + ", " + y + ", " + slope + ", " + Math.tan(diagonalize(x, y, false)) + ", " + curv + ", " + calc_F(x, y));
            //ynew = yold + 2*slope*delx;
            //yold = y;
            //y = ynew;
            y += slope*delx + 0.5*curv*delx*delx;
            x += delx;
        }
        System.out.println("                , " + x + ", " + y);
    }

    private static void solve_boundary()
    {
        final int MAXLOOP = 10;
        final int N = 10;                             // number of data points (must be even)
        double[] x = new double[N];
        double[] y = new double[N];
        double[] dydx = new double[N];
        double[][] arr = new double[N][N];
        double[] vec = new double[N];
        int i, loop = 0;
        double rms;

        java.text.DateFormat df = java.text.DateFormat.getInstance();
        System.out.println("boundary-value : " + df.format(new java.util.Date()));
        System.out.println("C20  C02  C11 = , " + C20 + ", " + C02 + ", " + C11);

        double A = C02*C30;
        double B = C02*C20 - C11*C11;
        double C = C02*C10 - C11*C01;
        double x0 = (-B - Math.sqrt(B*B - 4*A*C))/2/A;
        double y0 = (-C01 - C11*x0)/C02;
        double x1 = (-B + Math.sqrt(B*B - 4*A*C))/2/A;
        double y1 = (-C01 - C11*x1)/C02;
        double delx = 2*(x1 - x0)/(N - 1);                  // this is x[i + 1] - x[i - 1]
        double theta0 = diagonalize(x0, y0, true);
        double theta1 = diagonalize(x1, y1, true);
        System.out.println("(x0 y0) (x1 y1), " + x0 + ", " + y0 + ", " + x1 + ", " + y1 + ", " + theta0*180/Math.PI + ", " + theta1*180/Math.PI);

        // initiallize

        for (i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                arr[i][j] = 0;
        arr[0][0] = 1;
        arr[N - 1][N - 1] = 1;
        vec[0] = y0;
        vec[N - 1] = y1;
        double beta = C20*Math.tan(theta0)/(C20 + 1.5*C11*Math.tan(theta0) - 0.5*C02);
        System.out.println("alpha        = " + Math.tan(theta0));
        //System.out.println("guessed beta = " + 2*0.98*(y1 - x1*Math.tan(theta0))/x1/x1);
        System.out.println("calc    beta = " + beta);
        System.out.println("\ninit");
        System.out.println(" , x, y, dydx, calc_dFdx, calc_dFdy");
        for (i = 0; i < N; i++)
        {
            x[i] = x0 + i*delx/2;
            //y[i] = y0 + i*(y1 - y0)/(N - 1);                            // interpolate y
            //y[i] = y0 + Math.tan(theta0)*x[i] + 0.98*(y1 - x1*Math.tan(theta0))*x[i]*x[i]/x1/x1;
            y[i] = y0 + Math.tan(theta0)*x[i] + beta*x[i]*x[i]/2;
        }
        for (i = 1; i < N - 1; i++)
            dydx[i] = calc_dFdy(x[i], y[i])/calc_dFdx(x[i], y[i]);
        for (i = 0; i < N; i++)
            System.out.println(i + ", " + x[i] + ", " + y[i] + ", " + dydx[i]);

        // iterate

        do
        {
            loop++;
//            for (i = 2; i < N - 1; i += 2)
//            {
//                y[i] = y[i - 2] + dydx[i - 1]*(x[i] - x[i - 2]);
//                y[N - 1 - i] = y[N + 1 - i] + dydx[N - i]*(x[N - 1 - i] - x[N + 1 - i]);
//            }
            i = 1;
            arr[i][i - 1] = calc_ddydydx(x[i], y[i])/delx;
            arr[i][i] = 1/delx/delx;
            arr[i][i + 1] = -calc_ddydydx(x[i], y[i])/delx;
            arr[i][i + 2] = -1/delx/delx;
            vec[i] = -dydx[i]*calc_ddydydx(x[i], y[i]) - dydx[i + 1]/delx;
            for (i = 2; i < N - 2; i++)
            {
                arr[i][i - 2] = -1/delx/delx;
                arr[i][i - 1] = calc_ddydydx(x[i], y[i])/delx;
                arr[i][i] = 2/delx/delx;
                arr[i][i + 1] = -calc_ddydydx(x[i], y[i])/delx;
                arr[i][i + 2] = -1/delx/delx;
                vec[i] = -dydx[i]*calc_ddydydx(x[i], y[i]) + dydx[i - 1]/delx - dydx[i + 1]/delx;
            }
            i = N - 2;
            arr[i][i - 2] = -1/delx/delx;
            arr[i][i - 1] = calc_ddydydx(x[i], y[i])/delx;
            arr[i][i] = 1/delx/delx;
            arr[i][i + 1] = -calc_ddydydx(x[i], y[i])/delx;
            vec[i] = -dydx[i]*calc_ddydydx(x[i], y[i]) + dydx[i - 1]/delx;
            y = BSpline5.gaussj(arr, vec);
            //for (i = 0; i < N; i++)
            //{
            //    for (int j = 0; j < N; j++)
            //        System.out.print(arr[i][j] + ", ");
            //    System.out.println();
            //}
            //for (i = 0; i < N; i++)
            //    System.out.println(vec[i]);

            rms = 0;
            for (i = 1; i < N - 1; i++)
            {
                dydx[i] = calc_dFdy(x[i], y[i])/calc_dFdx(x[i], y[i]);
                rms += ((y[i + 1] - y[i - 1])/delx - dydx[i])*((y[i + 1] - y[i - 1])/delx - dydx[i]);
            }
            rms = Math.sqrt(rms/(N - 2));
            System.out.println("\niter " + loop + ": rms = " + (float) rms);
            System.out.println(" , x, y, dydx, calc_dFdx, calc_dFdy");
            for (i = 0; i < N; i++)
                System.out.println(i + ", " + x[i] + ", " + y[i] + ", " + dydx[i] + ", " + calc_dFdx(x[i], y[i]) + ", " + calc_dFdy(x[i], y[i]));
        } while (loop < MAXLOOP);
    }

    private static void solve_runge_kutta(double x, double y)
    {
        int Nstr = 200;
        double delx = 0.005;

        java.text.DateFormat df = java.text.DateFormat.getInstance();
        diagonalize(x, y, true);
        System.out.println("Runge-Kutta  - " + df.format(new java.util.Date()));
        System.out.println("C20  C02  C11 Nstr delx = , " + C20 + ", " + C02 + ", " + C11 + ", " + Nstr + ", " + delx);
        System.out.println(" , x, y, slope, tan theta0, curv, F");
        System.out.println("initial value = , " + x + ", " + y + ", " + Math.tan(diagonalize(x, y, false)) + ", " + Math.tan(diagonalize(x, y, false)) + ", " + C20*Math.tan(diagonalize(x, y, false))/(C20 + 1.5*Math.tan(diagonalize(x, y, false))*C11 - 0.5*C02) + ", " + calc_F(x, y));
        for (int i = 0; i < Nstr; i++)
        {
            y += runge_kutta(x, y, delx);
            x += delx;
            System.out.println("                , " + x + ", " + y + ", " + calc_dFdy(x, y)/calc_dFdx(x, y) + ", " + Math.tan(diagonalize(x, y, false)) + ", " + calc_d2ydx2(x, y) + ", " + calc_F(x, y));
        }
    }

    private static double runge_kutta(double x0, double y0, double h)
    {
        // perform one Runge-Kutta iteration
        // see Froberg, p. 268
        // h = del x, return k = del y

        double k1;
        if (x0 == 0 && y0 == 0)
            k1 = h*Math.tan(diagonalize(x0, y0, false));
        else
            k1 = h*calc_dFdy(x0, y0)/calc_dFdx(x0, y0);
        double k2 = h*calc_dFdy(x0 + h/2, y0 + k1/2)/calc_dFdx(x0 + h/2, y0 + k1/2);
        double k3 = h*calc_dFdy(x0 + h/2, y0 + k2/2)/calc_dFdx(x0 + h/2, y0 + k2/2);
        double k4 = h*calc_dFdy(x0 + h, y0 + k3)/calc_dFdx(x0 + h, y0 + k3);
        //System.out.println(k1 + ", " + k2 + ", " + k3 + ", " + k4);
        return (k1 + 2*k2 + 2*k3 + k4)/6;
    }

    private static double diagonalize(double x, double y, boolean print)
    {
        //System.out.println("\nC20 ,C02, C11 = " + C20 + ", " + C02 + ", " + C11);
        if (print) System.out.println("\ncritical point (" + x + ", " + y + ") F = " + calc_F(x, y));
        if (print) System.out.println("dFdx, dFdy = (" + calc_dFdx(x, y) + ", " + calc_dFdy(x, y) + ")");
        double M00 = C20 + 2*C30*x;
        double M01 = C11;
        double M11 = C02;
        double eig0 = (M00 + M11 - Math.sqrt((M00 - M11)*(M00 - M11) + 4*M01*M01))/2;
        double eig1 = (M00 + M11 + Math.sqrt((M00 - M11)*(M00 - M11) + 4*M01*M01))/2;
        if (print) System.out.println("eig0, eig1 = (" + eig0 + ", " + eig1 + ")");
        double theta0 = Math.atan(-M01/(M11 - eig0));
        if (print) System.out.println("theta0 =  " + theta0); //*180/Math.PI);
        return theta0;
    }

    private static void solve_theo_y()
    {
        // theoretical y based on following the eigenvector angle
        // see Spiro2SVG Book 8, page 60

        int Nstr = 440;
        double delx= 0.0025;
        double x, y;

        java.text.DateFormat df = java.text.DateFormat.getInstance();
        System.out.println("analytical y - " + df.format(new java.util.Date()));
        System.out.println("C20  C02  C11 Nstr delx = , " + C20 + ", " + C02 + ", " + C11 + ", " + Nstr + ", " + delx);
        System.out.println(" , x, y, F");
        for (int i = 0; i <= Nstr; i++)
        {
            x = delx*i;
            y = -C11*(integrate_z((2*C30*x + C20 - C02)/2/C11) - integrate_z((C20 - C02)/2/C11))/2/C30;
            System.out.println(" , " + x + ", " + y + ", " + calc_F(x, y));
        }
    }

    private static double integrate_z(double z)
    {
        return z*z + z*Math.sqrt(z*z + 1) + Math.log(z + Math.sqrt(z*z + 1));
    }

    private static double calc_F(double x, double y)
    {
        return C00 + C10*x + C20*x*x/2 + C30*x*x*x/3
                   + C01*y + C02*y*y/2
                           + C11*x*y;
    }

    private static double calc_dFdx(double x, double y)
    {
        return C10 + C20*x + C30*x*x
                   + C11*y;
    }

    private static double calc_dFdy(double x, double y)
    {
        return C01 + C02*y
                   + C11*x;
    }

    private static double calc_d2Fdx2(double x)
    {
        return C20 + 2*C30*x;
    }

    private static double calc_d2Fdxdy()
    {
        return C11;
    }

    private static double calc_d2Fdy2()
    {
        return C02;
    }

    private static double calc_ddydydx(double x, double y)
    {
        // d/dy(dFdy/dFdx)
        return (calc_dFdx(x, y)*calc_d2Fdy2() - calc_dFdy(x, y)*calc_d2Fdxdy())/calc_dFdx(x, y)/calc_dFdx(x, y);
    }

    private static double calc_d2ydx2(double x, double y)
    {
        return (calc_dFdx(x, y)*(calc_d2Fdxdy() + calc_d2Fdy2()*calc_dFdy(x, y)/calc_dFdx(x, y))
              - calc_dFdy(x, y)*(calc_d2Fdx2(x) + calc_d2Fdxdy()*calc_dFdy(x, y)/calc_dFdx(x, y)))/calc_dFdx(x, y)/calc_dFdx(x, y);
    }
}
