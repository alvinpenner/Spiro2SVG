
package spiro;

import java.awt.geom.Point2D;
import java.awt.geom.CubicCurve2D;
import java.util.Arrays;

//  fit Bezier to SpiroJ shape by matching slope and curvature at endpoints
//  Define slope and curvature using “Calculus” by James Stewart, page 902
//  Slope     : m = y′/x′, where all derivatives (′) are with respect to 't'
//  Curvature : κ = (x′y″ - y′x″)/((x′)^2 + (y′)^2)^(3/2)

public final class SpiroJCalc
{
    private static final double TOL = 0.0000001;
    private static double rx1, ry1, wx1, wy1, rx2, ry2, wx2, wy2;       // SpiroJ parameters
    private static boolean isCCW = true;                                // change sign at stationary point

    protected static CubicCurve2D.Float getBezier(double t1, double t2)
    {
        Point2D.Double[][] ptSpiro = new Point2D.Double[3][2];          // Point[derivative (0-2)][t = (0,1)]

        ptSpiro[0][0] = new Point2D.Double(getX(t1), getY(t1));
        ptSpiro[0][1] = new Point2D.Double(getX(t2), getY(t2));
        ptSpiro[1][0] = new Point2D.Double(getdX(t1), getdY(t1));
        ptSpiro[1][1] = new Point2D.Double(getdX(t2), getdY(t2));
        ptSpiro[2][0] = new Point2D.Double(getd2X(t1), getd2Y(t1));
        ptSpiro[2][1] = new Point2D.Double(getd2X(t2), getd2Y(t2));

        if (main.IS_DEBUG)
            System.out.println("");
        if (t1 == 0)
        {
            if ((Math.abs(ptSpiro[1][0].x) > TOL) || (Math.abs(ptSpiro[1][0].y) > TOL))
                main.theta[0] = Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x);
            else
            {
                if (main.IS_DEBUG)
                    System.out.println("spiroJ motion is stationary at t = " + t1);
                main.theta[0] = Math.atan2(ptSpiro[2][0].y, ptSpiro[2][0].x);   // use x″ & y″
            }
        }
        else
        {
            main.theta[0] = main.theta[1];
            if ((Math.abs(ptSpiro[1][0].x) < TOL) && (Math.abs(ptSpiro[1][0].y) < TOL))
            {
                if (main.IS_DEBUG)
                    System.out.println("spiroJ motion is stationary at t = " + t1);
                main.theta[0] -= isCCW ? Math.PI : -Math.PI;                // reverse direction
                isCCW = !isCCW;
            }
        }
        if ((Math.abs(ptSpiro[1][1].x) > TOL) || (Math.abs(ptSpiro[1][1].y) > TOL))
            main.theta[1] = Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x);
        else                                                                // point is stationary
        {
            if (main.IS_DEBUG)
                System.out.println("spiroJ motion is stationary at t = " + t2);
            main.theta[1] = Math.atan2(ptSpiro[2][1].y, ptSpiro[2][1].x);   // use x″ & y″
            main.theta[1] += isCCW ? Math.PI : -Math.PI;
        }

        if (isCCW)
            main.theta[1] = main.theta[0] + Math.PI + (main.theta[1] - main.theta[0] - Math.PI) % (-2*Math.PI);
        else
            main.theta[1] = main.theta[0] - Math.PI + (main.theta[1] - main.theta[0] + Math.PI) % (2*Math.PI);
        if (Math.abs(Math.abs(main.theta[0] - main.theta[1]) - Math.PI) < TOL)
            main.theta[1] = main.theta[0];                                  // anti-symmetric case

        for (int i = 0; i < ptSpiro[0].length; i++)                         // standard curvature
            if ((Math.abs(ptSpiro[1][i].x) > TOL) || (Math.abs(ptSpiro[1][i].y) > TOL))
                main.Cu[i] = (ptSpiro[1][i].x*ptSpiro[2][i].y - ptSpiro[1][i].y*ptSpiro[2][i].x)
                           /  Math.pow((ptSpiro[1][i].x*ptSpiro[1][i].x + ptSpiro[1][i].y*ptSpiro[1][i].y), 1.5);
            else
                main.Cu[i] = rx1*ry1*wx1*wx1*wy1*wy1*(wx1*wx1 - wy1*wy1)/3  // stationary point
                           / Math.pow(rx1*rx1*wx1*wx1*wx1*wx1 + ry1*ry1*wy1*wy1*wy1*wy1, 1.5);
        return main.calcBezier(ptSpiro, t1, t2, 1);
    }

    protected static int get_t_values(double[] t, double m_rx1, double m_ry1, double m_wx1, double m_wy1, double m_rx2, double m_ry2, double m_wx2, double m_wy2)
    {
        int i, N = 0;                                   // N = fit points per object
        rx1 = m_rx1; ry1 = m_ry1; wx1 = m_wx1; wy1 = m_wy1;
        rx2 = m_rx2; ry2 = m_ry2; wx2 = m_wx2; wy2 = m_wy2;

        double[][] rotors = new double[][] {{rx1*wx1, wx1, -Math.PI/2},
                                            {rx2*wx2, wx2, -Math.PI/2}};
        for (i = 0; i < rotors.length; i++)                                         // solve x′ = 0
            if ((Math.abs(rotors[i][0]) > TOL) && (Math.abs(rotors[i][1]) > TOL))   // check amplitude and frequency
                for (int j = 0; j < Math.round(Math.abs(2*rotors[i][1])); j++)
                    N = main.insert_t_value(N, N, t, solve_cos_t(rotors, Math.PI*j/Math.abs(rotors[i][1])));

        rotors = new double[][] {{ry1*wy1, wy1, 0},
                                 {ry2*wy2, wy2, 0}};
        for (i = 0; i < rotors.length; i++)                                         // solve y′ = 0
            if ((Math.abs(rotors[i][0]) > TOL) && (Math.abs(rotors[i][1]) > TOL))   // check amplitude and frequency
                for (int j = 0; j < Math.round(Math.abs(2*rotors[i][1])); j++)
                    N = main.insert_t_value(N, N, t, solve_cos_t(rotors, Math.PI*(j + 0.5)/Math.abs(rotors[i][1])));

        N = SpiroJCalc.sort_t_values(t, N);
        if (main.IS_DEBUG)
        {
            System.out.println("sorted solutions for x/y extrema N = " + N);
            for (i = 0; i < N; i++)
                System.out.println(i + ", " + t[i]);
        }

        rotors = new double[][] {{rx1*ry1*wx1*wy1*(wx1 + wy1)/2, wx1 - wy1, 0},
                                 {rx1*ry1*wx1*wy1*(wx1 - wy1)/2, wx1 + wy1, 0},
                                 {rx1*ry2*wx1*wy2*(wx1 + wy2)/2, wx1 - wy2, 0},
                                 {rx1*ry2*wx1*wy2*(wx1 - wy2)/2, wx1 + wy2, 0},
                                 {rx2*ry1*wx2*wy1*(wx2 + wy1)/2, wx2 - wy1, 0},
                                 {rx2*ry1*wx2*wy1*(wx2 - wy1)/2, wx2 + wy1, 0},
                                 {rx2*ry2*wx2*wy2*(wx2 + wy2)/2, wx2 - wy2, 0},
                                 {rx2*ry2*wx2*wy2*(wx2 - wy2)/2, wx2 + wy2, 0}};
        int N_old = N;                              // previous number of points (cusps or inflections)
        boolean tooclose;
        double temp_t;                              // temporary t value
        for (i = 0; i < rotors.length; i++)                                         // get inflection points
            if ((Math.abs(rotors[i][0]) > TOL) && (Math.abs(rotors[i][1]) > TOL))   // check amplitude and frequency
                for (int j = 0; j < Math.round(Math.abs(2*rotors[i][1])); j++)
                {
                    temp_t = solve_cos_t(rotors, Math.PI*(j + 0.5)/Math.abs(rotors[i][1]));
                    tooclose = false;
                    for (int k = 0; k < N_old; k++)
                        if ((Math.abs(temp_t - t[k]) < 100*TOL)                     // compare inflections to extrema
                        ||  (Math.abs(temp_t - t[k] - 2*Math.PI) < 100*TOL))        // with a very loose tolerance
                            tooclose = true;
                    if (tooclose)
                    {
                        if (main.IS_DEBUG)
                            System.out.println("inflection at " + temp_t + " is too close to x/y extremum");
                    }
                    else
                        N = main.insert_t_value(N, N, t, temp_t);
                }
//        for (i = 0; i < rotors.length; i++)
//            System.out.println(i + ", " + rotors[i][0] + ", " + rotors[i][1] + ", " + rotors[i][2]);
        if (main.IS_DEBUG)
            System.out.println("initial t array N = " + N);
        return sort_t_values(t, N);
    }

    protected static int sort_t_values(double[] t, int n)
    {
        int i;
        Arrays.sort(t, 0, n);
        if (n > 1)                                                      // check for duplicates
            for (i = 0; i < n - 1; i++)
                if (t[i + 1] < t[i] + TOL)
                    t[i] = 999999;
        Arrays.sort(t, 0, n);
        for (i = 0; i < n; i++)                                         // remove duplicates
            if (t[i] == 999999)
                break;
        return i;
    }

    protected static double solve_cos_t(double r[][], double t0)
    {
        // solve the equation : r[0] + r[1] + ... + r[n] = 0
        // where r[] is a sequence of rotors:
        // function r[i] = A[i]*cos(w[i]*t + phi[i])
        // elements r[i][] = {A[i], w[i], phi[i]}
        // t0 = initial estimate of t

        double f, fprime, del_t;
        int loop = 0;

        t0 = (t0 + 2*Math.PI) % (2*Math.PI);            // in case of negative phase shift
        if (Math.abs(t0) < TOL) t0 = 0;

        if (main.IS_DEBUG)
            System.out.println("solve_cos_t = " + r.length + ", " + t0*180/Math.PI + ", (" + t0 + ")");
        f = 0;
        fprime = 0;
        for (int i = 0; i < r.length; i++)              // test for multiple roots
            f += Math.abs(r[i][0]*Math.cos(r[i][1]*t0 + r[i][2]));
        if (f > TOL/10) do
        {
            f = 0;
            fprime = 0;
            for (int i = 0; i < r.length; i++)
            {
                f += r[i][0]*Math.cos(r[i][1]*t0 + r[i][2]);
                fprime -= r[i][0]*r[i][1]*Math.sin(r[i][1]*t0 + r[i][2]);
            }
            if (Math.abs(fprime) < 10*TOL)
            {
                if (main.IS_DEBUG)
                    System.out.println("too small slope t =, " + t0*180/Math.PI + ", " + f + ", " + fprime + " : Abort");
                return 999999;
            }
            if (loop > 100)
            {
                if (main.IS_DEBUG)
                    System.out.println("too many loops = " + loop + " : Abort");
                return 999999;
            }
            del_t = -f/fprime;
            t0 += del_t;
            loop++;
//            System.out.println("      t =, " + t0*180/Math.PI + ", " + f + ", " + fprime);
        } while (Math.abs(del_t) > TOL/4);

        t0 = (t0 + 2*Math.PI) % (2*Math.PI);        // in case of negative phase shift
        if (Math.abs(t0) < TOL) t0 = 0;
        if ((t0 < 0) || (t0 > 2*Math.PI - TOL))
            return 999999;
        if (main.IS_DEBUG)
            System.out.println("final t =, " + t0*180/Math.PI + ", (" + t0 + "), " + f + ", " + fprime);
        return t0;
    }

    private static double getX(double t)
    {
        return rx1*Math.cos(wx1*t) + rx2*Math.cos(wx2*t);
    }

    private static double getY(double t)
    {
        return ry1*Math.sin(wy1*t) + ry2*Math.sin(wy2*t);
    }

    private static double getdX(double t)
    {
        return -rx1*wx1*Math.sin(wx1*t) - rx2*wx2*Math.sin(wx2*t);
    }

    private static double getdY(double t)
    {
        return ry1*wy1*Math.cos(wy1*t) + ry2*wy2*Math.cos(wy2*t);
    }

    private static double getd2X(double t)
    {
        return -rx1*wx1*wx1*Math.cos(wx1*t) - rx2*wx2*wx2*Math.cos(wx2*t);
    }

    private static double getd2Y(double t)
    {
        return -ry1*wy1*wy1*Math.sin(wy1*t) - ry2*wy2*wy2*Math.sin(wy2*t);
    }
}
