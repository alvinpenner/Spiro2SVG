
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

    public static CubicCurve2D.Float getBezier(double t1, double t2)
    {
        Point2D.Double[][] ptSpiro = new Point2D.Double[3][2];          // Point[derivative (0-2)][t = (0,1)]

        ptSpiro[0][0] = new Point2D.Double(getX(t1), getY(t1));
        ptSpiro[0][1] = new Point2D.Double(getX(t2), getY(t2));
        ptSpiro[1][0] = new Point2D.Double(getdX(t1), getdY(t1));
        ptSpiro[1][1] = new Point2D.Double(getdX(t2), getdY(t2));
        ptSpiro[2][0] = new Point2D.Double(getd2X(t1), getd2Y(t1));
        ptSpiro[2][1] = new Point2D.Double(getd2X(t2), getd2Y(t2));

        System.out.println("");
        if (t1 == 0)
        {
            if ((Math.abs(ptSpiro[1][0].x) > TOL) || (Math.abs(ptSpiro[1][0].y) > TOL))
                main.theta[0] = Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x);
            else
            {
                System.out.println("spiro motion is stationary at t = " + t1);
                main.theta[0] = Math.atan2(ptSpiro[2][0].y, ptSpiro[2][0].x);       // use x″ & y″
            }
        }
        else
        {
            main.theta[0] = main.theta[1];
            if ((Math.abs(ptSpiro[1][0].x) < TOL) && (Math.abs(ptSpiro[1][0].y) < TOL))
            {
                System.out.println("spiro motion is stationary at t = " + t1);
                main.theta[0] -= isCCW ? Math.PI : -Math.PI;    // reverse direction
                isCCW = !isCCW;
            }
        }
        if ((Math.abs(ptSpiro[1][1].x) > TOL) || (Math.abs(ptSpiro[1][1].y) > TOL))
            main.theta[1] = Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x);
        else                                                    // point is stationary
        {
            System.out.println("spiro motion is stationary at t = " + t2);
            main.theta[1] = Math.atan2(ptSpiro[2][1].y, ptSpiro[2][1].x); // use x″ & y″
            main.theta[1] += isCCW ? Math.PI : -Math.PI;
        }

        if (isCCW)
            main.theta[1] = main.theta[0] + Math.PI + (main.theta[1] - main.theta[0] - Math.PI) % (-2*Math.PI);
        else
            main.theta[1] = main.theta[0] - Math.PI + (main.theta[1] - main.theta[0] + Math.PI) % (2*Math.PI);
        if (Math.abs(Math.abs(main.theta[0] - main.theta[1]) - Math.PI) < TOL)
            main.theta[1] = main.theta[0];                                          // anti-symmetric case

        for (int i = 0; i < ptSpiro[0].length; i++)                                 // standard curvature
            if ((Math.abs(ptSpiro[1][i].x) > TOL) || (Math.abs(ptSpiro[1][i].y) > TOL))
                main.Cu[i] = (ptSpiro[1][i].x*ptSpiro[2][i].y - ptSpiro[1][i].y*ptSpiro[2][i].x)
                           /  Math.pow((ptSpiro[1][i].x*ptSpiro[1][i].x + ptSpiro[1][i].y*ptSpiro[1][i].y), 1.5);
            else
                main.Cu[i] = rx1*ry1*wx1*wx1*wy1*wy1*(wx1*wx1 - wy1*wy1)/3          // stationary point
                           / Math.pow(rx1*rx1*wx1*wx1*wx1*wx1 + ry1*ry1*wy1*wy1*wy1*wy1, 1.5);
        return main.calcBezier(ptSpiro, t1, t2);
    }

    public static int get_t_values(double[] t, double m_rx1, double m_ry1, double m_wx1, double m_wy1, double m_rx2, double m_ry2, double m_wx2, double m_wy2)
    {
        int i, N = 0;                                               // points per object
        rx1 = m_rx1; ry1 = m_ry1; wx1 = m_wx1; wy1 = m_wy1;
        rx2 = m_rx2; ry2 = m_ry2; wx2 = m_wx2; wy2 = m_wy2;

        for (i = 0; i < Math.round(2*wx1); i++)                     // x′ = 0
            N = insert_t_value(N, N, t, Math.PI*i/wx1);
        for (i = 0; i < Math.round(2*wy1); i++)                     // y′ = 0
            N = insert_t_value(N, N, t, Math.PI*(i + 0.5)/wy1);
        if (Math.abs(wx1) != Math.abs(wy1))                         // reject ellipse
        {
            for (i = 0; i < Math.round(2*Math.abs(wx1 + wy1)); i++) // inflection using wx1+wy1
                N = insert_t_value(N, N, t, solve_cos_t((wx1 + wy1)/2, wx1 - wy1, (wx1 - wy1)/2, wx1 + wy1, Math.PI*(i + 0.5)/Math.abs(wx1 + wy1)));
            for (i = 0; i < Math.round(2*Math.abs(wx1 - wy1)); i++) // inflection using wx1-wy1
                N = insert_t_value(N, N, t, solve_cos_t((wx1 + wy1)/2, wx1 - wy1, (wx1 - wy1)/2, wx1 + wy1, Math.PI*(i + 0.5)/Math.abs(wx1 - wy1)));
        }
        Arrays.sort(t, 0, N);
        if (N > 1)                                                  // check for duplicates
            for (i = 0; i < N - 1; i++)
                if (t[i + 1] < t[i] + TOL)
                    t[i] = 999999;
        Arrays.sort(t, 0, N);
        for (i = 0; i < N; i++)                                     // remove duplicates
            if (t[i] == 999999)
                break;
        return i;
    }

    private static int insert_t_value(int N, int index, double[] t, double new_t)
    {
        // push t[index] up by one, insert new_t at location index
        if (N > index)
            for (int i = N; i > index; i--)
                t[i] = t[i - 1];
        t[index] = new_t;
        return N + 1;
    }

    private static double solve_cos_t(double a1, double w1, double a2, double w2, double t0)
    {
        // calculate the t value at inflection points
        // solve : 0 = a1*cos(w1*t) + a2*cos(w2*t)

        double f, fprime, del_t;

        System.out.println("\nsolve_cos_t = " + a1 + ", " + w1 + ", " + a2 + ", " + w2 + ", " + t0*180/Math.PI);
        fprime = -a1*w1*Math.sin(w1*t0) - a2*w2*Math.sin(w2*t0);
        if (Math.abs(fprime) < TOL)
        {
            System.out.println("initial slope = " + fprime + " : Abort");
            return 999999;
        }
        do
        {
            f = a1*Math.cos(w1*t0) + a2*Math.cos(w2*t0);
            fprime = -a1*w1*Math.sin(w1*t0) - a2*w2*Math.sin(w2*t0);
            del_t = -f/fprime;
            t0 += del_t;
            System.out.println("     t =, " + t0*180/Math.PI + ", " + fprime);
        } while ((Math.abs(del_t) > TOL/2) && (Math.abs(fprime) > 10*TOL));
        if ((Math.abs(fprime) <= 10*TOL)
        ||  (t0 < 0) || (t0 >= 2*Math.PI))
            return 999999;
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
