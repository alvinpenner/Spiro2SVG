
package spiro;

import java.awt.geom.Point2D;
import java.awt.geom.CubicCurve2D;
import java.util.Arrays;

//  fit Bezier to Farris shape (3 wheels) by matching slope and curvature at endpoints
//  Define slope and curvature using “Calculus” by James Stewart, page 902
//  Slope     : m = y′/x′, where all derivatives (′) are with respect to 't'
//  Curvature : κ = (x′y″ - y′x″)/((x′)^2 + (y′)^2)^(3/2)

public final class FarrisCalc
{
    private static final double TOL = 0.0000001;
    private static double r1, w1, phi1, r2, w2, phi2, r3, w3, phi3;         // Farris Wheel parameters
    private static boolean isCCW;                                           // change sign at stationary point
    private static boolean firstCall = true;

    protected static CubicCurve2D.Float getBezier(double t1, double t2)
    {
        Point2D.Double[][] ptSpiro = new Point2D.Double[3][2];              // Point[derivative (0-2)][t = (0,1)]

        ptSpiro[0][0] = new Point2D.Double(getX(t1), getY(t1));
        ptSpiro[0][1] = new Point2D.Double(getX(t2), getY(t2));
        ptSpiro[1][0] = new Point2D.Double(getdX(t1), getdY(t1));
        ptSpiro[1][1] = new Point2D.Double(getdX(t2), getdY(t2));
        ptSpiro[2][0] = new Point2D.Double(getd2X(t1), getd2Y(t1));
        ptSpiro[2][1] = new Point2D.Double(getd2X(t2), getd2Y(t2));

        System.out.println();
        if (firstCall)
        {
            firstCall = false;
            isCCW = w3 > 0;                     // assume w3 is the highest frequency
            if ((Math.abs(ptSpiro[1][0].x) > TOL) || (Math.abs(ptSpiro[1][0].y) > TOL))
                main.theta[0] = Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x);
            else
            {
                System.out.println("Farris motion is stationary at t = " + t1);
                main.theta[0] = Math.atan2(ptSpiro[2][0].y, ptSpiro[2][0].x);   // use x″ & y″
            }
        }
        else
        {
            main.theta[0] = main.theta[1];
            if ((Math.abs(ptSpiro[1][0].x) < TOL) && (Math.abs(ptSpiro[1][0].y) < TOL))
            {
                System.out.println("Farris motion is stationary at t = " + t1);
                main.theta[0] -= isCCW ? Math.PI : -Math.PI;                // reverse direction
                isCCW = !isCCW;
            }
        }
        if ((Math.abs(ptSpiro[1][1].x) > TOL) || (Math.abs(ptSpiro[1][1].y) > TOL))
            main.theta[1] = Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x);
        else                                                                // point is stationary
        {
            System.out.println("Farris motion is stationary at t = " + t2);
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
            {
                main.Cu[i] = 0;                                             // stationary point
                System.out.println("unexpected stationary point in Farris Wheel, needs fix!");
            }
        return main.calcBezier(ptSpiro, t1, t2);
    }

    protected static int get_t_values(double[] t, double m_r1, double m_w1, double m_phi1, double m_r2, double m_w2, double m_phi2, double m_r3, double m_w3, double m_phi3)
    {
        int i, N = 0;                                                       // points per object
        r1 = m_r1; w1 = m_w1; phi1 = m_phi1;
        r2 = m_r2; w2 = m_w2; phi2 = m_phi2;
        r3 = m_r3; w3 = m_w3; phi3 = m_phi3;
        System.out.println("Farris parms = " + r1 + ", " + w1 + ", " + phi1 + ", " + r2 + ", " + w2 + ", " + phi2);

        double[][] rotors = new double[][] {{r1*w1, w1, phi1 - Math.PI/2},
                                            {r2*w2, w2, phi2 - Math.PI/2},
                                            {r3*w3, w3, phi3 - Math.PI/2}};
        if ((r1 != 0) && (w1 != 0))
            for (i = 0; i < Math.round(Math.abs(2*w1)); i++)                // x′ = 0
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*i - phi1)/Math.abs(w1)));
        if ((r2 != 0) && (w2 != 0))
            for (i = 0; i < Math.round(Math.abs(2*w2)); i++)                // x′ = 0
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*i - phi2)/Math.abs(w2)));
        if ((r3 != 0) && (w3 != 0))
            for (i = 0; i < Math.round(Math.abs(2*w3)); i++)                // x′ = 0
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*i - phi3)/Math.abs(w3)));

        rotors = new double[][] {{r1*w1, w1, phi1},
                                 {r2*w2, w2, phi2},
                                 {r3*w3, w3, phi3}};
        if ((r1 != 0) && (w1 != 0))
            for (i = 0; i < Math.round(Math.abs(2*w1)); i++)                // y′ = 0
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(i + 0.5) - phi1)/Math.abs(w1)));
        if ((r2 != 0) && (w2 != 0))
            for (i = 0; i < Math.round(Math.abs(2*w2)); i++)                // y′ = 0
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(i + 0.5) - phi2)/Math.abs(w2)));
        if ((r3 != 0) && (w3 != 0))
            for (i = 0; i < Math.round(Math.abs(2*w3)); i++)                // y′ = 0
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(i + 0.5) - phi3)/Math.abs(w3)));

        rotors = new double[][] {{r1*r1*w1*w1*w1, 0, 0},
                                 {r2*r2*w2*w2*w2, 0, 0},
                                 {r3*r3*w3*w3*w3, 0, 0},
                                 {r1*r2*w1*w2*(w1 + w2), w1 - w2, phi1 - phi2},
                                 {r2*r3*w2*w3*(w2 + w3), w2 - w3, phi2 - phi3},
                                 {r3*r1*w3*w1*(w3 + w1), w3 - w1, phi3 - phi1}};
        if ((r1 != 0) && (r2 != 0) && (Math.abs(w1 - w2) > TOL))
            for (i = 0; i < Math.round(2*Math.abs(w1 - w2)); i++)           // inflection using w1-w2
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(i + 0.5) - phi1 + phi2)/Math.abs(w1 - w2)));
        if ((r2 != 0) && (r3 != 0) && (Math.abs(w2 - w3) > TOL))
            for (i = 0; i < Math.round(2*Math.abs(w2 - w3)); i++)           // inflection using w2-w3
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(i + 0.5) - phi2 + phi3)/Math.abs(w2 - w3)));
        if ((r3 != 0) && (r1 != 0) && (Math.abs(w3 - w1) > TOL))
            for (i = 0; i < Math.round(2*Math.abs(w3 - w1)); i++)           // inflection using w3-w1
                N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(i + 0.5) - phi3 + phi1)/Math.abs(w3 - w1)));

        N = main.insert_t_value(N, N, t, Math.PI/4);                    // fix fix temporary code
        System.out.println("initial t array N = " + N);
        Arrays.sort(t, 0, N);
        if (N > 1)                                                      // check for duplicates
            for (i = 0; i < N - 1; i++)
                if (t[i + 1] < t[i] + TOL)
                    t[i] = 999999;
        Arrays.sort(t, 0, N);
        for (i = 0; i < N; i++)                                         // remove duplicates
            if (t[i] == 999999)
                break;
        return i;
    }

    private static double getX(double t)
    {
        return r1*Math.cos(w1*t + phi1) + r2*Math.cos(w2*t + phi2) + r3*Math.cos(w3*t + phi3);
    }

    private static double getY(double t)
    {
        return r1*Math.sin(w1*t + phi1) + r2*Math.sin(w2*t + phi2) + r3*Math.sin(w3*t + phi3);
    }

    private static double getdX(double t)
    {
        return -r1*w1*Math.sin(w1*t + phi1) - r2*w2*Math.sin(w2*t + phi2) - r3*w3*Math.sin(w3*t + phi3);
    }

    private static double getdY(double t)
    {
        return r1*w1*Math.cos(w1*t + phi1) + r2*w2*Math.cos(w2*t + phi2) + r3*w3*Math.cos(w3*t + phi3);
    }

    private static double getd2X(double t)
    {
        return -r1*w1*w1*Math.cos(w1*t + phi1) - r2*w2*w2*Math.cos(w2*t + phi2) - r3*w3*w3*Math.cos(w3*t + phi3);
    }

    private static double getd2Y(double t)
    {
        return -r1*w1*w1*Math.sin(w1*t + phi1) - r2*w2*w2*Math.sin(w2*t + phi2) - r3*w3*w3*Math.sin(w3*t + phi3);
    }
}
