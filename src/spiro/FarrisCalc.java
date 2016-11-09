
package spiro;

import java.awt.geom.Point2D;
import java.awt.geom.CubicCurve2D;

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
        double A12 = r1*r2*w1*w2;
        double A23 = r2*r3*w2*w3;
        double A31 = r3*r1*w3*w1;
        double Sum22 = r1*r1*w1*w1 + r2*r2*w2*w2 + r3*r3*w3*w3;
        double Sum23 = r1*r1*w1*w1*w1 + r2*r2*w2*w2*w2 + r3*r3*w3*w3*w3;
        double[][] M = new double[][] {{3*(w1 + w2) - 2*(w1 + w2), 3*(w1 + w2) - 2*(w2 + w3), 3*(w1 + w2) - 2*(w3 + w1)},
                                       {3*(w2 + w3) - 2*(w1 + w2), 3*(w2 + w3) - 2*(w2 + w3), 3*(w2 + w3) - 2*(w3 + w1)},
                                       {3*(w3 + w1) - 2*(w1 + w2), 3*(w3 + w1) - 2*(w2 + w3), 3*(w3 + w1) - 2*(w3 + w1)}};
        double[] t_perp = new double[1000]; // t values associated with perpendiculat slope
        int N_perp = 0;                     // number of perpendicular solutions

        double[][] rotors = new double[][] {{Sum23, 0, 0},
                                            {A12*(w1 + w2), w1 - w2, phi1 - phi2},
                                            {A23*(w2 + w3), w2 - w3, phi2 - phi3},
                                            {A31*(w3 + w1), w3 - w1, phi3 - phi1}};
//        System.out.println("rotors for inflection points");
//        for (i = 0; i < rotors.length; i++)
//            System.out.println(i + ", " + rotors[i][0] + ", " + rotors[i][1] + ", " + rotors[i][2]);
        for (i = 0; i < rotors.length; i++)                                         // get inflection points
            if ((Math.abs(rotors[i][0]) > TOL) && (Math.abs(rotors[i][1]) > TOL))   // confirm amplitude and frequency
                for (int j = 0; j < Math.round(Math.abs(2*rotors[i][1])); j++)
                    N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(j + 0.5) - rotors[i][2])/Math.abs(rotors[i][1])));
//        System.out.println("raw solutions for inflection points N = " + N);
//        for (i = 0; i < N; i++)
//            System.out.println(i + ", " + t[i]);
        N = SpiroJCalc.sort_t_values(t, N);
        System.out.println("sorted solutions for inflection points N = " + N);
        for (i = 0; i < N; i++)
            System.out.println(i + ", " + t[i]);

        rotors = new double[][] {{A12*(w1 - w2)*(3*Sum23 - (w1 + w2)*Sum22), w1 - w2, phi1 - phi2 - Math.PI/2},
                                 {A23*(w2 - w3)*(3*Sum23 - (w2 + w3)*Sum22), w2 - w3, phi2 - phi3 - Math.PI/2},
                                 {A31*(w3 - w1)*(3*Sum23 - (w3 + w1)*Sum22), w3 - w1, phi3 - phi1 - Math.PI/2},
                                 {A12*A12*M[0][0]*(w1 - w2)/2, 2*(w1 - w2), 2*(phi1 - phi2) - Math.PI/2},
                                 {A23*A23*M[1][1]*(w2 - w3)/2, 2*(w2 - w3), 2*(phi2 - phi3) - Math.PI/2},
                                 {A31*A31*M[2][2]*(w3 - w1)/2, 2*(w3 - w1), 2*(phi3 - phi1) - Math.PI/2},
                                 {A12*A23*(M[1][0]*(w1 - w2) + M[0][1]*(w2 - w3))/2, w1 - w3, phi1 - phi3 - Math.PI/2},
                                 {A23*A31*(M[2][1]*(w2 - w3) + M[1][2]*(w3 - w1))/2, w2 - w1, phi2 - phi1 - Math.PI/2},
                                 {A31*A12*(M[0][2]*(w3 - w1) + M[2][0]*(w1 - w2))/2, w3 - w2, phi3 - phi2 - Math.PI/2},
                                 {A12*A23*(M[1][0]*(w1 - w2) - M[0][1]*(w2 - w3))/2, w1 - 2*w2 + w3, phi1 - 2*phi2 + phi3 - Math.PI/2},
                                 {A23*A31*(M[2][1]*(w2 - w3) - M[1][2]*(w3 - w1))/2, w2 - 2*w3 + w1, phi2 - 2*phi3 + phi1 - Math.PI/2},
                                 {A31*A12*(M[0][2]*(w3 - w1) - M[2][0]*(w1 - w2))/2, w3 - 2*w1 + w2, phi3 - 2*phi1 + phi2 - Math.PI/2}};
//        System.out.println("rotors for curvature extrema");
//        for (i = 0; i < rotors.length; i++)
//            System.out.println(i + ", " + rotors[i][0] + ", " + rotors[i][1] + ", " + rotors[i][2]);
/*
        System.out.println("dump of curvature and rotor fxn");
        double td;                                              // fix fix temporary code
        double f;
        for (i = 0; i < 1201; i++)
        {
            td = 2*Math.PI*i/1200.0/3.0;
            f = 0;
            for (int j = 0; j < rotors.length; j++)
                f += rotors[j][0]*Math.cos(rotors[j][1]*td + rotors[j][2]);
            System.out.println(i + ", " + (float)td + ", " + (float)getX(td) + ", " + (float)getY(td) + ", " + (float)(-(getdX(td)*getd2Y(td) - getdY(td)*getd2X(td))/Math.pow(getdX(td)*getdX(td) + getdY(td)*getdY(td), 1.5)) + ", " + (float)(-f/1E13));
        }
*/
        for (i = 0; i < rotors.length; i++)                                         // get extrema of curvature
            if ((Math.abs(rotors[i][0]) > TOL) && (Math.abs(rotors[i][1]) > TOL))   // check amplitude and frequency
                for (int j = 0; j < Math.round(Math.abs(2*rotors[i][1])); j++)
                    N = main.insert_t_value(N, N, t, SpiroJCalc.solve_cos_t(rotors, (Math.PI*(j - 0.5) - rotors[i][2])/Math.abs(rotors[i][1])));

//        N = main.insert_t_value(N, N, t, Math.PI/4);                    // fix fix temporary code
//        System.out.println("initial t array N = " + N);
        N = SpiroJCalc.sort_t_values(t, N);
        System.out.println("sorted inflection points plus curvature extrema N = " + N);
        for (i = 0; i < N; i++)
            System.out.println(i + ", " + (float)t[i] + ", " + (getdX(t[i])*getd2Y(t[i]) - getdY(t[i])*getd2X(t[i]))/Math.pow(getdX(t[i])*getdX(t[i]) + getdY(t[i])*getdY(t[i]), 1.5)
                                                      + ", " + Math.sqrt(getdX(t[i])*getdX(t[i]) + getdY(t[i])*getdY(t[i])));
        System.out.println("i, t,         x′,         y′,         v,         theta,     K");
        boolean bFirst = true;
        double old_t = 0, old_theta = 0, new_theta, old_K = -1, new_K;
        for (i = 0; i <= N; i++)                                                // check for perpendicular slope
        {
            new_K = (getdX(t[i])*getd2Y(t[i]) - getdY(t[i])*getd2X(t[i]))/Math.pow(getdX(t[i])*getdX(t[i]) + getdY(t[i])*getdY(t[i]), 1.5);
            new_theta = Math.atan2(getdY(t[i]), getdX(t[i]));
            System.out.printf("%d, %f, %f, %f, %f, %g\n", i, t[i], getdX(t[i]), getdY(t[i]),
                               180*new_theta/Math.PI, new_K);
            if (bFirst)
                t[N] = t[i] + 2*Math.PI;                                        // save for last loop
            else if (getdX(old_t)*getdX(t[i]) + getdY(old_t)*getdY(t[i]) < 0)   // greater than 90° arc (dot product)
            {
                if (Math.abs(new_K) > Math.abs(old_K))
                    old_theta = new_theta;
                rotors = new double[][] {{r1*w1, w1, phi1 - old_theta - Math.PI/2},
                                         {r2*w2, w2, phi2 - old_theta - Math.PI/2},
                                         {r3*w3, w3, phi3 - old_theta - Math.PI/2}};
                t_perp[N_perp] = SpiroJCalc.solve_cos_t(rotors, (old_t + t[i])/2);
                System.out.println("test t = ,,,,,,," + old_t + ", " + t_perp[N_perp] + ", " + old_theta*180/Math.PI + ", " + Math.atan2(getdY(t_perp[N_perp]), getdX(t_perp[N_perp]))*180/Math.PI);
                N_perp++;
            }
            old_K = new_K;
            old_theta = new_theta;
            old_t = t[i];
            bFirst = false;
        }
        System.out.println();
        System.out.println("perpendicular list = " + N_perp);                      // fix fix temporary code
        for (i = 0; i < N_perp; i++)
            System.out.println(i + ", " + t_perp[i]);
        System.arraycopy(t_perp, 0, t, N, N_perp);
        N += N_perp;
        System.out.println("unsorted appended t array N = " + N);
        return SpiroJCalc.sort_t_values(t, N);
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
