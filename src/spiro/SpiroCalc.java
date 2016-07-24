
package spiro;

import java.awt.geom.Point2D;
import java.awt.geom.CubicCurve2D;
import java.io.FileWriter;
import java.io.IOException;

//  fit Bezier to spirograph by matching slope and curvature at endpoints
//  Define slope and curvature using “Calculus” by James Stewart, page 902
//  Slope     : m = y′/x′, where all derivatives (′) are with respect to 't'
//  Curvature : κ = (x′y″ - y′x″)/((x′)^2 + (y′)^2)^(3/2)

public final class SpiroCalc
{
    private static final double TOL = 0.00001;
    private static double a, b, c;                                          // spiro parameters

    public static CubicCurve2D.Float getBezier(double t1, double t2)
    {
        Point2D.Double[][] ptSpiro = new Point2D.Double[3][2];              // Point[derivative (0-2)][t = (0,1)]

        ptSpiro[0][0] = new Point2D.Double(getX(t1), getY(t1));
        ptSpiro[0][1] = new Point2D.Double(getX(t2), getY(t2));
        ptSpiro[1][0] = new Point2D.Double(getdX(t1), getdY(t1));
        ptSpiro[1][1] = new Point2D.Double(getdX(t2), getdY(t2));
        ptSpiro[2][0] = new Point2D.Double(getd2X(t1), getd2Y(t1));
        ptSpiro[2][1] = new Point2D.Double(getd2X(t2), getd2Y(t2));

        if (t1 == 0)
        {
            if ((Math.abs(ptSpiro[1][0].x) > TOL) || (Math.abs(ptSpiro[1][0].y) > TOL))
                main.theta[0] = Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x);
            else
            {
                System.out.println("spiro motion is stationary at t = " + t1);
                main.theta[0] = Math.atan2(ptSpiro[2][0].y, ptSpiro[2][0].x);   // use x″ & y″
            }
        }
        else
            main.theta[0] = main.theta[1];
        if ((Math.abs(ptSpiro[1][1].x) > TOL) || (Math.abs(ptSpiro[1][1].y) > TOL))
            main.theta[1] = Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x);
        else                                                                // point is stationary
        {
            System.out.println("spiro motion is stationary at t = " + t2);
            main.theta[1] = Math.atan2(ptSpiro[2][1].y, ptSpiro[2][1].x);   // use x″ & y″
        }
        if ((b < 0) && (-c > -b))                   // hypo with extra loop
            main.theta[1] = main.theta[0] - Math.PI + (main.theta[1] - main.theta[0] + Math.PI) % (2*Math.PI);
        else                                        // hypo with no loop, plus all epi's
            main.theta[1] = main.theta[0] + Math.PI + (main.theta[1] - main.theta[0] - Math.PI) % (-2*Math.PI);
        if (Math.abs(Math.abs(main.theta[0] - main.theta[1]) - Math.PI) < TOL)
            main.theta[1] = main.theta[0];                                  // anti-symmetric case

        for (int i = 0; i < ptSpiro[0].length; i++)                         // standard curvature
            if ((Math.abs(ptSpiro[1][i].x) > TOL) || (Math.abs(ptSpiro[1][i].y) > TOL))
                main.Cu[i] = (ptSpiro[1][i].x*ptSpiro[2][i].y - ptSpiro[1][i].y*ptSpiro[2][i].x)
                           /  Math.pow((ptSpiro[1][i].x*ptSpiro[1][i].x + ptSpiro[1][i].y*ptSpiro[1][i].y), 1.5);
            else
                main.Cu[i] = 0;                                             // stationary point
        return main.calcBezier(ptSpiro, t1, t2);
    }

    public static void write_test_quadratic(FileWriter out)
    {
        // this will write a test quadratic bezier, just for testing purposes
        // this uses the standard parms from the last object rendered
        // it is assumed that : OriginX/Y = 0, InitialAngle = 0, Zoom = 1
        // fit only position and slope at t1, t2

        String strPath;
        double x0, y0, x2, y2;
        double dx0, dy0, dx2, dy2;
        double x1 = 0, y1 = 0;
        double org = main.PAGE_WIDTH*main.mm2px/2;
        if (main.t1 < 0 || main.t2 < 0) return;

        System.out.printf("quadratic test parms = %f, %f, %f, %f, %f\n", a, b, c, main.t1, main.t2);
        try
        {
            strPath = "    <path\n"
                    + "        d='";
            out.write(strPath);
            x0 = org + getX(main.t1);
            y0 = org + getY(main.t1);
            x2 = org + getX(main.t2);
            y2 = org + getY(main.t2);
            dx0 = getdX(main.t1);
            dy0 = getdY(main.t1);
            dx2 = getdX(main.t2);
            dy2 = getdY(main.t2);

            if (Math.abs(dy2*dx0 - dy0*dx2) < TOL)
                System.out.println("Error : parallel lines in write_test_quadratic, t = " + main.t1 + ", " + main.t2);
            else
            {
                x1 = (x2*dy2*dx0 - x0*dy0*dx2 - (y2 - y0)*dx0*dx2)/(dy2*dx0 - dy0*dx2);
                y1 = (y0*dy2*dx0 - y2*dy0*dx2 - (x0 - x2)*dy0*dy2)/(dy2*dx0 - dy0*dx2);
            }
            strPath = "M " + (float)x0 + ", " + (float)y0 + " Q " + (float)x1 + ", " + (float)y1 + " " + (float)x2 + ", " + (float)y2;
            out.write(strPath);
            strPath = "'\n"
                    + "        id='test_quadratic'\n"
                    + "        style='fill:none;stroke:#ff0000;stroke-width:1px' />\n";
            out.write(strPath);
        }
        catch (IOException e)
            {System.out.println("save test quadratic Bezier error = " + e);}
    }

    private static double calc_cos_t(double fAmp, double fFreq)
    {
        // calculate the t value at the point of maximum width
        // of a lobe of an epiTrochoid or hypoTrochoid
        // solve : cos(t) = fAmp*cos(fFreq*t)

        double t0 = 0, del_t;

        if (fAmp > 1)
            t0 = Math.sqrt(2*(fAmp - 1)/(fAmp*fFreq*fFreq - 1));
        else if (fAmp < 0)
            t0 = Math.PI/2/fFreq;
        else
            return 0;
        System.out.println("\ncalc_cos_t = " + fAmp + ", " + fFreq + ", " + t0);
        do
        {
            del_t = (fAmp*Math.cos(fFreq*t0) - Math.cos(t0))/(fFreq*fAmp*Math.sin(fFreq*t0) - Math.sin(t0));
            t0 += del_t;
            System.out.println("     t = " + t0);
        } while (Math.abs(del_t) > 0.0000001);
        return t0;
    }

    private static double calc_sin_t(double fAmp, double fFreq)
    {
        // calculate the t value at a point perpendicular to the maximum width
        // of a lobe of an epiTrochoid
        // solve : sin(t) = fAmp*sin(fFreq*t)

        double t0 = 0, del_t;

        if (fAmp*fFreq > 1)
            t0 = Math.PI/fFreq;
        else
            return 0;
        System.out.println("\ncalc_sin_t = " + fAmp + ", " + fFreq + ", " + t0);
        do
        {
            del_t = (fAmp*Math.sin(fFreq*t0) - Math.sin(t0))/(fFreq*fAmp*Math.cos(fFreq*t0) - Math.cos(t0));
            t0 -= del_t;
            System.out.println("     t = " + t0);
        } while (Math.abs(del_t) > 0.0000001);
        return t0;
    }

    public static int get_t_values(double[] t, double m_a, double m_b, double m_c)
    {
        int N = 0;                                                  // points per lobe
        a = m_a; b = m_b; c = m_c;

        if (b < 0)                                                  // hypo
        {
            N = main.insert_t_value(N, 0, t, 0);                         // at max R
            if (-c > -b)                                            // extra loop
            {
                N = main.insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));   // at max width of small lobe
                System.out.println("hypo max width : " + t[1]);
            }
            else if (c == b)                                        // cycloid
            {
                System.out.println("cycloid case : " + t[1]);
            }
            else if (-c > b/(1 + a/b))                              // no loop
            {
                N = main.insert_t_value(N, 1, t, -b/a*Math.acos((1 + c*c*(1 + a/b)/b/b)*b/c/(2 + a/b)));  // inflection point
                System.out.println("hypo inflection : " + t[1]);
            }
            else                                                    // no inflection
            {
                System.out.println("hypo no inflection");
            }
        }
        else if (b > 0)                                             // epi
        {
            N = main.insert_t_value(N, 0, t, 0);                         // at max R
            if (-c > b)                                             // extra loop
            {
                N = main.insert_t_value(N, 1, t, Math.PI*b/a - calc_sin_t(-c/b, 1 + a/b));   // perpendicular to max width of small lobe
                System.out.println("epi perpendicular to small width : " + t[1]);
                N = main.insert_t_value(N, 2, t, Math.PI*b/a - calc_cos_t(-c/b, 1 + a/b));   // at max width of small lobe
                System.out.println("epi small lobe width : " + t[2]);
                N = main.insert_t_value(N, 2, t, (t[1] + t[2])/2);       // interpolate
                if (a < 2*b)
                    N = main.insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
            }
            else if (-c == b)                                       // cycloid
            {
                N = main.insert_t_value(N, 1, t, Math.PI*b/a - calc_sin_t(-c/b, 1 + a/b));   // perpendicular to max width of small lobe
                System.out.println("cycloid case : " + t[1]);
                N = main.insert_t_value(N, 2, t, (t[1] + Math.PI*b/a)/2);// interpolate
                if (a < 2*b)
                    N = main.insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
            }
            else if (-c > b/(1 + a/b))                              // no loop
            {
                N = main.insert_t_value(N, 1, t, Math.PI*b/a - calc_sin_t(-c/b, 1 + a/b));   // perpendicular to max width of small lobe
                System.out.println("epi perpendicular to small width : " + t[1]);
                N = main.insert_t_value(N, 2, t, b/a*Math.acos((1 + c*c*(1 + a/b)/b/b)*b/c/(2 + a/b)));  // inflection point
                System.out.println("epi inflection : " + t[2]);
                N = main.insert_t_value(N, 2, t, (t[1] + t[2])/2);       // interpolate
                if (a < 2*b)
                {
                    N = main.insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
                    N = main.insert_t_value(N, 2, t, (t[1] + t[2])/2);   // interpolate
                }
                else
                    N = main.insert_t_value(N, 1, t, (t[0] + t[1])/2);   // interpolate
            }
            else                                                    // no inflection
            {
                if (a < 2*b)
                {
                    N = main.insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
                    N = main.insert_t_value(N, 2, t, (t[1] + Math.PI*b/a)/2);    // interpolate
                }
                else
                    N = main.insert_t_value(N, 1, t, Math.PI*b/a/2);     // interpolate
                System.out.println("epi no inflection");
            }
        }
        if (N > 0) t[N] = Math.PI*Math.abs(b)/a;                    // at min R
        if (N > 1)
            for (int i = 1; i < N; i++)
                t[N + i] = 2*Math.PI*Math.abs(b)/a - t[N - i];      // reflect first half
        return 2*N;
    }
/*
    private static int insert_t_value(int N, int index, double[] t, double new_t)
    {
        // push t[index] up by one, insert at location index
        if (N > index)
            for (int i = N; i > index; i--)
                t[i] = t[i - 1];
        t[index] = new_t;
        return N + 1;
    }
*/
    private static double getX(double t)
    {
        return (a + b)*Math.cos(t) - c*Math.cos(t*(a/b + 1));
    }

    private static double getY(double t)
    {
        return (a + b)*Math.sin(t) - c*Math.sin(t*(a/b + 1));
    }

    private static double getdX(double t)
    {
        return -(a + b)*Math.sin(t) + c*(a/b + 1)*Math.sin(t*(a/b + 1));
    }

    private static double getdY(double t)
    {
        return (a + b)*Math.cos(t) - c*(a/b + 1)*Math.cos(t*(a/b + 1));
    }

    private static double getd2X(double t)
    {
        return -(a + b)*Math.cos(t) + c*(a/b + 1)*(a/b + 1)*Math.cos(t*(a/b + 1));
    }

    private static double getd2Y(double t)
    {
        return -(a + b)*Math.sin(t) + c*(a/b + 1)*(a/b + 1)*Math.sin(t*(a/b + 1));
    }
}
