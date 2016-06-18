
package spiro;

import java.awt.geom.Point2D;
import java.awt.geom.CubicCurve2D;
import java.io.FileWriter;
import java.io.IOException;

//  fit Bezier to Spiro by matching slope and curvature at endpoints
//  Define slope and curvature using “Calculus” by James Stewart, page 902
//  Slope     : m = y′/x′, where all derivatives (′) are with respect to 't'
//  Curvature : κ = (x′y″ - y′x″)/((x′)^2 + (y′)^2)^(3/2)

public final class SpiroCalc
{
    private static final double TOL = 0.00001;
    private static Point2D.Double[][] ptSpiro = new Point2D.Double[3][2];   // Point[derivative (0-2)][t = (0,1)]
    private static Point2D.Double[] ptBez = new Point2D.Double[4];          // Point[point index (0-3)]
    private static double[] C = new double[2];                              // Spiro curvature[t = (0,1)]
    private static double[] theta = new double[2];                          // Spiro velocity angle[t = (0,1)]

    public static CubicCurve2D.Float getBezier(double a, double b, double c, double t1, double t2)
    {
        double delxrot0, delyrot0, delxrot3, delyrot3;

        ptSpiro[0][0] = new Point2D.Double(getX(a, b, c, t1), getY(a, b, c, t1));
        ptSpiro[0][1] = new Point2D.Double(getX(a, b, c, t2), getY(a, b, c, t2));
        ptSpiro[1][0] = new Point2D.Double(getdX(a, b, c, t1), getdY(a, b, c, t1));
        ptSpiro[1][1] = new Point2D.Double(getdX(a, b, c, t2), getdY(a, b, c, t2));
        ptSpiro[2][0] = new Point2D.Double(getd2X(a, b, c, t1), getd2Y(a, b, c, t1));
        ptSpiro[2][1] = new Point2D.Double(getd2X(a, b, c, t2), getd2Y(a, b, c, t2));
        ptBez[0] = new Point2D.Double(ptSpiro[0][0].x, ptSpiro[0][0].y);
        ptBez[1] = new Point2D.Double(ptBez[0].x, ptBez[0].y);
        ptBez[3] = new Point2D.Double(ptSpiro[0][1].x, ptSpiro[0][1].y);
        ptBez[2] = new Point2D.Double(ptBez[3].x, ptBez[3].y);

        if (t1 == 0)
        {
            if ((Math.abs(ptSpiro[1][0].x) > TOL) || (Math.abs(ptSpiro[1][0].y) > TOL))
                theta[0] = Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x);
            else
            {
                System.out.println("spiro motion is stationary at t = " + t1);
                theta[0] = Math.atan2(ptSpiro[2][0].y, ptSpiro[2][0].x);    // use x″ & y″
            }
        }
        else
            theta[0] = theta[1];
        if ((Math.abs(ptSpiro[1][1].x) > TOL) || (Math.abs(ptSpiro[1][1].y) > TOL))
            theta[1] = Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x);
        else                                            // point is stationary
        {
            System.out.println("spiro motion is stationary at t = " + t2);
            theta[1] = Math.atan2(ptSpiro[2][1].y, ptSpiro[2][1].x);    // use x″ & y″
        }
        if ((b < 0) && (-c > -b))                   // hypo with extra loop
            theta[1] = theta[0] - Math.PI + (theta[1] - theta[0] + Math.PI) % (2*Math.PI);
        else                                        // hypo with no loop, plus all epi's
            theta[1] = theta[0] + Math.PI + (theta[1] - theta[0] - Math.PI) % (-2*Math.PI);
        if (Math.abs(Math.abs(theta[0] - theta[1]) - Math.PI) < TOL)
            theta[1] = theta[0];                                        // anti-symmetric case
        double xrot0 = getrotX(ptBez[0].x, ptBez[0].y, (theta[0] + theta[1])/2);
        double yrot0 = getrotY(ptBez[0].x, ptBez[0].y, (theta[0] + theta[1])/2);
        double xrot3 = getrotX(ptBez[3].x, ptBez[3].y, (theta[0] + theta[1])/2);
        double yrot3 = getrotY(ptBez[3].x, ptBez[3].y, (theta[0] + theta[1])/2);
        double mrot0 = Math.tan((theta[0] - theta[1])/2);
        for (int i = 0; i < ptSpiro[0].length; i++)                     // effective 'rotated' curvature
            if ((Math.abs(ptSpiro[1][i].x) > TOL) || (Math.abs(ptSpiro[1][i].y) > TOL))
                C[i] = Math.signum(getrotX(ptSpiro[1][i].x, ptSpiro[1][i].y, (theta[0] + theta[1])/2))*
                       (ptSpiro[1][i].x*ptSpiro[2][i].y - ptSpiro[1][i].y*ptSpiro[2][i].x)*
                       Math.pow((1 + mrot0*mrot0)/
                       (ptSpiro[1][i].x*ptSpiro[1][i].x + ptSpiro[1][i].y*ptSpiro[1][i].y), 1.5);
            else
                C[i] = 0;                                       // stationary point
        System.out.println("angles = " + Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x)*180/Math.PI + ", " + Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x)*180/Math.PI + ", " +  mrot0);
        System.out.println("angles = " + (float)t1 + ", " + (float)ptSpiro[0][0].x + ", " + (float)ptSpiro[0][0].y + ", " + (float)(theta[0]*180/Math.PI) + ", " + (float)C[0]);
        System.out.println("angles = " + (float)t2 + ", " + (float)ptSpiro[0][1].x + ", " + (float)ptSpiro[0][1].y + ", " + (float)(theta[1]*180/Math.PI) + ", " + (float)C[1]);

        if (Math.abs(theta[1] - theta[0]) < TOL)                // parallel finite slopes
        {
            // the parallel case can be solved as two decoupled quadratic equations
            System.out.println("parallel finite : " + (float)(t1/2/Math.PI) + ", " + (float)(t2/2/Math.PI) + ", " + (float)mrot0 + ", " + (float)C[0] + ", " + (float)C[1] + ", " + (float)yrot0 + ", " + (float)yrot3);
            delxrot0 =  2*(yrot3 - yrot0)/C[0]/3;
            if (delxrot0 < 0)
            {
                System.out.println("parallel slope, wrong curvature at t = 0 : " + delxrot0 + " : abort");
                delxrot0 = Double.NaN;
                return new CubicCurve2D.Float();
            }
            else
                delxrot0 = Math.signum(getrotX(ptSpiro[1][0].x, ptSpiro[1][0].y, (theta[0] + theta[1])/2)*(t2 - t1))*Math.sqrt(delxrot0);
            delxrot3 = -2*(yrot3 - yrot0)/C[1]/3;
            if (delxrot3 < 0)
            {
                System.out.println("parallel slope, wrong curvature at t = 1 : " + delxrot3 + " : abort");
                delxrot3 = Double.NaN;
                return new CubicCurve2D.Float();
            }
            else
                delxrot3 = Math.signum(getrotX(ptSpiro[1][1].x, ptSpiro[1][1].y, (theta[0] + theta[1])/2)*(t2 - t1))*Math.sqrt(delxrot3);
        }
        else if (Math.abs(C[0]) < TOL)                          // zero curvature at t = 0
        {
            System.out.println("zero curvature at t = 0 : " + Math.signum(ptSpiro[1][0].y) + ", " + Math.signum(t2 - t1) + ", " + C[0]);
            delxrot3 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))/2/mrot0;
            if (-c == Math.abs(b))                              // for cycloid, clamp this so it cannot change sign
                delxrot0 = 0;
            else
                delxrot0 =  (yrot3 - yrot0 + mrot0*(xrot3 - xrot0) + 3*C[1]*delxrot3*delxrot3/2)/2/mrot0;
            System.out.println("delxrot0/3 = " + delxrot0 + ", " + delxrot3 + ", " + c + ", " + b);
        }
        else if (Math.abs(C[1]) < TOL)                          // zero curvature at t = 1
        {
            System.out.println("zero curvature at t = 1 : " + Math.signum(ptSpiro[1][0].y) + ", " + Math.signum(t2 - t1) + ", " + C[1]);
            delxrot0 = (yrot3 - yrot0 + mrot0*(xrot3 - xrot0))/2/mrot0;
            if (-c == Math.abs(b))                              // for cycloid, clamp this so it cannot change sign
                delxrot3 = 0;
            else
                delxrot3 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*C[0]*delxrot0*delxrot0/2)/2/mrot0;
            System.out.println("delxrot0/3 = " + delxrot0 + ", " + delxrot3 + ", " + c + ", " + b);
        }
        else                                                    // general slopes
        {
            System.out.println("general slope = " + t1 + ", " + t2 + ", " + theta[0]*180/Math.PI + ", " + theta[1]*180/Math.PI + ", " + (theta[0] + theta[1])*90/Math.PI + ", " + mrot0 + ", " + C[0] + ", " + C[1]);
            System.out.println("general slope = " + t1 + ", " + t2 + ", " + xrot0 + ", " + yrot0 + ", " + xrot3 + ", " + yrot3);
            delxrot0 = solve_quartic(-27*C[1]*C[0]*C[0]/8,
                                      0,
                                      9*C[1]*C[0]*(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))/2,
                                      8*mrot0*mrot0*mrot0,
                                     -(yrot3 - yrot0 + mrot0*(xrot3 - xrot0))*4*mrot0*mrot0
                                     -3*C[1]*(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))
                                            *(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))/2,
                                      Math.signum(getrotX(ptSpiro[1][0].x, ptSpiro[1][0].y, (theta[0] + theta[1])/2)*(t2 - t1)));
            delxrot3 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*C[0]*delxrot0*delxrot0/2)/2/mrot0;
        }
        delyrot0 = mrot0*delxrot0;
        delyrot3 = -mrot0*delxrot3;
        ptBez[1].x += getrotX(delxrot0, delyrot0, -(theta[0] + theta[1])/2);
        ptBez[1].y += getrotY(delxrot0, delyrot0, -(theta[0] + theta[1])/2);
        ptBez[2].x -= getrotX(delxrot3, delyrot3, -(theta[0] + theta[1])/2);
        ptBez[2].y -= getrotY(delxrot3, delyrot3, -(theta[0] + theta[1])/2);
        return new CubicCurve2D.Float((float)ptBez[0].x, (float)ptBez[0].y, (float)ptBez[1].x, (float)ptBez[1].y, (float)ptBez[2].x, (float)ptBez[2].y, (float)ptBez[3].x, (float)ptBez[3].y);
/*
        System.out.println("\nptSpiro[3][2]");
        for (int i = 0; i < ptSpiro.length; i++)
            for (int j = 0; j < ptSpiro[0].length; j++)
                System.out.println(i + " " + j + " : " + ptSpiro[i][j]);
        System.out.println("\nptBez[4]");
        for (int i = 0; i < ptBez.length; i++)
            System.out.println(i + " : " + ptBez[i]);
        System.out.println("\nm[2]");
        for (int i = 0; i < m.length; i++)
            System.out.println(i + " : " + m[i]);
        System.out.println("\nC[2]");
        for (int i = 0; i < C.length; i++)
            System.out.println(i + " : " + C[i]);
        System.out.println("\ndel[2]");
        for (int i = 0; i < del.length; i++)
            System.out.println(i + " : " + del[i]);
*/
    }

    public static void write_test_cubic(FileWriter out, double a, double b, double c)
    {
        // this will write a test cubic bezier, just for testing purposes

        getBezier(a, b, c, main.t1, main.t2);                   // refresh ptBez[]
        String strPath;
        for (int i = 0; i < ptBez.length; i++)                  // re-define origin
        {
            ptBez[i].x += main.PAGE_WIDTH*main.mm2px/2;
            ptBez[i].y += main.PAGE_WIDTH*main.mm2px/2;
        }
        try
        {
            strPath = "    <path\n"
                    + "        d='";
            out.write(strPath);
            strPath =  "M " + (float)ptBez[0].x + ", " + (float)ptBez[0].y +
                      " C " + (float)ptBez[1].x + ", " + (float)ptBez[1].y + " " + (float)ptBez[2].x + ", " + (float)ptBez[2].y + " " + (float)ptBez[3].x + ", " + (float)ptBez[3].y;
            out.write(strPath);
            strPath = "'\n"
                    + "        id='test_cubic'\n"
                    + "        style='fill:none;stroke:#0000ff;stroke-width:1px' />\n";
            out.write(strPath);
        }
        catch (IOException e)
            {System.out.println("save test cubic Bezier error = " + e);}
    }

    public static void write_test_quadratic(FileWriter out, double a, double b, double c)
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

        System.out.printf("test parms = %f, %f, %f, %f, %f\n", a, b, c, main.t1, main.t2);
        try
        {
            strPath = "    <path\n"
                    + "        d='";
            out.write(strPath);
            x0 = org + getX(a, b, c, main.t1);
            y0 = org + getY(a, b, c, main.t1);
            x2 = org + getX(a, b, c, main.t2);
            y2 = org + getY(a, b, c, main.t2);
            dx0 = getdX(a, b, c, main.t1);
            dy0 = getdY(a, b, c, main.t1);
            dx2 = getdX(a, b, c, main.t2);
            dy2 = getdY(a, b, c, main.t2);

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

    private static double solve_quartic(double lead, double a, double b, double c, double d, double sgn)
    {
        double sol, R, D, E;

        a /= lead;
        b /= lead;
        c /= lead;
        d /= lead;
        System.out.println("\nquartic      a,b,c,d = " + a + ", " + b + ", " + c + ", " + d + ", " + sgn);
        sol = solve_cubic(-b, a*c - 4*d, -a*a*d + 4*b*d - c*c);
        R = Math.sqrt(a*a/4 - b + sol);
        D = Math.sqrt(3*a*a/4 - R*R - 2*b + (4*a*b - 8*c - a*a*a)/4/R);
        E = Math.sqrt(3*a*a/4 - R*R - 2*b - (4*a*b - 8*c - a*a*a)/4/R);
        System.out.println("\ncubic sol = " + sol + ", " + R + ", " + D + ", " + E);
//        System.out.println(main.t1/2/Math.PI + ", " + main.t2/2/Math.PI);
        System.out.print("root 1 = "); check_quartic(-a/4 + R/2 + D/2);
        System.out.print("root 2 = "); check_quartic(-a/4 + R/2 - D/2);
        System.out.print("root 3 = "); check_quartic(-a/4 - R/2 + E/2);
        System.out.print("root 4 = "); check_quartic(-a/4 - R/2 - E/2);
        if (!Double.isNaN(D) && (sgn == Math.signum(-a/4 + R/2 + D/2)))
        {
            System.out.println("using root 1 = " + sgn + ", " + (-a/4 + R/2 + D/2));
            return (-a/4 + R/2 + D/2);
        }
        else if (!Double.isNaN(E) && (sgn == Math.signum(-a/4 - R/2 - E/2)))
        {
            System.out.println("using root 4 = " + sgn + ", " + (-a/4 - R/2 - E/2));
            return (-a/4 - R/2 - E/2);
        }
        else if (!Double.isNaN(E) && (sgn == Math.signum(-a/4 - R/2 + E/2)))
        {
            System.out.println("using root 3 = " + sgn + ", " + (-a/4 - R/2 + E/2));
            return (-a/4 - R/2 + E/2);
        }
        else if (!Double.isNaN(D) && (sgn == Math.signum(-a/4 + R/2 - D/2)))
        {
            System.out.println("using root 2 = " + sgn + ", " + (-a/4 + R/2 - D/2));
            return (-a/4 + R/2 - D/2);
        }
        System.out.println("general quartic : Bad solution = " + R + ", " + D + ", " + E);
        return Double.NaN;
    }

    private static void check_quartic(double delx0)
    {
        double xrot0 = getrotX(ptBez[0].x, ptBez[0].y, (theta[0] + theta[1])/2);
        double yrot0 = getrotY(ptBez[0].x, ptBez[0].y, (theta[0] + theta[1])/2);
        double xrot3 = getrotX(ptBez[3].x, ptBez[3].y, (theta[0] + theta[1])/2);
        double yrot3 = getrotY(ptBez[3].x, ptBez[3].y, (theta[0] + theta[1])/2);
        double mrot0 = Math.tan((theta[0] - theta[1])/2);

        // calculate other quartic root (delx1), and cross-check
        double delx1 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*C[0]*delx0*delx0/2)/2/mrot0;
        System.out.print(delx0 + ", " + delx1);
        System.out.println(", " + (-yrot3 + yrot0 - mrot0*(xrot3 - xrot0) + delx0*2*mrot0 - 3*C[1]*delx1*delx1/2));
    }

    private static double solve_cubic(double p, double q, double r)
    {
        // see Math CRC book, page 392

        double a = (3*q - p*p)/3;
        double b = (2*p*p*p - 9*p*q + 27*r)/27;
        double d = b*b/4 + a*a*a/27;

//        System.out.println("\ncubic p,q,r = " + p + ", " + q + ", " + r);
//        System.out.println("\ncubic a,b,d = " + a + ", " + b + ", " + d);
        if (d < 0)
        {
            double phi = Math.acos(-b/2/Math.sqrt(-a*a*a/27));
//            System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-a/3)*Math.cos(phi/3) - p/3) + ", " + (2*Math.sqrt(-a/3)*Math.cos(phi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-a/3)*Math.cos(phi/3 + 4*Math.PI/3) - p/3));
            return 2*Math.sqrt(-a/3)*Math.cos(phi/3 + 2*Math.PI/3) - p/3;
        }
        else
        {
//            System.out.println("1 cubic d > 0 : " + (Math.cbrt(-b/2 + Math.sqrt(d)) + Math.cbrt(-b/2 - Math.sqrt(d)) - p/3));
            return Math.cbrt(-b/2 + Math.sqrt(d)) + Math.cbrt(-b/2 - Math.sqrt(d)) - p/3;
        }
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

    private static int insert_t_value(int N, int index, double[] t, double new_t)
    {
        // push t[index] up by one, insert at location index
        if (N > index)
            for (int i = N; i > index; i--)
                t[i] = t[i - 1];
        t[index] = new_t;
        return N + 1;
    }

    public static int get_t_values(double[] t, double a, double b, double c)
    {
        int N = 0;                                              // points per lobe

        if (b < 0)                                              // hypo
        {
            N = insert_t_value(N, 0, t, 0);                     // at max R
            if (-c > -b)                                        // extra loop
            {
                N = insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));   // at max width of small lobe
                System.out.println("hypo max width : " + t[1] + ", " + getX(a, b, c, t[1]) + ", " + getY(a, b, c, t[1]) + ", " + getdX(a, b, c, t[1]) + ", " + getdY(a, b, c, t[1]));
            }
            else if (c == b)                                    // cycloid
            {
                System.out.println("cycloid case : " + t[1] + ", " + getX(a, b, c, t[1]) + ", " + getY(a, b, c, t[1]) + ", " + getdX(a, b, c, t[1]) + ", " + getdY(a, b, c, t[1]));
            }
            else if (-c > b/(1 + a/b))                          // no loop
            {
                N = insert_t_value(N, 1, t, -b/a*Math.acos((1 + c*c*(1 + a/b)/b/b)*b/c/(2 + a/b)));  // inflection point
                System.out.println("hypo inflection : " + t[1] + ", " + getX(a, b, c, t[1]) + ", " + getY(a, b, c, t[1]) + ", " + (getdX(a, b, c, t[1])*getd2Y(a, b, c, t[1]) - getdY(a, b, c, t[1])*getd2X(a, b, c, t[1])));
            }
            else                                                // no inflection
            {
                System.out.println("hypo no inflection");
            }
        }
        else if (b > 0)                                             // epi
        {
            N = insert_t_value(N, 0, t, 0);                         // at max R
            if (-c > b)                                             // extra loop
            {
                N = insert_t_value(N, 1, t, Math.PI*b/a - calc_sin_t(-c/b, 1 + a/b));   // perpendicular to max width of small lobe
                System.out.println("epi perp to small width : " + t[1] + ", " + getX(a, b, c, t[1]) + ", " + getY(a, b, c, t[1]) + ", " + getdX(a, b, c, t[1]) + ", " + getdY(a, b, c, t[1]));
                N = insert_t_value(N, 2, t, Math.PI*b/a - calc_cos_t(-c/b, 1 + a/b));   // at max width of small lobe
                System.out.println("epi small width : " + t[2] + ", " + getX(a, b, c, t[2]) + ", " + getY(a, b, c, t[2]) + ", " + getdX(a, b, c, t[2]) + ", " + getdY(a, b, c, t[2]));
                N = insert_t_value(N, 2, t, (t[1] + t[2])/2);       // interpolate
                if (a < 2*b)
                    N = insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
            }
            else if (-c == b)                                       // cycloid
            {
                N = insert_t_value(N, 1, t, Math.PI*b/a - calc_sin_t(-c/b, 1 + a/b));   // perpendicular to max width of small lobe
                System.out.println("cycloid case : " + t[1] + ", " + getX(a, b, c, t[1]) + ", " + getY(a, b, c, t[1]) + ", " + getdX(a, b, c, t[1]) + ", " + getdY(a, b, c, t[1]));
                N = insert_t_value(N, 2, t, (t[1] + Math.PI*b/a)/2);   // interpolate
                if (a < 2*b)
                    N = insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
            }
            else if (-c > b/(1 + a/b))                              // no loop
            {
                N = insert_t_value(N, 1, t, Math.PI*b/a - calc_sin_t(-c/b, 1 + a/b));   // perpendicular to max width of small lobe
                System.out.println("epi perp to small width : " + t[1] + ", " + getX(a, b, c, t[1]) + ", " + getY(a, b, c, t[1]) + ", " + getdX(a, b, c, t[1]) + ", " + getdY(a, b, c, t[1]));
                N = insert_t_value(N, 2, t, b/a*Math.acos((1 + c*c*(1 + a/b)/b/b)*b/c/(2 + a/b)));  // inflection point
                System.out.println("epi inflection : " + t[2] + ", " + getX(a, b, c, t[2]) + ", " + getY(a, b, c, t[2]) + ", " + (getdX(a, b, c, t[2])*getd2Y(a, b, c, t[2]) - getdY(a, b, c, t[2])*getd2X(a, b, c, t[2])));
                N = insert_t_value(N, 2, t, (t[1] + t[2])/2);       // interpolate
                if (a < 2*b)
                {
                    N = insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
                    N = insert_t_value(N, 2, t, (t[1] + t[2])/2);   // interpolate
                }
                else
                    N = insert_t_value(N, 1, t, (t[0] + t[1])/2);   // interpolate
            }
            else                                                    // no inflection
            {
                if (a < 2*b)
                {
                    N = insert_t_value(N, 1, t, calc_cos_t(c/b, 1 + a/b));  // at max width of large lobe
                    N = insert_t_value(N, 2, t, (t[1] + Math.PI*b/a)/2);    // interpolate
                }
                else
                    N = insert_t_value(N, 1, t, Math.PI*b/a/2);     // interpolate
                System.out.println("epi no inflection");
            }
        }
        if (N > 0) t[N] = Math.PI*Math.abs(b)/a;                    // at min R
        if (N > 1)
            for (int i = 1; i < N; i++)
                t[N + i] = 2*Math.PI*Math.abs(b)/a - t[N - i];      // reflect first half
        return 2*N;
    }

    private static double getrotX(double m_x, double m_y, double m_theta)
    {
        return m_x*Math.cos(m_theta) + m_y*Math.sin(m_theta);
    }

    private static double getrotY(double m_x, double m_y, double m_theta)
    {
        return -m_x*Math.sin(m_theta) + m_y*Math.cos(m_theta);
    }

    private static double getX(double a, double b, double c, double t)
    {
        return (a + b)*Math.cos(t) - c*Math.cos(t*(a/b + 1));
    }

    private static double getY(double a, double b, double c, double t)
    {
        return (a + b)*Math.sin(t) - c*Math.sin(t*(a/b + 1));
    }

    private static double getdX(double a, double b, double c, double t)
    {
        return -(a + b)*Math.sin(t) + c*(a/b + 1)*Math.sin(t*(a/b + 1));
    }

    private static double getdY(double a, double b, double c, double t)
    {
        return (a + b)*Math.cos(t) - c*(a/b + 1)*Math.cos(t*(a/b + 1));
    }

    private static double getd2X(double a, double b, double c, double t)
    {
        return -(a + b)*Math.cos(t) + c*(a/b + 1)*(a/b + 1)*Math.cos(t*(a/b + 1));
    }

    private static double getd2Y(double a, double b, double c, double t)
    {
        return -(a + b)*Math.sin(t) + c*(a/b + 1)*(a/b + 1)*Math.sin(t*(a/b + 1));
    }
}
