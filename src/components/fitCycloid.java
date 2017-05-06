
package components;

import java.awt.geom.Point2D;

// calculate arm length of a Bezier curve (quartic in d1)
// by fitting the curvature of a cycloid at the endpoints
// see middle of Spiro2SVG Book 2, Feb 2017
// x = t - c*sin(t)
// y = 1 - c*cos(t)

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\fitCycloid.java

public class fitCycloid
{
    private static final double TOL = 1E-6;
    private static final Point2D.Double origin = new Point2D.Double(100, 500);
//    private static final Point2D.Double origin = new Point2D.Double(0, 0);
    private static final double size = 200; // 1;
    private static double c;                                // spiro 'c', dimensionless

    public static void main (String[] args)
    {
        c = 0.9;
        double t1 = Math.acos(c);                           // zero curvature (c < 1)
        double t2 = Math.acos((2*c*c - 1)/c);               // max negative curvature (.5 < c < 1)
        t1 = Math.acos(1/c);                                // vertical slope (c > 1)

        System.out.println("fitCycloid");
        System.out.printf("t c x y dydx dxdy d2ydx2 d2xdy2 =, %.12f, %f, %.10f, %.10f, %g, %g, %g, %g\n", 0F, c, C_x(0), C_y(0), C_dydx(0), C_dxdy(0), C_d2ydx2(0), C_d2xdy2(0));
//        System.out.printf("t c x y dydx dxdy d2ydx2 d2xdy2 =, %.12f, %f, %.10f, %.10f, %g, %g, %g, %g\n", t1, c, C_x(t1), C_y(t1), C_dydx(t1), C_dxdy(t1), C_d2ydx2(t1), C_d2xdy2(t1));
        System.out.printf("t c x y dydx dxdy d2ydx2 d2xdy2 =, %.12f, %f, %.10f, %.10f, %g, %g, %g, %g\n", t2, c, C_x(t2), C_y(t2), C_dydx(t2), C_dxdy(t2), C_d2ydx2(t2), C_d2xdy2(t2));
        System.out.printf("t c x y dydx dxdy d2ydx2 d2xdy2 =, %.12f, %f, %.10f, %.10f, %g, %g, %g, %g\n", Math.PI, c, C_x(Math.PI), C_y(Math.PI), C_dydx(Math.PI), C_dxdy(Math.PI), C_d2ydx2(Math.PI), C_d2xdy2(Math.PI));
        gen_points(0, Math.PI, 90);
//        fit_inflect_to_d2ydx2(t1, t2);
//        fit_inflect_to_d2ydx2(t1, 0);
        fit_d2ydx2_to_d2ydx2(t2, Math.PI);
//        fit_d2ydx2_to_d2ydx2(Math.PI, t2);
//        fit_vertical_C0_to_horizontal_C1(t1, Math.PI);
//        fit_vertical_C0_to_horizontal_C1(t1, 0);
    }

    public static void set_c (double m_c)
    {
        c = m_c;
    }

    public static Point2D.Double[] fit_inflect_to_d2ydx2(double t1, double t2)
    {
        // assume t1 is an inflection point
        // assume t2 has finite d2ydx2

        double delx0, delx3;                            // horizontal Bezier arm length
        Point2D.Double[] ptBez = new Point2D.Double[4]; // Point[0-3] (origin @ bottom-left)

        if (Math.abs(C_d2ydx2(t1)) > TOL)
        {
            System.out.println("t1 not an inflection point, abort.");
            return null;
        }
        if (Math.abs(C_dxdy(t2)) < TOL)
        {
            System.out.println("t2 slope too high, abort.");
            return null;
        }
        delx3 =  (C_y(t2) - C_y(t1) - C_dydx(t1)*(C_x(t2) - C_x(t1)))/(C_dydx(t2) - C_dydx(t1));
        delx0 = -(C_y(t2) - C_y(t1) - C_dydx(t2)*(C_x(t2) - C_x(t1)) + 3*C_d2ydx2(t2)*delx3*delx3/2)/(C_dydx(t2) - C_dydx(t1));
        System.out.printf("M %f, %f C %f, %f %f, %f %f, %f\n",
                          origin.x + size*C_x(t1), origin.y - size*C_y(t1),
                          origin.x + size*(C_x(t1) + delx0), origin.y - size*(C_y(t1) + C_dydx(t1)*delx0),
                          origin.x + size*(C_x(t2) - delx3), origin.y - size*(C_y(t2) - C_dydx(t2)*delx3),
                          origin.x + size*C_x(t2), origin.y - size*C_y(t2));
        ptBez[0] = new Point2D.Double(C_x(t1), C_y(t1));
        ptBez[1] = new Point2D.Double(C_x(t1) + delx0, C_y(t1) + C_dydx(t1)*delx0);
        ptBez[2] = new Point2D.Double(C_x(t2) - delx3, C_y(t2) - C_dydx(t2)*delx3);
        ptBez[3] = new Point2D.Double(C_x(t2), C_y(t2));
        return ptBez;
    }

    public static Point2D.Double[] fit_d2ydx2_to_d2ydx2(double t1, double t2)
    {
        // assume t1 and t2 do not have infinite slope

        double delx0, delx3;                            // horizontal Bezier arm length
        Point2D.Double[] ptBez = new Point2D.Double[4]; // Point[0-3] (origin @ bottom-left)

        if (Math.abs(C_dxdy(t1)) < TOL)
        {
            System.out.println("t1 slope too high, abort.");
            return null;
        }
        if (Math.abs(C_dxdy(t2)) < TOL)
        {
            System.out.println("t2 slope too high, abort.");
            return null;
        }
        delx0 = fitymoment.solve_quartic(-27*C_d2ydx2(t2)*C_d2ydx2(t1)*C_d2ydx2(t1)/8,
                                          0,
                                          9*C_d2ydx2(t2)*C_d2ydx2(t1)*(C_y(t2) - C_y(t1) - C_dydx(t1)*(C_x(t2) - C_x(t1)))/2,
                                         -(C_dydx(t2) - C_dydx(t1))*(C_dydx(t2) - C_dydx(t1))*(C_dydx(t2) - C_dydx(t1)),
                                         -(C_y(t2) - C_y(t1) - C_dydx(t2)*(C_x(t2) - C_x(t1)))*(C_dydx(t2) - C_dydx(t1))*(C_dydx(t2) - C_dydx(t1))
                                         -3*C_d2ydx2(t2)*(C_y(t2) - C_y(t1) - C_dydx(t1)*(C_x(t2) - C_x(t1)))
                                                        *(C_y(t2) - C_y(t1) - C_dydx(t1)*(C_x(t2) - C_x(t1)))/2,
                                          true);
        delx3 = (C_y(t2) - C_y(t1) - C_dydx(t1)*(C_x(t2) - C_x(t1)) - 3*C_d2ydx2(t1)*delx0*delx0/2)/(C_dydx(t2) - C_dydx(t1));
        System.out.println("delx0/3 = ," + c + ", " + delx0 + ", " + delx0*Math.sqrt(1 + C_dydx(t1)*C_dydx(t1)) + ", " + delx3);
        System.out.printf("M %f, %f C %f, %f %f, %f %f, %f\n",
                          origin.x + size*C_x(t1), origin.y - size*C_y(t1),
                          origin.x + size*(C_x(t1) + delx0), origin.y - size*(C_y(t1) + C_dydx(t1)*delx0),
                          origin.x + size*(C_x(t2) - delx3), origin.y - size*(C_y(t2) - C_dydx(t2)*delx3),
                          origin.x + size*C_x(t2), origin.y - size*C_y(t2));
        ptBez[0] = new Point2D.Double(C_x(t1), C_y(t1));
        ptBez[1] = new Point2D.Double(C_x(t1) + delx0, C_y(t1) + C_dydx(t1)*delx0);
        ptBez[2] = new Point2D.Double(C_x(t2) - delx3, C_y(t2) - C_dydx(t2)*delx3);
        ptBez[3] = new Point2D.Double(C_x(t2), C_y(t2));
        return ptBez;
    }

    public static Point2D.Double[] fit_vertical_C0_to_horizontal_C1(double t1, double t2)
    {
        // assume t1 has infinite slope
        // assume t2 has zero slope

        double dely0, delx3;                            // vertical/horizontal Bezier arm length
        Point2D.Double[] ptBez = new Point2D.Double[4]; // Point[0-3] (origin @ bottom-left)

        if (Math.abs(C_dxdy(t1)) > TOL)
        {
            System.out.println("t1 slope too low, abort.");
            return null;
        }
        if (Math.abs(C_dydx(t2)) > TOL)
        {
            System.out.println("t2 slope too high, abort.");
            return null;
        }
        dely0 = fitymoment.solve_quartic(-27*C_d2ydx2(t2)*C_d2xdy2(t1)*C_d2xdy2(t1)/8,
                                          0,
                                          9*C_d2ydx2(t2)*C_d2xdy2(t1)*(C_x(t2) - C_x(t1))/2,
                                          1,
                                         -(C_y(t2) - C_y(t1))
                                         -3*C_d2ydx2(t2)*(C_x(t2) - C_x(t1))*(C_x(t2) - C_x(t1))/2,
                                         t1 > t2);
        delx3 = C_x(t2) - C_x(t1) - 3*C_d2xdy2(t1)*dely0*dely0/2;
        System.out.println("dely0/x3 = ," + c + ", " + dely0 + ", " + delx3);
        System.out.printf("M %f, %f C %f, %f %f, %f %f, %f\n",
                          origin.x + size*C_x(t1), origin.y - size*C_y(t1),
                          origin.x + size*C_x(t1), origin.y - size*(C_y(t1) + dely0),
                          origin.x + size*(C_x(t2) - delx3), origin.y - size*C_y(t2),
                          origin.x + size*C_x(t2), origin.y - size*C_y(t2));
        ptBez[0] = new Point2D.Double(C_x(t1), C_y(t1));
        ptBez[1] = new Point2D.Double(C_x(t1), C_y(t1) + dely0);
        ptBez[2] = new Point2D.Double(C_x(t2) - delx3, C_y(t2));
        ptBez[3] = new Point2D.Double(C_x(t2), C_y(t2));
        return ptBez;
    }

    public static Point2D.Double get_center(double t1, double t2)
    {
        // t1, t2 are cycloid segment endpoints

        return new Point2D.Double((C_dydx(t2)*C_x(t1) - C_dydx(t1)*C_x(t2) - C_dydx(t1)*C_dydx(t2)*(C_y(t2) - C_y(t1)))/(C_dydx(t2) - C_dydx(t1)),
                                  (C_x(t2) - C_x(t1) + C_dydx(t2)*C_y(t2) - C_dydx(t1)*C_y(t1))/(C_dydx(t2) - C_dydx(t1)));
    }

    private static void gen_points(double t1, double t2, int N)
    {
        System.out.printf("M");
        for (int i = 0; i <= N; i++)
            System.out.printf(" %f, %f", origin.x + size*C_x(t1 + i*(t2 - t1)/N), origin.y - size*C_y(t1 + i*(t2 - t1)/N));
        System.out.println();
    }

    public static double C_x(double t)
    {
        return t - c*Math.sin(t);
    }

    public static double C_y(double t)
    {
        return 1 - c*Math.cos(t);
    }

    private static double C_dydx(double t)
    {
        return c*Math.sin(t)/(1 - c*Math.cos(t));
    }

    private static double C_dxdy(double t)
    {
        return 1/C_dydx(t);
    }

    private static double C_d2ydx2(double t)
    {
        return -c*(c - Math.cos(t))/(1 - c*Math.cos(t))/(1 - c*Math.cos(t))/(1 - c*Math.cos(t));
    }

    private static double C_d2xdy2(double t)
    {
        return c*(c - Math.cos(t))/(c*Math.sin(t))/(c*Math.sin(t))/(c*Math.sin(t));
    }
}
