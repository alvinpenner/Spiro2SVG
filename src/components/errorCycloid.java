
package components;

import java.awt.geom.Point2D;
import java.io.FileWriter;
import java.io.IOException;

// for a cycloid, calculate the two (or three) points that are perpendicular to the tangents at the endpoints
// where the endpoints are points of maximum or zero curvature.
// use these points as the 'center' of an 'arc' and calculate the distance to the cycloid and the Bezier.
// first calculate cycloid x, y, radius, theta as a function of t
// then calculate Bezier t, x, y, r, error as a function of theta.
// then calculate the rms error over half a cycle
// see Spiro2SVG Book2, Dec 2016
// use FIT_CURV to switch between a curvature fit and a center-of-mass fit

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\errorCycloid.java
// see also     : \Documents\NetBeansProjects\MyDemo\src\components\fitCycloid.java

public class errorCycloid
{
    private static final boolean FIT_CURV = true;         // fit curvature at endpoints
    private static final double c = .5;                    // cycloid 'c'

    public static void main (String[] args)
    {
        double[] t;                                     // t values at zero or maximum curvature
        Point2D.Double[][] Bezs;                        // array of cubic Beziers
        Point2D.Double[] centers;                       // array of 'center' points

        // define Beziers for each segment

        fitCycloid.set_c(c);                            // initiallize fitCycloid
        if (c > 1)                                      // extra loop
        {
            t = new double[] {Math.acos(1/c), Math.PI}; // vertical slope
            if (FIT_CURV)
                Bezs = new Point2D.Double[][] {fitCycloid.fit_vertical_C0_to_horizontal_C1(t[0], 0),
                                               fitCycloid.fit_vertical_C0_to_horizontal_C1(t[0], Math.PI)};
            else
                Bezs = new Point2D.Double[][] {setup_general_quartic_cofmx_cofmy(t[0], 0, 1),
                                               setup_general_quartic_cofmx_cofmy(t[0], Math.PI, 1)};
            centers = new Point2D.Double[] {new Point2D.Double(0, 0), new Point2D.Double(Math.PI, 0)};
        }
        else if (c == 1)                                // cusp
        {
            t = new double[] {Math.PI};
            if (FIT_CURV)
                Bezs = new Point2D.Double[][] {{new Point2D.Double(0, 0), new Point2D.Double(0, 0),
                                                new Point2D.Double(Math.PI - Math.sqrt(16.0/3.0), 2), new Point2D.Double(Math.PI, 2)}};
            else
                Bezs = new Point2D.Double[][] {setup_general_quartic_cofmx_cofmy(0, Math.PI, 1)};
            centers = new Point2D.Double[] {new Point2D.Double(Math.PI, -Math.PI*(Math.PI - Math.sqrt(16.0/3.0))/2)};
            fitCycloid.gen_Bezier(Bezs[0]);
        }
        else if (c > 0.5)                               // inflection plus max curvature
        {
            t = new double[] {Math.acos(c), Math.acos((2*c*c - 1)/c), Math.PI};  // zero curvature, max negative curvature
            if (FIT_CURV)
                Bezs = new Point2D.Double[][] {fitCycloid.fit_inflect_to_d2ydx2(t[0], 0),
                                               fitCycloid.fit_inflect_to_d2ydx2(t[0], t[1]),
                                               fitCycloid.fit_d2ydx2_to_d2ydx2(t[1], Math.PI)};
            else
                Bezs = new Point2D.Double[][] {setup_general_quartic_cofmx_cofmy(t[0], 0, -1),
                                               setup_general_quartic_cofmx_cofmy(t[0], t[1], 1),
                                               setup_general_quartic_cofmx_cofmy(t[1], Math.PI, 1)};
            centers = new Point2D.Double[] {fitCycloid.get_center(0, t[0]), fitCycloid.get_center(t[0], t[1]), fitCycloid.get_center(t[1], Math.PI)};
        }
        else if (c > 0)                                 // inflection point
        {
            t = new double[] {Math.acos(c), Math.PI};   // zero curvature
            if (FIT_CURV)
                Bezs = new Point2D.Double[][] {fitCycloid.fit_inflect_to_d2ydx2(t[0], 0),
                                               fitCycloid.fit_inflect_to_d2ydx2(t[0], Math.PI)};
            else
                Bezs = new Point2D.Double[][] {setup_general_quartic_cofmx_cofmy(t[0], 0, -1),
                                               setup_general_quartic_cofmx_cofmy(t[0], Math.PI, 1)};
            centers = new Point2D.Double[] {fitCycloid.get_center(0, t[0]), fitCycloid.get_center(t[0], Math.PI)};
        }
        else
        {
            System.out.println("c is negative, abort.");
            return;
        }
        System.out.printf("errorCycloid: c t[] =, %f", c);
        for (double tout : t)
            System.out.printf(" %f", tout);
        System.out.println();
        for (Point2D.Double[] Bez : Bezs)
        {
            System.out.println();
            for (Point2D.Double pt : Bez)
                System.out.println(pt.x + ", " + pt.y);
        }
        System.out.println();
        for (Point2D.Double center : centers)
            System.out.println(center.x + ", " + center.y);

        // scan the cycloid t values

        double t_cycloid, x_cycloid, y_cycloid, r_cycloid;
        double t_bez, x_bez, y_bez, r_bez;
        double theta_cycloid, theta_bez;                // should be the same
        double x0, x1, x2, x3;
        double y0, y1, y2, y3;
        double tempa, tempb, tempc, tempd;
        double err = 0;
        int N = 200;

        try {
            FileWriter out = new FileWriter("\\Windows\\Temp\\cycloidout" + c + ".csv");
            out.write("cycloid : c = " + c + " (FIT_CURV = " + FIT_CURV + ")\r\n");
            out.write("t_cycloid, index,  x_cycloid,  y_cycloid,  r_cycloid, theta_cycloid,    t_bez,    x_bez,    y_bez,    r_bez,    theta_bez,    error\r\n");
            int index = 0;                                          // segment id
            for (int i = 0; i <= N; i++)
            {
                t_cycloid = i*Math.PI/N;
                if (t_cycloid > t[index]) index++;
                x_cycloid = fitCycloid.C_x(t_cycloid);
                y_cycloid = fitCycloid.C_y(t_cycloid);
                out.write(t_cycloid + ", " + index + ", " +  x_cycloid + ", " + y_cycloid + ", ");
                x_cycloid -= centers[index].x;                      // re-define relative to center
                y_cycloid -= centers[index].y;
                r_cycloid = Math.sqrt(x_cycloid*x_cycloid + y_cycloid*y_cycloid);
                theta_cycloid = Math.atan2(y_cycloid, x_cycloid);   // common to both cycloid and bezier
                out.write(r_cycloid + ", " + theta_cycloid*180/Math.PI + ", ");
                x0 = Bezs[index][0].x - centers[index].x;
                x1 = Bezs[index][1].x - centers[index].x;
                x2 = Bezs[index][2].x - centers[index].x;
                x3 = Bezs[index][3].x - centers[index].x;
                y0 = Bezs[index][0].y - centers[index].y;
                y1 = Bezs[index][1].y - centers[index].y;
                y2 = Bezs[index][2].y - centers[index].y;
                y3 = Bezs[index][3].y - centers[index].y;
                tempa = (-x0 + 3*x1 - 3*x2 + x3)*y_cycloid - (-y0 + 3*y1 - 3*y2 + y3)*x_cycloid;
                tempb = (3*x0 - 6*x1 + 3*x2)*y_cycloid - (3*y0 - 6*y1 + 3*y2)*x_cycloid;
                tempc = (-3*x0 + 3*x1)*y_cycloid - (-3*y0 + 3*y1)*x_cycloid;
                tempd = x0*y_cycloid - y0*x_cycloid;
                //System.out.println("i  = " + i + ", " + t_cycloid + ", " + x_cycloid + ", " + y_cycloid + ", " + tempa + ", " + tempb + ", " + tempc + ", " + tempd);
                //System.out.println("x0 = " + x0 + ", " + x1 + ", " + x2 + ", " + x3 + ", " + y0 + ", " + y1 + ", " + y2 + ", " + y3);
                t_bez = Beziererror.solve_cubic(tempb/tempa, tempc/tempa, tempd/tempa); // t as a function of theta_cycloid
                x_bez = x0*(1-t_bez)*(1-t_bez)*(1-t_bez) + 3*x1*t_bez*(1-t_bez)*(1-t_bez) + 3*x2*t_bez*t_bez*(1-t_bez) + x3*t_bez*t_bez*t_bez;
                y_bez = y0*(1-t_bez)*(1-t_bez)*(1-t_bez) + 3*y1*t_bez*(1-t_bez)*(1-t_bez) + 3*y2*t_bez*t_bez*(1-t_bez) + y3*t_bez*t_bez*t_bez;
                r_bez = Math.sqrt(x_bez*x_bez + y_bez*y_bez);
                theta_bez = Math.atan2(y_bez, x_bez);           // just a double check
                x_bez += centers[index].x;
                y_bez += centers[index].y;
                err += (r_bez - r_cycloid)*(r_bez - r_cycloid);
                out.write(t_bez + ", " +  x_bez + ", " + y_bez + ", " + r_bez + ", " + theta_bez*180/Math.PI + ", " + (r_bez - r_cycloid) + "\r\n");
            }
            out.close(); 
        } catch (IOException e) {System.err.println("Unable to save file: " + e);}
        System.out.println("c err = ," + c + ", " + Math.sqrt(err/N));

        //setup_general_quartic_cofmx_cofmy(Math.acos((2*c*c - 1)/c), Math.PI);
        //setup_general_quartic_cofmx_cofmy(Math.acos((2*c*c - 1)/c), Math.acos(c));
        //setup_general_quartic_cofmx_cofmy(Math.acos(c), Math.PI, 1);
        //setup_general_quartic_cofmx_cofmy(Math.acos(c), Math.acos((2*c*c - 1)/c));
        //setup_general_quartic_cofmx_cofmy(Math.PI, Math.acos(c));
        //setup_general_quartic_cofmx_cofmy(Math.acos(1/c), 0);
        //setup_general_quartic_cofmx_cofmy(Math.PI, Math.PI - 1);
        //setup_general_quartic_cofmx_cofmy(0, Math.PI, 1);
        //variational_coeffs();
    }

    private static Point2D.Double[] setup_general_quartic_cofmx_cofmy(double t1, double t2, double sgn)
    {
        // this will fit the tilted <x>/<1> center-of-mass at two different angles, May 9, 2017
        // used only for the cycloid fit (see Spiro2SVG Book 2)

        Point2D.Double[] ptBez = new Point2D.Double[4]; // Point[0-3] (origin @ bottom-left)
        final Point2D.Double center = fitCycloid.get_center(t1, t2);
        final double l1 = fitCycloid.get_l1(t1, t2);
        final double l2 = fitCycloid.get_l1(t2, t1);
        // for convex upwards, use +PI/2; for convex downwards, use -PI/2
        double theta1 = Math.atan2(c*Math.sin(t1), 1 - c*Math.cos(t1)) + sgn*Math.PI/2; // use d/dt
        double theta2 = Math.atan2(c*Math.sin(t2), 1 - c*Math.cos(t2)) + sgn*Math.PI/2; // use d/dt

        if ((1 == c*Math.cos(t1)) && (0 == c*Math.sin(t1)))     // true only at cusp
            theta1 = Math.PI/2 + sgn*Math.PI/2;
        if ((1 == c*Math.cos(t2)) && (0 == c*Math.sin(t2)))     // true only at cusp
            theta2 = Math.PI/2 + sgn*Math.PI/2;

        // express 20<1>/3 = a0 = -a1*d1 + a2*d2 + a3*d1*d2
        final double a0 = 20*cycloid_area(t1, t2)/3;
        final double a1 = -2*(l1 - l2*Math.cos(theta2 - theta1));
        final double a2 =  2*(l2 - l1*Math.cos(theta2 - theta1));
        final double a3 = -Math.sin(theta2 - theta1);

        // express 280<x>@θ1 = r1*d1 + r3*d2 + r4*d1*d2 + r6*d2**2 + r7*d1*d2**2
        final double r1 = (50*l1 + 34*l2*Math.cos(theta2 - theta1))*(l1 - l2*Math.cos(theta2 - theta1));
        final double r3 = (34*l1 + 50*l2*Math.cos(theta2 - theta1))*(l2 - l1*Math.cos(theta2 - theta1));
        final double r4 = (-9*l1 - 33*l2*Math.cos(theta2 - theta1))*Math.sin(theta2 - theta1);
        final double r6 = 15*(l2 - l1*Math.cos(theta2 - theta1))*Math.sin(theta2 - theta1);
        final double r7 = -9*Math.sin(theta2 - theta1)*Math.sin(theta2 - theta1);

        // express 280<x>@θ2 = s1*d1 + s2*d1**2 + s3*d2 + s4*d1*d2 + s5*d1**2*d2
        final double s3 = (50*l2 + 34*l1*Math.cos(theta2 - theta1))*(l2 - l1*Math.cos(theta2 - theta1));
        final double s1 = (34*l2 + 50*l1*Math.cos(theta2 - theta1))*(l1 - l2*Math.cos(theta2 - theta1));
        final double s4 = (-9*l2 - 33*l1*Math.cos(theta2 - theta1))*Math.sin(theta2 - theta1);
        final double s2 = 15*(l1 - l2*Math.cos(theta2 - theta1))*Math.sin(theta2 - theta1);
        final double s5 = -9*Math.sin(theta2 - theta1)*Math.sin(theta2 - theta1);

        // express d1   = (u1*d2 + u2*d2*d2)/(u3 + u4*d2 + u5*d2*d2) - transformed to angle θ1
        final double u0 = 280*transform_x_moment(t1, t2, center, theta1);   // transformed <x> moment at θ1
        final double u1 = -r3 + a2*u0/a0;
        final double u2 = -r6;
        final double u3 =  r1 + a1*u0/a0;
        final double u4 =  r4 - a3*u0/a0;
        final double u5 =  r7;

        // express d2   = (b1*d1 + b2*d1*d1)/(b3 + b4*d1 + b5*d1*d1) - transformed to angle θ2
        final double b0 = 280*transform_x_moment(t1, t2, center, theta2);   // transformed <x> moment at θ2
        final double b1 = -s1 - a1*b0/a0;
        final double b2 = -s2;
        final double b3 =  s3 - a2*b0/a0;
        final double b4 =  s4 - a3*b0/a0;
        final double b5 =  s5;

        System.out.println("\nsetup_general_quartic_cofmx_cofmy : " + c + ", " + center.x + ", " + center.y);
        System.out.println("t1 t2 l1 l2 theta1 theta2         : ," + t1 + ", " + t2 + ", " + l1 + ", " + l2 + ", " + 180*theta1/Math.PI + ", " + 180*theta2/Math.PI);
        System.out.println("<x> <y> tr1<x> tr2<x>             : ," + cycloid_x_moment(t1, t2) + ", " + cycloid_y_moment(t1, t2) + ", " + transform_x_moment(t1, t2, center, theta1) + ", " + transform_x_moment(t1, t2, center, theta2));
        fitCycloid.gen_point(new Point2D.Double(cycloid_x_moment(t1, t2)/cycloid_area(t1, t2), cycloid_y_moment(t1, t2)/cycloid_area(t1, t2)));
        fitCycloid.gen_point(center);

        double d1, d2;
        d1 = fitymoment.solve_quartic(u3*b5*b5 + u4*b2*b5 + u5*b2*b2,
                                      2*u3*b4*b5 + u4*(b1*b5 + b2*b4) - u1*b2*b5 + 2*u5*b1*b2 - u2*b2*b2,
                                      u3*(2*b3*b5 + b4*b4) + u4*(b1*b4 + b2*b3)
                                    - u1*(b1*b5 + b2*b4) + u5*b1*b1 - 2*u2*b1*b2,
                                      2*u3*b3*b4 + u4*b1*b3
                                    - u1*(b1*b4 + b2*b3) - u2*b1*b1,
                                      u3*b3*b3 - u1*b1*b3,
                                      t1 > t2);
        d2 = (b1*d1 + b2*d1*d1)/(b3 + b4*d1 + b5*d1*d1);
        System.out.println("solve_quartic c,  d1,  d2 = " + c + ", " + d1 + ", " + d2);
        System.out.println("cycloid     <1>, <x>, <y> = " + 3*a0/20 + ", " + u0/280 + ", " + b0/280);
        System.out.println("Bezier      <1>, <x>, <y> = " + 3*(-a1*d1 + a2*d2 + a3*d1*d2)/20 + ", " +
                                                            (r3*d2 + r6*d2*d2 + (r1 + r4*d2 + r7*d2*d2)*d1)/280 + ", " +
                                                            (s1*d1 + s2*d1*d1 + (s3 + s4*d1 + s5*d1*d1)*d2)/280);

        ptBez[0] = new Point2D.Double(fitCycloid.C_x(t1), fitCycloid.C_y(t1));
        ptBez[1] = new Point2D.Double(fitCycloid.C_x(t1) - d1*Math.sin(theta1), fitCycloid.C_y(t1) + d1*Math.cos(theta1));
        ptBez[2] = new Point2D.Double(fitCycloid.C_x(t2) + d2*Math.sin(theta2), fitCycloid.C_y(t2) - d2*Math.cos(theta2));
        ptBez[3] = new Point2D.Double(fitCycloid.C_x(t2), fitCycloid.C_y(t2));
        fitCycloid.gen_Bezier(ptBez);
        return ptBez;
    }

    private static double transform_x_moment(double t1, double t2, Point2D.Double c, double theta)
    {
        return Math.cos(theta)*(cycloid_x_moment(t1, t2) - c.x*cycloid_area(t1, t2))
             + Math.sin(theta)*(cycloid_y_moment(t1, t2) - c.y*cycloid_area(t1, t2));
    }

    private static double cycloid_area(double t1, double t2)
    {
        return cycloid_area_fxn(t2) - cycloid_area_fxn(t1) - (fitCycloid.C_x(t2) - fitCycloid.C_x(t1))*(fitCycloid.C_y(t2) + fitCycloid.C_y(t1))/2;
    }

    private static double cycloid_area_fxn(double t)
    {
        return t - 2*c*Math.sin(t) + c*c*(t + Math.sin(t)*Math.cos(t))/2;
    }

    private static double cycloid_x_moment(double t1, double t2)
    {
        return cycloid_x_fxn(t2) - cycloid_x_fxn(t1) - (fitCycloid.C_x(t2) - fitCycloid.C_x(t1))
             * (fitCycloid.C_y(t1)*(2*fitCycloid.C_x(t1) + fitCycloid.C_x(t2))
             +  fitCycloid.C_y(t2)*(fitCycloid.C_x(t1) + 2*fitCycloid.C_x(t2)))/6;
    }

    private static double cycloid_x_fxn(double t)
    {
        return t*t/2 - c*(2*t*Math.sin(t) + Math.cos(t))
             + c*c*(t*t + 2*t*Math.sin(t)*Math.cos(t) + 3*Math.sin(t)*Math.sin(t))/4
             + c*c*c*Math.cos(t)*Math.cos(t)*Math.cos(t)/3;
    }

    private static double cycloid_y_moment(double t1, double t2)
    {
        return cycloid_y_fxn(t2) - cycloid_y_fxn(t1) - (fitCycloid.C_x(t2) - fitCycloid.C_x(t1))
             * (fitCycloid.C_y(t1)*fitCycloid.C_y(t1) + fitCycloid.C_y(t1)*fitCycloid.C_y(t2) + fitCycloid.C_y(t2)*fitCycloid.C_y(t2))/6;
    }

    private static double cycloid_y_fxn(double t)
    {
        return t/2 - 3*c*Math.sin(t)/2
             + 3*c*c*(t + Math.sin(t)*Math.cos(t))/4
             - c*c*c*Math.sin(t)*(1 - Math.sin(t)*Math.sin(t)/3)/2;
    }

    private static void variational_coeffs()
    {
        // calculate cubic coeffs for d from a variational calc
        // see Spiro2SVG Book2, page 46

        double c0 = 324*(I(4, 8) + I(6, 6));
        double c1 = 324*(I(4, 6) + 4*I(5, 7));
        double c2 = -36*I(2, 4)
                  +  18*(I(4, 6) + 4*I(5, 6) + 4*I(6, 6))
                  +  18*(I(2, 8) + 4*I(3, 8) + 4*I(4, 8))
                  +  18*(9*I(2, 8) - 12*I(2, 9) + 4*I(2, 10))
                  +  18*(9*I(4, 6) - 12*I(4, 7) + 4*I(4, 8))
                  +  36*(I(4, 4) + 8*I(4, 5) + 8*I(4, 6) - 32*I(4, 7) + 16*I(4, 8));
        double c3 =  -6*(I(2, 2) + 4*I(2, 3) - 4*I(2, 4))
                  +   6*(I(2, 6) + 8*I(3, 6) + 16*I(4, 6) - 16*I(6, 6))
                  +   6*(9*I(2, 6) + 24*I(2, 7) - 80*I(2, 8) + 64*I(2, 9) - 16*I(2, 10));
        System.out.println("\nvariational_coeffs = ," + c0 + ", " + c1 + ", " + c2 + ", " + c3);
        System.out.println("cubic root = " + Beziererror.solve_cubic(c1/c0, c2/c0, c3/c0));
        System.out.println("cubic root = " + Beziererror.solve_cubic(3996./1161., 1314./1161., -2138./1161.));
    }

    private static double I(int i, int j)
    {
        return factorial(i)*factorial(j)/factorial(i + j + 1);
    }

    private static double factorial(int N)
    {
        double fact = 1;
        if (N < 1) return 0;
        for (int i = 1; i <= N; i++)
            fact *= i;
        return fact;
    }
}
