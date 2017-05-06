
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

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\errorCycloid.java
// see also     : \Documents\NetBeansProjects\MyDemo\src\components\fitCycloid.java

public class errorCycloid
{
    private static final double c = 1.2;         // cycloid 'c'

    public static void main (String[] args)
    {
        double[] t;                             // t values at zero or maximum curvature
        Point2D.Double[][] Bezs;                // array of cubic Beziers
        Point2D.Double[] centers;               // array of 'center' points

        // define Beziers for each segment

        fitCycloid.set_c(c);                    // initiallize fitCycloid
        if (c > 1)                              // extra loop
        {
            t = new double[] {Math.acos(1/c), Math.PI};  // vertical slope
            Bezs = new Point2D.Double[][] {fitCycloid.fit_vertical_C0_to_horizontal_C1(t[0], 0),
                                           fitCycloid.fit_vertical_C0_to_horizontal_C1(t[0], Math.PI)};
            centers = new Point2D.Double[] {new Point2D.Double(0, 0), new Point2D.Double(Math.PI, 0)};
        }
        else if (c == 1)                        // cusp
        {
            t = new double[] {Math.PI};
            Bezs = new Point2D.Double[][] {{new Point2D.Double(0, 0), new Point2D.Double(0, 0),
                                            new Point2D.Double(Math.PI - Math.sqrt(16.0/3.0), 2), new Point2D.Double(Math.PI, 2)}};
            centers = new Point2D.Double[] {new Point2D.Double(Math.PI, -Math.PI*(Math.PI - Math.sqrt(16.0/3.0))/2)};
        }
        else if (c > 0.5)                       // inflection plus max curvature
        {
            t = new double[] {Math.acos(c), Math.acos((2*c*c - 1)/c), Math.PI};  // zero curvature, max negative curvature
            Bezs = new Point2D.Double[][] {fitCycloid.fit_inflect_to_d2ydx2(t[0], 0),
                                           fitCycloid.fit_inflect_to_d2ydx2(t[0], t[1]),
                                           fitCycloid.fit_d2ydx2_to_d2ydx2(t[1], Math.PI)};
            centers = new Point2D.Double[] {fitCycloid.get_center(0, t[0]), fitCycloid.get_center(t[0], t[1]), fitCycloid.get_center(t[1], Math.PI)};
        }
        else if (c > 0)                         // inflection point
        {
            t = new double[] {Math.acos(c), Math.PI};    // zero curvature
            Bezs = new Point2D.Double[][] {fitCycloid.fit_inflect_to_d2ydx2(t[0], 0),
                                           fitCycloid.fit_inflect_to_d2ydx2(t[0], Math.PI)};
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
        double theta_cycloid, theta_bez;                            // should be the same
        double x0, x1, x2, x3;
        double y0, y1, y2, y3;
        double tempa, tempb, tempc, tempd;
        double err = 0;
        int N = 200;

        try {
            FileWriter out = new FileWriter("\\Windows\\Temp\\cycloidout" + c + ".csv");
            out.write("cycloid : c = " + c + "\r\n");
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
    }
}
