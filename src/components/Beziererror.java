
package components;

// Sept 2016, 'Bezier fit' paper
// calculate spiro x, y, radius, theta as a function of t
// then calculate Bezier t, x, y, r, error as a function of theta
// see middle of CofM book, Sep 2016
// use the asymmetric Bezier construction for the asymmetric spiro, with arc 45º
// scan d1, calculate d2 as a function of d1, keeping bezier area fixed
// then plot bezier curvatures, moments, error as function of d1

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\Beziererror.java
// see also     : \Documents\NetBeansProjects\MyDemo\src\components\fitymoment.java

import java.io.FileWriter;

public class Beziererror
{
    private static final double TOL = 1E-9;
    private static final double a_b = 180;          // spiro 'a - b'
    private static final double c = -2; //3.371577; //3.5297135;      // spiro 'c'
    private static final double l1 = a_b + c;       // distance to start point (l1, 0)
    private static final double l2 = a_b - c;       // distance to end   point (l2/√2, l2/√2)

    private static final double theta = Math.PI/8;                         // 22.5º
    // coeff of Bezier area: <1>
    private static final double a1 = -2*(l1 - l2*Math.cos(2*theta));
    private static final double a2 =  2*(l2 - l1*Math.cos(2*theta));
    private static final double a3 = -Math.sin(2*theta);
    // coeff of amti-symmetric Bezier moment: <y>
    private static final double e1 = (l1 - l2*Math.cos(2*theta))*(-50*l1 + 34*l2)*Math.sin(theta);
    private static final double e2 = 15*(l1 - l2*Math.cos(2*theta))*Math.cos(theta);
    private static final double e3 = (l2 - l1*Math.cos(2*theta))*(-34*l1 + 50*l2)*Math.sin(theta);
    private static final double e4 = (l1 - l2)*(21*Math.sin(2*theta)*Math.sin(theta) - 12*(1 + Math.cos(2*theta))*Math.cos(theta));
    private static final double e5 = -9*Math.sin(2*theta)*Math.cos(theta);
    private static final double e6 = -15*(l2 - l1*Math.cos(2*theta))*Math.cos(theta);
    private static final double e7 = 9*Math.sin(2*theta)*Math.cos(theta);
    // coeff of symmetric Bezier moment: <x>
    private static final double f1 = (l1 - l2*Math.cos(2*theta))*(50*l1 + 34*l2)*Math.cos(theta);
    private static final double f2 = 15*(l1 - l2*Math.cos(2*theta))*Math.sin(theta);
    private static final double f3 = (l2 - l1*Math.cos(2*theta))*(34*l1 + 50*l2)*Math.cos(theta);
    private static final double f4 = (l1 + l2)*(-21*Math.sin(2*theta)*Math.cos(theta) + 12*(1 - Math.cos(2*theta))*Math.sin(theta));
    private static final double f5 = -9*Math.sin(2*theta)*Math.sin(theta);
    private static final double f6 = 15*(l2 - l1*Math.cos(2*theta))*Math.sin(theta);
    private static final double f7 = -9*Math.sin(2*theta)*Math.sin(theta);

    public static void main (String[] args)
    {
        double d1, d2;
        double start = 0;                           // d1 range
        double incr = 1;
        int steps = 0;
        double err0 = 0, err1 = 0, err2;            // previous error levels
        double d2old = 0;
        String strout = "";                         // record extrema

        System.out.printf("a-b c =, %f, %f\n", a_b, c);
        System.out.println("    ,         ,       C0,         C1,        area,        <x>,        <y>");
        System.out.printf( "    ,         ,       %f, %f, %f, %f, %f\n\n", spiro_curvature(c), spiro_curvature(-c), spiro_area(), spiro_moment("x"), spiro_moment("y"));

        if (steps > 0)
        {
            System.out.println("  d1,       d2,       C0,         C1,        area,        <x>,        <y>,        err");
            for (d1 = start; d1 <= start + steps*incr; d1 += incr)
            {
//                d2 = -l2 - (3*spiro_curvature(c)*d1*d1/2 - l1)*Math.sqrt(2);    // satisfy C0 constraint
                d2 = (20*spiro_area()/3 - 2*d1*(l1 - l2/Math.sqrt(2)))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2)); // satisfy area constraint
                err2 = calc_error(d1, d2);                  // search for extrema in rms error
                System.out.printf("%f, %f, %f, %f, %f, %f, %f, %.8f\n", d1, d2, bezier_C0(d1, d2), bezier_C1(d1, d2), bezier_area(d1, d2), bezier_moment_x(d1, d2), bezier_moment_y(d1, d2), err2);
//                err2 = bezier_moment_y(d1, d2);             // search for extrema in <y>
                if ( d1 > start + incr)
                {
                    if (err2 > err1 && err0 > err1) strout += "\nlocal min at, " + (d1 - incr) + ", " + d2old + ", " + err1;
                    if (err2 < err1 && err0 < err1) strout += "\nlocal max at, " + (d1 - incr) + ", " + d2old + ", " + err1;
                }
                err0 = err1;
                err1 = err2;
                d2old = d2;
            }
            if (!strout.isEmpty())
                System.out.println("\nList of extrema:" + strout);
        }
//        check_neighbours(50, 39.5);
//        d2 = (20*spiro_area()/3 - 2*start*(l1 - l2/Math.sqrt(2)))/(2*(l2 - l1/Math.sqrt(2)) - start/Math.sqrt(2)); // satisfy area constraint
//        calc_one(start, d2);
        calc_one(43.64, 50.72977);
//        scan_moment_y(54, 0, 0.03, 100);
//        scan_moment_y(-350, 0, 1, 1000);
//        scan_area(50, 0, 0.1, 100);
//        scan_moment_y(-350, 0, 1, 1000, e1, e2, e3, e4, e5, e6, e7, false);
//        calc_array();       // 2D array = calc_error(d1, d2) plus four 1D scans
    }

    private static void calc_one(double d1, double d2)
    {
//        d1 /= Math.cos(Math.PI/8);            // use only if data comes from Spiro2SVG
//        d2 /= Math.cos(Math.PI/8);
//        System.out.printf("c         ,d1        ,d2        ,spiro_area  ,spiro <x>     ,spiro <y>     ,bezier_area ,bezier <x>    ,bezier <y>    ,err\n");
//        System.out.printf("%f ,%f ,%f ,%f ,%f ,%f ,%f ,%f ,%f ,%.7f\n", c, d1, d2, spiro_area(), spiro_moment("x"), spiro_moment("y"), bezier_area(d1, d2), bezier_moment_x(d1, d2), bezier_moment_y(d1, d2), calc_error(d1, d2));
        System.out.printf("c         ,d1        ,d2        ,err\n");
        System.out.printf("%f ,%f ,%f ,%.8f\n", c, d1, d2, calc_error(d1, d2));
    }

    private static void calc_array()
    {
        // generate data suitable for matplotlib
        // rms error as function of d1,d2

        double a0 = 20*spiro_area()/3;
        double e0 = -280*spiro_moment("y");
        double f0 = -280*spiro_moment("x");
        double d1, d2;
//        double startd1 = 0, startd2 = 0;
//        double incr = 1;
        double startd1 = 56.5, startd2 = 28;
        double incr = .02;
//        double startd1 = 56, startd2 = 28;
//        double incr = .04;
        int steps = 100;

        try
        {
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\scan_error.txt");
            java.text.DateFormat df = java.text.DateFormat.getInstance();
            out.write(df.format(new java.util.Date()) + " : a-b, c = " + a_b + ", " + c + "\r\n");
            out.write("extent = (, " + startd1 + ", " + (startd1 + steps*incr) + ", " + startd2 + ", " + (startd2 + steps*incr) + ",) step size = " + incr + "\r\n");
            for (d1 = startd1; d1 < startd1 + steps*incr + TOL; d1 += incr)
                out.write(String.format("     %.1f,", d1));
            out.write("\r\n");
            for (d2 = startd2; d2 < startd2 + steps*incr + TOL; d2 += incr)
            {
                for (d1 = startd1; d1 < startd1 + steps*incr + TOL; d1 += incr)
                    out.write(String.format(" %.7f", calc_error(d1, d2)));
                out.write("\r\n");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("calc_array() save error = " + e);}
        scan_C0(startd1, startd2, incr, steps);
        scan_C1(startd1, startd2, incr, steps);
        scan_area(startd1, startd2, incr, steps);

        // scan <y> : requires branching because there is a crossover of the two solutions
//        if (c > 3.37158)    // this is where the discriminant has a tangential zero
//            scan_moment_y(startd1, startd2, incr, steps, e0, e1, e2, e3, e4, e5, e6, e7, false);
//        else
//            scan_moment_y(startd2, startd1, incr, steps, e0, e3, e6, e1, e4, e7, e2, e5, true); // flip d1 and d2

        // scan <y>/<1> : no branching required
//        scan_moment_y(startd1, startd2, incr, steps, 0, e1 - a1*e0/a0, e2, e3 + a2*e0/a0, e4 + a3*e0/a0, e5, e6, e7, false);

        // scan <x> : this re-uses scan_moment_y() to scan <x> instead
//        scan_moment_y(startd1, startd2, incr, steps, f0, f1, f2, f3, f4, f5, f6, f7, false);

        // scan <x>/<1> : no branching required
        scan_moment_y(startd1, startd2, incr, steps, 0, f1 - a1*f0/a0, f2, f3 + a2*f0/a0, f4 + a3*f0/a0, f5, f6, f7, false);
    }

    private static void scan_C0(double startd1, double startd2, double incr, double steps)
    {
        // generate data suitable for matplotlib
        // d2 as a fxn of d1, with C0 constant

        try
        {
            double d1, d2;
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\scan_C0.csv");
            out.write("d1, d2\r\n");
            for (d1 = startd1; d1 < startd1 + steps*incr + TOL; d1 += incr)
            {
                d2 = -l2 - (3*spiro_curvature(c)*d1*d1/2 - l1)*Math.sqrt(2);    // satisfy C0 constraint
                if (d2 >= startd2 && d2 < startd2 + steps*incr + TOL)
                    out.write(d1 + ", " + d2 + "\r\n");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("scan_C0() save error = " + e);}
    }

    private static void scan_C1(double startd1, double startd2, double incr, double steps)
    {
        // generate data suitable for matplotlib
        // d1 as a fxn of d2, with C1 constant

        try
        {
            double d1, d2;
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\scan_C1.csv");
            out.write("d1, d2\r\n");
            for (d2 = startd2; d2 < startd2 + steps*incr + TOL; d2 += incr)
            {
                d1 = -l1 - (3*spiro_curvature(-c)*d2*d2/2 - l2)*Math.sqrt(2);    // satisfy C1 constraint
                if (d1 >= startd1 && d1 < startd1 + steps*incr + TOL)
                    out.write(d1 + ", " + d2 + "\r\n");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("scan_C0() save error = " + e);}
    }

    private static void scan_area(double startd1, double startd2, double incr, double steps)
    {
        // generate data suitable for matplotlib
        // d2 as a fxn of d1, with area constant

        try
        {
            double d1, d2;
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\scan_area.csv");
            out.write("d1, d2\r\n");
            for (d1 = startd1; d1 < startd1 + steps*incr + TOL; d1 += incr)
            {
                d2 = (20*spiro_area()/3 - 2*d1*(l1 - l2/Math.sqrt(2)))/(2*(l2 - l1/Math.sqrt(2)) - d1/Math.sqrt(2)); // satisfy area constraint
                if (d2 >= startd2 && d2 < startd2 + steps*incr + TOL)
                    out.write(d1 + ", " + d2 + "\r\n");
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("scan_C0() save error = " + e);}
    }

    private static void scan_moment_y(double startd1, double startd2, double incr, double steps, double e0_m, double e1_m, double e2_m, double e3_m, double e4_m, double e5_m, double e6_m, double e7_m, boolean flip)
    {
        // generate data suitable for matplotlib
        // d2 as a fxn of d1, with anti-symmetric Bezier moment <y> constant

        try
        {
            double d1, d2, d1final = 0;
            double quada, quadb, quadc;                 // a*d2^2 + b*d2 + c = 0
            System.out.println("e1_m, e2_m, e3_m, e4_m, e5_m, e6_m, e7_m, 280<y>");
            System.out.println(e1_m + ", " + e2_m + ", " + e3_m + ", " + e4_m + ", " + e5_m + ", " + e6_m + ", " + e7_m + ", " + (-e0_m));
            FileWriter out = new FileWriter("\\APP\\MATPLOTLIB\\scan_moment_y.csv");
            out.write("d1, d2\r\n");
            for (d1 = startd1; d1 < startd1 + steps*incr + TOL; d1 += incr)
            {
                quada = e6_m + e7_m*d1;
                quadb = e3_m + e4_m*d1 + e5_m*d1*d1;
                quadc = e0_m + e1_m*d1 + e2_m*d1*d1;
//                System.out.printf("%f, %f, %f, %f, %f\n", d1, quada, quadb, quadc, quadb*quadb - 4*quada*quadc);
                if (quadb*quadb - 4*quada*quadc >= 0)   // satisfy anti-symmetric <y> constraint
                {
                    d2 = (-quadb + Math.sqrt(quadb*quadb - 4*quada*quadc))/2/quada;
                    if (d2 >= startd2 && d2 < startd2 + steps*incr + TOL)
                    {
                        if (flip)
                        {
                            d1final = d2;
                            out.write(d2 + ", " + d1 + "\r\n");
                        }
                        else
                        {
                            d1final = d1;
                            out.write(d1 + ", " + d2 + "\r\n");
                        }
                    }
                }
            }
            if (flip)
            {
                out.write(d1final + ", " + (startd1 + steps*incr) + "\r\n");
                out.write(startd2 + ", " + (startd1 + steps*incr) + "\r\n");
            }
            else
            {
                out.write(d1final + ", " + startd2 + "\r\n");
                out.write(startd1 + ", " + startd2 + "\r\n");
            }
            d1final = startd2 + steps*incr + 1;
            for (d1 = startd1; d1 < startd1 + steps*incr + TOL; d1 += incr)
            {
                quada = e6_m + e7_m*d1;
                quadb = e3_m + e4_m*d1 + e5_m*d1*d1;
                quadc = e0_m + e1_m*d1 + e2_m*d1*d1;
                if (quadb*quadb - 4*quada*quadc >= 0)   // satisfy anti-symmetric <y> constraint
                {
                    d2 = (-quadb - Math.sqrt(quadb*quadb - 4*quada*quadc))/2/quada;
                    if (d2 >= startd2 && d2 < startd2 + steps*incr + TOL)
                    {
                        if (d1final > startd2 + steps*incr)
                        {
                            if (flip)
                                out.write(startd2 + ", " + d1 + "\r\n");
                            else
                                out.write(startd1 + ", " + d2 + "\r\n");
                            d1final = d2;
                        }
                        if (flip)
                            out.write(d2 + ", " + d1 + "\r\n");
                        else
                            out.write(d1 + ", " + d2 + "\r\n");
                    }
                }
            }
            out.close();
        }
        catch (java.io.IOException e)
            {System.out.println("scan_moment_y() save error = " + e);}
    }

    private static double calc_error(double d1, double d2)
    {
        // calculate radial error as a function of theta, given Bezier arms d1, d2
        // assume a = 4b, spiro has 4-fold symmetry
        // perform scan over 45º, half a cycle

        double t_spiro, x_spiro, y_spiro, r_spiro;
        double t_bez, x_bez, y_bez, r_bez;
        double theta_spiro, theta_bez;                            // should be the same
        double x0 = l1, y0 = 0;
        double x1 = l1, y1 = d1;
        double x2 = l2/Math.sqrt(2) + d2/Math.sqrt(2), y2 = l2/Math.sqrt(2) - d2/Math.sqrt(2);
        double x3 = l2/Math.sqrt(2), y3 = l2/Math.sqrt(2);
        double tempa, tempb, tempc, tempd;
        double err = 0;

        System.out.printf("t_spiro,  x_spiro,  y_spiro,  r_spiro, theta_spiro,    t_bez,    x_bez,    y_bez,    r_bez,    theta_bez,    error\n");
        for (int i = 0; i <= 100; i++)
        {
            t_spiro = i/100.0;
            x_spiro = a_b*Math.cos(Math.PI*t_spiro/4) + c*Math.cos(-3*Math.PI*t_spiro/4);
            y_spiro = a_b*Math.sin(Math.PI*t_spiro/4) + c*Math.sin(-3*Math.PI*t_spiro/4);
            r_spiro = Math.sqrt(x_spiro*x_spiro + y_spiro*y_spiro);
            theta_spiro = Math.atan2(y_spiro, x_spiro);           // common to both spiro and bezier
            System.out.printf("%f, %f, %f, %f, %f, ", t_spiro, x_spiro, y_spiro, r_spiro, theta_spiro*180/Math.PI);
            tempa = (-x0 + 3*x1 - 3*x2 + x3)*y_spiro - (-y0 + 3*y1 - 3*y2 + y3)*x_spiro;
            tempb = (3*x0 - 6*x1 + 3*x2)*y_spiro - (3*y0 - 6*y1 + 3*y2)*x_spiro;
            tempc = (-3*x0 + 3*x1)*y_spiro - (-3*y0 + 3*y1)*x_spiro;
            tempd = x0*y_spiro - y0*x_spiro;
            t_bez = solve_cubic(tempb/tempa, tempc/tempa, tempd/tempa); // t as a function of theta_spiro
//            t_bez = t_spiro;                              // for testing only
            x_bez = x0*(1-t_bez)*(1-t_bez)*(1-t_bez) + 3*x1*t_bez*(1-t_bez)*(1-t_bez) + 3*x2*t_bez*t_bez*(1-t_bez) + x3*t_bez*t_bez*t_bez;
            y_bez = y0*(1-t_bez)*(1-t_bez)*(1-t_bez) + 3*y1*t_bez*(1-t_bez)*(1-t_bez) + 3*y2*t_bez*t_bez*(1-t_bez) + y3*t_bez*t_bez*t_bez;
            r_bez = Math.sqrt(x_bez*x_bez + y_bez*y_bez);
            err += (r_bez - r_spiro)*(r_bez - r_spiro);
            theta_bez = Math.atan2(y_bez, x_bez);           // just a double check
            System.out.printf("%f, %f, %f, %f, %f, %.7f\n", t_bez, x_bez, y_bez, r_bez, theta_bez*180/Math.PI, r_bez/r_spiro - 1);
        }
        return Math.sqrt(err/100)/a_b;
    }

    private static double bezier_C0(double d1, double d2)       // Bezier curvature at start
    {
        return 2*(l1 - l2/Math.sqrt(2) - d2/Math.sqrt(2))/3/d1/d1;
    }

    private static double bezier_C1(double d1, double d2)       // Bezier curvature at end
    {
        return 2*(2*l2/Math.sqrt(2) - d1 - l1)/Math.sqrt(2)/3/d2/d2;
    }

    private static double bezier_area(double d1, double d2)     // Bezier area: <1>
    {
        return 3*(2*d1*(l1 - l2/Math.sqrt(2)) + 2*d2*(l2 - l1/Math.sqrt(2)) - d1*d2/Math.sqrt(2))/20;
    }

    private static double bezier_moment_x(double d1, double d2) // Bezier moment: <x>
    {
        // express 280<x> =  f1*d1 + f2*d1*d1
        //                + (f3 + f4*d1 + f5*d1*d1)*d2
        //                + (f6 + f7*d1)*d2*d2

        return (f1*d1 + f2*d1*d1 + (f3 + f4*d1 + f5*d1*d1)*d2 + (f6 + f7*d1)*d2*d2)/280;
    }

    private static double bezier_moment_y(double d1, double d2) // Bezier moment: <y>
    {
        // express 280<y> =  e1*d1 + e2*d1*d1
        //                + (e3 + e4*d1 + e5*d1*d1)*d2
        //                + (e6 + e7*d1)*d2*d2

        return (e1*d1 + e2*d1*d1 + (e3 + e4*d1 + e5*d1*d1)*d2 + (e6 + e7*d1)*d2*d2)/280;
    }

    private static double spiro_curvature(double myc)
    {
        return (1 + 9*myc/a_b)/(1 - 3*myc/a_b)/(1 - 3*myc/a_b)/a_b;
    }

    private static double spiro_area()
    {
        return a_b*a_b*(Math.PI/2 - Math.sqrt(2))/4                     // <1>
               - c*c*(3*Math.PI/2 - Math.sqrt(2))/4;
    }

    private static double spiro_moment(String dir)
    {
        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;
        final double h = Math.sqrt(4 - 2*Math.sqrt(2));

        if (dir.equals("x"))                                            // <x>
            return h*(a_b*a_b*a_b*A + a_b*c*c*C)/12;
        else                                                            // <y>
            return -h*(Math.sqrt(2) + 1)*(a_b*a_b*c*B + c*c*c*D)/12;
    }

    public static double solve_cubic(double p, double q, double r)
    {
        // see Math CRC book, page 392

        double cua = (3*q - p*p)/3;
        double cub = (2*p*p*p - 9*p*q + 27*r)/27;
        double cud = cub*cub/4 + cua*cua*cua/27;
        double test;

        //System.out.println("\ncubic p,q,r = " + p + ", " + q + ", " + r);
        //System.out.println("\ncubic a,b,d = " + cua + ", " + cub + ", " + cud);
        if (cud < TOL)
        {
            double myphi = Math.acos(-Math.signum(cub));                // coalesced roots
            if (cud < 0)
                myphi = Math.acos(-cub/2/Math.sqrt(-cua*cua*cua/27));   // distinct roots
//            System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3));
            test = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3) - p/3;
            if (test > -TOL && test < 1 + TOL) return test;
            test = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3;
            if (test > -TOL && test < 1 + TOL) return test;
            test = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3;
            if (test > -TOL && test < 1 + TOL) return test;
            return Double.NaN;
        }
        else
        {
//            System.out.println("1 cubic d > 0 : " + (Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3));
            return Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3;
        }
    }

    private static void check_neighbours(double d1, double d2)
    {
        // check for the existence of a local minimum by incrementing d1 and d2 by +/- 1
        // choose d1, d2 to be a local minimum at a fixed area

        System.out.printf("\nchecking neighbours at %f, %f, %f, %f, %f, %f, %f, %.7f\n", d1, d2, bezier_C0(d1, d2), bezier_C1(d1, d2), bezier_area(d1, d2), bezier_moment_x(d1, d2), bezier_moment_y(d1, d2), calc_error(d1, d2));
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
                System.out.printf("%f, %f, %.7f\n", d1 + i, d2 + j, calc_error(d1 + i, d2 + j));
    }
}
