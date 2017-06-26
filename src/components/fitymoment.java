
package components;

// calculate arm length of a Bezier curve (quartic in d1)
// given that we know the spiro area and the first moment <y>
// see middle of CofM book, Sep 2016
// use the symmetric Bezier construction for the asymmetric spiro, with arc 45º

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\fitymoment.java

public class fitymoment
{
    private static final double a_b = 180;          // spiro 'a - b'
    private static final double theta = Math.PI/8;  // 22.5º
    private static double c;                        // spiro 'c'
    private static double l1;                       // distance to start point (-22.5º)
    private static double l2;                       // distance to end   point ( 22.5º)
    private static double phi;                      // angle perpendicular to line joining p0 and p3
    private static double a0, a1, a2, a3;
    private static double e1, e2, e3, e4, e5, e6, e7;
    private static double f1, f2, f3, f4, f5, f6, f7;

    public static void main (String[] args)
    {
        c = -3.5298;         // 3.5297135;
        l1 = a_b + c;
        l2 = a_b - c;
//        phi = Math.atan(l1*Math.sqrt(2)/l2 - 1);    // general transform that rotates with the object baseline
        phi = theta;                                // simplified transform at a fixed 22.5º

        // express d2 = (a0 + a1*d1)/(a2 + a3*d1)
        a0 = 20*spiro_area()/3;
        a1 = -2*(l1 - l2*Math.cos(2*theta));
        a2 =  2*(l2 - l1*Math.cos(2*theta));
        a3 = -Math.sin(2*theta);

        // express 280<x> =  f1*d1 + f2*d1*d1
        //                + (f3 + f4*d1 + f5*d1*d1)*d2
        //                + (f6 + f7*d1)*d2*d2
        f1 = (l1 - l2*Math.cos(2*theta))*(50*l1 + 34*l2)*Math.cos(theta);
        f2 = 15*(l1 - l2*Math.cos(2*theta))*Math.sin(theta);
        f3 = (l2 - l1*Math.cos(2*theta))*(34*l1 + 50*l2)*Math.cos(theta);
        f4 = (l1 + l2)*(-21*Math.sin(2*theta)*Math.cos(theta) + 12*(1 - Math.cos(2*theta))*Math.sin(theta));
        f5 = -9*Math.sin(2*theta)*Math.sin(theta);
        f6 = 15*(l2 - l1*Math.cos(2*theta))*Math.sin(theta);
        f7 = -9*Math.sin(2*theta)*Math.sin(theta);

        // express 280<y> =  e1*d1 + e2*d1*d1
        //                + (e3 + e4*d1 + e5*d1*d1)*d2
        //                + (e6 + e7*d1)*d2*d2
        e1 = (l1 - l2*Math.cos(2*theta))*(-50*l1 + 34*l2)*Math.sin(theta);
        e2 = 15*(l1 - l2*Math.cos(2*theta))*Math.cos(theta);
        e3 = (l2 - l1*Math.cos(2*theta))*(-34*l1 + 50*l2)*Math.sin(theta);
        e4 = (l1 - l2)*(21*Math.sin(2*theta)*Math.sin(theta) - 12*(1 + Math.cos(2*theta))*Math.cos(theta));
        e5 = -9*Math.sin(2*theta)*Math.cos(theta);
        e6 = -15*(l2 - l1*Math.cos(2*theta))*Math.cos(theta);
        e7 = 9*Math.sin(2*theta)*Math.cos(theta);

        System.out.printf("a-b, c, phi, <x>, <y> = %f, %f, %f, %f, %f\n", a_b, c, phi*180/Math.PI, spiro_moment_x(), spiro_moment_y());
//        setup_quartic_area_y();
//        setup_quartic_area_x();
//        setup_quartic_extremum_y();
//        setup_quintic_x_y();
//        setup_quartic_cofmx_cofmy();
//        scan_spiro_moment_y();              // for testing only fix fix
//        scan_discriminant_vs_c();
        scan_d();
    }

    private static void setup_quartic_area_y()
    {
        // this will fit the area and the tilted <y> moment (leads to S-shaped curve)

        final double b1 = -Math.sin(phi - theta)*f1 + Math.cos(phi - theta)*e1;
        final double b2 = -Math.sin(phi - theta)*f2 + Math.cos(phi - theta)*e2;
        final double b3 = -Math.sin(phi - theta)*f3 + Math.cos(phi - theta)*e3;
        final double b4 = -Math.sin(phi - theta)*f4 + Math.cos(phi - theta)*e4;
        final double b5 = -Math.sin(phi - theta)*f5 + Math.cos(phi - theta)*e5;
        final double b6 = -Math.sin(phi - theta)*f6 + Math.cos(phi - theta)*e6;
        final double b7 = -Math.sin(phi - theta)*f7 + Math.cos(phi - theta)*e7;

        double d1, d2, beziermomenty;
//        System.out.printf("spiro  <y'> = %f\n", spiro_moment_y());
//        double testd1 = (0.265202 + .01)*a_b;
//        double testd2 = (0.265202 - .01)*a_b;
//        System.out.printf("bezier <y'> = %f\n", (b1*testd1 + b2*testd1*testd1 + (b3 + b4*testd1 + b5*testd1*testd1)*testd2 + (b6 + b7*testd1)*testd2*testd2)/280);    // test code only
        d1 = solve_quartic(a1*a3*b5 + a3*a3*b2,
             a1*a1*b7 + (a1*a2 + a0*a3)*b5 + a1*a3*b4 + 2*a2*a3*b2 + a3*a3*b1,
             a1*a1*b6 + 2*a0*a1*b7 + a0*a2*b5 + (a1*a2 + a0*a3)*b4
           + a1*a3*b3 + a2*a2*b2 + 2*a2*a3*b1 - a3*a3*280*spiro_moment_y(),
             2*a0*a1*b6 + a0*a0*b7 + a0*a2*b4 + (a1*a2 + a0*a3)*b3
           + a2*a2*b1 - 2*a2*a3*280*spiro_moment_y(),
             a0*a0*b6 + a0*a2*b3 - a2*a2*280*spiro_moment_y(),
             true);
        d2 = (a0 + a1*d1)/(a2 + a3*d1);
        beziermomenty = b1*d1 + b2*d1*d1 + (b3 + b4*d1 + b5*d1*d1)*d2 + (b6 + b7*d1)*d2*d2;
        System.out.println("d1, d2 = " + d1 + ", " + d2);
        System.out.println("area   = " + spiro_area() + ", " + bezier_area(d1, d2));
        System.out.println("<y>    = " + spiro_moment_y() + ", " + beziermomenty/280);
    }

    private static void setup_quartic_area_x()
    {
        // this will fit the area and the tilted <x> moment (result is not interesting, too localized)

        final double b1 = f1;   // use only the contribution from <x> tilted by 22.5º
        final double b2 = f2;
        final double b3 = f3;
        final double b4 = f4;
        final double b5 = f5;
        final double b6 = f6;
        final double b7 = f7;

        double d1, d2, beziermomentx;
//        System.out.printf("spiro  <y'> = %f\n", spiro_moment_y());
//        double testd1 = (0.265202 + .01)*a_b;
//        double testd2 = (0.265202 - .01)*a_b;
//        System.out.printf("bezier <y'> = %f\n", (b1*testd1 + b2*testd1*testd1 + (b3 + b4*testd1 + b5*testd1*testd1)*testd2 + (b6 + b7*testd1)*testd2*testd2)/280);    // test code only
        d1 = solve_quartic(a1*a3*b5 + a3*a3*b2,
             a1*a1*b7 + (a1*a2 + a0*a3)*b5 + a1*a3*b4 + 2*a2*a3*b2 + a3*a3*b1,
             a1*a1*b6 + 2*a0*a1*b7 + a0*a2*b5 + (a1*a2 + a0*a3)*b4
           + a1*a3*b3 + a2*a2*b2 + 2*a2*a3*b1 - a3*a3*280*spiro_moment_x(),
             2*a0*a1*b6 + a0*a0*b7 + a0*a2*b4 + (a1*a2 + a0*a3)*b3
           + a2*a2*b1 - 2*a2*a3*280*spiro_moment_x(),
             a0*a0*b6 + a0*a2*b3 - a2*a2*280*spiro_moment_x(),
             true);
        d2 = (a0 + a1*d1)/(a2 + a3*d1);
        beziermomentx = b1*d1 + b2*d1*d1 + (b3 + b4*d1 + b5*d1*d1)*d2 + (b6 + b7*d1)*d2*d2;
        System.out.println("d1, d2 = " + d1 + ", " + d2);
        System.out.println("area   = " + spiro_area() + ", " + bezier_area(d1, d2));
        System.out.println("<x>    = " + spiro_moment_x() + ", " + beziermomentx/280);
    }

    private static void setup_quartic_extremum_y()
    {
        // this will calculate extrema of the tilted <y> moment, holding area fixed
        // this should produce a closed oblong shape about the center where c = 0

        double d1, d2;
//        System.out.printf("spiro  <y'> = %f\n", spiro_moment_y());
        d1 = solve_quartic(2*a3*a3*a3*e2 + 2*a1*a3*a3*e5,
             a3*a3*a3*e1 + 6*a2*a3*a3*e2 + a1*a3*a3*e4
           + 2*a0*a3*a3*e5 + 4*a1*a2*a3*e5 + a1*a1*a3*e7
           + (a1*a2 - a0*a3)*a3*e5,
             3*a2*a3*a3*e1 + 6*a2*a2*a3*e2 + a0*a3*a3*e4
           + 2*a1*a2*a3*e4 + 4*a0*a2*a3*e5 + 2*a1*a2*a2*e5
           + a1*a1*a2*e7 + 2*a0*a1*a3*e7
           + (a1*a2 - a0*a3)*(a3*e4 + 2*a1*e7 + a2*e5),
             3*a2*a2*a3*e1 + 2*a2*a2*a2*e2 + 2*a0*a2*a3*e4
           + a1*a2*a2*e4 + 2*a0*a2*a2*e5 + 2*a0*a1*a2*e7
           + a0*a0*a3*e7 + (a1*a2 - a0*a3)*(a3*e3 + 2*a1*e6 + a2*e4 + 2*a0*e7),
             a2*a2*a2*e1 + a0*a2*a2*e4 + a0*a0*a2*e7
           + (a1*a2 - a0*a3)*(a2*e3 + 2*a0*e6),
           true);
        d2 = (a0 + a1*d1)/(a2 + a3*d1);
        System.out.println("quartic_extremum_y c d1 d2 = ," + c + ", " + d1 + ", " + d2);
        System.out.println("area   = " + spiro_area() + ", " + bezier_area(d1, d2));
    }

    private static void setup_quintic_x_y()
    {
        // this will fit the tilted <x> and the tilted <y> moment
        // results are not interesting, essentially the same as setup_quartic_area_y()

        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;
        final double momxprime = (a_b*a_b*a_b*A + a_b*c*c*C)/12/Math.cos(theta);    // expressed at angle 22.5º
        final double momyprime = -(Math.sqrt(2) + 1)*(a_b*a_b*c*B + c*c*c*D)/12/Math.cos(theta);

        // express d1   = (u0 + u1*d2 + u2*d2*d2)/(u3 + u4*d2 + u5*d2*d2) - transformed to angle 0º
        final double u0 = 280*(momxprime*Math.cos(theta) - momyprime*Math.sin(theta));    // same as momx
        final double u1 = -(f3*Math.cos(theta) - e3*Math.sin(theta));
        final double u2 = -(f6*Math.cos(theta) - e6*Math.sin(theta));
        final double u3 = f1*Math.cos(theta) - e1*Math.sin(theta);
        final double u4 = f4*Math.cos(theta) - e4*Math.sin(theta);
        final double u5 = f7*Math.cos(theta) - e7*Math.sin(theta);
//        System.out.println("trans1 f2 = " + (f2*Math.cos(theta) - e2*Math.sin(theta)));
//        System.out.println("trans1 f5 = " + (f5*Math.cos(theta) - e5*Math.sin(theta)));

        // express d2   = (b0 + b1*d1 + b2*d1*d1)/(b3 + b4*d1 + b5*d1*d1) - transformed to angle 45º
        final double b0 = 280*(momxprime*Math.cos(theta) + momyprime*Math.sin(theta));
        final double b1 = -(f1*Math.cos(theta) + e1*Math.sin(theta));
        final double b2 = -(f2*Math.cos(theta) + e2*Math.sin(theta));
        final double b3 = f3*Math.cos(theta) + e3*Math.sin(theta);
        final double b4 = f4*Math.cos(theta) + e4*Math.sin(theta);
        final double b5 = f5*Math.cos(theta) + e5*Math.sin(theta);
        double d1, d2;
//        System.out.println("trans2 f6 = " + (f6*Math.cos(theta) + e6*Math.sin(theta)));
//        System.out.println("trans2 f7 = " + (f7*Math.cos(theta) + e7*Math.sin(theta)));
        d1 = solve_quintic(u3*b5*b5 + u4*b2*b5 + u5*b2*b2,
                           2*u3*b4*b5 - u0*b5*b5 + u4*(b1*b5 + b2*b4) - u1*b2*b5 + 2*u5*b1*b2 - u2*b2*b2,
                           u3*(2*b3*b5 + b4*b4) - 2*u0*b4*b5 + u4*(b0*b5 + b1*b4 + b2*b3)
                         - u1*(b1*b5 + b2*b4) + u5*(2*b0*b2 + b1*b1) - 2*u2*b1*b2,
                           2*u3*b3*b4 - u0*(2*b3*b5 + b4*b4) + u4*(b0*b4 + b1*b3)
                         - u1*(b0*b5 + b1*b4 + b2*b3) + 2*u5*b0*b1 - u2*(2*b0*b2 + b1*b1),
                           u3*b3*b3 - 2*u0*b3*b4 + u4*b0*b3 - u1*(b0*b4 + b1*b3) + u5*b0*b0 - 2*u2*b0*b1,
                           -u0*b3*b3 - u1*b0*b3 - u2*b0*b0);
        d2 = (b0 + b1*d1 + b2*d1*d1)/(b3 + b4*d1 + b5*d1*d1);
        System.out.println("solve_quintic c,  d1,  d2 = " + c + ", " + d1 + ", " + d2);
        System.out.println("spiro       <1>, <x>, <y> = " + spiro_area() + ", " + spiro_moment_x() + ", " + spiro_moment_y());
        System.out.println("Bezier      <1>, <x>, <y> = " + bezier_area(d1, d2) + ", " + 
                                                            (f1*d1 + f2*d1*d1 + (f3 + f4*d1 + f5*d1*d1)*d2 + (f6 + f7*d1)*d2*d2)/280 + ", " +
                                                            (e1*d1 + e2*d1*d1 + (e3 + e4*d1 + e5*d1*d1)*d2 + (e6 + e7*d1)*d2*d2)/280);
    }

    private static double solve_quintic(double lead, double qua, double qub, double quc, double qud, double que)
    {
        // this routine is not really a solution of the quintic
        // it is a customizable scan followed by a linear interpolation
        qua /= lead;
        qub /= lead;
        quc /= lead;
        qud /= lead;
        que /= lead;
        System.out.println("\nquintic a,b,c,d,e = " + (float)qua + ", " + (float)qub + ", " + (float)quc + ", " + (float)qud + ", " + (float)que);
        double start = 3;
        double incr = .01;
        double d, res0 = 0, res1;
        double retval = 0;                              // record interpolated root
        for (double i = 0; i <= 100; i++)
        {
            d = start + i*incr;
            res1 = d*d*d*d*d + qua*d*d*d*d + qub*d*d*d + quc*d*d + qud*d + que;
            System.out.println(d + ", " + res1);
            if ((res1 >= 0 && res0 < 0) || (res1 <= 0 && res0 > 0))
                retval = d - incr*res1/(res1 - res0);   // linear interpolation
            res0 = res1;
        }
        return retval;
    }

    private static void setup_quartic_cofmx_cofmy()
    {
        // this will fit the tilted <x>/<1> and the tilted <y>/<1> center of masses, Oct 22
        // results are much better than the setup_quartic_area_y fit, this is the best yet

        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;
        final double momxprime = (a_b*a_b*a_b*A + a_b*c*c*C)/12/Math.cos(theta);    // expressed at angle 22.5º
        final double momyprime = -(Math.sqrt(2) + 1)*(a_b*a_b*c*B + c*c*c*D)/12/Math.cos(theta);

        // express d1   = (u1*d2 + u2*d2*d2)/(u3 + u4*d2 + u5*d2*d2) - transformed to angle 0º
        final double u0 = 280*(momxprime*Math.cos(theta) - momyprime*Math.sin(theta));    // same as momx
        final double u1 = -(f3*Math.cos(theta) - e3*Math.sin(theta) - a2*u0/a0);
        final double u2 = -(f6*Math.cos(theta) - e6*Math.sin(theta));
        final double u3 = f1*Math.cos(theta) - e1*Math.sin(theta) + a1*u0/a0;
        final double u4 = f4*Math.cos(theta) - e4*Math.sin(theta) - a3*u0/a0;
        final double u5 = f7*Math.cos(theta) - e7*Math.sin(theta);
//        System.out.println("setup_quartic_cofmx_cofmy = " + momxprime + ", " + momyprime);
//        System.out.println("trans1 f2 = " + (f2*Math.cos(theta) - e2*Math.sin(theta)));
//        System.out.println("trans1 f5 = " + (f5*Math.cos(theta) - e5*Math.sin(theta)));

        // express d2   = (b1*d1 + b2*d1*d1)/(b3 + b4*d1 + b5*d1*d1) - transformed to angle 45º
        final double b0 = 280*(momxprime*Math.cos(theta) + momyprime*Math.sin(theta));
        final double b1 = -(f1*Math.cos(theta) + e1*Math.sin(theta) + a1*b0/a0);
        final double b2 = -(f2*Math.cos(theta) + e2*Math.sin(theta));
        final double b3 = f3*Math.cos(theta) + e3*Math.sin(theta) - a2*b0/a0;
        final double b4 = f4*Math.cos(theta) + e4*Math.sin(theta) - a3*b0/a0;
        final double b5 = f5*Math.cos(theta) + e5*Math.sin(theta);
        double d1, d2;
//        System.out.println("trans2 f6 = " + (f6*Math.cos(theta) + e6*Math.sin(theta)));
//        System.out.println("trans2 f7 = " + (f7*Math.cos(theta) + e7*Math.sin(theta)));

        d1 = solve_quartic(u3*b5*b5 + u4*b2*b5 + u5*b2*b2,
                           2*u3*b4*b5 + u4*(b1*b5 + b2*b4) - u1*b2*b5 + 2*u5*b1*b2 - u2*b2*b2,
                           u3*(2*b3*b5 + b4*b4) + u4*(b1*b4 + b2*b3)
                         - u1*(b1*b5 + b2*b4) + u5*b1*b1 - 2*u2*b1*b2,
                           2*u3*b3*b4 + u4*b1*b3
                         - u1*(b1*b4 + b2*b3) - u2*b1*b1,
                           u3*b3*b3 - u1*b1*b3,
                           true);
        d2 = (b1*d1 + b2*d1*d1)/(b3 + b4*d1 + b5*d1*d1);
        System.out.println("solve_quartic c,  d1,  d2 = " + c + ", " + d1 + ", " + d2 + ", " + bezier_area(d1, d2)/spiro_area());
        System.out.println("spiro       <1>, <x>, <y> = " + spiro_area() + ", " + spiro_moment_x() + ", " + spiro_moment_y());
        System.out.println("Bezier      <1>, <x>, <y> = " + bezier_area(d1, d2) + ", " +
                                                            (f1*d1 + f2*d1*d1 + (f3 + f4*d1 + f5*d1*d1)*d2 + (f6 + f7*d1)*d2*d2)/280 + ", " +
                                                            (e1*d1 + e2*d1*d1 + (e3 + e4*d1 + e5*d1*d1)*d2 + (e6 + e7*d1)*d2*d2)/280);
    }

    private static double spiro_area()
    {
        return a_b*a_b*(Math.PI/2 - Math.sqrt(2))/4
               - c*c*(3*Math.PI/2 - Math.sqrt(2))/4;
    }

    private static double spiro_moment_y()
    {
        // this is a general transform to get the projection onto y using the angle phi
        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;

        final double momx = (a_b*a_b*a_b*A + a_b*a_b*c*B + a_b*c*c*C + c*c*c*D)/12;
        final double momy = (a_b*a_b*a_b*A*(Math.sqrt(2) - 1) - a_b*a_b*c*B/(Math.sqrt(2) - 1) + a_b*c*c*C*(Math.sqrt(2) - 1) - c*c*c*D/(Math.sqrt(2) - 1))/12;
//        System.out.printf("moment_y A, B, C, D, momx, momy = %f, %f, %f, %f, %f, %f\n", A, B, C, D, momx, momy);
        return -momx*Math.sin(phi) + momy*Math.cos(phi);
    }

    private static double spiro_moment_x()
    {
        // this is a specific transform to get the projection onto x using only 22.5º
        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;

        final double momx = (a_b*a_b*a_b*A + a_b*a_b*c*B + a_b*c*c*C + c*c*c*D)/12;
        final double momy = (a_b*a_b*a_b*A*(Math.sqrt(2) - 1) - a_b*a_b*c*B/(Math.sqrt(2) - 1) + a_b*c*c*C*(Math.sqrt(2) - 1) - c*c*c*D/(Math.sqrt(2) - 1))/12;
//        System.out.printf("moment_x A, B, C, D, momx, momy = %f, %f, %f, %f, %f, %f\n", A, B, C, D, momx, momy);
        return momx*Math.cos(theta) + momy*Math.sin(theta);
    }

    private static double bezier_area(double d1, double d2)
    {
        return 3*(2*d1*(l1 - l2*Math.cos(2*theta)) + 2*d2*(l2 - l1*Math.cos(2*theta)) - d1*d2*Math.sin(2*theta))/20;
    }

    public static double solve_quartic(double lead, double qua, double qub, double quc, double qud, boolean sgn)
    {
        double sol, R, D, E;

        qua /= lead;
        qub /= lead;
        quc /= lead;
        qud /= lead;
        System.out.println("\nquartic      a,b,c,d = " + (float)qua + ", " + (float)qub + ", " + (float)quc + ", " + (float)qud);
        sol = solve_cubic(-qub, qua*quc - 4*qud, -qua*qua*qud + 4*qub*qud - quc*quc);
        R = Math.sqrt(qua*qua/4 - qub + sol);
        D = Math.sqrt(3*qua*qua/4 - R*R - 2*qub + (4*qua*qub - 8*quc - qua*qua*qua)/4/R);
        E = Math.sqrt(3*qua*qua/4 - R*R - 2*qub - (4*qua*qub - 8*quc - qua*qua*qua)/4/R);
        System.out.println("cubic sol = " + sol + ", " + R + ", " + D + ", " + E);
        System.out.print("roots = , " + (-qua/4 + R/2 + D/2));
        System.out.print(", " + (-qua/4 + R/2 - D/2));
        System.out.print(", " + (-qua/4 - R/2 + E/2));
        System.out.println(", " + (-qua/4 - R/2 - E/2));
        if (!Double.isNaN(E) && sgn)
        {
            System.out.println("using root 4 = " + (-qua/4 - R/2 - E/2));
            return (-qua/4 - R/2 - E/2);
        }
        if (!Double.isNaN(D) && sgn && false)
        {
            System.out.println("using root 1 = " + (-qua/4 + R/2 + D/2));
            return (-qua/4 + R/2 + D/2);
        }
        if (!Double.isNaN(D) && sgn)
        {
            System.out.println("using root 2 = " + (-qua/4 + R/2 - D/2));
            return (-qua/4 + R/2 - D/2);
        }
        if (!Double.isNaN(E))
        {
            System.out.println("using root 3 = " + (-qua/4 - R/2 + E/2));
            return (-qua/4 - R/2 + E/2);
        }
        System.out.println("general quartic : Bad solution = " + R + ", " + D + ", " + E);
        return Double.NaN;
    }

    private static double solve_cubic(double p, double q, double r)
    {
        // see Math CRC book, page 392

        double cua = (3*q - p*p)/3;
        double cub = (2*p*p*p - 9*p*q + 27*r)/27;
        double cud = cub*cub/4 + cua*cua*cua/27;

//        System.out.println("\ncubic p,q,r = " + p + ", " + q + ", " + r);
//        System.out.println("\ncubic a,b,d = " + cua + ", " + cub + ", " + cud);
        if (cud < 0)
        {
            double myphi = Math.acos(-cub/2/Math.sqrt(-cua*cua*cua/27));
//            System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3));
            return 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3;
        }
        else
        {
//            System.out.println("1 cubic d > 0 : " + (Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3));
            return Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3;
        }
    }

    private static void scan_spiro_moment_y()
    {
        // calculate moments as a function of c, compare to midpoint of baseline

        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;
        double myc;
        double myphi;
        double area, momx, momy, midx, midy;

        System.out.println("\nscan moments");
        System.out.println("c, theta, momx, momy, area, midx, midy");
//        System.out.printf("A, B, C, D = %f, %f, %f, %f\n", A, B, C, D);

        for (myc = -2; myc <= 20; myc++)
        {
            myphi = Math.atan((a_b + myc)*Math.sqrt(2)/(a_b - myc) - 1);
            momx = (a_b*a_b*a_b*A + a_b*a_b*myc*B + a_b*myc*myc*C + myc*myc*myc*D)/12;
            momy = (a_b*a_b*a_b*A*(Math.sqrt(2) - 1) - a_b*a_b*myc*B/(Math.sqrt(2) - 1) + a_b*myc*myc*C*(Math.sqrt(2) - 1) - myc*myc*myc*D/(Math.sqrt(2) - 1))/12;
            area = a_b*a_b*(Math.PI/2 - Math.sqrt(2))/4 - myc*myc*(3*Math.PI/2 - Math.sqrt(2))/4;
            midx = (a_b + myc + (a_b - myc)/Math.sqrt(2))/2;
            midy = (a_b - myc)/Math.sqrt(2)/2;
            System.out.printf("%f, %f, %f, %f, %f, %f, %f\n", myc, myphi*180/Math.PI,
                    momx*Math.cos(myphi) + momy*Math.sin(myphi),
                   -momx*Math.sin(myphi) + momy*Math.cos(myphi),
                    area,
                    midx*Math.cos(myphi) + midy*Math.sin(myphi),
                   -midx*Math.sin(myphi) + midy*Math.cos(myphi));
        }
    }

    private static void scan_discriminant_vs_c()
    {
        // try to determine where two solutions of d2 vs d1 coalesce
        // when fitting the anti-symmetric <y> moment (Spiro2SVG Book2, p.10)
        // by scanning the discriminant of the quartic equation for d1

        double startc = -3.37158;
        double incr = .000001;
        int steps = 10;
        double lead, qua, qub, quc, qud;

        for (c = startc; c < startc + steps*incr + .00001; c += incr)
        {
            l1 = a_b + c;
            l2 = a_b - c;

            e1 = (l1 - l2*Math.cos(2*theta))*(-50*l1 + 34*l2)*Math.sin(theta);
            e2 = 15*(l1 - l2*Math.cos(2*theta))*Math.cos(theta);
            e3 = (l2 - l1*Math.cos(2*theta))*(-34*l1 + 50*l2)*Math.sin(theta);
            e4 = (l1 - l2)*(21*Math.sin(2*theta)*Math.sin(theta) - 12*(1 + Math.cos(2*theta))*Math.cos(theta));
            e5 = -9*Math.sin(2*theta)*Math.cos(theta);
            e6 = -15*(l2 - l1*Math.cos(2*theta))*Math.cos(theta);
            e7 = 9*Math.sin(2*theta)*Math.cos(theta);

            lead = e5*e5;
            qua = 2*e4*e5 - 4*e2*e7;
            qub = 2*e3*e5 + e4*e4 - 4*e1*e7 - 4*e2*e6;
            quc = 2*e3*e4 + 1120*spiro_moment_y()*e7 - 4*e1*e6;
            qud = e3*e3 + 1120*spiro_moment_y()*e6;
            System.out.printf("a-b, c, <y> = %f, %f, %f, ", a_b, c, 280*spiro_moment_y());
            qua /= lead;
            qub /= lead;
            quc /= lead;
            qud /= lead;
            //System.out.println("quartic      a,b,c,d =, " + (float)qua + ", " + (float)qub + ", " + (float)quc + ", " + (float)qud);
            solve_cubic(-qub, qua*quc - 4*qud, -qua*qua*qud + 4*qub*qud - quc*quc);
        }
    }

    private static void scan_d()
    {
        double d1 = 0.5519;
        double d2 = 0.5520;
        System.out.println("\nscan_d from " + d1 + " to " + d2);
        for (int i = 0; i <= 20; i++)
            scan_r2_vs_t(d1 + i*(d2 - d1)/20);
    }

    private static void scan_r2_vs_t(double d)
    {
        // calculate Bezier r**2 for t = (0,1) for a quarter circle
        // as function of arm length d
        // see Spiro2SVG Book 2, p.45

        double r2;
        double sum = 0;
//        System.out.println("\nscan_r2_vs_t");
        for (double t = 0; t < 1.001; t += 0.005)
        {
            r2 = (1 - t)*(1 - t)*(1 - t)*(1 - t)*(1 + 2*t)*(1 + 2*t)
               + t*t*t*t*(3 - 2*t)*(3 - 2*t)
               + 6*d*t*t*(1 - t)*(1 - t)*(1 - t)*(1 + 2*t)
               + 6*d*t*t*t*(1 - t)*(1 - t)*(3 - 2*t)
               + 9*d*d*t*t*t*t*(1 - t)*(1 - t)
               + 9*d*d*t*t*(1 - t)*(1 - t)*(1 - t)*(1 - t);
            sum += (1 - r2)*(1 - r2);
//            System.out.println(t + ", " + r2);
        }
        System.out.printf("d sum = , %f, %.7f\n", d, Math.sqrt(sum/200));
    }
}
