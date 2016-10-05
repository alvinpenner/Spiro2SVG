
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

    public static void main (String[] args)
    {
        c = 7;              // 3.53;
        l1 = a_b + c;
        l2 = a_b - c;
//        phi = Math.atan(l1*Math.sqrt(2)/l2 - 1);    // general transform that rotates with the object baseline
        phi = theta;                                // simplified transform at a fixed 22.5º
        System.out.printf("a-b, c, phi = %f, %f, %f\n", a_b, c, phi*180/Math.PI);
        setup_quartic();
//        scan_spiro_moment_y();              // for testing only fix fix
    }

    private static void setup_quartic()
    {
        // express d2     = (a0 + a1*d1)/(a2 + a3*d1)
        // express 280<x> =  f1*d1 + f2*d1*d1
        //                + (f3 + f4*d1 + f5*d1*d1)*d2
        //                + (f6 + f7*d1)*d2*d2
        // express 280<y> =  e1*d1 + e2*d1*d1
        //                + (e3 + e4*d1 + e5*d1*d1)*d2
        //                + (e6 + e7*d1)*d2*d2

        final double a0 = 20*spiro_area()/3;
        final double a1 = -2*(l1 - l2*Math.cos(2*theta));
        final double a2 =  2*(l2 - l1*Math.cos(2*theta));
        final double a3 = -Math.sin(2*theta);

        final double f1 = (l1 - l2*Math.cos(2*theta))*(50*l1 + 34*l2)*Math.cos(theta);
        final double f2 = 15*(l1 - l2*Math.cos(2*theta))*Math.sin(theta);
        final double f3 = (l2 - l1*Math.cos(2*theta))*(34*l1 + 50*l2)*Math.cos(theta);
        final double f4 = (l1 + l2)*(-21*Math.sin(2*theta)*Math.cos(theta) + 12*(1 - Math.cos(2*theta))*Math.sin(theta));
        final double f5 = -9*Math.sin(2*theta)*Math.sin(theta);
        final double f6 = 15*(l2 - l1*Math.cos(2*theta))*Math.sin(theta);
        final double f7 = -9*Math.sin(2*theta)*Math.sin(theta);

        final double e1 = (l1 - l2*Math.cos(2*theta))*(-50*l1 + 34*l2)*Math.sin(theta);
        final double e2 = 15*(l1 - l2*Math.cos(2*theta))*Math.cos(theta);
        final double e3 = (l2 - l1*Math.cos(2*theta))*(-34*l1 + 50*l2)*Math.sin(theta);
        final double e4 = (l1 - l2)*(21*Math.sin(2*theta)*Math.sin(theta) - 12*(1 + Math.cos(2*theta))*Math.cos(theta));
        final double e5 = -9*Math.sin(2*theta)*Math.cos(theta);
        final double e6 = -15*(l2 - l1*Math.cos(2*theta))*Math.cos(theta);
        final double e7 = 9*Math.sin(2*theta)*Math.cos(theta);

        final double b1 = -Math.sin(phi - theta)*f1 + Math.cos(phi - theta)*e1;
        final double b2 = -Math.sin(phi - theta)*f2 + Math.cos(phi - theta)*e2;
        final double b3 = -Math.sin(phi - theta)*f3 + Math.cos(phi - theta)*e3;
        final double b4 = -Math.sin(phi - theta)*f4 + Math.cos(phi - theta)*e4;
        final double b5 = -Math.sin(phi - theta)*f5 + Math.cos(phi - theta)*e5;
        final double b6 = -Math.sin(phi - theta)*f6 + Math.cos(phi - theta)*e6;
        final double b7 = -Math.sin(phi - theta)*f7 + Math.cos(phi - theta)*e7;

        double d1, d2, beziermoment;
//        System.out.printf("spiro  <y'> = %f\n", spiro_moment_y());             // test code only
//        double testd1 = (0.265202 + .01)*a_b;
//        double testd2 = (0.265202 - .01)*a_b;
//        System.out.printf("bezier <y'> = %f\n", (b1*testd1 + b2*testd1*testd1 + (b3 + b4*testd1 + b5*testd1*testd1)*testd2 + (b6 + b7*testd1)*testd2*testd2)/280);    // test code only
        d1 = solve_quartic(a1*a3*b5 + a3*a3*b2,
             a1*a1*b7 + (a1*a2 + a0*a3)*b5 + a1*a3*b4 + 2*a2*a3*b2 + a3*a3*b1,
             a1*a1*b6 + 2*a0*a1*b7 + a0*a2*b5 + (a1*a2 + a0*a3)*b4
           + a1*a3*b3 + a2*a2*b2 + 2*a2*a3*b1 - a3*a3*280*spiro_moment_y(),
             2*a0*a1*b6 + a0*a0*b7 + a0*a2*b4 + (a1*a2 + a0*a3)*b3
           + a2*a2*b1 - 2*a2*a3*280*spiro_moment_y(),
             a0*a0*b6 + a0*a2*b3 - a2*a2*280*spiro_moment_y());
        d2 = (a0 + a1*d1)/(a2 + a3*d1);
        beziermoment = b1*d1 + b2*d1*d1 + (b3 + b4*d1 + b5*d1*d1)*d2 + (b6 + b7*d1)*d2*d2;
        System.out.println("d1, d2 = " + d1 + ", " + d2);
        System.out.println("area   = " + spiro_area() + ", " + bezier_area(d1, d2));
        System.out.println("<y>    = " + spiro_moment_y() + ", " + beziermoment/280);
    }

    private static double spiro_area()
    {
        return a_b*a_b*(Math.PI/2 - Math.sqrt(2))/4
               - c*c*(3*Math.PI/2 - Math.sqrt(2))/4;
    }

    private static double spiro_moment_y()
    {
        final double A = Math.sqrt(2) - 1;
        final double B = -3*Math.sqrt(2)/5 + 1;
        final double C = -47*Math.sqrt(2)/7 + 1;
        final double D = -Math.sqrt(2) - 1;

        final double momx = (a_b*a_b*a_b*A + a_b*a_b*c*B + a_b*c*c*C + c*c*c*D)/12;
        final double momy = (a_b*a_b*a_b*A*(Math.sqrt(2) - 1) - a_b*a_b*c*B/(Math.sqrt(2) - 1) + a_b*c*c*C*(Math.sqrt(2) - 1) - c*c*c*D/(Math.sqrt(2) - 1))/12;
        System.out.printf("A, B, C, D, momx, momy = %f, %f, %f, %f, %f, %f\n", A, B, C, D, momx, momy);
        return -momx*Math.sin(phi) + momy*Math.cos(phi);
    }

    private static double bezier_area(double d1, double d2)
    {
        return 3*(2*d1*(l1 - l2*Math.cos(2*theta)) + 2*d2*(l2 - l1*Math.cos(2*theta)) - d1*d2*Math.sin(2*theta))/20;
    }

    private static double solve_quartic(double lead, double qua, double qub, double quc, double qud)
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
        System.out.print("roots = ," + c + ", " + (-qua/4 + R/2 + D/2));
        System.out.print(", " + (-qua/4 + R/2 - D/2));
        System.out.print(", " + (-qua/4 - R/2 + E/2));
        System.out.println(", " + (-qua/4 - R/2 - E/2));
        if (!Double.isNaN(E) && false)
        {
            System.out.println("using root 4 = " + (-qua/4 - R/2 - E/2));
            return (-qua/4 - R/2 - E/2);
        }
        if (!Double.isNaN(D) && false)
        {
            System.out.println("using root 1 = " + (-qua/4 + R/2 + D/2));
            return (-qua/4 + R/2 + D/2);
        }
        if (!Double.isNaN(D))
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
//            System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-cua/3)*Math.cos(phi/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(phi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(phi/3 + 4*Math.PI/3) - p/3));
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
}
