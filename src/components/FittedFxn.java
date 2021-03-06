
package components;

// abstract class to represent a function we are attempting to fit
// using Orthogonal Distance Fitting in class t2_vs_t1.java
// the argument 't' in this class will normally be 't1' in t2_vs_t1.java

public abstract class FittedFxn
{
    private final double origin_x = 250;             // just for svg output
    private final double origin_y = 250; // 500; // 700;            // just for svg output
    private final double scale = 1; // 1.5; // 2; // 200;               // just for svg output

    protected abstract double getc();
    protected abstract double getx(double t);
    protected abstract double gety(double t);
    protected abstract double gettheta(double t);   // angle of the velocity vector
    protected abstract double getkappa(double t);   // curvature

    protected void gen_points(double t1, double t2, int N)
    {
        System.out.printf("M");
        for (int i = 0; i <= N; i++)
            System.out.printf(" %f, %f", origin_x + scale*getx(t1 + i*(t2 - t1)/N), origin_y - scale*gety(t1 + i*(t2 - t1)/N));
        System.out.println("\n");
    }

    protected void gen_Bezier(double[] pts)
    {
        System.out.printf("M %f, %f C %f, %f %f, %f %f, %f\n",
                          origin_x + scale*pts[0], origin_y - scale*pts[1],
                          origin_x + scale*pts[2], origin_y - scale*pts[3],
                          origin_x + scale*pts[4], origin_y - scale*pts[5],
                          origin_x + scale*pts[6], origin_y - scale*pts[7]);
    }

    protected void gen_Bezier2(double[][] ptx, double[][] pty)
    {
        System.out.printf("M %f, %f C %f, %f %f, %f %f, %f C %f, %f %f, %f %f, %f\n",
                          origin_x + scale*ptx[0][0], origin_y - scale*pty[0][0],
                          origin_x + scale*ptx[0][1], origin_y - scale*pty[0][1],
                          origin_x + scale*ptx[0][2], origin_y - scale*pty[0][2],
                          origin_x + scale*ptx[0][3], origin_y - scale*pty[0][3],
                          origin_x + scale*ptx[1][1], origin_y - scale*pty[1][1],
                          origin_x + scale*ptx[1][2], origin_y - scale*pty[1][2],
                          origin_x + scale*ptx[1][3], origin_y - scale*pty[1][3]);
    }

    protected void gen_Bezier3(double[][] ptx, double[][] pty)
    {
        System.out.printf("M %f, %f C %f, %f %f, %f %f, %f C %f, %f %f, %f %f, %f C %f, %f %f, %f %f, %f\n",
                          origin_x + scale*ptx[0][0], origin_y - scale*pty[0][0],
                          origin_x + scale*ptx[0][1], origin_y - scale*pty[0][1],
                          origin_x + scale*ptx[0][2], origin_y - scale*pty[0][2],
                          origin_x + scale*ptx[0][3], origin_y - scale*pty[0][3],
                          origin_x + scale*ptx[1][1], origin_y - scale*pty[1][1],
                          origin_x + scale*ptx[1][2], origin_y - scale*pty[1][2],
                          origin_x + scale*ptx[1][3], origin_y - scale*pty[1][3],
                          origin_x + scale*ptx[2][1], origin_y - scale*pty[2][1],
                          origin_x + scale*ptx[2][2], origin_y - scale*pty[2][2],
                          origin_x + scale*ptx[2][3], origin_y - scale*pty[2][3]);
    }

    protected void gen_BezierN(double[] ptx, double[] pty)
    {
        // ptx[], pty[] are N-point uniform cubic B-Splines (closed)
        double[] ptx1 = new double[ptx.length];
        double[] ptx2 = new double[ptx.length];
        double[] pty1 = new double[pty.length];
        double[] pty2 = new double[pty.length];

        for (int i = 0; i < ptx.length; i++)
        {
            ptx1[i] = (2*ptx[i] + ptx[(i + 1) % ptx.length])/3;
            ptx2[i] = (ptx[i] + 2*ptx[(i + 1) % ptx.length])/3;
            pty1[i] = (2*pty[i] + pty[(i + 1) % pty.length])/3;
            pty2[i] = (pty[i] + 2*pty[(i + 1) % pty.length])/3;
            //System.out.println("gen_BezierN = ," + ptx[i] + ", " + ptx1[i] + ", " + ptx2[i] + ", " + pty[i] + ", " + pty1[i] + ", " + pty2[i]);
        }
        System.out.printf("M %f, %f ", origin_x + scale*(ptx2[ptx.length - 1] + ptx1[0])/2, origin_y - scale*(pty2[pty.length - 1] + pty1[0])/2);
        for (int i = 0; i < ptx.length; i++)
            System.out.printf("C %f, %f %f, %f %f, %f ",
                                        origin_x + scale*ptx1[i], origin_y - scale*pty1[i],
                                        origin_x + scale*ptx2[i], origin_y - scale*pty2[i],
                                        origin_x + scale*(ptx2[i] + ptx1[(i + 1)%ptx.length])/2, origin_y - scale*(pty2[i] + pty1[(i + 1)%pty.length])/2);
        System.out.printf("\n");
    }
}

class CycloidFxn extends FittedFxn
{
    private double c;

    public CycloidFxn(double m_c)
    {
        c = m_c;
        //gen_points(0, Math.PI, 100);
    }

    protected double getc()
    {
        return c;
    }

    protected double getx(double t)
    {
        return t - c*Math.sin(t);
    }

    protected double gety(double t)
    {
        return 1 - c*Math.cos(t);
    }

    protected double gettheta(double t)
    {
        if (c == 1 && t == 0)
            return Math.PI/2;
        return Math.atan2(c*Math.sin(t), 1 - c*Math.cos(t));
    }

    protected double getkappa(double t)
    {
        double xdot = 1 - c*Math.cos(t);
        double ydot = c*Math.sin(t);
        double x2dot = c*Math.sin(t);
        double y2dot = c*Math.cos(t);
        double v2 = xdot*xdot + ydot*ydot;
        if (v2 > 0)
            return (xdot*y2dot - ydot*x2dot)/Math.pow(v2, 1.5);
        return Double.POSITIVE_INFINITY;
    }
}

class epiTrochoidFxn extends FittedFxn
{
    //  see : http://turnbull.mcs.st-and.ac.uk/~history/Curves/Epitrochoid.html
    //  a = stator, b = rotor, c = pen distance
    //  x = (a + b) cos(t) + c cos((a/b + 1)t)
    //  y = (a + b) sin(t) + c sin((a/b + 1)t)

    protected final double a = 240; // 120;     // 90;  // 240;
    protected final double b = -60; // 60;    // 90;  // -60;
    private static double c;

    public epiTrochoidFxn(double m_c)
    {
        c = m_c;
        //gen_points(0, 2*Math.PI, 800);
    }

    protected double getc()
    {
        return c;
    }

    protected double getdxdc(double t)
    {
        return Math.cos(t*(a/b + 1));
    }

    protected double getdydc(double t)
    {
        return Math.sin(t*(a/b + 1));
    }

    protected double getx(double t)
    {
        return (a + b)*Math.cos(t) + c*Math.cos(t*(a/b + 1));
    }

    protected double gety(double t)
    {
        return (a + b)*Math.sin(t) + c*Math.sin(t*(a/b + 1));
    }

    protected double gettheta(double t)
    {
        return Math.atan2(Math.cos(t) + c/b*Math.cos(t*(a/b + 1)), -Math.sin(t) - c/b*Math.sin(t*(a/b + 1)));
    }

    protected double getkappa(double t)
    {
        double xdot = -(a + b)*Math.sin(t) - c*(a/b + 1)*Math.sin(t*(a/b + 1));
        double ydot =  (a + b)*Math.cos(t) + c*(a/b + 1)*Math.cos(t*(a/b + 1));
        double x2dot = -(a + b)*Math.cos(t) - c*(a/b + 1)*(a/b + 1)*Math.cos(t*(a/b + 1));
        double y2dot = -(a + b)*Math.sin(t) - c*(a/b + 1)*(a/b + 1)*Math.sin(t*(a/b + 1));
        double v2 = xdot*xdot + ydot*ydot;
        if (v2 > 0)
            return (xdot*y2dot - ydot*x2dot)/Math.pow(v2, 1.5);
        return Double.POSITIVE_INFINITY;
    }
}

class CircleFxn extends FittedFxn
{
    private double c;

    public CircleFxn(double m_c)
    {
        c = m_c;
        //gen_points(0, Math.PI, 100);
    }

    protected double getc()
    {
        return c;
    }

    protected double getx(double t)
    {
        return c*Math.cos(t);
    }

    protected double gety(double t)
    {
        return c*Math.sin(t);
    }

    protected double gettheta(double t)
    {
        if (t == 0)
            return Math.PI/2;
        return Math.atan2(Math.cos(t), -Math.sin(t));
    }

    protected double getkappa(double t)
    {
        return 1.0/c;
    }
}

class HippopedeFxn extends FittedFxn
{
    // see : J. D. Lawrence, p. 145
    // our c = b/a: 0 < c < 0.5
    // define 4*a*b = 180*180
    // x = 2 cos(t) sqrt(a*b - b*b*sin(t)*sin(t))
    // y = 2 sin(t) sqrt(a*b - b*b*sin(t)*sin(t))

    protected double a;
    protected double b;
    private double c;

    public HippopedeFxn(double m_c)
    {
        c = m_c;
        a = 90/Math.sqrt(c);
        b = 90*Math.sqrt(c);
        System.out.println("HippopedeFxn: " + a + ", " + b + ", " + c);
        //gen_points(0, 2*Math.PI, 800);
        //for (int i = 0; i <= 800; i++)
        //    System.out.printf("%d, %f, %f, %f, %f, %f\n", i, getx(i*2*Math.PI/800), gety(i*2*Math.PI/800), getdxdc(i*2*Math.PI/800), getdydc(i*2*Math.PI/800), gettheta(i*2*Math.PI/800));
        //for (int i = 0; i <= 800; i++)
        //    System.out.println(i + ", " + getx(i*2*Math.PI/800) + ", " + gety(i*2*Math.PI/800) + ", " + getdxdc(i*2*Math.PI/800) + ", " + getdydc(i*2*Math.PI/800) + ", " + getd2xdc2(i*2*Math.PI/800) + ", " + getd2ydc2(i*2*Math.PI/800) + ", " + gettheta(i*2*Math.PI/800));
    }

    protected double getc()
    {
        return c;
    }

    protected double getx(double t)
    {
        return 2*Math.cos(t)*Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t));
    }

    protected double gety(double t)
    {
        return 2*Math.sin(t)*Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t));
    }

    protected double getdxdc(double t)
    {
        return -Math.cos(t)*a*b*Math.sin(t)*Math.sin(t)/Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t));
    }

    protected double getdydc(double t)
    {
        return -Math.sin(t)*a*b*Math.sin(t)*Math.sin(t)/Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t));
    }

    protected double getd2xdc2(double t)
    {
        return a*b*Math.sin(t)*Math.sin(t)/2/(a*b - b*b*Math.sin(t)*Math.sin(t))*getdxdc(t);
    }

    protected double getd2ydc2(double t)
    {
        return a*b*Math.sin(t)*Math.sin(t)/2/(a*b - b*b*Math.sin(t)*Math.sin(t))*getdydc(t);
    }

    protected double getdxdt(double t)
    {
        return -2*Math.sin(t)*Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t))
               - Math.cos(t)*b*b*Math.sin(2*t)/Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t));
    }

    protected double getdydt(double t)
    {
        return 2*Math.cos(t)*Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t))
               - Math.sin(t)*b*b*Math.sin(2*t)/Math.sqrt(a*b - b*b*Math.sin(t)*Math.sin(t));
    }

    protected double gettheta(double t)
    {
        return Math.atan2(getdydt(t), getdxdt(t));
    }

    protected double getkappa(double t)
    {
        return Double.NaN;
    }
}

class SuperEllipse extends FittedFxn
{
    // see : S. J. Ahn - Orthogonal distance fitting of implicit curves and surfaces
    // https://doi.org/10.1109/34.1000237
    // (x/a)^(2/c) + (y/a)^(2/c) = 1
    // x/a = cos(t)^c
    // y/a = sin(t)^c

    protected final double a = 180;
    protected final double b = 0;
    private double c;

    public SuperEllipse(double m_c)
    {
        c = m_c;
        System.out.println("SuperEllipse: ," + a + ", " + c);
        //gen_points(0, 2*Math.PI, 800);
        //System.out.println("i, getx, gety, getdxdt, getdydt, gettheta");
        //for (int i = 0; i <= 800; i++)
        //    System.out.printf("%d, %f, %f, %f, %f, %f, %f\n", i, getx(i*2*Math.PI/800), gety(i*2*Math.PI/800), getdxdc(i*2*Math.PI/800), getdydc(i*2*Math.PI/800), getd2xdc2(i*2*Math.PI/800), getd2ydc2(i*2*Math.PI/800));
        //    System.out.printf("%d, %f, %f, %f, %f, %f\n", i, getx(i*2*Math.PI/800), gety(i*2*Math.PI/800), getdxdt(i*2*Math.PI/800), getdydt(i*2*Math.PI/800), gettheta(i*2*Math.PI/800));
        //for (int i = 0; i <= 800; i++)
        //    System.out.println(i + ", " + getx(i*2*Math.PI/800) + ", " + gety(i*2*Math.PI/800) + ", " + getdxdc(i*2*Math.PI/800) + ", " + getdydc(i*2*Math.PI/800) + ", " + getd2xdc2(i*2*Math.PI/800) + ", " + getd2ydc2(i*2*Math.PI/800) + ", " + gettheta(i*2*Math.PI/800));
        //for (int i = -5; i <= 5; i++)
        //{
        //    double th = i*Math.PI/180;
        //    System.out.println(th + ", " + getdxdt(th) + ", " + getdydt(th) + ", " + gettheta(th));
        //    th += 2*Math.PI;
        //    System.out.println(th + ", " + getdxdt(th) + ", " + getdydt(th) + ", " + gettheta(th));
        //}
    }

    protected double getc()
    {
        return c;
    }

    protected double getx(double t)
    {
        return Math.signum(Math.cos(t))*a*Math.pow(Math.abs(Math.cos(t)), c);
    }

    protected double gety(double t)
    {
        return Math.signum(Math.sin(t))*a*Math.pow(Math.abs(Math.sin(t)), c);
    }

    protected double getdxdc(double t)
    {
        if (Math.cos(t) == 0)
            return 0;
        else
            return getx(t)*Math.log(Math.abs(Math.cos(t)));
    }

    protected double getdydc(double t)
    {
        if (Math.sin(t) == 0)
            return 0;
        else
            return gety(t)*Math.log(Math.abs(Math.sin(t)));
    }

    protected double getd2xdc2(double t)
    {
        if (Math.cos(t) == 0)
            return 0;
        else
            return getx(t)*Math.log(Math.abs(Math.cos(t)))*Math.log(Math.abs(Math.cos(t)));
    }

    protected double getd2ydc2(double t)
    {
        if (Math.sin(t) == 0)
            return 0;
        else
            return gety(t)*Math.log(Math.abs(Math.sin(t)))*Math.log(Math.abs(Math.sin(t)));
    }

    protected double getdxdt(double t)
    {
        return -a*c*Math.pow(Math.abs(Math.cos(t)), c - 1)*Math.sin(t);
    }

    protected double getdydt(double t)
    {
        return a*c*Math.pow(Math.abs(Math.sin(t)), c - 1)*Math.cos(t);
    }

    protected double gettheta(double t)
    {
        //System.out.println("gettheta = " + t + ", " + getdxdt(t) + ", " + getdydt(t));
        if (t == 0)
            return Math.PI/2;
        else if (t == 2*Math.PI)
            return Math.PI/2;
        else
            return Math.atan2(getdydt(t), getdxdt(t));
    }

    protected double getkappa(double t)                 // fix fix bug bug
    {
        double xdot  = -a*c*Math.pow(Math.abs(Math.cos(t)), c - 1)*Math.sin(t);
        double ydot  =  a*c*Math.pow(Math.abs(Math.sin(t)), c - 1)*Math.cos(t);
        double x2dot =  a*c*(c - 1)*Math.pow(Math.abs(Math.cos(t)), c - 2)*Math.sin(t)*Math.sin(t)
                     -  a*c*Math.pow(Math.abs(Math.cos(t)), c);
        double y2dot =  a*c*(c - 1)*Math.pow(Math.abs(Math.sin(t)), c - 2)*Math.cos(t)*Math.cos(t)
                     -  a*c*Math.pow(Math.abs(Math.sin(t)), c);
        double v2 = xdot*xdot + ydot*ydot;
        //System.out.print(a*Math.pow(Math.abs(Math.cos(t)), c) + ", " + xdot + ", " + x2dot + ", ");
        //System.out.print(a*Math.pow(Math.abs(Math.sin(t)), c) + ", " + ydot + ", " + y2dot + ", ");
        if (v2 > 0)
            return (xdot*y2dot - ydot*x2dot)/Math.pow(v2, 1.5);
        return Double.POSITIVE_INFINITY;
    }
}
