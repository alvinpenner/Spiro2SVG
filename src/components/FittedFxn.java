
package components;

// abstract class to represent a function we are attempting to fit
// using Orthogonal Distance Fitting in class t2_vs_t1.java
// the argument 't' in this class will normally be 't1' in t2_vs_t1.java

public abstract class FittedFxn
{
    private final double origin_x = 80;             // just for svg output
    private final double origin_y = 500;            // just for svg output
    private final double scale = 200;               // just for svg output

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

    private static final double a = 240;
    private static final double b = -60;
    private static double c;

    public epiTrochoidFxn(double m_c)
    {
        c = m_c;
        //gen_points(0, Math.PI/4, 100);
    }

    protected double getc()
    {
        return c;
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
