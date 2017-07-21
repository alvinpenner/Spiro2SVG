
package components;

// abstract class to represent a function we are attempting to fit
// using Orthogonal Distance Fitting in class t2_vs_t1.java
// the argument 't' in this class will normally be 't1' in t2_vs_t1.java

public abstract class FittedFxn
{
    private final double origin_x = -100;        // just for svg output
    private final double origin_y = 400;         // just for svg output
    private final double scale = 2;              // just for svg output

    protected abstract double getc();
    protected abstract double getx(double t);
    protected abstract double gety(double t);
    protected abstract double gettheta(double t);

    protected void gen_points(double t1, double t2, int N)
    {
        System.out.printf("M");
        for (int i = 0; i <= N; i++)
            System.out.printf(" %f, %f", origin_x + scale*getx(t1 + i*(t2 - t1)/N), origin_y - scale*gety(t1 + i*(t2 - t1)/N));
            //System.out.printf("%f, %f, %f\n", getx(t1 + i*(t2 - t1)/N), gety(t1 + i*(t2 - t1)/N), getdydx(t1 + i*(t2 - t1)/N));
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
}

class CycloidFxn extends FittedFxn
{
    private double c;

    public CycloidFxn(double m_c)
    {
        c = m_c;
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
        // angle of the velocity vector
        return Math.atan2(c*Math.sin(t), 1 - c*Math.cos(t));
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
        // angle of the velocity vector
        return Math.atan2(Math.cos(t) + c/b*Math.cos(t*(a/b + 1)), -Math.sin(t) - c/b*Math.sin(t*(a/b + 1)));
    }
}
