
package components;

// abstract class to represent a function we are attempting to fit
// using Orthogonal Distance Fitting in class t2_vs_t1.java
// the argument 't' in this class will normally be 't1' in t2_vs_t1.java

public abstract class FittedFxn
{
    protected abstract double getc();
    protected abstract double getx(double t);
    protected abstract double gety(double t);
    protected abstract double getdydx(double t);
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

    protected double getdydx(double t)
    {
        if (c*Math.cos(t) == 1)
            return Double.POSITIVE_INFINITY;
        return c*Math.sin(t)/(1 - c*Math.cos(t));
    }
}

class CircleFxn extends FittedFxn
{
    private double c;

    public CircleFxn(double m_c)
    {
        c = m_c;
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

    protected double getdydx(double t)
    {
        if (Math.sin(t) == 0)
            return Double.NEGATIVE_INFINITY;
        return -Math.cos(t)/Math.sin(t);
    }
}
