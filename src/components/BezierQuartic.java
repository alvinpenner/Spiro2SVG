
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 5-point quartic Bezier (P0 - P4) to it, using parameter 0 < t2 < 1.
// constrain only the slopes at the endpoints and keep P2 arbitrary
// linearize the equations wrt (d1, d2, x2, y2) and solve a 4x4 system of equations
// Bezier[5] = f(x0, x1, x2, x3, x4, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book2, Dec 2016, page 65

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BezierQuartic.java

public class BezierQuartic
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI; // Math.PI/4;
    public static final int N = 100;
    public static double[] Bezx;                // quartic Bezier, 5 points, x component
    public static double[] Bezy;                // quartic Bezier, 5 points, y component
    //private static CircleFxn fitted;
    private static CycloidFxn fitted;
    //private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[] t2dd1 = new double[N+1];            // partial wrt d1
    private static double[] t2dd2 = new double[N+1];            // partial wrt d2
    private static double[] t2dx2 = new double[N+1];            // partial wrt x2
    private static double[] t2dy2 = new double[N+1];            // partial wrt y2
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 80;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //fitted = new epiTrochoidFxn(20);
        paste_data(80, 0.9886277018143461, 0.2624690492231421, 0.7640885237135162, 1.8275481762035848);
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        // convert BSpline5 (d1,d2) to Quartic Bezier (d1, d2)
        //System.out.println("solve_at_P2 = " + solve_at_P2(1.5*0.3820442618567581, 1.5*0.9137740881017924, 0.7263404994170831, 1.393169269182819, true));
        //System.out.println("solve_at_P2 = " + solve_at_P2(1.5*0.382, 1.5*0.914, 0.726, 1.393, true) + "\n");
        //grid_search_at_P2(.571, 1.369, 0.724, 1.394);
    }

    private static void paste_data(double phi, double tempc, double dummy, double d1, double d2)
    {
        // paste cubic Bezier (d1, d2) data from a file
        // field 1 = phi for a Cycloid / root number for an epiTrochoid
        // field 2 = c for an epiTrochoid
        // convert to d1, d2, and P2 for a 5 point quartic Bezier, and initiallize 'solve_at_P2()'

        tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        fitted = new CycloidFxn(tempc);
        t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        System.out.println("cubic Bezier data at phi c t d1 d2   = ," + phi + ", , " + tempc + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + fitted.getkappa(t1_start) + ", " + fitted.getkappa(t1_end));

        // calculate P2 for a quartic Bezier

        theta_start = fitted.gettheta(t1_start);                    // possibly redundant
        theta_end = fitted.gettheta(t1_end);
        double[] Bez = new double[] {fitted.getx(t1_start),
                                     fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                     fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                     fitted.getx(t1_end)};
        double P2x = (Bez[1] + Bez[2])/2;       // convert from Bezier to B-Spline
        Bez = new double[] {fitted.gety(t1_start),
                            fitted.gety(t1_start) + d1*Math.sin(theta_start),
                            fitted.gety(t1_end) - d2*Math.sin(theta_end),
                            fitted.gety(t1_end)};
        double P2y = (Bez[1] + Bez[2])/2;
        System.out.println("quartic Bezier solve_at_P2 = " + solve_at_P2(0.75*d1, 0.75*d2, P2x, P2y, true) + "\n");
    }

    private static void grid_search_at_P2(double d1, double d2, double x2, double y2)
    {
        double del = 0.01;

        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
                for (int k = -1; k < 2; k++)
                    for (int l = -1; l < 2; l++)
                    {
                        System.out.println("grid =, " + i + ", " + j + ", " + k + ", " + l);
                        solve_at_P2(d1 + i*del, d2 + j*del, x2 + k*del, y2 + l*del, false);
                    }
    }

    private static double solve_at_P2(double d1, double d2, double x2, double y2, boolean print)
    {
        // 5-point quartic Bezier curve
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        Bezx = new double[] {fitted.getx(t1_start),
                             fitted.getx(t1_start) + d1*Math.cos(theta_start),
                             x2,
                             fitted.getx(t1_end) - d2*Math.cos(theta_end),
                             fitted.getx(t1_end)};
        Bezy = new double[] {fitted.gety(t1_start),
                             fitted.gety(t1_start) + d1*Math.sin(theta_start),
                             y2,
                             fitted.gety(t1_end) - d2*Math.sin(theta_end),
                             fitted.gety(t1_end)};

        if (t2[N] == 0)
            System.out.println("__start at theta c t d1 d2 = , " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + calc_error());
        //gen_BezierQuad(Bezx, Bezy);
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);
        //System.out.println(Bezx[4] + "\t " + Bezy[4]);

        if (print) System.out.println("\n   , t1, t2, t2dd1, t2dd2, t2dx2, t2dy2");
        for (int i = 0; i <= N; i++)
        {
            solve_septic_for_t2(i);
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(1 - t2[i]) > TOL)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //scan_septic_near_t2(i, t2[i]);
                return Double.NaN;
            }
            t2dd1[i] = calc_t2dxy(i, t2[i], "d1");
            t2dd2[i] = calc_t2dxy(i, t2[i], "d2");
            t2dx2[i] = calc_t2dxy(i, t2[i], "x2");
            t2dy2[i] = calc_t2dxy(i, t2[i], "y2");
            if (print)
                System.out.println("quartic Bez, " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i] + ", " + t2dx2[i] + ", " + t2dy2[i]);
        }
        double retVal = calc_error();
        System.out.println("__new t2[] at theta c t d1 d2 rms = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + retVal);
        return retVal;
    }

    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

        //double a_b = 180;         // scale factor to make rms error dimensionless
        double a_b = 1;             // Cycloid only
        double t1 = t1_start;
        double[] trap_in = new double[N+1];

        if ((Math.abs(1 - t2[N]) > TOL) || (Math.abs(t2[0]) > TOL))
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i <= N; i++)
        {
            trap_in[i] = (fn(Bezx, t2[i]) - fitted.getx(t1))*(fn(Bezx, t2[i]) - fitted.getx(t1))
                       + (fn(Bezy, t2[i]) - fitted.gety(t1))*(fn(Bezy, t2[i]) - fitted.gety(t1));
            //System.out.println(i + ", " + ", " + (fn(Bezx, t2[i]) - fitted.getx(t1)) + ", " + (fn(Bezy, t2[i]) - fitted.gety(t1)) + ", " + Math.sqrt(trap_in[i]));
            t1 += (t1_end - t1_start)/N;
        }
        return Math.sqrt(integrate(trap_in))/a_b;
    }

    private static double integrate(double[] trap)
    {
        // trapezoidal rule integration of a fxn of t1 (N+1 points)

        //System.out.println("trap length = " + trap.length);
        double ret = (trap[0] + trap[trap.length - 1])/2;
        for (int i = 1; i < trap.length - 1; i++)
            ret += trap[i];
        return ret/(trap.length - 1);
    }

    private static void solve_septic_for_t2(int i)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the quartic Bezier t-value
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double f, fprime, f2prime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double t;
        int loop = 0;

        // initial estimate using quadratic approximation

        if (i == 0) t = 0;
        else t = t2[i-1];
        f = (fn(Bezx, t) - X)*dfn(Bezx, t) + (fn(Bezy, t) - Y)*dfn(Bezy, t);
        fprime = dfn(Bezx, t)*dfn(Bezx, t) + (fn(Bezx, t) - X)*d2fn(Bezx, t) + dfn(Bezy, t)*dfn(Bezy, t) + (fn(Bezy, t) - Y)*d2fn(Bezy, t);
        f2prime = 3*dfn(Bezx, t)*d2fn(Bezx, t) + (fn(Bezx, t) - X)*d3fn(Bezx, t) + 3*dfn(Bezy, t)*d2fn(Bezy, t) + (fn(Bezy, t) - Y)*d3fn(Bezy, t);
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
            del_t = -fprime/f2prime;
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            if (del_t < 0 || del_t > 1)
            {
                del_t = (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
                if (del_t < 0 || del_t > 1)
                {
                    System.out.println("\nBad init  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
                    t2[i] = Double.NaN;
                    return;
                }
            }
        }
        t += del_t;

        //System.out.println("roots = " + (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime + ", " + (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime);
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
        do
        {
            f = (fn(Bezx, t) - X)*dfn(Bezx, t) + (fn(Bezy, t) - Y)*dfn(Bezy, t);
            fprime = dfn(Bezx, t)*dfn(Bezx, t) + (fn(Bezx, t) - X)*d2fn(Bezx, t) + dfn(Bezy, t)*dfn(Bezy, t) + (fn(Bezy, t) - Y)*d2fn(Bezy, t);
            if (loop > 100)
            {
                t2[i] = Double.NaN;
                return;
            }
            if (f == 0 && fprime == 0)
                del_t = 0;
            else
                del_t = -f/fprime;
            t += del_t;
            loop++;
            //System.out.println("         t2 =, " + t + ", " + f + ", " + fprime);
        } while (Math.abs(del_t) > TOL);
        t2[i] = t;
    }

    private static double calc_t2dxy(int i, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double denom = dfn(Bezx, t2)*dfn(Bezx, t2) + (fn(Bezx, t2) - X)*d2fn(Bezx, t2) + dfn(Bezy, t2)*dfn(Bezy, t2) + (fn(Bezy, t2) - Y)*d2fn(Bezy, t2);
        double numer = Double.NaN;

        if (type.equals("x2"))
            numer = dfn(Bezx, t2)*N44(t2)[2] + (fn(Bezx, t2) - X)*dN44(t2)[2];
        else if(type.equals("y2"))
            numer = dfn(Bezy, t2)*N44(t2)[2] + (fn(Bezy, t2) - Y)*dN44(t2)[2];
        else if(type.equals("d1"))
            numer = (Math.cos(theta_start)*dfn(Bezx, t2) + Math.sin(theta_start)*dfn(Bezy, t2))*N44(t2)[1]
                  + (Math.cos(theta_start)*(fn(Bezx, t2) - X) + Math.sin(theta_start)*(fn(Bezy, t2) - Y))*dN44(t2)[1];
        else if(type.equals("d2"))
            numer = -(Math.cos(theta_end)*dfn(Bezx, t2) + Math.sin(theta_end)*dfn(Bezy, t2))*N44(t2)[3]
                  -  (Math.cos(theta_end)*(fn(Bezx, t2) - X) + Math.sin(theta_end)*(fn(Bezy, t2) - Y))*dN44(t2)[3];
        //System.out.println("calc_t2dx2 = " + t1 + ", " + t2 + ", " + rhs1 + ", " + rhs2 + ", " + fprime + ", " + ((-3*t2*(1 - t2)*(1 - t2)*rhs1 - 3*(1 - t2)*(1 - 3*t2)*rhs2)/fprime));
        return -numer/denom;
    }

    private static double[] N44(double u)
    {
        if (u < -TOL)
            return null;
        else
            return new double[] {(1 - u)*(1 - u)*(1 - u)*(1 - u),
                                 4*u*(1 - u)*(1 - u)*(1 - u),
                                 6*u*u*(1 - u)*(1 - u),
                                 4*u*u*u*(1 - u),
                                 u*u*u*u};
    }

    private static double[] dN44(double u)
    {
        if (u < -TOL)
            return null;
        else
            return new double[] {-4*(1 - u)*(1 - u)*(1 - u),
                                  4*(1 - 4*u)*(1 - u)*(1 - u),
                                 12*(1 - 2*u)*u*(1 - u),
                                  4*(3 - 4*u)*u*u,
                                  4*u*u*u};
    }

    private static double[] d2N44(double u)
    {
        if (u < -TOL)
            return null;
        else
            return new double[] { 12*(1 - u)*(1 - u),
                                 -24*(1 - 2*u)*(1 - u),
                                  12*(1 - 6*u + 6*u*u),
                                  24*(1 - 2*u)*u,
                                  12*u*u};
    }

    private static double[] d3N44(double u)
    {
        if (u < -TOL)
            return null;
        else
            return new double[] {-24*(1 - u),
                                  24*(3 - 4*u),
                                 -72*(1 - 2*u),
                                  24*(1 - 4*u),
                                  24*u};
    }

    private static void gen_BezierQuad(double[] ptx, double[] pty)
    {
        double origin_x = 80;                           // just for svg output
        double origin_y = 500;                          // just for svg output
        double scale = 200;                             // just for svg output

        System.out.printf("M");
        for (int i = 0; i <= N; i++)
            System.out.printf(" %f, %f", origin_x + scale*fn(ptx, (double) i/N), origin_y - scale*fn(pty, (double) i/N));
        System.out.println("\n");
    }

    private static double fn(double[] Bez, double t)
    {
        return BSpline5.multvv(Bez, N44(t));
    }

    private static double dfn(double[] Bez, double t)
    {
        return BSpline5.multvv(Bez, dN44(t));
    }

    private static double d2fn(double[] Bez, double t)
    {
        return BSpline5.multvv(Bez, d2N44(t));
    }

    private static double d3fn(double[] Bez, double t)
    {
        return BSpline5.multvv(Bez, d3N44(t));
    }
}
