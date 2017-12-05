
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1.
// fit a 5-point Beta2-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// constrain only the slopes at the endpoints and keep P2 (Q3/Q4) arbitrary.
// the Beta2 spline has an adjustable arm length, d, at the Bezier splice at Q3/Q4.
// psi, ψ, is the angle at the Bezier splice (Q3-Q2/Q5-Q4) = angle of (Q6-Q1).
// linearize the equations wrt (d1, d2, x2, y2, d) and solve a 5x5 system of equations.
// decompose the Beta2-Spline into 2 Beziers, range (0,1) and (1,2).
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2).
// t2 must be chosen to minimize the distance to g(t1).
// see Spiro2SVG Book4, Dec 2017, page 13.

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\Beta2Spline.java

public class Beta2Spline
{
    public static double t1_start = 0;
    public static final double t1_end = Math.PI/4; // Math.PI/4;
    public static final int N = 100;
    public static double[][] Bezx;              // 2 Beziers, 4 points each, x component
    public static double[][] Bezy;              // 2 Beziers, 4 points each, y component
    //private static CircleFxn fitted;
    //private static CycloidFxn fitted;
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[] t2dd1 = new double[N+1];            // partial wrt d1
    private static double[] t2dd2 = new double[N+1];            // partial wrt d2
    private static double[] t2dx2 = new double[N+1];            // partial wrt x2
    private static double[] t2dy2 = new double[N+1];            // partial wrt y2
    private static double[] t2dd = new double[N+1];             // partial wrt symmetric d
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 40;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        fitted = new epiTrochoidFxn(16.2);
        //System.out.println("Beta2-Spline solve_at_P2 = " + convert_at_P2(0.4653789140116872, 0.4270959656087598, 1.8406794857301048, 1.6140740428609932, true) + "\n");
        System.out.println("Beta2-Spline solve_at_P2 = " + convert_at_P2(8.504261491368203, 51.59554117494055, 192.45848525164885, 32.91614151837636, true) + "\n");
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
    }

    private static double convert_at_P2(double d1, double d2, double x2, double y2, boolean print)
    {
        // convert a 5-point B-Spline to a pair of spliced cubic Beziers
        System.out.println("convert_at_P2: B-Spline theta c t d1 d2 = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        double x1 = fitted.getx(t1_start) + d1*Math.cos(theta_start);
        double y1 = fitted.gety(t1_start) + d1*Math.sin(theta_start);
        double x3 = fitted.getx(t1_end) - d2*Math.cos(theta_end);
        double y3 = fitted.gety(t1_end) - d2*Math.sin(theta_end);
        double d = 0.25*Math.sqrt((x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1));

        return solve_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, d, print);
    }

    private static double solve_at_P2(double d1, double d2, double x2, double y2, double d, boolean print)
    {
        // 5-point cubic Beta2-Spline curve, decomposed as two Beziers
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        Bezx = new double[][] {{fitted.getx(t1_start),
                                fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                x2,
                                x2},
                               {x2,
                                x2,
                                fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                fitted.getx(t1_end)}};
        Bezy = new double[][] {{fitted.gety(t1_start),
                                fitted.gety(t1_start) + d1*Math.sin(theta_start),
                                y2,
                                y2},
                               {y2,
                                y2,
                                fitted.gety(t1_end) - d2*Math.sin(theta_end),
                                fitted.gety(t1_end)}};


        double psi = Math.atan2(Bezy[1][2] - Bezy[0][1], Bezx[1][2] - Bezx[0][1]);
        Bezx[0][2] -= d*Math.cos(psi);
        Bezx[1][1] += d*Math.cos(psi);
        Bezy[0][2] -= d*Math.sin(psi);
        Bezy[1][1] += d*Math.sin(psi);

        if (t2[N] == 0)
            System.out.println("__start Beta2-Spline at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d + ", " + calc_error());
        //fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2, t2dd1, t2dd2, t2dx2, t2dy2, t2dd");
        int seg = 0;                // Bezier segment, before or after the splice
        for (int i = 0; i <= N; i++)
        {
            solve_quintic_for_t2(i, seg);
            if (seg == 0 && t2[i] > 1)
            {
                seg++;
                solve_quintic_for_t2(i, seg);   // re-calculate
            }
            if (seg == 1 && t2[i] < 1)
                return Double.NaN;
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //scan_quintic_near_t2(i, seg, t2[i]);
                return Double.NaN;
            }
            t2dd1[i] = calc_t2dxy(i, t2[i], "d1");
            t2dd2[i] = calc_t2dxy(i, t2[i], "d2");
            t2dx2[i] = calc_t2dxy(i, t2[i], "x2");
            t2dy2[i] = calc_t2dxy(i, t2[i], "y2");
            if (print)
                System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i] + ", " + t2dx2[i] + ", " + t2dy2[i]);
        }
        //System.out.println("new t2[] profile rms   = ," + d1 + ", " + d2 + ", " + calc_error());
        double retVal = calc_error();
        System.out.println("__new t2[] at theta c t d1 d2 rms = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + d + ", " + retVal);
        return retVal;
    }

    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

        double a_b = 180;         // scale factor to make rms error dimensionless
        //double a_b = 1;             // Cycloid only
        double t1 = t1_start;
        double[] trap_in = new double[N+1];
        int seg = 0;                // Bezier segment, before or after the splice

        if ((Math.abs(2 - t2[N]) > TOL) || (Math.abs(t2[0]) > TOL))
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort = " + t2[0] + ", " + t2[N]);
            return Double.NaN;
        }
        for (int i = 0; i <= N; i++)
        {
            if (seg == 0 && t2[i] > 1)
                seg++;
            trap_in[i] = (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))*(t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1))
                       + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1))*(t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1));
            //System.out.println(i + ", " + seg + ", " + (t2_vs_t1.fn(Bezx[seg], t2[i] - seg) - fitted.getx(t1)) + ", " + (t2_vs_t1.fn(Bezy[seg], t2[i] - seg) - fitted.gety(t1)) + ", " + Math.sqrt(trap_in[i]));
            t1 += (t1_end - t1_start)/N;
        }
        return Math.sqrt(t2_vs_t1.integrate(trap_in))/a_b;
    }

    private static void solve_quintic_for_t2(int i, int seg)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the cubic Bezier t-value

        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double f, fprime, f2prime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double t;
        int loop = 0;

        // initial estimate using quadratic approximation

        if (i == 0) t = 0;
        else t = t2[i-1];
        t -= seg;               // compensate for Bezier segment offset
        f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
        fprime = t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d2fn(Bezx[seg], t) + t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.dfn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d2fn(Bezy[seg], t);
        f2prime = 3*t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.d2fn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d3fn(Bezx[seg], t) + 3*t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.d2fn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d3fn(Bezy[seg], t);
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
        if (f == 0)
            del_t = 0;
        else if (fprime * fprime < 2 * f * f2prime)
            del_t = -fprime/f2prime;
        else
        {
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
            if (del_t < 0 || del_t > 2)
            {
                del_t = (-fprime - Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
                if (del_t < 0 || del_t > 2)
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
            f = (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.dfn(Bezy[seg], t);
            fprime = t2_vs_t1.dfn(Bezx[seg], t)*t2_vs_t1.dfn(Bezx[seg], t) + (t2_vs_t1.fn(Bezx[seg], t) - X)*t2_vs_t1.d2fn(Bezx[seg], t) + t2_vs_t1.dfn(Bezy[seg], t)*t2_vs_t1.dfn(Bezy[seg], t) + (t2_vs_t1.fn(Bezy[seg], t) - Y)*t2_vs_t1.d2fn(Bezy[seg], t);
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
        t2[i] = t + seg;               // compensate for Bezier segment offset
    }

    private static double calc_t2dxy(int i, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
/*        double fnx = multvv(Splinex, N43(t2));
        double fny = multvv(Spliney, N43(t2));
        double dfnx = multvv(Splinex, dN43(t2));
        double dfny = multvv(Spliney, dN43(t2));
        double d2fnx = multvv(Splinex, d2N43(t2));
        double d2fny = multvv(Spliney, d2N43(t2));
        double denom = dfnx*dfnx + (fnx - X)*d2fnx + dfny*dfny + (fny - Y)*d2fny;
        double numer = Double.NaN;
        if (type.equals("x2"))
            numer = dfnx*N43(t2)[2] + (fnx - X)*dN43(t2)[2];
        else if(type.equals("y2"))
            numer = dfny*N43(t2)[2] + (fny - Y)*dN43(t2)[2];
        else if(type.equals("d1"))
            numer = (Math.cos(theta_start)*dfnx + Math.sin(theta_start)*dfny)*N43(t2)[1]
                  + (Math.cos(theta_start)*(fnx - X) + Math.sin(theta_start)*(fny - Y))*dN43(t2)[1];
        else if(type.equals("d2"))
            numer = -(Math.cos(theta_end)*dfnx + Math.sin(theta_end)*dfny)*N43(t2)[3]
                  -  (Math.cos(theta_end)*(fnx - X) + Math.sin(theta_end)*(fny - Y))*dN43(t2)[3];
        //System.out.println("calc_t2dx2 = " + t1 + ", " + t2 + ", " + rhs1 + ", " + rhs2 + ", " + fprime + ", " + ((-3*t2*(1 - t2)*(1 - t2)*rhs1 - 3*(1 - t2)*(1 - 3*t2)*rhs2)/fprime));
        return -numer/denom;*/
        return Double.NaN;
    }
}
