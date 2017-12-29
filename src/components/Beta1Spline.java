
package components;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1.
// fit a 5-point Beta1-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// constrain only the slopes at the endpoints and keep P2 (Q3/Q4) arbitrary.
// the Beta1 spline has an asymmetric arm length, d- and d+, at the Bezier splice at Q3/Q4.
// linearize the equations wrt (d1, d2, x2, y2, beta1) and solve a 5x5 system of equations.
// decompose the Beta1-Spline into 2 Beziers, range (0,1) and (1,2).
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2).
// t2 must be chosen to minimize the distance to g(t1).
// see Spiro2SVG Book 4, Dec 2017, page 24.

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\Beta1Spline.java

public class Beta1Spline
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
//    private static double[] t2dd1 = new double[N+1];            // partial wrt d1
//    private static double[] t2dd2 = new double[N+1];            // partial wrt d2
//    private static double[] t2dx2 = new double[N+1];            // partial wrt x2
//    private static double[] t2dy2 = new double[N+1];            // partial wrt y2
//    private static double[] t2dd = new double[N+1];             // partial wrt symmetric d
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 60;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        fitted = new epiTrochoidFxn(10);
        //System.out.println("Beta1-Spline solve_at_P2 = " + convert_at_P2(0.5027981574023782, 0.6554626509923203, 1.2475458906992178, 1.5807488348561327, true) + "\n");
        System.out.println("Beta1-Spline solve_at_P2 = " + convert_at_P2(19.983314292966483, 26.42763336958588, 175.47633731103565, 59.05668195284478, true) + "\n");
        //System.out.println("Beta1-Spline solve_at_P2 = " + convert_at_P2(15, 20, 170, 59, true) + "\n");
        //System.out.println("Beta1-Spline convert_at_P2 = " + convert_at_P2(23.264222028261724, 23.34619627704507, 171.41612180193727, 67.26310370327987, true) + "\n");
        //System.out.println("Beta1-Spline iterate_at_P2 = " + iterate_at_P2( 22.8770976550121,24.818420825802733, 165.94185238942046, 69.55892810479992, 23.390620200745754) + "\n");
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

        //return iterate_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, 1);
        return solve_at_P2(d1, d2, (x1 + 2*x2 + x3)/4, (y1 + 2*y2 + y3)/4, 1, print);
    }

    private static double solve_at_P2(double d1, double d2, double x2, double y2, double beta1, boolean print)
    {
        // 5-point cubic Beta1-Spline curve, decomposed as two Beziers
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);

        beta1 += 0.01;
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

        double delx = (beta1*beta1*(Bezx[0][3] - Bezx[0][1]) + Bezx[1][2] - Bezx[1][0])/2/(beta1 + 1);
        double dely = (beta1*beta1*(Bezy[0][3] - Bezy[0][1]) + Bezy[1][2] - Bezy[1][0])/2/(beta1 + 1);
        Bezx[0][2] -= delx/beta1;
        Bezx[1][1] += delx;
        Bezy[0][2] -= dely/beta1;
        Bezy[1][1] += dely;
        //double dmins = Math.sqrt((Bezx[0][3] - Bezx[0][2])*(Bezx[0][3] - Bezx[0][2]) + (Bezy[0][3] - Bezy[0][2])*(Bezy[0][3] - Bezy[0][2]));
        //double dplus = Math.sqrt((Bezx[1][1] - Bezx[1][0])*(Bezx[1][1] - Bezx[1][0]) + (Bezy[1][1] - Bezy[1][0])*(Bezy[1][1] - Bezy[1][0]));
        //System.out.println("beta1 = ," + beta1 + ", " + dmins + ", " + dplus);

        if (t2[N] == 0)
            System.out.println("__start Beta1-Spline at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + (float) t1_start + ", " + (float) t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1 + ", " + calc_error());
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
            //t2dd1[i] = calc_t2dxy(i, seg, t2[i] - seg, "d1");
            //t2dd2[i] = calc_t2dxy(i, seg, t2[i] - seg, "d2");
            //t2dx2[i] = calc_t2dxy(i, seg, t2[i] - seg, "x2");
            //t2dy2[i] = calc_t2dxy(i, seg, t2[i] - seg, "y2");
            //t2dd[i] = calc_t2dxy(i, seg, t2[i] - seg, "d");
            if (print)
                System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i] + ", " + t2dx2[i] + ", " + t2dy2[i] + ", " + t2dd[i]);
        }
        double retVal = calc_error();
        System.out.println("__new t2[] at theta c t d1 d2 rms = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + beta1 + ", " + retVal);
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
/*
    private static double calc_t2dxy(int i, int seg, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        double denom = t2_vs_t1.dfn(Bezx[seg], t2)*t2_vs_t1.dfn(Bezx[seg], t2) + (t2_vs_t1.fn(Bezx[seg], t2) - X)*t2_vs_t1.d2fn(Bezx[seg], t2) + t2_vs_t1.dfn(Bezy[seg], t2)*t2_vs_t1.dfn(Bezy[seg], t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*t2_vs_t1.d2fn(Bezy[seg], t2);
        double numer = Double.NaN;
        if (type.equals("x2"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdx2(seg, t2) + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudx2(seg, t2);
        else if (type.equals("y2"))
            numer = t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydy2(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudy2(seg, t2);
        else if (type.equals("d1"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdd1(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydd1(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudd1(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudd1(seg, t2);
        else if (type.equals("d2"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdd2(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydd2(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudd2(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudd2(seg, t2);
        else if (type.equals("d"))
            numer = t2_vs_t1.dfn(Bezx[seg], t2)*calc_dfxdd(seg, t2) + t2_vs_t1.dfn(Bezy[seg], t2)*calc_dfydd(seg, t2)
                  + (t2_vs_t1.fn(Bezx[seg], t2) - X)*calc_d2fxdudd(seg, t2) + (t2_vs_t1.fn(Bezy[seg], t2) - Y)*calc_d2fydudd(seg, t2);
        return -numer/denom;
    }

    private static double calc_dfxdx2(int seg, double t2)
    {
        return N33(t2)[2 - 2*seg] + N33(t2)[3 - 2*seg];
    }

    private static double calc_d2fxdudx2(int seg, double t2)
    {
        return dN33(t2)[2 - 2*seg] + dN33(t2)[3 - 2*seg];
    }

    private static double calc_dfydy2(int seg, double t2)
    {
        return N33(t2)[2 - 2*seg] + N33(t2)[3 - 2*seg];
    }

    private static double calc_d2fydudy2(int seg, double t2)
    {
        return dN33(t2)[2 - 2*seg] + dN33(t2)[3 - 2*seg];
    }

    private static double calc_dfxdd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd1 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_start) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_start))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.cos(theta_start)*N33(t2)[1];
        return retVal + sgn*dpsidd1*d*Math.sin(psi)*N33(t2)[2 - seg];
    }

    private static double calc_dfydd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd1 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_start) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_start))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.sin(theta_start)*N33(t2)[1];
        return retVal - sgn*dpsidd1*d*Math.cos(psi)*N33(t2)[2 - seg];
    }

    private static double calc_d2fxdudd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd1 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_start) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_start))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.cos(theta_start)*dN33(t2)[1];
        return retVal + sgn*dpsidd1*d*Math.sin(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fydudd1(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd1 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_start) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_start))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg == 0)
            retVal = Math.sin(theta_start)*dN33(t2)[1];
        return retVal - sgn*dpsidd1*d*Math.cos(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_dfxdd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd2 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_end) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_end))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.cos(theta_end)*N33(t2)[2];
        return retVal + sgn*dpsidd2*d*Math.sin(psi)*N33(t2)[2 - seg];
    }

    private static double calc_dfydd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd2 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_end) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_end))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.sin(theta_end)*N33(t2)[2];
        return retVal - sgn*dpsidd2*d*Math.cos(psi)*N33(t2)[2 - seg];
    }

    private static double calc_d2fxdudd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd2 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_end) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_end))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.cos(theta_end)*dN33(t2)[2];
        return retVal + sgn*dpsidd2*d*Math.sin(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fydudd2(int seg, double t2)
    {
        double retVal = 0;
        double d = Math.sqrt((Bezx[0][2] - Bezx[0][3])*(Bezx[0][2] - Bezx[0][3]) + (Bezy[0][2] - Bezy[0][3])*(Bezy[0][2] - Bezy[0][3]));
        double dpsidd2 = (-(Bezx[1][2] - Bezx[0][1])*Math.sin(theta_end) + (Bezy[1][2] - Bezy[0][1])*Math.cos(theta_end))/((Bezx[1][2] - Bezx[0][1])*(Bezx[1][2] - Bezx[0][1]) + (Bezy[1][2] - Bezy[0][1])*(Bezy[1][2] - Bezy[0][1]));
        double sgn = 1;

        if (seg != 0) sgn = -1;
        if (seg != 0)
            retVal = -Math.sin(theta_end)*dN33(t2)[2];
        return retVal - sgn*dpsidd2*d*Math.cos(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_dfxdd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.cos(psi)*N33(t2)[2 - seg];
    }

    private static double calc_dfydd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.sin(psi)*N33(t2)[2 - seg];
    }

    private static double calc_d2fxdudd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.cos(psi)*dN33(t2)[2 - seg];
    }

    private static double calc_d2fydudd(int seg, double t2)
    {
        double sgn = 1;

        if (seg != 0) sgn = -1;
        return -sgn*Math.sin(psi)*dN33(t2)[2 - seg];
    }

    private static double[] N33(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:   N33 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:   N33 too large u value = " + u);
        return new double[] {(1 - u)*(1 - u)*(1 - u),
                             3*u*(1 - u)*(1 - u),
                             3*u*u*(1 - u),
                             u*u*u};
    }

    private static double[] dN33(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:  dN33 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:  dN33 too large u value = " + u);
        return new double[] {-3*(1 - u)*(1 - u),
                              3*(1 - u)*(1 - 3*u),
                              3*u*(2 - 3*u),
                              3*u*u};
    }

    private static double[] d2N33(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING: d2N33 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING: d2N33 too large u value = " + u);
        return new double[] { 6*(1 - u),
                              6*(-2 + 3*u),
                              6*(1 - 3*u),
                              6*u};
    }
*/
 }
