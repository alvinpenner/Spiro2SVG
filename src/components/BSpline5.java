
package components;

import java.io.*;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a 5-point B-Spline (P0 - P4) to it, using parameter 0 < t2 < 2.
// fit the curvature at the endpoints and keep P2 arbitrary
// decompose the B-Spline into 2 Beziers, range (0,1) and (0,2).
// Bezier[2][4] = f[2](x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book2, Dec 2016, page 50

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\BSpline5.java

public class BSpline5
{
    public static double t1_start = 0;           // Math.PI/3;
    public static final double t1_end = Math.PI; // Math.PI/4;
    public static final int N = 100;
    public static double[] Splinex, Spliney;    // 5 point spline
    public static double[][] Bezx;              // 2 Beziers, 4 points each, x component
    public static double[][] Bezy;              // 2 Beziers, 4 points each, y component
    private static CycloidFxn fitted;       // = new CycloidFxn(.5);           // set c value
    //private static epiTrochoidFxn fitted; // = new epiTrochoidFxn(-2);       // set c value
    private static double[] t2 = new double[N+1];
    //private static double[] t2dd1 = new double[N+1];            // partial wrt d1
    //private static double[] t2dd2 = new double[N+1];            // partial wrt d2
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 80;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        read_data(80, -1);                          // cycloid data at angle theta
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //System.out.println("BSpline phi c t1_start t1_end  = ," + phi + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end);
        //System.out.println("solve_at_P2 = " + solve_at_P2(.8, 1.5, true));
    }

    private static void read_data(int theta, double tempc)
    {
        // read Bezier (d1, d2) data from a file
        // match field 1 = theta for a Cycloid / root number for an epiTrochoid
        // match field 2 = c for an epiTrochoid
        // convert to P2 for a 5 point B-Spline, and initiallize 'solve_at_P2()'

        String str, firstline = "";
        double d1 = 0;                  // Bezier arm length
        double d2 = 0;

        try
        {
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig5_hypocofmd1d2.csv"));
            BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig2_Cycloidcurvd1d2.csv"));
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig7_hypoODFd1d2.csv"));
            try
            {
                if (istr.ready())
                    firstline = istr.readLine();
                while (istr.ready())
                {
                    str = istr.readLine();
                    //System.out.println(str);
                    if (!str.isEmpty())
                        if ((firstline.startsWith("theta") && (Integer.parseInt(str.split(",")[0].trim()) == theta))
                        ||  (firstline.startsWith("root") && (Integer.parseInt(str.split(",")[0].trim()) == theta) && (Double.parseDouble(str.split(",")[1]) == tempc)))
                        {
                            d1 = Double.parseDouble(str.split(",")[3]);
                            d2 = Double.parseDouble(str.split(",")[4]);
                            break;
                        }
                }
            }
            catch (IOException e)
            {
                System.out.println("read error : " + e.getMessage());
                return;
            }
        }
        catch (FileNotFoundException e)
        {
            System.out.println("file not found : " + e.getMessage());
            return;
        }
        if (d1 == 0 && d2 == 0)
        {
            System.out.println("file data match not found, abort");
            return;
        }
        if (firstline.startsWith("theta"))
        {
            // calculate the point of maximum curvature of a cycloid from a tangent angle
            tempc = Math.sqrt(1 - .75*Math.cos(theta*Math.PI/180)*Math.cos(theta*Math.PI/180)); // override for cycloid only
            fitted = new CycloidFxn(tempc);                             // set c value
            t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        }
        //if (firstline.startsWith("root"))
        //    fitted = new epiTrochoidFxn(tempc);                           // set c value
        //d1 = 0.7641;
        //d2 = 1.8275;
        System.out.println("file data at theta c t d1 d2   = ," + theta + ", " + tempc + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2);

        // calculate P2 for a B-Spline

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
        P2x = 0.618;
        P2y = 1.342;
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x, P2y, true));
/*        double del = 0.001;
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x - del, P2y - del, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x - del, P2y, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x - del, P2y + del, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x, P2y - del, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x, P2y, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x, P2y + del, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x + del, P2y - del, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x + del, P2y, false));
        System.out.println("solve_at_P2 = " + solve_at_P2(P2x + del, P2y + del, false));
*/
        //iterate_at_d1_d2(d1, d2);     // default
    }

    private static boolean solve_at_P2(double x2, double y2, boolean print)
    {
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error
        // (calculate d1, d2 to satisfy the curvature at the endpoints)

        double d1, d2;
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);
        d1 = Math.sqrt(((y2 - fitted.gety(t1_start))*Math.cos(theta_start) - (x2 - fitted.getx(t1_start))*Math.sin(theta_start))/3/fitted.getkappa(t1_start));
        d2 = Math.sqrt(((y2 - fitted.gety(t1_end))*Math.cos(theta_end) - (x2 - fitted.getx(t1_end))*Math.sin(theta_end))/3/fitted.getkappa(t1_end));
        System.out.println("solve_at_P2                = " + x2 + ", " + y2 + ", " + fitted.getkappa(t1_start) + ", " + fitted.getkappa(t1_end) + ", " + d1 + ", " + d2);
        Splinex = new double[] {fitted.getx(t1_start),
                                fitted.getx(t1_start) + d1*Math.cos(theta_start),
                                x2,
                                fitted.getx(t1_end) - d2*Math.cos(theta_end),
                                fitted.getx(t1_end)};
        Spliney = new double[] {fitted.gety(t1_start),
                                fitted.gety(t1_start) + d1*Math.sin(theta_start),
                                y2,
                                fitted.gety(t1_end) - d2*Math.sin(theta_end),
                                fitted.gety(t1_end)};
        Bezx = new double[][] {{Splinex[0],
                                Splinex[1],
                                (Splinex[1] + Splinex[2])/2,
                                (Splinex[1] + 2*Splinex[2] + Splinex[3])/4},
                               {(Splinex[1] + 2*Splinex[2] + Splinex[3])/4,
                                (Splinex[2] + Splinex[3])/2,
                                Splinex[3],
                                Splinex[4]}};
        Bezy = new double[][] {{Spliney[0],
                                Spliney[1],
                                (Spliney[1] + Spliney[2])/2,
                                (Spliney[1] + 2*Spliney[2] + Spliney[3])/4},
                               {(Spliney[1] + 2*Spliney[2] + Spliney[3])/4,
                                (Spliney[2] + Spliney[3])/2,
                                Spliney[3],
                                Spliney[4]}};
        if (t2[N] == 0)
            System.out.println("__start at theta c t d1 d2 = , " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + calc_error());
        fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2");
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
                return false;
//            t2dd1[i] = calc_t2dd1(t1, t2[i]);
//            t2dd2[i] = calc_t2dd2(t1, t2[i]);
            if (print) System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || (t2[i] == Double.NaN))
            {
                System.out.println("abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
//                scan_quintic_near_t2(t1, t2[i]);
                return false;
            }
        }
        //System.out.println("new t2[] profile rms   = ," + d1 + ", " + d2 + ", " + calc_error());
        System.out.println("__new t2[] at theta c t d1 d2 rms = , " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + calc_error());
        return true;
    }

    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

        //double a_b = 180;         // scale factor to make rms error dimensionless
        double a_b = 1;             // Cycloid only
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
        if (fprime*fprime < 2*f*f2prime)
            del_t = -fprime/f2prime;
        else
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
        t += del_t;

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
            del_t = -f/fprime;
            t += del_t;
            loop++;
            //System.out.println("         t2 =, " + t + ", " + f + ", " + fprime);
        } while (Math.abs(del_t) > TOL);
        t2[i] = t + seg;               // compensate for Bezier segment offset
    }
}
