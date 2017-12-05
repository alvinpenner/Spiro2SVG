
package components;

import java.io.*;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a cubic Bezier to it, using parameter t2.
// Bezier = f(x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book2, Dec 2016, page 50

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\t2_vs_t1.java

public class t2_vs_t1
{
    public static double t1_start = 0; //Math.PI/3;         // Math.PI/3;
    public static final double t1_end = Math.PI/4; // 25*Math.PI/180; // Math.PI/4;
    public static final int N = 100;
    public static double[] Bezx;
    public static double[] Bezy;
    //private static CycloidFxn fitted;       // = new CycloidFxn(.5);           // set c value
    private static epiTrochoidFxn fitted; // = new epiTrochoidFxn(-2);       // set c value
    private static double[] t2 = new double[N+1];
    private static double[] t2dd1 = new double[N+1];            // partial wrt d1
    private static double[] t2dd2 = new double[N+1];            // partial wrt d2
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 80;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        //read_data(2, -3.54);                                       // initiallize the routine solve_at_d1_d2()
        //read_data(80, 0);
        fitted = new epiTrochoidFxn(0);
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //System.out.println("c t1_start t1_end      = ," + fitted.getc() + ", " + t1_start + ", " + t1_end);
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.460221195, 1.260856660, true));   // cofm, c = 0.5
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.494071922, 1.228369703, true));   // curv, c = 0.5
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.544262981, 2.050759223, true));   // cofm, c = 1
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.000000000, 2.309401076, true));   // curv, c = 1
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(-0.552, -0.552, true));             // circle, c = 1
        System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(47.73733159182497, 47.737331591589715, true));
        //iterate_at_d1_d2(0.494071922, 1.228369703);         // curv c = 0.5 init
        //iterate_at_d1_d2(0.460221195, 1.260856660);           // cofm c = 0.5 init
        //iterate_at_d1_d2(0.54426, 2.05076);                   // cofm, c = 1
        //iterate_at_d1_d2(0.5000, 2.30940);                    // curv, c = 1
        //iterate_at_d1_d2(-fitted.getc()*t1_end/3, -fitted.getc()*t1_end/3); // circle, c = 4.5
        //iterate_at_d1_d2(47, 47);
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(47.73733159182497, 47.737331591589715, false));
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(49.18949608004587, 46.032318702281785, true));
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.46858825099, 1.2567247398, true));   // cofm, c = 0.5 test test
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(9.000000 ,87.986060, true));   // epiTrochoid c = 5
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(10.325965305049323, 89.81941244540867, false));
        //grid_search_at_d1_d2(50.28654971690186, 38.97605103554598);
    }

    private static void iterate_at_d1_d2(double d1, double d2)
    {
        // calculate a new estimate of (d1, d2) by setting dF/dd1 = dF/dd2 = 0
        // see Spiro2SVG Book 3, pages 6-8
        // setup two elliptical equations for d1, d2

        final int MAXLOOP = 200;
        double A0, B0, C0, D0, E0, F0;
        double A1, B1, C1, D1, E1, F1;
        double A2, C2, D2, E2, F2;
        double t1, old_d1, old_d2;
        double[] h_gx = new double[N+1];
        double[] h_gy = new double[N+1];
        double[] dhx = new double[N+1];
        double[] dhy = new double[N+1];
        double[] ux = new double[N+1];
        double[] uy = new double[N+1];
        double[] vx = new double[N+1];
        double[] vy = new double[N+1];
        double[] dux = new double[N+1];
        double[] duy = new double[N+1];
        double[] dvx = new double[N+1];
        double[] dvy = new double[N+1];
        double[] trap_in = new double[N+1];
        int i, loop = 0;

        do
        {
            loop++;
            if (!solve_at_d1_d2(d1, d2, false))                    // initiallize at (d1, d2)
            {
                System.out.println("fail at " + d1 + ", " + d2);
                return;
            }

            for (i = 0; i <= N; i++)
            {
                t1 = t1_start + i*(t1_end - t1_start)/N;
                h_gx[i] = fnh(Bezx, t2[i]) - fitted.getx(t1);
                h_gy[i] = fnh(Bezy, t2[i]) - fitted.gety(t1);
                dhx[i] = dfnh(Bezx, t2[i]);
                dhy[i] = dfnh(Bezy, t2[i]);
                ux[i] = 3*Math.cos(theta_start)*t2[i]*(1 - t2[i])*(1 - t2[i]);
                uy[i] = 3*Math.sin(theta_start)*t2[i]*(1 - t2[i])*(1 - t2[i]);
                vx[i] = -3*Math.cos(theta_end)*t2[i]*t2[i]*(1 - t2[i]);
                vy[i] = -3*Math.sin(theta_end)*t2[i]*t2[i]*(1 - t2[i]);
                dux[i] = 3*Math.cos(theta_start)*(1 - t2[i])*(1 - 3*t2[i]);
                duy[i] = 3*Math.sin(theta_start)*(1 - t2[i])*(1 - 3*t2[i]);
                dvx[i] = -3*Math.cos(theta_end)*t2[i]*(2 - 3*t2[i]);
                dvy[i] = -3*Math.sin(theta_end)*t2[i]*(2 - 3*t2[i]);
            }

            // set dFdd1 = 0

            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dux[i] + uy[i]*duy[i])*t2dd1[i];
            A0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (vx[i]*dvx[i] + vy[i]*dvy[i])*t2dd1[i];
            B0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dvx[i] + dux[i]*vx[i] + uy[i]*dvy[i] + duy[i]*vy[i])*t2dd1[i];
            C0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = ux[i]*ux[i] + (ux[i]*dhx[i] + dux[i]*h_gx[i])*t2dd1[i]
                           + uy[i]*uy[i] + (uy[i]*dhy[i] + duy[i]*h_gy[i])*t2dd1[i];
            D0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = vx[i]*ux[i] + (vx[i]*dhx[i] + dvx[i]*h_gx[i])*t2dd1[i]
                           + vy[i]*uy[i] + (vy[i]*dhy[i] + dvy[i]*h_gy[i])*t2dd1[i];
            E0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = h_gx[i]*(ux[i] + dhx[i]*t2dd1[i])
                           + h_gy[i]*(uy[i] + dhy[i]*t2dd1[i]);
            F0 = integrate(trap_in);

            // set dFdd2 = 0

            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dux[i] + uy[i]*duy[i])*t2dd2[i];
            A1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (vx[i]*dvx[i] + vy[i]*dvy[i])*t2dd2[i];
            B1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dvx[i] + dux[i]*vx[i] + uy[i]*dvy[i] + duy[i]*vy[i])*t2dd2[i];
            C1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = ux[i]*vx[i] + (ux[i]*dhx[i] + dux[i]*h_gx[i])*t2dd2[i]
                           + uy[i]*vy[i] + (uy[i]*dhy[i] + duy[i]*h_gy[i])*t2dd2[i];
            D1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = vx[i]*vx[i] + (vx[i]*dhx[i] + dvx[i]*h_gx[i])*t2dd2[i]
                           + vy[i]*vy[i] + (vy[i]*dhy[i] + dvy[i]*h_gy[i])*t2dd2[i];
            E1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = h_gx[i]*(vx[i] + dhx[i]*t2dd2[i])
                           + h_gy[i]*(vy[i] + dhy[i]*t2dd2[i]);
            F1 = integrate(trap_in);

            // see: https://math.stackexchange.com/questions/1767225/algorithm-intersection-of-two-conics

            A2 = A0*B1 - A1*B0;
            C2 = C0*B1 - C1*B0;
            D2 = D0*B1 - D1*B0;
            E2 = E0*B1 - E1*B0;
            F2 = F0*B1 - F1*B0;
            old_d1 = d1;
            old_d2 = d2;
            d1 = fitymoment.solve_quartic(A0*C2*C2 + B0*A2*A2 - C0*C2*A2,
                                          2*A0*C2*E2 + D0*C2*C2 + 2*B0*A2*D2 - C0*C2*D2 - (C0*E2 + C2*E0)*A2,
                                          A0*E2*E2 + 2*D0*C2*E2 + F0*C2*C2 + 2*B0*A2*F2 + B0*D2*D2
                                          - C0*C2*F2 - (C0*E2 + C2*E0)*D2 - E0*E2*A2,
                                          D0*E2*E2 + 2*F0*C2*E2 + 2*B0*D2*F2 - (C0*E2 + C2*E0)*F2 - E0*E2*D2,
                                          F0*E2*E2 + B0*F2*F2 - E0*E2*F2,
                                          true);
            d2 = -(A2*d1*d1 + D2*d1 + F2)/(C2*d1 + E2);
            //System.out.println("new           d1 d2    = ," + d1 + ", " + d2
            //                 + ", " + (A0*d1*d1 + B0*d2*d2 + C0*d1*d2 + D0*d1 + E0*d2 + F0)
            //                 + ", " + (A1*d1*d1 + B1*d2*d2 + C1*d1*d2 + D1*d1 + E1*d2 + F1));

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                t2[i] += t2dd1[i]*(d1 - old_d1) + t2dd2[i]*(d2 - old_d2);   // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(old_d1 - d1) < TOL) && (Math.abs(old_d2 - d2) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 = , , , , , , " + d1 + ", " + d2 + ", " + Math.abs(old_d1 - d1) + ", " + Math.abs(old_d2 - d2));
            solve_at_d1_d2(d1, d2, true);                  // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops!");
    }

    private static void grid_search_at_d1_d2(double d1, double d2)
    {
        // calculate rms error in neighbouring d1, d2 region
        double del = 0.0001;

        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
            {
                System.out.println("grid =, " + i + ", " + j);
                solve_at_d1_d2(d1 + i*del, d2 + j*del, false);
            }
    }

    private static void grid_at_d1_d2(double d1, double d2)
    {
        // calculate rms error in neighbouring d1, d2 region
        // obsolete code, response to d1/d2 change is not optimized!
        double del_d = 0.0001; // 0.000002;
        double tempd1, tempd2;

        System.out.println("\nd1\\d2, " + (d2 - del_d) + ", " + d2 + ", " + (d2 + del_d));
        for (int i = -1; i < 2; i++)
        {
            tempd1 = d1 + i*del_d;
            System.out.print(tempd1 + ", ");
            for (int j = -1; j < 2; j++)
            {
                tempd2 = d2 + j*del_d;
                Bezx = new double[] {fitted.getx(t1_start),
                                     fitted.getx(t1_start) + tempd1*Math.cos(theta_start),
                                     fitted.getx(t1_end) - tempd2*Math.cos(theta_end),
                                     fitted.getx(t1_end)};
                Bezy = new double[] {fitted.gety(t1_start),
                                     fitted.gety(t1_start) + tempd1*Math.sin(theta_start),
                                     fitted.gety(t1_end) - tempd2*Math.sin(theta_end),
                                     fitted.gety(t1_end)};
                System.out.print(calc_error() + ", ");
            }
            System.out.println();
        }
    }

    private static void read_data(int ln, double tempc)
    {
        // match field 1==ln from a data file, and initiallize 'solve_at_d1_d2()'
        // match field 2 = c for an epiTrochoid

        String str;
        double d1 = 0;
        double d2 = 0;

        try
        {
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig5_hypocofmd1d2.csv"));
            BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig2_Cycloidcofmd1d2.csv"));
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\ODFPaper\\Fig7_hypoODFd1d2.csv"));
            try
            {
                while (istr.ready())
                {
                    str = istr.readLine();
                    //System.out.println(str);
                    if (!str.isEmpty() && !str.startsWith("theta") && !str.startsWith("root"))
                        if ((Integer.parseInt(str.split(",")[0].trim()) == ln))
                        //&&  (Double.parseDouble(str.split(",")[1]) == tempc))
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
        // extract the point of maximum curvature of a cycloid from a tangent angle
        //tempc = Math.sqrt(1 - .75*Math.cos(ln*Math.PI/180)*Math.cos(ln*Math.PI/180)); // override for cycloid only
        //fitted = new epiTrochoidFxn(tempc);                           // set c value
        //fitted = new epiTrochoidFxn(-1.5);       // override c value
        //fitted = new CycloidFxn(tempc);                             // set c value
        t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        System.out.println("file data at theta c t d1 d2   = ," + ln + ", " + tempc + ", " + t1_start + ", " + d1 + ", " + d2);
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(d1, d2, true));
        //iterate_at_d1_d2(9, 75);       // override default
        //iterate_at_d1_d2(78.9, 7.99);       // override default
        iterate_at_d1_d2(d1, d2);     // default
    }

    private static boolean solve_at_d1_d2(double d1, double d2, boolean print)
    {
        // perform a single calculation of a complete t2[] profile
        // at a given (d1, d2), and calculate the rms error

        double t1;
        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);
        Bezx = new double[] {fitted.getx(t1_start),
                             fitted.getx(t1_start) + d1*Math.cos(theta_start),
                             fitted.getx(t1_end) - d2*Math.cos(theta_end),
                             fitted.getx(t1_end)};
        Bezy = new double[] {fitted.gety(t1_start),
                             fitted.gety(t1_start) + d1*Math.sin(theta_start),
                             fitted.gety(t1_end) - d2*Math.sin(theta_end),
                             fitted.gety(t1_end)};
        if (t2[N] == 0)
            System.out.println("__start cubic Bezier at theta c t d1 d2        = , " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + calc_error());
        //fitted.gen_Bezier(new double[] {Bezx[0], Bezy[0], Bezx[1], Bezy[1], Bezx[2], Bezy[2], Bezx[3], Bezy[3]});
        //System.out.println("theta = " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI);
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);

        if (print) System.out.println("\n t1, t2, t2dd1, t2dd2");
        for (int i = 0; i <= N; i++)
        {
            t1 = t1_start + i*(t1_end - t1_start)/N;
            if (i == 0)
                t2[i] = solve_quintic_for_t2(t1, 0);
            else
                t2[i] = solve_quintic_for_t2(t1, t2[i-1]);
            t2dd1[i] = calc_t2dd1(t1, t2[i]);
            t2dd2[i] = calc_t2dd2(t1, t2[i]);
            if (print) System.out.println(t1 + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i]);
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(1 - t2[i]) > TOL)
            ||  (t2[i] < -TOL) || (t2[i] == Double.NaN))
            {
                System.out.println("abort at " + i + " : " + t1 + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i]);
                scan_quintic_near_t2(t1, t2[i]);
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

        double a_b = 180;                   // scale factor to make rms error dimensionless
        //double a_b = 1;                     // Cycloid only
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
            System.out.println(i + ", " + (fn(Bezx, t2[i]) - fitted.getx(t1)) + ", " + (fn(Bezy, t2[i]) - fitted.gety(t1)) + ", " + Math.sqrt(trap_in[i]));
            t1 += (t1_end - t1_start)/N;
        }
        return Math.sqrt(integrate(trap_in))/a_b;
    }

    public static double integrate(double[] trap)
    {
        // trapezoidal rule integration of a fxn of t1 (N+1 points)

        //System.out.println("trap length = " + trap.length);
        double ret = (trap[0] + trap[trap.length - 1])/2;
        for (int i = 1; i < trap.length - 1; i++)
            ret += trap[i];
        return ret/(trap.length - 1);
    }

    private static double solve_quintic_for_t2(double t1, double t2_init)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the cubic Bezier t-value

        double f, fprime, f2prime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        //double t = t2_init + 1.0/N + .01; // previous t2 value + linear increment as an initial estimate
        double t = t2_init;
        int loop = 0;

        // initial estimate using quadratic approximation

        f = (fn(Bezx, t) - X)*dfn(Bezx, t) + (fn(Bezy, t) - Y)*dfn(Bezy, t);
        fprime = dfn(Bezx, t)*dfn(Bezx, t) + (fn(Bezx, t) - X)*d2fn(Bezx, t) + dfn(Bezy, t)*dfn(Bezy, t) + (fn(Bezy, t) - Y)*d2fn(Bezy, t);
        f2prime = 3*dfn(Bezx, t)*d2fn(Bezx, t) + (fn(Bezx, t) - X)*d3fn(Bezx, t) + 3*dfn(Bezy, t)*d2fn(Bezy, t) + (fn(Bezy, t) - Y)*d3fn(Bezy, t);
        if (fprime*fprime < 2*f*f2prime)
            del_t = -fprime/f2prime;
        else
            del_t = (-fprime + Math.sqrt(fprime*fprime - 2*f*f2prime))/f2prime;
        t += del_t;

        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
        do
        {
            f = (fn(Bezx, t) - X)*dfn(Bezx, t) + (fn(Bezy, t) - Y)*dfn(Bezy, t);
            fprime = dfn(Bezx, t)*dfn(Bezx, t) + (fn(Bezx, t) - X)*d2fn(Bezx, t) + dfn(Bezy, t)*dfn(Bezy, t) + (fn(Bezy, t) - Y)*d2fn(Bezy, t);
            if (loop > 100)
                return Double.NaN;
            del_t = -f/fprime;
            t += del_t;
            loop++;
            //System.out.println("         t2 =, " + t + ", " + f + ", " + fprime);
        } while (Math.abs(del_t) > TOL);
        return t;
    }

    private static void scan_quintic_near_t2(double t1, double t2_bad)
    {
        // if solve_quintic_for_t2 fails, scan the area for other roots
        double f, t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        System.out.println("\nscanning at " + t1 + ", " + t2_bad);
        for (int i = -20; i <= 200; i++)
        {
            t = -i*t2_bad/10;
            f = (fn(Bezx, t) - X)*dfn(Bezx, t) + (fn(Bezy, t) - Y)*dfn(Bezy, t);
            System.out.println(t + ", " + f);
        }
        System.out.println();
    }

    private static double calc_t2dd1(double t1, double t2)
    {
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double rhs1 = Math.cos(theta_start)*dfn(Bezx, t2) + Math.sin(theta_start)*dfn(Bezy, t2);
        double rhs2 = Math.cos(theta_start)*(fn(Bezx, t2) - X) + Math.sin(theta_start)*(fn(Bezy, t2) - Y);
        double fprime = dfn(Bezx, t2)*dfn(Bezx, t2) + (fn(Bezx, t2) - X)*d2fn(Bezx, t2) + dfn(Bezy, t2)*dfn(Bezy, t2) + (fn(Bezy, t2) - Y)*d2fn(Bezy, t2);

        //System.out.println("calc_t2dd1 = " + t1 + ", " + t2 + ", " + rhs1 + ", " + rhs2 + ", " + fprime + ", " + ((-3*t2*(1 - t2)*(1 - t2)*rhs1 - 3*(1 - t2)*(1 - 3*t2)*rhs2)/fprime));
        return (-3*t2*(1 - t2)*(1 - t2)*rhs1 - 3*(1 - t2)*(1 - 3*t2)*rhs2)/fprime;
    }

    private static double calc_t2dd2(double t1, double t2)
    {
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double rhs1 = Math.cos(theta_end)*dfn(Bezx, t2) + Math.sin(theta_end)*dfn(Bezy, t2);
        double rhs2 = Math.cos(theta_end)*(fn(Bezx, t2) - X) + Math.sin(theta_end)*(fn(Bezy, t2) - Y);
        double fprime = dfn(Bezx, t2)*dfn(Bezx, t2) + (fn(Bezx, t2) - X)*d2fn(Bezx, t2) + dfn(Bezy, t2)*dfn(Bezy, t2) + (fn(Bezy, t2) - Y)*d2fn(Bezy, t2);

        return (3*t2*t2*(1 - t2)*rhs1 + 3*t2*(2 - 3*t2)*rhs2)/fprime;
    }

    public static double fn(double[] Bez, double t)
    {
        return Bez[0]*(1 - t)*(1 - t)*(1 - t) + 3*Bez[1]*t*(1 - t)*(1 - t) + 3*Bez[2]*t*t*(1 - t) + Bez[3]*t*t*t;
    }

    public static double dfn(double[] Bez, double t)
    {
        return -3*Bez[0]*(1 - t)*(1 - t) + 3*Bez[1]*(1 - t)*(1 - 3*t) + 3*Bez[2]*t*(2 - 3*t) + 3*Bez[3]*t*t;
    }

    public static double d2fn(double[] Bez, double t)
    {
        return 6*Bez[0]*(1 - t) + 6*Bez[1]*(-2 + 3*t) + 6*Bez[2]*(1 - 3*t) + 6*Bez[3]*t;
    }

    public static double d3fn(double[] Bez, double t)
    {
        return -6*Bez[0] + 18*Bez[1] - 18*Bez[2] + 6*Bez[3];
    }

    private static double fnh(double[] Bez, double t)
    {
        // this is fn with assumed degenerate end-points
        return Bez[0]*(1 - t)*(1 - t)*(1 - t) + 3*Bez[0]*t*(1 - t)*(1 - t) + 3*Bez[3]*t*t*(1 - t) + Bez[3]*t*t*t;
    }

    private static double dfnh(double[] Bez, double t)
    {
        // this is dfn with assumed degenerate end-points
        return -3*Bez[0]*(1 - t)*(1 - t) + 3*Bez[0]*(1 - t)*(1 - 3*t) + 3*Bez[3]*t*(2 - 3*t) + 3*Bez[3]*t*t;
    }
}
