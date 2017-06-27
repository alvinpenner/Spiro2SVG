
package components;

import java.awt.geom.Point2D;

// consider a parametric curve, g, either cycloid or trochoid, with parameter t1
// fit a cubic Bezier to it, using parameter t2.
// Bezier = f(x0, x1, x2, x3, t2)
// t2 must be chosen to minimize the distance to g(t1)
// see Spiro2SVG Book2, Dec 2016, page 50

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\t2_vs_t1.java

public class t2_vs_t1
{
    public static final double t1_start = 0; // Math.PI/3;        // Math.PI/3;
    public static final double t1_end = Math.PI/2;    // Math.PI;
    public static final int N = 100;
    public static double[] Bezx;
    public static double[] Bezy;
    //private static CycloidFxn fitted = new CycloidFxn(1);       // set c value
    private static CircleFxn fitted = new CircleFxn(1);           // set c value
    private static double[] t2 = new double[N+1];
    private static double[] t2dd1 = new double[N+1];            // partial wrt d1
    private static double[] t2dd2 = new double[N+1];            // partial wrt d2
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        System.out.println("c t1_start t1_end      = ," + fitted.getc() + ", " + t1_start + ", " + t1_end);
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.460221195, 1.260856660, true));   // cofm, c = 0.5
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.494071922, 1.228369703, true));   // curv, c = 0.5
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.544262981, 2.050759223, true));   // cofm, c = 1
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.000000000, 2.309401076, true));     // curv, c = 1
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(-0.552, -0.552, true));                       // circle, c = 1
        //iterate_at_d1_d2(0.494071922, 1.228369703);         // curv c = 0.5 init
        //iterate_at_d1_d2(0.460221195, 1.260856660);           // cofm c = 0.5 init
        //iterate_at_d1_d2(0.54426, 2.05076);                   // cofm, c = 1
        //iterate_at_d1_d2(0.5000, 2.30940);                    // curv, c = 1
        iterate_at_d1_d2(-0.5, -0.5);                           // circle, c = 1
        //System.out.println("solve_at_d1_d2 = " + solve_at_d1_d2(0.46858825099, 1.2567247398, true));   // cofm, c = 0.5 test test
    }

    private static void iterate_at_d1_d2(double d1, double d2)
    {
        // calculate a new estimate of (d1, d2) by setting dF/dd1 = dF/dd2 = 0
        // see Spiro2SVG Book 3, pages 6-8
        // setup two elliptical equations for d1, d2

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
            System.out.println("new           d1 d2    = ," + d1 + ", " + d2
                             + ", " + (A0*d1*d1 + B0*d2*d2 + C0*d1*d2 + D0*d1 + E0*d2 + F0)
                             + ", " + (A1*d1*d1 + B1*d2*d2 + C1*d1*d2 + D1*d1 + E1*d2 + F1));
        } while ((loop < 10) && !((Math.abs(old_d1 - d1) < TOL) && (Math.abs(old_d2 - d2) < TOL)));
        if (loop < 10)
            System.out.println("\nconverged new d1 d2    = ," + d1 + ", " + d2);
        else
            System.out.println("\nNOT converged !");
    }

    private static void grid_at_d1_d2(double d1, double d2)
    {
        // calculate rms error in neighbouring d1, d2 region
        double del_d = 0.000002;
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

    private static boolean solve_at_d1_d2(double d1, double d2, boolean print)
    {
        // perform a single calculation of a complete t2[] profile
        // at a given (d1, d2), and calculate the rms error

        double t1;
        theta_start = Math.atan(fitted.getdydx(t1_start));
        theta_end = Math.atan(fitted.getdydx(t1_end));
        Bezx = new double[] {fitted.getx(t1_start),
                             fitted.getx(t1_start) + d1*Math.cos(theta_start),
                             fitted.getx(t1_end) - d2*Math.cos(theta_end),
                             fitted.getx(t1_end)};
        Bezy = new double[] {fitted.gety(t1_start),
                             fitted.gety(t1_start) + d1*Math.sin(theta_start),
                             fitted.gety(t1_end) - d2*Math.sin(theta_end),
                             fitted.gety(t1_end)};
        if (t2[N] == 0)
            System.out.println("start at     d1 d2     = ," + d1 + ", " + d2);
        else
            System.out.println("solve at new d1 d2 rms = ," + d1 + ", " + d2 + ", " + calc_error());
        Point2D.Double[] ptBez = new Point2D.Double[] {new Point2D.Double(Bezx[0], Bezy[0]), new Point2D.Double(Bezx[1], Bezy[1]), new Point2D.Double(Bezx[2], Bezy[2]), new Point2D.Double(Bezx[3], Bezy[3])};
        //fitCycloid.gen_Bezier(ptBez);
        //System.out.println("theta = " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI);
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);

        if (print) System.out.println("\n t1 , t2");
        for (int i = 0; i <= N; i++)
        {
            t1 = t1_start + i*(t1_end - t1_start)/N;
            t2[i] = solve_quintic_for_t2(t1);
            t2dd1[i] = calc_t2dd1(t1, t2[i]);
            t2dd2[i] = calc_t2dd2(t1, t2[i]);
            if (print) System.out.println(t1 + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i]);
            if (t2[i] == Double.NaN) return false;
        }
        System.out.println("new t2[] profile rms   = ," + d1 + ", " + d2 + ", " + calc_error());

/*        if (shoot_t2(clst2))    // generate profile of t2 versus t1 with boundary conditions
        {
            System.out.println("rms err at d1 d2 = ," + d1 + ", " + d2 + ", " + calc_error());
            if (iloop == 200)
            {
                System.out.println("shooting profile rms err at d1 d2 = ," + d1 + ", " + d2 + ", " + calc_error());
                for (int i = 0; i <= N; i++)
                {
                    double t1 = t1_start + i*(t1_end - t1_start)/N;
                    System.out.println(t1 + ", " + clst2.data[i] + ", " +
                    Math.sqrt((fn(Bezx, clst2.data[i]) - fitted.getx(t1))*(fn(Bezx, clst2.data[i]) - fitted.getx(t1))
                           +  (fn(Bezy, clst2.data[i]) - fitted.gety(t1))*(fn(Bezy, clst2.data[i]) - fitted.gety(t1))));
                }
                // generate quintic solution for t2 profile
                //for (int i = 0; i <= N; i++)        // fix fix temporary usage of clst2.data[]
                //    clst2.data[i] = solve_quintic_for_t2(t1_start + i*(t1_end - t1_start)/N, 0.5); // init used to be t2, didn't work
                //System.out.println("quintic profile rms err at d1 d2 = ," + d1 + ", " + d2 + ", " + calc_error(iloop));
                for (int i = 0; i <= N; i++)
                {
                    double t1 = t1_start + i*(t1_end - t1_start)/N;
                    System.out.println(t1 + ", " + clst2.data[i] + ", " +
                    Math.sqrt((fn(Bezx, clst2.data[i]) - fitted.getx(t1))*(fn(Bezx, clst2.data[i]) - fitted.getx(t1))
                           +  (fn(Bezy, clst2.data[i]) - fitted.gety(t1))*(fn(Bezy, clst2.data[i]) - fitted.gety(t1))));
                }
            }
            //System.out.println("............................................");
            //for (int i = 0; i <= N; i++)
            //    System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + clst2.data[i] + ", " + clst2dd1.fxn(t1_start + i*(t1_end - t1_start)/N, clst2.data[i]));
            //System.out.println("--------------------------------------------");
            //shoot_t2(clst2dd1);
            return true;
        } */
        return true;
    }
/*
    private static boolean shoot_t2(BoundaryFxn cls)
    {
        // satisfy a boundary condition on a function t2.
        // iterate to calculate an initial condition t2_init at index 1
        // in order to satisfy a final condition at the endpoint index N

        double t2_init_old, t2_init_new, t2_temp;
        double t2_final_old, t2_final_new;
        int loop = 0;

        t2_init_old = (t1_end - t1_start)/N*cls.fxn(t1_start, 0);
        //t2_init_old = 0.0084;     // curvature fit at c = 1 (override default)
        //t2_init_old = 0.011104627543711588;     // fix fix test code
        //t2_init_old = 0;
        t2_final_old = scan_t2_vs_t1(cls, t2_init_old, false, true);
        System.out.println("\nt2_init_old = " + t2_init_old + ", " + t2_final_old);
        t2_init_new = t2_init_old + .0001;
        t2_final_new = scan_t2_vs_t1(cls, t2_init_new, false, true);
        System.out.println("t2_init_new = " + t2_init_new + ", " + t2_final_new);
        do
        {
            loop++;
            t2_temp = t2_init_new + (cls.t2_final_set - t2_final_new)*(t2_init_new - t2_init_old)/(t2_final_new - t2_final_old);
            t2_final_old = t2_final_new;
            t2_final_new = scan_t2_vs_t1(cls, t2_temp, false, true);
            t2_init_old = t2_init_new;
            t2_init_new = t2_temp;
            System.out.println("t2_init_new = " + t2_init_new + ", " + t2_final_new);
        } while ((loop < 10) && (Math.abs(cls.t2_final_set - t2_final_new) > TOL));
        if (loop == 10)
        {
            System.out.println("too many loops, aborting!");
            return false;
        }
        else if (Double.isNaN(t2_final_new))
        {
            System.out.println("Bad data, aborting!");
            return false;
        }
        scan_t2_vs_t1(cls, t2_init_new, false, true);
        return true;
    }

    private static double scan_t2_vs_t1(BoundaryFxn cls, double t2_init, boolean print, boolean save)
    {
        // calculate t2 from a non-linear differential equation
        // assuming a given initial condition t2_init at index 1

        double t1 = t1_start;
        double t = t2_init;             // t2 at index 1, to be incremented and returned
        double t2_old = 0;
        double temp;                    // temporary storage of new t2

        if (print) System.out.println("\n t1 , " + cls.getClass().getSimpleName());
        if (print) System.out.println(t1 + ", 0");

        t1 += (t1_end - t1_start)/N;
        if (print) System.out.println(t1 + ", " + t);
        if (save) cls.data[1] = t;
        for (int i = 1; i < N; i++)
        {
            //temp = t2_old + 2*(t1_end - t1_start)/N*cls.fxn(t1, t);
            temp = t2_old + 2*(t1_end - t1_start)/N*cls.fxn(t1, clst2.data[i]);
            t2_old = t;
            t = temp;
            t1 += (t1_end - t1_start)/N;
            if (print) System.out.println(t1 + ", " + t);   // + ", " + clst2.data[i]);
            if (save) cls.data[i + 1] = t;      // save function data
        }
        return t;
    }
*/
    private static double calc_error()
    {
        // calculate rms error function assuming the error is zero at the endpoints
        // and assuming t2[i] is known

        double t1 = t1_start;
        double[] trap_in = new double[N+1];

        if ((Math.abs(1 - t2[N]) > TOL) || (Math.abs(t2[0]) > TOL))
        {
            System.out.println("calc_error() t2 sequence not bounded correctly, abort");
            return Double.NaN;
        }
        for (int i = 0; i <= N; i++)
        {
            trap_in[i] = (fn(Bezx, t2[i]) - fitted.getx(t1))*(fn(Bezx, t2[i]) - fitted.getx(t1))
                       + (fn(Bezy, t2[i]) - fitted.gety(t1))*(fn(Bezy, t2[i]) - fitted.gety(t1));
            t1 += (t1_end - t1_start)/N;
        }
        return Math.sqrt(integrate(trap_in));
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

    private static double solve_quintic_for_t2(double t1)
    {
        // calculate t2 at a known, fixed value of t1 (Newton-Raphson)
        // t = initial estimate of t2, the cubic Bezier t-value

        double f, fprime, del_t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double t = 0.5;                             // initial estimate of t2
        int loop = 0;

        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t);
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
