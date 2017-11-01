
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
        read_data(10, -1);                          // cycloid data at angle theta
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        //System.out.println("BSpline phi c t1_start t1_end  = ," + phi + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end);
        //System.out.println("solve_at_P2 = " + solve_at_P2(.8, 1.5, true));
        //for (int i = 0; i <= 200; i++)
        //    System.out.println(i + ", " + d2N43(i/100.0)[0] + ", " + d2N43(i/100.0)[1] + ", " + d2N43(i/100.0)[2] + ", " + d2N43(i/100.0)[3] + ", " + d2N43(i/100.0)[4]);
        //System.out.println("test mmult = " + mmult(Spliney, N43(1.7)));

/*        // exercise matrix functions, to be deleted
        double[] v1 = new double[] {1.1, 3.2, 5.7, 6.8, 9.3};
        double[] v2 = new double[] {11.1, 3.7, 7.2, 9.1, 10.0};
        double[][] m = new double[][] {{1.7, 7.2, 9.7, 12.8, -9.3},
                                       {2.1, 4.2, 5.7,  5.1, 21.4},
                                       {4.1, 5.3, 7.0,  6.8, -3.2}};
        //System.out.println(multvv(v1, v2));
        double[] ret = multmv(m, v2);
        System.out.println(ret[0] + ", " + ret[1] + ", " + ret[2]);
        m = new double[][] {{1.7, 7.2, 9.7, 12.8},
                            {2.1, 4.2, 5.7,  5.1},
                            {4.1, 5.3, 7.0,  6.8},
                            {8.1, 7.3, 9.0,  2.8}};
        System.out.println("detm   = " + (float) detm(m, 2, 1));
        double[][] minv = invertm(m);
        System.out.println("invert = " + minv.toString());
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
                System.out.print((float)minv[i][j] + ", ");
            System.out.println();
        }
*/
    }

    private static void iterate_at_P2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        // calculate a new estimate of (x2, y2) by setting dF/dx2 = dF/dy2 = 0
        // fit curvature at endpoints to determine d1, d2 as fxns of P2
        // include only first-order responses
        // see Spiro2SVG Book 3, page 52
        // setup two elliptical equations for x2, y2

        final int MAXLOOP = 200;
        double A0, B0, C0, D0, E0, F0;
        double A1, B1, C1, D1, E1, F1;
        double A2, C2, D2, E2, F2;
        double t1, old_x2, old_y2;
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
            if (Double.isNaN(solve_at_P2(x2, y2, false)))               // initiallize at (x2, y2)
            {
                System.out.println("fail at " + x2 + ", " + y2 + ", " + calc_d1(x2, y2) + ", " + calc_d2(x2, y2));
                return;
            }

            for (i = 0; i <= N; i++)
            {
                t1 = t1_start + i*(t1_end - t1_start)/N;
                h_gx[i] = multvv(Splinex, N43(t2[i])) - fitted.getx(t1);
                h_gy[i] = multvv(Spliney, N43(t2[i])) - fitted.gety(t1);
                dhx[i] = multvv(Splinex, dN43(t2[i]));
                dhy[i] = multvv(Spliney, dN43(t2[i]));
                ux[i] = N43(t2[i])[2] + Math.cos(theta_start)*N43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.cos(theta_end)*N43(t2[i])[3]*calc_d2dx2(x2, y2);
                vx[i] =                 Math.cos(theta_start)*N43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.cos(theta_end)*N43(t2[i])[3]*calc_d2dy2(x2, y2);
                uy[i] =                 Math.sin(theta_start)*N43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.sin(theta_end)*N43(t2[i])[3]*calc_d2dx2(x2, y2);
                vy[i] = N43(t2[i])[2] + Math.sin(theta_start)*N43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.sin(theta_end)*N43(t2[i])[3]*calc_d2dy2(x2, y2);
                dux[i] = dN43(t2[i])[2] + Math.cos(theta_start)*dN43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.cos(theta_end)*dN43(t2[i])[3]*calc_d2dx2(x2, y2);
                dvx[i] =                  Math.cos(theta_start)*dN43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.cos(theta_end)*dN43(t2[i])[3]*calc_d2dy2(x2, y2);
                duy[i] =                  Math.sin(theta_start)*dN43(t2[i])[1]*calc_d1dx2(x2, y2) - Math.sin(theta_end)*dN43(t2[i])[3]*calc_d2dx2(x2, y2);
                dvy[i] = dN43(t2[i])[2] + Math.sin(theta_start)*dN43(t2[i])[1]*calc_d1dy2(x2, y2) - Math.sin(theta_end)*dN43(t2[i])[3]*calc_d2dy2(x2, y2);
                h_gx[i] -= ux[i]*x2 + vx[i]*y2;
                h_gy[i] -= uy[i]*x2 + vy[i]*y2;
                dhx[i] -= dux[i]*x2 + dvx[i]*y2;
                dhy[i] -= duy[i]*x2 + dvy[i]*y2;
            }

            // set dFdx2 = 0

            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dux[i] + uy[i]*duy[i])*t2dx2[i];
            A0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (vx[i]*dvx[i] + vy[i]*dvy[i])*t2dx2[i];
            B0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dvx[i] + dux[i]*vx[i] + uy[i]*dvy[i] + duy[i]*vy[i])*t2dx2[i];
            C0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = ux[i]*ux[i] + (ux[i]*dhx[i] + dux[i]*h_gx[i])*t2dx2[i]
                           + uy[i]*uy[i] + (uy[i]*dhy[i] + duy[i]*h_gy[i])*t2dx2[i];
            D0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = vx[i]*ux[i] + (vx[i]*dhx[i] + dvx[i]*h_gx[i])*t2dx2[i]
                           + vy[i]*uy[i] + (vy[i]*dhy[i] + dvy[i]*h_gy[i])*t2dx2[i];
            E0 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = h_gx[i]*(ux[i] + dhx[i]*t2dx2[i])
                           + h_gy[i]*(uy[i] + dhy[i]*t2dx2[i]);
            F0 = integrate(trap_in);

            // set dFdy2 = 0

            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dux[i] + uy[i]*duy[i])*t2dy2[i];
            A1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (vx[i]*dvx[i] + vy[i]*dvy[i])*t2dy2[i];
            B1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = (ux[i]*dvx[i] + dux[i]*vx[i] + uy[i]*dvy[i] + duy[i]*vy[i])*t2dy2[i];
            C1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = ux[i]*vx[i] + (ux[i]*dhx[i] + dux[i]*h_gx[i])*t2dy2[i]
                           + uy[i]*vy[i] + (uy[i]*dhy[i] + duy[i]*h_gy[i])*t2dy2[i];
            D1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = vx[i]*vx[i] + (vx[i]*dhx[i] + dvx[i]*h_gx[i])*t2dy2[i]
                           + vy[i]*vy[i] + (vy[i]*dhy[i] + dvy[i]*h_gy[i])*t2dy2[i];
            E1 = integrate(trap_in);
            for (i = 0; i <= N; i++)
                trap_in[i] = h_gx[i]*(vx[i] + dhx[i]*t2dy2[i])
                           + h_gy[i]*(vy[i] + dhy[i]*t2dy2[i]);
            F1 = integrate(trap_in);

            // see: https://math.stackexchange.com/questions/1767225/algorithm-intersection-of-two-conics

            A2 = A0*B1 - A1*B0;
            C2 = C0*B1 - C1*B0;
            D2 = D0*B1 - D1*B0;
            E2 = E0*B1 - E1*B0;
            F2 = F0*B1 - F1*B0;
            old_x2 = x2;
            old_y2 = y2;
            x2 = fitymoment.solve_quartic(A0*C2*C2 + B0*A2*A2 - C0*C2*A2,
                                          2*A0*C2*E2 + D0*C2*C2 + 2*B0*A2*D2 - C0*C2*D2 - (C0*E2 + C2*E0)*A2,
                                          A0*E2*E2 + 2*D0*C2*E2 + F0*C2*C2 + 2*B0*A2*F2 + B0*D2*D2
                                          - C0*C2*F2 - (C0*E2 + C2*E0)*D2 - E0*E2*A2,
                                          D0*E2*E2 + 2*F0*C2*E2 + 2*B0*D2*F2 - (C0*E2 + C2*E0)*F2 - E0*E2*D2,
                                          F0*E2*E2 + B0*F2*F2 - E0*E2*F2,
                                          true);
            y2 = -(A2*x2*x2 + D2*x2 + F2)/(C2*x2 + E2);
            //System.out.println("new           d1 d2    = ," + d1 + ", " + d2
            //                 + ", " + (A0*d1*d1 + B0*d2*d2 + C0*d1*d2 + D0*d1 + E0*d2 + F0)
            //                 + ", " + (A1*d1*d1 + B1*d2*d2 + C1*d1*d2 + D1*d1 + E1*d2 + F1));

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)
            {
                t2[i] += t2dx2[i]*(x2 - old_x2) + t2dy2[i]*(y2 - old_y2);   // first-order response
                //System.out.println((t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(old_x2 - x2) < TOL) && (Math.abs(old_y2 - y2) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new x2 y2 = , , , , , , " + x2 + ", " + y2 + ", " + Math.abs(old_x2 - x2) + ", " + Math.abs(old_y2 - y2));
            solve_at_P2(x2, y2, true);                  // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops!");
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
        if (firstline.startsWith("theta"))                      // cycloid, calculate c
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
        System.out.println("file data at theta c t d1 d2   = ," + theta + ", " + tempc + ", " + t1_start + ", " + t1_end + ", " + fitted.getkappa(t1_start) + ", " + fitted.getkappa(t1_end) + ", " + d1 + ", " + d2);

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
        //P2x = 0.67; // 1.2545; // 2.408125;
        //P2y = 1.39; // 1.5819; // 1.540;
        System.out.println("return code solve_at_P2 = " + solve_at_P2(P2x, P2y, true));
        //grid_search_at_P2(P2x, P2y);
        //iterate_at_P2(P2x, P2y);     // default
    }

    private static void grid_search_at_P2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        final int MAXCOUNT = 1000;
        double del = 0.000001;
        double err, errmin = 999999;
        int iold = 0, jold = 0;
        int imin = 0, jmin = 0;
        int count = 0;

        do
        {
            count++;
            errmin = 999999;
            for (int i = -1; i < 2; i++)
                for (int j = -1; j < 2; j++)
                {
                    err = solve_at_P2(x2 + (i + iold)*del, y2 + (j + jold)*del, false);
                    if (err < errmin)
                    {
                        imin = i;
                        jmin = j;
                        errmin = err;
                    }
                }
            iold += imin;
            jold += jmin;
            System.out.println(",,,,,,,,,,, " + (float)(x2 + iold*del) + ", " + (float)(y2 + jold*del) + ", " + errmin);
        } while ((count < MAXCOUNT) && !(imin == 0 && jmin == 0));
        if (count == MAXCOUNT)
            System.out.println("Not converged!");
        else
            System.out.println("Converged in " + count);
    }

    private static double solve_at_P2(double x2, double y2, boolean print)
    {
        // perform a single calculation of a complete t2[] profile
        // at a given P2, and calculate the rms error
        // (calculate d1, d2 to satisfy the curvature at the endpoints)

        theta_start = fitted.gettheta(t1_start);
        theta_end = fitted.gettheta(t1_end);
        double d1 = calc_d1(x2, y2);
        double d2 = calc_d2(x2, y2);
        if (Double.isNaN(d1) || Double.isNaN(d2))
            return Double.NaN;

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
            System.out.println("__start at theta c t d1 d2 = , " + theta_start*180/Math.PI + ", " + theta_end*180/Math.PI + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + calc_error());
        //fitted.gen_Bezier2(Bezx, Bezy);
        //System.out.println(Bezx[0][0] + "\t " + Bezy[0][0]);
        //System.out.println(Bezx[0][1] + "\t " + Bezy[0][1]);
        //System.out.println(Bezx[0][2] + "\t " + Bezy[0][2]);
        //System.out.println(Bezx[0][3] + "\t " + Bezy[0][3]);
        //System.out.println(Bezx[1][0] + "\t " + Bezy[1][0]);
        //System.out.println(Bezx[1][1] + "\t " + Bezy[1][1]);
        //System.out.println(Bezx[1][2] + "\t " + Bezy[1][2]);
        //System.out.println(Bezx[1][3] + "\t " + Bezy[1][3]);

        if (print) System.out.println("\nseg, t1, t2, t2dx2, t2dy2");
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
            t2dd1[i] = calc_t2dxy(i, t2[i], "d1");
            t2dd2[i] = calc_t2dxy(i, t2[i], "d2");
            t2dx2[i] = calc_t2dxy(i, t2[i], "x2");
            t2dy2[i] = calc_t2dxy(i, t2[i], "y2");
            if (print) System.out.println(seg + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd1[i] + ", " + t2dd2[i] + ", " + t2dx2[i] + ", " + t2dy2[i]);
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(2 - t2[i]) > TOL)
            ||  (i == N && seg != 1)
            ||  (t2[i] < -TOL) || (t2[i] == Double.NaN))
            {
                System.out.println("abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
//                scan_quintic_near_t2(t1, t2[i]);
                return Double.NaN;
            }
        }
        //System.out.println("new t2[] profile rms   = ," + d1 + ", " + d2 + ", " + calc_error());
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
        //System.out.println("\ninit  t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime);
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

    private static double calc_t2dxy(int i, double t2, String type)
    {
        double t1 = t1_start + i*(t1_end - t1_start)/N;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);
        double fnx = multvv(Splinex, N43(t2));
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
        return -numer/denom;
    }

    private static double calc_d1(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        return Math.sqrt(((y2 - fitted.gety(t1_start))*Math.cos(theta_start) - (x2 - fitted.getx(t1_start))*Math.sin(theta_start))/3/fitted.getkappa(t1_start));
    }

    private static double calc_d2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        return Math.sqrt(((y2 - fitted.gety(t1_end))*Math.cos(theta_end) - (x2 - fitted.getx(t1_end))*Math.sin(theta_end))/3/fitted.getkappa(t1_end));
    }

    private static double calc_d1dx2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        return -Math.sin(theta_start)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start);
    }

    private static double calc_d1dy2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        return Math.cos(theta_start)/calc_d1(x2, y2)/6/fitted.getkappa(t1_start);
    }

    private static double calc_d2dx2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        return -Math.sin(theta_end)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end);
    }

    private static double calc_d2dy2(double x2, double y2) // obsolete, used only if constraining the curvature
    {
        return Math.cos(theta_end)/calc_d2(x2, y2)/6/fitted.getkappa(t1_end);
    }

    private static double[] N43(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {(1 - u)*(1 - u)*(1 - u),
                                 u*(12 - 18*u + 7*u*u)/4,
                                 u*u*(3 - 2*u)/2,
                                 u*u*u/4,
                                 0};
        else if (u < 2 + TOL)
            return new double[] {0,
                                 (2 - u)*(2 - u)*(2 - u)/4,
                                 (2 - u)*(2 - u)*(2*u - 1)/2,
                                 (2 - u)*(7*u*u - 10*u + 4)/4,
                                 (u - 1)*(u - 1)*(u - 1)};
        else
            return null;
    }

    private static double[] dN43(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {-3*(1 - u)*(1 - u),
                                 (12 - 36*u + 21*u*u)/4,
                                 3*u*(1 - u),
                                 3*u*u/4,
                                 0};
        else if (u < 2 + TOL)
            return new double[] {0,
                                 -3*(2 - u)*(2 - u)/4,
                                 -3*(2 - u)*(u - 1),
                                 -(21*u*u - 48*u + 24)/4,
                                 3*(u - 1)*(u - 1)};
        else
            return null;
    }

    private static double[] d2N43(double u)
    {
        if (u < -TOL)
            return null;
        else if (u < 1)
            return new double[] {6*(1 - u),
                                 (-36 + 42*u)/4,
                                 3*(1 - 2*u),
                                 3*u/2,
                                 0};
        else if (u < 2 + TOL)
            return new double[] {0,
                                 3*(2 - u)/2,
                                 3*(2*u - 3),
                                 (-42*u + 48)/4,
                                 6*(u - 1)};
        else
            return null;
    }

    private static double multvv(double[] v1, double[] v2)
    {
        if (v1.length != v2.length) return Double.NaN;
        double retVal = 0;

        for (int i = 0; i < v1.length; i++)
            retVal += v1[i]*v2[i];
        return retVal;
    }

    private static double[] multmv(double[][] m, double[] v)
    {
        if (m[0].length != v.length) return null;
        double[] retVal = new double[m.length];

        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < v.length; j++)
                retVal[i] += m[i][j]*v[j];
        return retVal;
    }

    private static double[][] invertm(double[][] m)
    {
        if (m.length != 4) return null;
        if (m[0].length != 4) return null;
        double[][] retVal = new double[4][4];
        double det = 0;
        int sgn = 1;

        for (int i = 0; i < 4; i++)
        {
            det += sgn*m[0][i]*detm(m, 0, i);
            sgn *= -1;
        }

        sgn = 1;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                retVal[i][j] = sgn*detm(m, j, i)/det;
                if (j < 3) sgn *= -1;
            }
        return retVal;
    }

    private static double detm(double[][] m, int row, int col)
    {
        if (row < 0 || row > 3) return Double.NaN;
        if (col < 0 || col > 3) return Double.NaN;
        if (m.length != 4) return Double.NaN;
        if (m[0].length != 4) return Double.NaN;
        int[] arri = new int[3];
        int[] arrj = new int[3];
        int index = -1;

        for (int i = 0; i < 4; i++)
            if (i != row)
            {
                index++;
                arri[index] = i;
            }
        index = -1;
        for (int i = 0; i < 4; i++)
            if (i != col)
            {
                index++;
                arrj[index] = i;
            }
        return m[arri[0]][arrj[0]]*m[arri[1]][arrj[1]]*m[arri[2]][arrj[2]]
             + m[arri[0]][arrj[1]]*m[arri[1]][arrj[2]]*m[arri[2]][arrj[0]]
             + m[arri[0]][arrj[2]]*m[arri[1]][arrj[0]]*m[arri[2]][arrj[1]]
             - m[arri[0]][arrj[2]]*m[arri[1]][arrj[1]]*m[arri[2]][arrj[0]]
             - m[arri[0]][arrj[0]]*m[arri[1]][arrj[2]]*m[arri[2]][arrj[1]]
             - m[arri[0]][arrj[1]]*m[arri[1]][arrj[0]]*m[arri[2]][arrj[2]];
    }
}
