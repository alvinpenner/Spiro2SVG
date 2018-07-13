
package components;

import java.io.*;

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
    public static final double t1_end = Math.PI/4; // Math.PI;
    public static final int N = 100;
    public static double[] Bezx;                // quartic Bezier, 5 points, x component
    public static double[] Bezy;                // quartic Bezier, 5 points, y component
    //private static CircleFxn fitted;
    //private static CycloidFxn fitted;
    private static epiTrochoidFxn fitted;
    private static double[] t2 = new double[N+1];
    private static double[][] t2dd = new double[4][N+1];    // partial wrt (d1, d2, x2, y2, x3, y3)
    private static double Jacdet = Double.NaN;
    private static double[] eig = new double[] {0,0,0,0};
    public static double theta_start, theta_end;
    private static final double TOL = 0.000000001;

    public static void main (String[] args)
    {
        // extract the point of maximum curvature of a cycloid from a tangent angle phi
        //double phi = 53;
        //double tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //fitted = new CycloidFxn(tempc);
        fitted = new epiTrochoidFxn(0.25592);
        //paste_data(80, 0.9886277018143461, 0.2624690492231421, 0.7640885237135162, 1.8275481762035848);
        //read_Quartic_Bezier_data(230);
        iterate_at_P2(32.02406276623338, 38.62191169612965, 172.49594071073963, 66.40511369595288); // 19.25 minimum
        //iterate_at_P2(31.464595462598876, 122.8023528021591, 161.91830044656243, 57.58332982440786); // 19.25 saddle
        //iterate_at_P2(40.923936413081165, 174.54348410425095, 135.45126668176417, 77.10113250461616); // 19.25
        //iterate_at_P2(37.81048776975601, 142.1606207834995, 151.90527738346833, 63.000956232673424); // 18.0
        //iterate_at_P2(40.38173870009711, 171.26075879775118, 137.5341366809092, 75.39842833631714); // 19.2
        //iterate_at_P2(40.923936413081165, 174.54348410425095, 135.45126668176417, 77.10113250461616); // 19.26
        //System.out.println("quartic Bezier solve_at_P2 = " + solve_at_P2(35.564623789595046, 35.564623789653915, 170.63083487797272, 70.67760596545668, true) + "\n");
        //iterate_at_P2(30, 36, 173, 60);         // test code at c = 10
        if (fitted == null)
        {
            System.out.println("class 'fitted' is not defined, abort");
            return;
        }
        // convert BSpline5 (d1,d2) to Quartic Bezier (d1, d2)
        //System.out.println("solve_at_P2 = " + solve_at_P2(1.5*0.3820442618567581, 1.5*0.9137740881017924, 0.7263404994170831, 1.393169269182819, true));
        //System.out.println("solve_at_P2 = " + solve_at_P2(1.5*0.382, 1.5*0.914, 0.726, 1.393, true) + "\n");
        //grid_search_at_P2(.571, 1.369, 0.724, 1.394);
        //System.out.println(fitted.getClass().getCanonicalName());
        //scan_septic_simple(43);
        //System.out.printf("M");
        //for (int i = 0; i <= N; i++)
        //    System.out.printf(" %f, %f", 80 + 2*BSpline5.multvv(Bezx, N44(t2[i])), 500 - 2*BSpline5.multvv(Bezy, N44(t2[i])));
        //System.out.println("\n");
    }

    private static void iterate_at_P2(double d1, double d2, double x2, double y2)
    {
        // calculate a new estimate of (d1, d2, x2, y2) by setting dF = 0
        // include only first-order responses
        // see Spiro2SVG Book 3, page 54 (applied to quartic Bezier)
        // setup 4-variable Newton-Raphson iteration

        final int MAXLOOP = 1000;
        double[] f_gx = new double[N+1];
        double[] f_gy = new double[N+1];
        double[] dfxdu = new double[N+1];
        double[] d2fxdudu = new double[N+1];
        double[] dfydu = new double[N+1];
        double[] d2fydudu = new double[N+1];
        double[][] dfxdd = new double[4][N+1];
        double[][] dfydd = new double[4][N+1];

        double[][] Jac = new double[4][4];
        double[] dFdd = new double[4];
        double[] trap_in = new double[N+1];
        double[] deld;                                              // (-Δd1, -Δd2, -Δx2, -Δy2)
        int i, j, k, loop = 0;
        double t1;

        do
        {
            loop++;
            if (Double.isNaN(solve_at_P2(d1, d2, x2, y2, false)))   // initiallize at (x2, y2)
            {
                System.out.println("fail at (d1, d2, x2, y2): " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
                return;
            }

            for (i = 0; i <= N; i++)
            {
                t1 = t1_start + i*(t1_end - t1_start)/N;
                f_gx[i] = BSpline5.multvv(Bezx, N44(t2[i])) - fitted.getx(t1);
                f_gy[i] = BSpline5.multvv(Bezy, N44(t2[i])) - fitted.gety(t1);
                dfxdu[i] = BSpline5.multvv(Bezx, dN44(t2[i]));
                d2fxdudu[i] = BSpline5.multvv(Bezx, d2N44(t2[i]));
                dfydu[i] = BSpline5.multvv(Bezy, dN44(t2[i]));
                d2fydudu[i] = BSpline5.multvv(Bezy, d2N44(t2[i]));
                dfxdd[0][i] = Math.cos(theta_start)*N44(t2[i])[1];
                dfydd[0][i] = Math.sin(theta_start)*N44(t2[i])[1];
                dfxdd[1][i] = -Math.cos(theta_end)*N44(t2[i])[3];
                dfydd[1][i] = -Math.sin(theta_end)*N44(t2[i])[3];
                dfxdd[2][i] = N44(t2[i])[2];
                dfydd[2][i] = 0;
                dfxdd[3][i] = 0;
                dfydd[3][i] = N44(t2[i])[2];
            }

            // calc dFdd[j] at current (d1, d2, x2, y2)

            for (i = 0; i < 4; i++)
            {
                for (k = 0; k <= N; k++)
                    trap_in[k] = f_gx[k]*(dfxdd[i][k] + dfxdu[k]*t2dd[i][k]) + f_gy[k]*(dfydd[i][k] + dfydu[k]*t2dd[i][k]); // original code
                    //trap_in[k] = f_gx[k]*dfxdd[i][k] + f_gy[k]*dfydd[i][k];         // new code
                dFdd[i] = t2_vs_t1.integrate(trap_in);
            }

            // calc d2Fdd[i]dd[j] (Jacobean matrix)

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                {
                    //System.out.println(i + ", "+ j);
                    for (k = 0; k <= N; k++)
                    {
                        trap_in[k] = dfxdd[i][k]*dfxdd[j][k] + dfydd[i][k]*dfydd[j][k]
                                   - (dfxdu[k]*dfxdu[k] + dfydu[k]*dfydu[k] + f_gx[k]*d2fxdudu[k] + f_gy[k]*d2fydudu[k])*t2dd[i][k]*t2dd[j][k];
                        //System.out.println(k + ", " + trap_in[k]);
                    }
                    Jac[i][j] = t2_vs_t1.integrate(trap_in);
                }

            //deld = BSpline5.multmv(BSpline5.invertm(Jac), dFdd);  // this is actually the negative of Δd
            deld = BSpline5.gaussj(Jac, dFdd);                      // this is actually the negative of Δd
            d1 -= deld[0];
            d2 -= deld[1];
            x2 -= deld[2];
            y2 -= deld[3];

            Jacdet = BSpline5.detm(Jac);
            // calculate four eigenvalues
            double qua = -Jac[0][0] - Jac[1][1] - Jac[2][2] - Jac[3][3];
            double qub =  Jac[0][0]*Jac[1][1] + Jac[0][0]*Jac[2][2] + Jac[0][0]*Jac[3][3] + Jac[1][1]*Jac[2][2] + Jac[1][1]*Jac[3][3] + Jac[2][2]*Jac[3][3]
                       -  Jac[0][1]*Jac[0][1] - Jac[0][2]*Jac[0][2] - Jac[0][3]*Jac[0][3] - Jac[1][2]*Jac[1][2] - Jac[1][3]*Jac[1][3] - Jac[2][3]*Jac[2][3];
            double quc = -Jac[0][0]*Jac[1][1]*Jac[2][2] - Jac[0][0]*Jac[1][1]*Jac[3][3] - Jac[0][0]*Jac[2][2]*Jac[3][3] - Jac[1][1]*Jac[2][2]*Jac[3][3]
                       +  Jac[0][1]*Jac[0][1]*(Jac[2][2] + Jac[3][3]) + Jac[0][2]*Jac[0][2]*(Jac[1][1] + Jac[3][3]) + Jac[0][3]*Jac[0][3]*(Jac[1][1] + Jac[2][2])
                       +  Jac[1][2]*Jac[1][2]*(Jac[0][0] + Jac[3][3]) + Jac[1][3]*Jac[1][3]*(Jac[0][0] + Jac[2][2]) + Jac[2][3]*Jac[2][3]*(Jac[0][0] + Jac[1][1])
                       - 2*Jac[0][1]*Jac[1][2]*Jac[0][2] - 2*Jac[0][2]*Jac[2][3]*Jac[0][3] - 2*Jac[0][1]*Jac[1][3]*Jac[0][3] - 2*Jac[1][2]*Jac[2][3]*Jac[1][3];
            double qud = Jac[0][0]*Jac[1][1]*Jac[2][2]*Jac[3][3]
                       + Jac[0][1]*Jac[0][1]*Jac[2][3]*Jac[2][3] + Jac[0][2]*Jac[0][2]*Jac[1][3]*Jac[1][3] + Jac[0][3]*Jac[0][3]*Jac[1][2]*Jac[1][2]
                       - Jac[0][1]*Jac[0][1]*Jac[2][2]*Jac[3][3] - Jac[0][2]*Jac[0][2]*Jac[1][1]*Jac[3][3] - Jac[0][3]*Jac[0][3]*Jac[1][1]*Jac[2][2]
                       - Jac[1][2]*Jac[1][2]*Jac[0][0]*Jac[3][3] - Jac[1][3]*Jac[1][3]*Jac[0][0]*Jac[2][2] - Jac[2][3]*Jac[2][3]*Jac[0][0]*Jac[1][1]
                       + 2*Jac[0][0]*Jac[1][2]*Jac[2][3]*Jac[1][3] + 2*Jac[1][1]*Jac[0][2]*Jac[2][3]*Jac[0][3] + 2*Jac[2][2]*Jac[0][1]*Jac[1][3]*Jac[0][3] + 2*Jac[3][3]*Jac[0][1]*Jac[1][2]*Jac[0][2]
                       - 2*Jac[0][1]*Jac[1][2]*Jac[2][3]*Jac[0][3] - 2*Jac[0][1]*Jac[0][2]*Jac[2][3]*Jac[1][3] - 2*Jac[0][2]*Jac[0][3]*Jac[1][3]*Jac[1][2];
            eig = fitymoment.solve_quartic_all(1, qua, qub, quc, qud);
            //System.out.println("eigenvalue = " + eig[1] + ", " + eig[3] + ", " + eig[2] + ", " + eig[0]);
            System.out.println("dFdd = " + dFdd[0] + ", " + dFdd[1] + ", " + dFdd[2] + ", " + dFdd[3] + ", " + Jacdet);
            System.out.println("deld = " + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3]);
            BSpline5.dump_Jac(Jac);

            // perform a preliminary first-order recalculation of t2[i]
            // just for the purpose of improving the calc_error() result
            //System.out.println("\npreliminary recalc of t2[i]\n t1, t2");
            for (i = 0; i <= N; i++)            // disabled due to crashes
            {
                t2[i] -= t2dd[0][i]*deld[0] + t2dd[1][i]*deld[1] + t2dd[2][i]*deld[2] + t2dd[3][i]*deld[3];   // first-order response
                // System.out.println("incrementing t2 array : " + i + ", " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
            }
        } while ((loop < MAXLOOP) && !((Math.abs(deld[0]) < TOL) && (Math.abs(deld[1]) < TOL) && (Math.abs(deld[2]) < TOL) && (Math.abs(deld[3]) < TOL)));
        if (loop < MAXLOOP)
        {
            System.out.println("\n__converged in " + loop + " at new d1 d2 x2 y2 = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
            solve_at_P2(d1, d2, x2, y2, true);                          // final run just for good measure
        }
        else
            System.out.println("\nNOT converged after " + loop + " loops! (" + deld[0] + ", " + deld[1] + ", " + deld[2] + ", " + deld[3] + ")");
    }

    private static void read_Quartic_Bezier_data(int ln)
    {
        // read quartic Bezier data (c, d1, d2, x2, y2) data from a file
        // read only one line (ln), and assume c is positive
        // convert to negative c parameters, and initiallize 'iterate_at_P2()'

        String str = "";
        double c = 0, d1 = 0, d2 = 0, x2 = 0, y2 = 0, x2new = 0, y2new = 0;
        int row = 0;

        try
        {
            //BufferedReader istr = new BufferedReader(new FileReader("C:\\APP\\Java\\SpiroGraph\\BezierQuartic\\hypoQuarticBezd1d2rms.csv"));
            BufferedReader istr = new BufferedReader(new FileReader("\\Windows\\Temp\\hypoQuarticBezd1d2rms.csv"));
            try
            {
                while (istr.ready() && row < ln)
                {
                    str = istr.readLine();
                    row++;
                    //System.out.println(str);
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

        System.out.println(str);
        if (str.startsWith("__new"))
        {
            c = Double.parseDouble(str.split(",")[3]);
            d1 = Double.parseDouble(str.split(",")[6]);
            d2 = Double.parseDouble(str.split(",")[7]);
            x2 = Double.parseDouble(str.split(",")[8]);
            y2 = Double.parseDouble(str.split(",")[9]);
        }
        if (d1 == 0 && d2 == 0)
        {
            System.out.println("file data match not found, abort");
            return;
        }
        if (c < 0)
        {
            System.out.println("c must be positive : " + c);
            return;
        }
        System.out.println("file data at c d1 d2   = ," + c + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2);

        // calculate P2 for a Quartic Bezier with a negative c

        fitted = new epiTrochoidFxn(-c);
        theta_start = fitted.gettheta(t1_start);                        // possibly redundant
        theta_end = fitted.gettheta(t1_end);
        x2new =  Math.cos(Math.PI/8)*x2 + Math.sin(Math.PI/8)*y2;       // transform to symmetric frame
        y2new = -Math.sin(Math.PI/8)*x2 + Math.cos(Math.PI/8)*y2;
        x2 =  Math.cos(Math.PI/8)*x2new + Math.sin(Math.PI/8)*y2new;    // reflect and back-transform
        y2 =  Math.sin(Math.PI/8)*x2new - Math.cos(Math.PI/8)*y2new;
        //System.out.println("return code solve_at_P2 = " + solve_at_P2(d1, d2, x2, y2, true));
        iterate_at_P2(d2, d1, x2, y2);
    }

    private static void paste_data(double phi, double tempc, double dummy, double d1, double d2)
    {
        // paste cubic Bezier (d1, d2) data from a file
        // field 1 = phi for a Cycloid / root number for an epiTrochoid
        // field 2 = c for an epiTrochoid
        // convert to d1, d2, and P2 for a 5 point quartic Bezier, and initiallize 'solve_at_P2()'

        //tempc = Math.sqrt(1 - .75*Math.cos(phi*Math.PI/180)*Math.cos(phi*Math.PI/180));
        //fitted = new CycloidFxn(tempc);
        //t1_start = Math.acos((2*tempc*tempc - 1)/tempc);
        //tempc = 4.5;                                        // override for circle
        //fitted = new CircleFxn(tempc);                      // set c value (radius)
        //t1_start = Math.PI - phi*Math.PI/180;               // set arc angle
        fitted = new epiTrochoidFxn(tempc);
        
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

        // convert from cubic Bezier (d1, d2) to quartic Bezier (d1, d2)
        //System.out.println("quartic Bezier solve_at_P2 = " + solve_at_P2(0.75*d1, 0.75*d2, P2x, P2y, true) + "\n");
        iterate_at_P2(0.75*d1, 0.75*d2, P2x, P2y);
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
            System.out.println("__start quartic Bezier at theta c t d1 d2 = , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + fitted.getc() + ", " + t1_start + ", " + t1_end + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2);
        else
            System.out.println("__solve at new d1 d2 rms = , , , , , , " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + calc_error());
        if (d1 < 0 || d2 < 0)
            System.out.println("WARNING: negative arm length = " + d1 + ", " + d2);
        //gen_BezierQuartic(Bezx, Bezy);
        //System.out.println(Bezx[0] + "\t " + Bezy[0]);
        //System.out.println(Bezx[1] + "\t " + Bezy[1]);
        //System.out.println(Bezx[2] + "\t " + Bezy[2]);
        //System.out.println(Bezx[3] + "\t " + Bezy[3]);
        //System.out.println(Bezx[4] + "\t " + Bezy[4]);

        if (print) System.out.println("\n           , t1, t2, t2dd1, t2dd2, t2dx2, t2dy2");
        for (int i = 0; i <= N; i++)
        {
            solve_septic_for_t2(i);
            if ((i == 0 && Math.abs(t2[i]) > TOL)
            ||  (i == N && Math.abs(1 - t2[i]) > TOL)
            ||  (t2[i] < -TOL) || Double.isNaN(t2[i]))
            {
                System.out.println("t2[i] abort at " + i + " : " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i]);
                //for (int ii = 0; ii <= i; ii++)
                //    System.out.println(ii + ", " + t2[ii]);
                //scan_septic_simple(i - 1);
                //scan_septic_simple(i);
                return Double.NaN;
            }
            t2dd[0][i] = calc_t2dxy(i, t2[i], "d1");
            t2dd[1][i] = calc_t2dxy(i, t2[i], "d2");
            t2dd[2][i] = calc_t2dxy(i, t2[i], "x2");
            t2dd[3][i] = calc_t2dxy(i, t2[i], "y2");
            if (print)
                System.out.println("quartic Bez, " + (t1_start + i*(t1_end - t1_start)/N) + ", " + t2[i] + ", " + t2dd[0][i] + ", " + t2dd[1][i] + ", " + t2dd[2][i] + ", " + t2dd[3][i]);
        }
        double retVal = calc_error();
        //System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + (float) retVal + ", " + (float) Jacdet + ", " + (float) eig[1] + ", " + (float) eig[3] + ", " + (float) eig[2] + ", " + (float) eig[0]);
        System.out.println("gauss t2[] @ , " + (float) (theta_start*180/Math.PI) + ", " + (float) (theta_end*180/Math.PI) + ", " + (float) fitted.getc() + ", " + ", " + ", " + d1 + ", " + d2 + ", " + x2 + ", " + y2 + ", " + (float) retVal + ", " + (float) Jacdet + ", ");
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
        return Math.sqrt(t2_vs_t1.integrate(trap_in))/a_b;
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
            if (loop > 200)
            {
                System.out.println("\ntoo many loops t1 t2 =, " + t1 + ", " + t + ", " + f + ", " + fprime + ", " + f2prime + ", " + del_t);
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

    private static void scan_septic_simple(int index)
    {
        // if solve_septic_for_t2 fails, scan the area for other roots
        double t1 = t1_start + index*(t1_end - t1_start)/N;
        double f, t;
        double X = fitted.getx(t1);
        double Y = fitted.gety(t1);

        System.out.println("\nscanning at " + index + ", " + t1);
        for (int i = 0; i <= 100; i++)
        {
            t = i/100.0;
            f = (fn(Bezx, t) - X)*dfn(Bezx, t) + (fn(Bezy, t) - Y)*dfn(Bezy, t);
            System.out.println(t + ", " + f);
        }
        System.out.println();
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
            System.out.println("WARNING:   N44 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:   N44 too large u value = " + u);
        return new double[] {(1 - u)*(1 - u)*(1 - u)*(1 - u),
                             4*u*(1 - u)*(1 - u)*(1 - u),
                             6*u*u*(1 - u)*(1 - u),
                             4*u*u*u*(1 - u),
                             u*u*u*u};
    }

    private static double[] dN44(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING:  dN44 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING:  dN44 too large u value = " + u);
        return new double[] {-4*(1 - u)*(1 - u)*(1 - u),
                              4*(1 - 4*u)*(1 - u)*(1 - u),
                             12*(1 - 2*u)*u*(1 - u),
                              4*(3 - 4*u)*u*u,
                              4*u*u*u};
    }

    private static double[] d2N44(double u)
    {
        if (u < -TOL)
            System.out.println("WARNING: d2N44 negative u value = " + u);
        if (u > 1 + 1000*TOL)
            System.out.println("WARNING: d2N44 too large u value = " + u);
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

    private static void gen_BezierQuartic(double[] ptx, double[] pty)
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
