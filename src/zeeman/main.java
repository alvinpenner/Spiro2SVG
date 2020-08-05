
// ZeemanMachine - May 2020 - Alvin Penner - penner@vaxxine.com

// calculate critical points of a potential energy function
// for the Zeeman Catastrophe Machine (Poston and Stewart p. 75)
// radius of wheel = 1
// position of point A = (0, -A), e = distance to wheel
// position of point B = (x,  y), e = distance to wheel
// default length of elastics = A/2
// y in range 2*(1.404, 2.455) for A = 4

// this is file : \Documents\NetBeansProjects\ZeemanMachine\src\zeeman\main.java

package zeeman;

import java.awt.geom.Point2D;
import javax.swing.*;
import java.io.*;
import java.util.Properties;

public class main
{
    private static Properties pgmProp = new Properties();
    protected static double[][] xbound = new double[4][];
    protected static double[][] ybound = new double[4][];
    protected static double A, x, y, keyincr;                       // variables from staticComponent
    protected static double ystart, yend;                           // variables from Bifurcate
    protected static double theta0, w0, x0, xa, y0, c, Tx, phi0;    // variables from PhaseSpace
    protected static double dtheta0dy, dw0dy;
    protected static int NLimit;                                    // # of Tx per limit cycle (to be determined)
    protected static final String VERSION_NO = "0.1";
    private static final double TOL = 0.0000000001;

    public static void main(String[] args)
    {
        staticFrame frame = new staticFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
        for (int i = 0; i < xbound.length; i++)
            get_boundary(i);
        //double theta = solve_for_critical(Math.PI, -0.1, 4);
        //System.out.println("critical theta = ," + 180/Math.PI*theta + ", " + calc_d2Fdth2(theta, -.1, 4));

        //double th = 5.04;
        //for (int i = 5780; i <= 5800; i += 4)
        //    th = solve_for_boundary(th, 6, 1.36, i*0.001);
        //double tmpx = 0.0;
        //double tmpy = 5.5;
        //for (int i = 0; i <= 360; i ++)
            //System.out.println(i + ", " + calc_F(i*Math.PI/180, A, tmpx, tmpy) + ", " + calc_dFdth(i*Math.PI/180, A, tmpx, tmpy) + ", " + calc_d2Fdth2(i*Math.PI/180, A, tmpx, tmpy) + ", " + calc_d2Fdthdy(i*Math.PI/180, A, tmpx, tmpy));
            //System.out.println(i + ", " + dedth(i*Math.PI/180, tmpx, tmpy) + ", " + d2edthdx(i*Math.PI/180, tmpx, tmpy) + ", " + d2edthdy(i*Math.PI/180, tmpx, tmpy));
            //System.out.println(i + ", " + e(i*Math.PI/180, tmpx, tmpy) + ", " + dedx(i*Math.PI/180, tmpx, tmpy) + ", " + dedy(i*Math.PI/180, tmpx, tmpy));
   }

    protected static void load_prefs()
    {
        try                                         // recall program properties
        {
            if (new File(System.getProperty("user.home"), "ZCMPrefs.ini").exists())
            {
                pgmProp.load(new FileInputStream(new File(System.getProperty("user.home"), "ZCMPrefs.ini")));
                A = Double.parseDouble(pgmProp.getProperty("initA", "4"));
                x = Double.parseDouble(pgmProp.getProperty("initx", "200"));
                y = Double.parseDouble(pgmProp.getProperty("inity", "400"));
                keyincr = Double.parseDouble(pgmProp.getProperty("keyincr", "1"));
                theta0 = Double.parseDouble(pgmProp.getProperty("theta0", "0"));
                w0 = Double.parseDouble(pgmProp.getProperty("w0", "0"));
                x0 = Double.parseDouble(pgmProp.getProperty("x0", "0"));
                xa = Double.parseDouble(pgmProp.getProperty("xa", "0.6"));
                y0 = Double.parseDouble(pgmProp.getProperty("y0", "6"));
                c = Double.parseDouble(pgmProp.getProperty("c", "1"));
                Tx = Double.parseDouble(pgmProp.getProperty("Tx", "1"));
                phi0 = Double.parseDouble(pgmProp.getProperty("phi0", "0"));
                NLimit = Integer.parseInt(pgmProp.getProperty("NLimit", "0"));
                ystart = Double.parseDouble(pgmProp.getProperty("ystart", "6"));
                yend = Double.parseDouble(pgmProp.getProperty("yend", "6.5"));
                dtheta0dy = Double.parseDouble(pgmProp.getProperty("dtheta0dy", "0"));
                dw0dy = Double.parseDouble(pgmProp.getProperty("dw0dy", "0"));
            }
            else                                    // factory default
            {
                A = 4;
                x = 200;
                y = 400;
                keyincr = 1;
                theta0 = 0;
                w0 = 0;
                x0 = 0;
                xa = 0.6;
                y0 = 6;
                c = 1;
                Tx = 1;
                phi0 = 0;
                NLimit = 0;
                ystart = 6;
                yend = 6.5;
                dtheta0dy = 0;
                dw0dy = 0;
            }
        }
        catch (IOException e)
            {System.out.println("error reading ZCMPrefs.ini : " + e);}
    }

    protected static void save_prefs()
    {
        pgmProp.setProperty("initA", "" + A);
        pgmProp.setProperty("initx", "" + x);
        pgmProp.setProperty("inity", "" + y);
        pgmProp.setProperty("keyincr", "" + keyincr);
        pgmProp.setProperty("theta0", "" + theta0);
        pgmProp.setProperty("w0", "" + w0);
        pgmProp.setProperty("x0", "" + x0);
        pgmProp.setProperty("xa", "" + xa);
        pgmProp.setProperty("y0", "" + y0);
        pgmProp.setProperty("c", "" + c);
        pgmProp.setProperty("Tx", "" + Tx);
        pgmProp.setProperty("phi0", "" + phi0);
        pgmProp.setProperty("NLimit", "" + NLimit);
        pgmProp.setProperty("ystart", "" + ystart);
        pgmProp.setProperty("yend", "" + yend);
        pgmProp.setProperty("dtheta0dy", "" + dtheta0dy);
        pgmProp.setProperty("dw0dy", "" + dw0dy);
        try
            {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "ZCMPrefs.ini"), "Zeeman Catastrophe Machine v" + VERSION_NO + " Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }

    protected static String getInfo()
    {
	return "java.version = " + System.getProperty("java.version") + "\n" +
               "java.home = " + System.getProperty("java.home") + "  \n" +
               "user.dir = " + System.getProperty("user.dir") + "  \n" +
               "user.home = " + System.getProperty("user.home") + "\n" +
               "os.name = " + System.getProperty("os.name") + " - " +
               System.getProperty("os.arch") + " - " +
               System.getProperty("os.version") + "\n\n" +
               "Preferences File = " + System.getProperty("user.home") + System.getProperty("file.separator") + "ZCMPrefs.ini\n" +
               "Output Data File = " + System.getProperty("user.home") + System.getProperty("file.separator") + "ZCM_Output.csv\n";
    }

    protected static double solve_for_critical(double th, double A, double x, double y)
    {
        // calculate critical theta at a given (x, y) (Newton-Raphson)
        // th = initial estimate of theta, f = dF/dtheta

        double f, fprime, del_th;
        int tries = 0;
        do
        {
        int loop = 0;
        //System.out.println("iter =, " + loop);
        do
        {
            f = calc_dFdth(th, A, x, y);             // extremum
            fprime = calc_d2Fdth2(th, A, x, y);
            //f = calc_d2Fdth2(th, x, y);             // inflection point
            //fprime = calc_d3Fdth3(th, x, y);
            del_th = -f/fprime;
            th += del_th;
            loop++;
            //System.out.println("iter =, " + loop + ", " + th + ", " + x + ", " + y + ", " + f + ", " + fprime);
        } while (Math.abs(del_th) > TOL && loop < 100);
        if (loop == 100)
            System.out.println("BAD many loops   : " + tries + ", " + loop + ", " + f + ", " + fprime);
        if (fprime < 0)
            System.out.println("BAD data (saddle): " + tries + ", " + loop + ", " + f + ", " + fprime);
        if (loop == 100 || fprime < 0)
        {
            tries++;
            //for (int i = 0; i <= 360; i++)
            //    System.out.println(i + ", " + calc_F(i*Math.PI/180, x, y) + ", " + calc_dFdth(i*Math.PI/180, x, y) + ", " + calc_d2Fdth2(i*Math.PI/180, x, y) + ", " + calc_d3Fdth3(i*Math.PI/180, x, y));
            //th += -2*fprime/calc_d3Fdth3(th, x, y);
            int index = 0;                      // scan for minimum F
            double minF = 1000000;
            for (int i = 0; i < 36; i++)
                if (calc_F((0.5 + i)*Math.PI/18, A, x, y) < minF)
                {
                    index = i;
                    minF = calc_F((0.5 + i)*Math.PI/18, A, x, y);
                }
            th = (0.5 + index)*Math.PI/18;
            //System.out.println("restart: " + tries + ", " + loop + ", " + f + ", " + fprime + ", " + th);
        }
        else
            tries--;
        } while (tries == 1);
        if (tries > 1) return Double.NaN;
        th = th % (2.0*Math.PI);
        if (th < 0) th += 2.0*Math.PI;
        //System.out.println("theta =, " + 180*th/Math.PI + ", " + x + ", " + y + ", " + f + ", " + fprime);
        return th;
    }

    private static void get_boundary(int i)
    {
        // read pre-calculated boundary from a resource file
        String[] vals;
        Properties fileProp = new Properties();
        fileProp.clear();
        try
        {
            fileProp.load(new BufferedInputStream(main.class.getResourceAsStream("images/boundary" + (i + 3) + ".rsc")));
            vals = fileProp.getProperty("xbound").split(",");
            xbound[i] = new double[vals.length];
            for (int iLoop = 0; iLoop < vals.length; iLoop++)
                xbound[i][iLoop] = Double.parseDouble(vals[iLoop].trim());
            vals = fileProp.getProperty("ybound").split(",");
            ybound[i] = new double[vals.length];
            for (int iLoop = 0; iLoop < vals.length; iLoop++)
                ybound[i][iLoop] = Double.parseDouble(vals[iLoop].trim());
            if (xbound[i].length != ybound[i].length)
                JOptionPane.showMessageDialog(null, "length of xbound and ybound data is not the same", " Warning: bad data ", JOptionPane.WARNING_MESSAGE);
        }
        catch (IOException e)
        {
            JOptionPane.showMessageDialog(null, "images/boundary.rsc : " + e.getMessage(), " Could not load resource ", JOptionPane.WARNING_MESSAGE);
        }
    }

    private static double solve_for_boundary(double th, double A, double x, double y)
    {
        // calculate theta and x such that both d2Fdtheta2 = 0 and dFdtheta = 0
        // th, x = initial estimates, subsequently recalulated
        // fprime*del = - f

        double f0, f1;                  // constraint i = (d2Fdth2, dFdth)
        double fp00, fp01, fp10, fp11;  // variable   j = (theta, x)
        double del_th, del_x;
        int loop = 0;
        //System.out.println("init , " + 180*th/Math.PI + ", " + x + ", " + y + ", " + calc_d2Fdth2(th, x, y) + ", " + calc_dFdth(th, x, y));

        do
        {
            //System.out.println(loop + "   , " + th + ", " + x + ", " + y + ", " + calc_d2Fdth2(th, x, y) + ", " + calc_dFdth(th, x, y));
            f0 = calc_d2Fdth2(th, A, x, y);
            f1 = calc_dFdth(th, A, x, y);
            fp00 = calc_d3Fdth3(th, A, x, y);
            fp01 = calc_d3Fdth2dx(th, A, x, y);
            fp10 = calc_d2Fdth2(th, A, x, y);
            fp11 = calc_d2Fdthdx(th, A, x, y);
            if (loop > 100)
            {
                System.out.println("too many loops = " + loop);
                return Double.NaN;
            }
            del_th = -( fp11*f0 - fp01*f1)/(fp00*fp11 - fp01*fp10);
            del_x  = -(-fp10*f0 + fp00*f1)/(fp00*fp11 - fp01*fp10);
            th += del_th;
            x  += del_x;
            loop++;
            //System.out.println("iter           =, " + th + ", " + f + ", " + fprime);
        } while (Math.abs(del_th) > TOL || Math.abs(del_x) > TOL);
        System.out.println(loop + "   , " + 180*th/Math.PI + ", " + A + ", " + x + ", " + y + ", " + calc_d2Fdth2(th, A, x, y) + ", " + calc_dFdth(th, A, x, y));
        return th;
    }

    protected static Point2D.Double runge_kutta(int itime, double th, double w, double delt, double y)
    {
        double xt;
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;

        xt = x0 + xa*Math.cos(2*Math.PI*itime*delt/Tx + phi0);
        k1 = delt*w;
        l1 = delt*(-c*calc_dFdth(th, A, xt, y) -w);

        xt = x0 + xa*Math.cos(2*Math.PI*(itime*delt + delt/2)/Tx + phi0);
        k2 = delt*(w + l1/2);
        l2 = delt*(-c*calc_dFdth(th + k1/2, A, xt, y) -(w + l1/2));

        k3 = delt*(w + l2/2);
        l3 = delt*(-c*calc_dFdth(th + k2/2, A, xt, y) -(w + l2/2));

        xt = x0 + xa*Math.cos(2*Math.PI*(itime*delt + delt)/Tx + phi0);
        k4 = delt*(w + l3);
        l4 = delt*(-c*calc_dFdth(th + k3, A, xt, y) -(w + l3));
        return new Point2D.Double(th + (k1 + 2*k2 + 2*k3 + k4)/6, w + (l1 + 2*l2 + 2*l3 + l4)/6);
    }

    protected static Point2D.Double runge_kutta_4(int itime, double[] pt4, double delt, double y)
    {
        double xt;
        double th = pt4[0];
        double w = pt4[1];
        double dthdy = pt4[2];
        double dwdy = pt4[3];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;
        double n1, n2, n3, n4;

        xt = x0 + xa*Math.cos(2*Math.PI*itime*delt/Tx + phi0);
        k1 = delt*w;
        l1 = delt*(-c*calc_dFdth(th, A, xt, y) -w);
        m1 = delt*dwdy;
        n1 = delt*(-c*calc_d2Fdthdy(th, A, xt, y)
                   -c*calc_d2Fdth2(th, A, xt, y)*dthdy - dwdy);

        xt = x0 + xa*Math.cos(2*Math.PI*(itime*delt + delt/2)/Tx + phi0);
        k2 = delt*(w + l1/2);
        l2 = delt*(-c*calc_dFdth(th + k1/2, A, xt, y) - (w + l1/2));
        m2 = delt*(dwdy + n1/2);
        n2 = delt*(-c*calc_d2Fdthdy(th + k1/2, A, xt, y)
                   -c*calc_d2Fdth2(th + k1/2, A, xt, y)*(dthdy + m1/2) - (dwdy + n1/2));

        k3 = delt*(w + l2/2);
        l3 = delt*(-c*calc_dFdth(th + k2/2, A, xt, y) - (w + l2/2));
        m3 = delt*(dwdy + n2/2);
        n3 = delt*(-c*calc_d2Fdthdy(th + k2/2, A, xt, y)
                   -c*calc_d2Fdth2(th + k2/2, A, xt, y)*(dthdy + m2/2) - (dwdy + n2/2));

        xt = x0 + xa*Math.cos(2*Math.PI*(itime*delt + delt)/Tx + phi0);
        k4 = delt*(w + l3);
        l4 = delt*(-c*calc_dFdth(th + k3, A, xt, y) - (w + l3));
        m4 = delt*(dwdy + n3);
        n4 = delt*(-c*calc_d2Fdthdy(th + k3, A, xt, y)
                   -c*calc_d2Fdth2(th + k3, A, xt, y)*(dthdy + m3) - (dwdy + n3));
        pt4[0] = th + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt4[1] = w + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt4[2] = dthdy + (m1 + 2*m2 + 2*m3 + m4)/6;
        pt4[3] = dwdy + (n1 + 2*n2 + 2*n3 + n4)/6;
        //System.out.println(itime + ", " + (th + (k1 + 2*k2 + 2*k3 + k4)/6) + ", " + (w + (l1 + 2*l2 + 2*l3 + l4)/6));
        return new Point2D.Double(pt4[2], pt4[3]);
    }

    private static double e(double th, double x, double y)
    {
        return Math.sqrt((x - Math.sin(th))*(x - Math.sin(th)) + (y + Math.cos(th))*(y + Math.cos(th)));
    }

    private static double dedth(double th, double x, double y)
    {
        return -(x*Math.cos(th) + y*Math.sin(th))/e(th, x, y);
    }

    private static double d2edth2(double th, double x, double y)
    {
        return (x*Math.sin(th) - y*Math.cos(th) - dedth(th, x, y)*dedth(th, x, y))/e(th, x, y);
    }

    private static double d3edth3(double th, double x, double y)
    {
        return -dedth(th, x, y)*(e(th, x, y) + 3*d2edth2(th, x, y))/e(th, x, y);
    }

    private static double dedx(double th, double x, double y)
    {
        return (x - Math.sin(th))/e(th, x, y);
    }

    private static double dedy(double th, double x, double y)
    {
        return (y + Math.cos(th))/e(th, x, y);
    }

    private static double d2edthdx(double th, double x, double y)
    {
        return (-Math.cos(th) - dedth(th, x, y)*dedx(th, x, y))/e(th, x, y);
    }

    private static double d2edthdy(double th, double x, double y)
    {
        return (-Math.sin(th) - dedth(th, x, y)*dedy(th, x, y))/e(th, x, y);
    }

    private static double d3edth2dx(double th, double x, double y)
    {
        return (Math.sin(th) - d2edth2(th, x, y)*dedx(th, x, y) - 2*dedth(th, x, y)*d2edthdx(th, x, y))/e(th, x, y);
    }

    protected static double calc_F(double th, double A, double x, double y)
    {
        return 0.5*(e(th, 0, -A) - A/2)*(e(th, 0, -A) - A/2) + 0.5*(e(th, x, y) - A/2)*(e(th, x, y) - A/2);
    }

    protected static double calc_dFdth(double th, double A, double x, double y)
    {
        return (e(th, 0, -A) - A/2)*dedth(th, 0, -A) + (e(th, x, y) - A/2)*dedth(th, x, y);
    }

    private static double calc_d2Fdthdy(double th, double A, double x, double y)
    {
        return dedy(th, x, y)*dedth(th, x, y)   + (e(th, x, y) - A/2)*d2edthdy(th, x, y);
    }

    private static double calc_dFdx(double th, double A, double x, double y)
    {
        return (e(th, x, y) - A/2)*dedx(th, x, y);
    }

    private static double calc_d2Fdthdx(double th, double A, double x, double y)
    {
        return dedth(th, x, y)*dedx(th, x, y) + (e(th, x, y) - A/2)*d2edthdx(th, x, y);
    }

    private static double calc_d2Fdth2(double th, double A, double x, double y)
    {
        return dedth(th, 0, -A)*dedth(th, 0, -A) + (e(th, 0, -A) - A/2)*d2edth2(th, 0, -A)
             + dedth(th, x,  y)*dedth(th, x,  y) + (e(th, x,  y) - A/2)*d2edth2(th, x,  y);
    }

    private static double calc_d3Fdth3(double th, double A, double x, double y)
    {
        return 3*dedth(th, 0, -A)*d2edth2(th, 0, -A) + (e(th, 0, -A) - A/2)*d3edth3(th, 0, -A)
             + 3*dedth(th, x,  y)*d2edth2(th, x,  y) + (e(th, x,  y) - A/2)*d3edth3(th, x,  y);
    }

    private static double calc_d3Fdth2dx(double th, double A, double x, double y)
    {
        return d2edth2(th, x, y)*dedx(th, x, y) + 2*dedth(th, x, y)*d2edthdx(th, x, y) + (e(th, x, y) - A/2)*d3edth2dx(th, x, y);
    }
}
