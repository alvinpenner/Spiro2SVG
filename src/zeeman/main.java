
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

import javax.swing.*;
import java.io.*;
import java.util.Properties;

public class main
{
    protected static double[][] xbound = new double[4][];
    protected static double[][] ybound = new double[4][];
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

    private static double d2edthdx(double th, double x, double y)
    {
        return (-Math.cos(th) - dedth(th, x, y)*dedx(th, x, y))/e(th, x, y);
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
