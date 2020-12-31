
// RosslerSystem - Dec 2020 - Alvin Penner - penner@vaxxine.com

// simulate Rossler system of chaos
// see Strogatz, p. 377
// OLSEN AND DEGN, Chaos in biological systems
// xdot = -y - z
// ydot = x + a*y
// zdot = b + z*(x - c)

// this is file : \Documents\NetBeansProjects\RosslerSystem\src\rossler\Main.java

package rossler;

import java.awt.*;
import javax.swing.*;
import java.io.*;
import java.util.Properties;

public class Main
{
    private static Properties pgmProp = new Properties();
    protected static double a, b, c;                        // parameters
    protected static double cstart, cend;                   // bifurcate range
    protected static double x0, y0, z0;                     // initial values
    protected static double dx0dc, dy0dc, dz0dc;            // initial values
    protected static String type;
    protected static final String VERSION_NO = "0.1";
    protected static Rossler_y_vs_x plot_y_vs_x;
    protected static Rossler_z_bifurcate z_bifurcate;
    private static double[][] arrinv = new double[][] {{ 2,-16, 16,-2},
                                                       {-1, 16, 16,-1},
                                                       {-2,  4, -4, 2},
                                                       { 1, -4, -4, 1}};
    public static void main(String[] args)
    {
        load_prefs();
        if (type.equals("phase"))
            plot_y_vs_x = new Rossler_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        else
            z_bifurcate = new Rossler_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        //double qa = 2, qb = 3.7654, qc = 1.276;
        //double qx0 = 1.23, qx1 = 2.15, qx2 = 7.21;
        //System.out.println(parabola(qx0, qx1, qx2, qa + qb*qx0 + qc*qx0*qx0, qa + qb*qx1 + qc*qx1*qx1, qa + qb*qx2 + qc*qx2*qx2, true));
        //System.out.println(parabola(qx0, qx1, qx2, qa + qb*qx0 + qc*qx0*qx0, qa + qb*qx1 + qc*qx1*qx1, qa + qb*qx2 + qc*qx2*qx2, false));
    }

    protected static void runge_kutta_rossler3(double[] pt3, double delt, double tempc)
    {
        double x = pt3[0];
        double y = pt3[1];
        double z = pt3[2];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;

        k1 = delt*(-y - z);
        l1 = delt*( x + a*y);
        m1 = delt*(b + z*(x - tempc));

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + a*(y + l1/2));
        m2 = delt*(b + (z + m1/2)*(x + k1/2 - tempc));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + a*(y + l2/2));
        m3 = delt*(b + (z + m2/2)*(x + k2/2 - tempc));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + a*(y + l3));
        m4 = delt*(b + (z + m3)*(x + k3 - tempc));

        pt3[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt3[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt3[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
    }

    protected static void runge_kutta_rossler6(double[] pt6, double delt, double tempc)
    {
        double x = pt6[0];
        double y = pt6[1];
        double z = pt6[2];
        double dx = pt6[3];
        double dy = pt6[4];
        double dz = pt6[5];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;
        double dk1, dk2, dk3, dk4;
        double dl1, dl2, dl3, dl4;
        double dm1, dm2, dm3, dm4;

        k1 = delt*(-y - z);
        l1 = delt*( x + a*y);
        m1 = delt*(b + z*(x - tempc));
        dk1 = delt*(-dy - dz);
        dl1 = delt*(dx + a*dy);
        dm1 = delt*(dz*(x - tempc) + z*(dx - 1));

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + a*(y + l1/2));
        m2 = delt*(b + (z + m1/2)*(x + k1/2 - tempc));
        dk2 = delt*(-dy - dl1/2 - dz - dm1/2);
        dl2 = delt*(dx + dk1/2 + a*(dy + dl1/2));
        dm2 = delt*((dz + dm1/2)*(x + k1/2 - tempc) + (z + m1/2)*(dx + dk1/2 - 1));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + a*(y + l2/2));
        m3 = delt*(b + (z + m2/2)*(x + k2/2 - tempc));
        dk3 = delt*(-dy - dl2/2 - dz - dm2/2);
        dl3 = delt*(dx + dk2/2 + a*(dy + dl2/2));
        dm3 = delt*((dz + dm2/2)*(x + k2/2 - tempc) + (z + m2/2)*(dx + dk2/2 - 1));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + a*(y + l3));
        m4 = delt*(b + (z + m3)*(x + k3 - tempc));
        dk4 = delt*(-dy - dl3 - dz - dm3);
        dl4 = delt*(dx + dk3 + a*(dy + dl3));
        dm4 = delt*((dz + dm3)*(x + k3 - tempc) + (z + m3)*(dx + dk3 - 1));

        pt6[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt6[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt6[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
        pt6[3] = dx + (dk1 + 2*dk2 + 2*dk3 + dk4)/6;
        pt6[4] = dy + (dl1 + 2*dl2 + 2*dl3 + dl4)/6;
        pt6[5] = dz + (dm1 + 2*dm2 + 2*dm3 + dm4)/6;
    }

    protected static double parabola(double x0, double x1, double x2, double y0, double y1, double y2, boolean ymax)
    {
        // y = A + B*x + C*x*x
        // calculate xmax, ymax at maximum
        double denom = (x2 - x1)*(x1 - x0)*(x0 - x2);
        double A = (y0*x1*x2*(x1 - x2) + y1*x2*x0*(x2 - x0) + y2*x0*x1*(x0 - x1))/denom;
        double B = (y0*(x2*x2 - x1*x1) + y1*(x0*x0 - x2*x2) + y2*(x1*x1 - x0*x0))/denom;
        double C = (y0*(x1 - x2) + y1*(x2 - x0) + y2*(x0 - x1))/denom;
        //System.out.println("A, B, C = " + A + ", " + B + ", " + C);
        if (ymax)                                   // ymax
            return A - B*B/4/C;
        else                                        // xmax
            return -B/2/C;
    }

    protected static double quarticT(double dely0, double dely1, double dely3, double dely4)
    {
        // calculate the x position of the maximum of a quartic function
        // assume x = (-2, -1, 0, 1, 2), y2 = max y (approx)
        // see Chaos II book, page 18

        double[] dely = new double[] {dely0, dely1, dely3, dely4};
        double[] coef = new double[] {0, 0, 0, 0};

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                coef[i] += arrinv[i][j]*dely[j];
        //System.out.println("coef, " + coef[0]/24 + ", " + coef[1]/24 + ", " + coef[2]/24 + ", " + coef[3]/24);
        return solve_cubic(3*coef[2]/4/coef[3], 2*coef[1]/4/coef[3], coef[0]/4/coef[3]);
    }

    private static double solve_cubic(double p, double q, double r)
    {
        // see Math CRC book, page 392

        double cua = (3*q - p*p)/3;
        double cub = (2*p*p*p - 9*p*q + 27*r)/27;
        double cud = cub*cub/4 + cua*cua*cua/27;

//        System.out.println("\ncubic p,q,r = " + p + ", " + q + ", " + r);
//        System.out.println("\ncubic a,b,d = " + cua + ", " + cub + ", " + cud);
        if (cud < 0)
        {
            double myphi = Math.acos(-cub/2/Math.sqrt(-cua*cua*cua/27));
//            System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3));
            return 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3;
        }
        else
        {
//            System.out.println("1 cubic d > 0 : " + (Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3));
            return Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3;
        }
    }

    protected static void load_prefs()
    {
        try                                         // recall program properties
        {
            if (new File(System.getProperty("user.home"), "RosslerPrefs.ini").exists())
            {
                pgmProp.load(new FileInputStream(new File(System.getProperty("user.home"), "RosslerPrefs.ini")));
                type = pgmProp.getProperty("type", "phase");
                a = Double.parseDouble(pgmProp.getProperty("a", "0.2"));
                b = Double.parseDouble(pgmProp.getProperty("b", "0.2"));
                c = Double.parseDouble(pgmProp.getProperty("c", "2.5"));
                cstart = Double.parseDouble(pgmProp.getProperty("cstart", "2"));
                cend = Double.parseDouble(pgmProp.getProperty("cend", "6"));
                x0 = Double.parseDouble(pgmProp.getProperty("x0", "0"));
                y0 = Double.parseDouble(pgmProp.getProperty("y0", "0"));
                z0 = Double.parseDouble(pgmProp.getProperty("z0", "0"));
                dx0dc = Double.parseDouble(pgmProp.getProperty("dx0dc", "0"));
                dy0dc = Double.parseDouble(pgmProp.getProperty("dy0dc", "0"));
                dz0dc = Double.parseDouble(pgmProp.getProperty("dz0dc", "0"));
            }
            else                                    // factory default
            {
                type = "phase";
                a = 0.2;
                b = 0.2;
                c = 2.5;
                cstart = 2;
                cend = 6;
                x0 = 0;
                y0 = 0;
                z0 = 0;
                dx0dc = 0;
                dy0dc = 0;
                dz0dc = 0;
            }
        }
        catch (IOException e)
            {System.out.println("error reading RosslerPrefs.ini : " + e);}
    }

    protected static void save_prefs()
    {
        pgmProp.setProperty("type", type);
        pgmProp.setProperty("a", "" + a);
        pgmProp.setProperty("b", "" + b);
        pgmProp.setProperty("c", "" + c);
        pgmProp.setProperty("cstart", "" + cstart);
        pgmProp.setProperty("cend", "" + cend);
        pgmProp.setProperty("x0", "" + x0);
        pgmProp.setProperty("y0", "" + y0);
        pgmProp.setProperty("z0", "" + z0);
        pgmProp.setProperty("dx0dc", "" + dx0dc);
        pgmProp.setProperty("dy0dc", "" + dy0dc);
        pgmProp.setProperty("dz0dc", "" + dz0dc);
        try
            {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "RosslerPrefs.ini"), "Rossler System v" + VERSION_NO + " Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }
}
