
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
    protected static String type;
    protected static final String VERSION_NO = "0.1";
    protected static Rossler_y_vs_x plot_y_vs_x;
    protected static Rossler_z_bifurcate z_bifurcate;

    public static void main(String[] args)
    {
        load_prefs();
        if (type.equals("phase"))
            plot_y_vs_x = new Rossler_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        else
            z_bifurcate = new Rossler_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
    }

    protected static void runge_kutta_rossler(double[] pt3, double delt)
    {
        double x = pt3[0];
        double y = pt3[1];
        double z = pt3[2];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;

        k1 = delt*(-y - z);
        l1 = delt*( x + a*y);
        m1 = delt*(b + z*(x - c));

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + a*(y + l1/2));
        m2 = delt*(b + (z + m1/2)*(x + k1/2 - c));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + a*(y + l2/2));
        m3 = delt*(b + (z + m2/2)*(x + k2/2 - c));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + a*(y + l3));
        m4 = delt*(b + (z + m3)*(x + k3 - c));

        pt3[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt3[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt3[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
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
        try
            {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "RosslerPrefs.ini"), "Rossler System v" + VERSION_NO + " Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }
}
