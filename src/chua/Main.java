
// ChuaOscillator - Dec 2022 - Alvin Penner - penner@vaxxine.com

// simulate Chua Oscillator in chaos
// see Huang, Eq. 6, with k = 1
// see book Chaos IV, p. 50

// xdot = alpha*y - alpha*(a*x*x*x + c*x)
// ydot = x - y + z
// zdot = -beta*y - gamma*z

// this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Main.java

package chua;

import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.awt.geom.Point2D;
import javax.imageio.ImageIO;
import javax.swing.*;
import java.io.*;
import java.util.Properties;

public class Main
{
    private static Properties pgmProp = new Properties();
    protected static boolean skew_transform = false;
    protected static final double a = -1;               // parameters
    protected static double alpha, beta, gamma, c;      // parameters
    protected static double alpha_s, alpha_e;           // bifurcate range
    protected static double beta_s, beta_e;             // bifurcate range
    protected static double gamma_s, gamma_e;           // bifurcate range
    protected static double c_s, c_e;                   // bifurcate range
    protected static double x0, y0, z0;                 // initial values
    protected static double dx0du, dy0du, dz0du;        // initial values
    protected static double xc, yc, zc;                 // x', y', z' crossover value for x-y scatter (transformed)
    protected static double xmin, xmax, ymin, ymax;     // (x', y') plot range for x_y_scatter (transformed)
    protected static double project_phi, project_theta; // Euler angles
    protected static double project_psi;
    protected static double final_x, final_y, final_z, final_delt;
    protected static int final_Period;
    protected static int iTinit;                        // initiallize the line number in 'fit_linear_response'
    protected static double final_Re_V21, final_Im_V21; // skew transform first-order response to 'normal' form
    protected static String type;
    protected static final String VERSION_NO = "0.1";
    protected static Chua_y_vs_x plot_y_vs_x;
    protected static Chua_z_bifurcate z_bifurcate;
    protected static Chua_x_y_scatter x_y_scatter;
    protected static Chua_Euler_Slider chua_euler_slider;
    private static double[][] arrinv = new double[][] {{ 2,-16, 16,-2},
                                                       {-1, 16, 16,-1},
                                                       {-2,  4, -4, 2},
                                                       { 1, -4, -4, 1}};

    public static void main(String[] args)
    {
        load_prefs();
        if (type.equals("phase"))
        {
            plot_y_vs_x = new Chua_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
            chua_euler_slider = new Chua_Euler_Slider(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        }
        else if (type.equals("bifurcate"))
            z_bifurcate = new Chua_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        else
        {
            x_y_scatter = new Chua_x_y_scatter(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
            chua_euler_slider = new Chua_Euler_Slider(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        }
        //double qa = 2, qb = 3.7654, qc = 1.276;
        //double qx0 = 1.23, qx1 = 2.15, qx2 = 7.21;
        //System.out.println(parabola(qx0, qx1, qx2, qa + qb*qx0 + qc*qx0*qx0, qa + qb*qx1 + qc*qx1*qx1, qa + qb*qx2 + qc*qx2*qx2, true));
        //System.out.println(parabola(qx0, qx1, qx2, qa + qb*qx0 + qc*qx0*qx0, qa + qb*qx1 + qc*qx1*qx1, qa + qb*qx2 + qc*qx2*qx2, false));
    }

    protected static void runge_kutta_chua3(double[] pt3, double delt)
    {
        double x = pt3[0];
        double y = pt3[1];
        double z = pt3[2];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;

        k1 = delt*alpha*(y - a*x*x*x - c*x);
        l1 = delt*(x - y + z);
        m1 = delt*(-beta*y - gamma*z);

        k2 = delt*alpha*(y + l1/2 - a*(x + k1/2)*(x + k1/2)*(x + k1/2) - c*(x + k1/2));
        l2 = delt*(x + k1/2 - y - l1/2 + z + m1/2);
        m2 = delt*(-beta*(y + l1/2) - gamma*(z + m1/2));

        k3 = delt*alpha*(y + l2/2 - a*(x + k2/2)*(x + k2/2)*(x + k2/2) - c*(x + k2/2));
        l3 = delt*(x + k2/2 - y - l2/2 + z + m2/2);
        m3 = delt*(-beta*(y + l2/2) - gamma*(z + m2/2));

        k4 = delt*alpha*(y + l3 - a*(x + k3)*(x + k3)*(x + k3) - c*(x + k3));
        l4 = delt*(x + k3 - y - l3 + z + m3);
        m4 = delt*(-beta*(y + l3) - gamma*(z + m3));

        pt3[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt3[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt3[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
    }

    protected static void runge_kutta_chua6_ddu3(double[] pt6, double delt)
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

        k1 = delt*alpha*(y - a*x*x*x - c*x);
        l1 = delt*(x - y + z);
        m1 = delt*(-beta*y - gamma*z);
        dk1 = delt*alpha*(dy - 3*a*x*x*dx - c*dx);
        dl1 = delt*(dx - dy + dz);
        dm1 = delt*(-beta*dy - gamma*dz);

        k2 = delt*alpha*(y + l1/2 - a*(x + k1/2)*(x + k1/2)*(x + k1/2) - c*(x + k1/2));
        l2 = delt*(x + k1/2 - y - l1/2 + z + m1/2);
        m2 = delt*(-beta*(y + l1/2) - gamma*(z + m1/2));
        dk2 = delt*alpha*(dy + dl1/2 - 3*a*(x + k1/2)*(x + k1/2)*(dx + dk1/2) - c*(dx + dk1/2));
        dl2 = delt*(dx + dk1/2 - dy - dl1/2 + dz + dm1/2);
        dm2 = delt*(-beta*(dy + dl1/2) - gamma*(dz + dm1/2));

        k3 = delt*alpha*(y + l2/2 - a*(x + k2/2)*(x + k2/2)*(x + k2/2) - c*(x + k2/2));
        l3 = delt*(x + k2/2 - y - l2/2 + z + m2/2);
        m3 = delt*(-beta*(y + l2/2) - gamma*(z + m2/2));
        dk3 = delt*alpha*(dy + dl2/2 - 3*a*(x + k2/2)*(x + k2/2)*(dx + dk2/2) - c*(dx + dk2/2));
        dl3 = delt*(dx + dk2/2 - dy - dl2/2 + dz + dm2/2);
        dm3 = delt*(-beta*(dy + dl2/2) - gamma*(dz + dm2/2));

        k4 = delt*alpha*(y + l3 - a*(x + k3)*(x + k3)*(x + k3) - c*(x + k3));
        l4 = delt*(x + k3 - y - l3 + z + m3);
        m4 = delt*(-beta*(y + l3) - gamma*(z + m3));
        dk4 = delt*alpha*(dy + dl3 - 3*a*(x + k3)*(x + k3)*(dx + dk3) - c*(dx + dk3));
        dl4 = delt*(dx + dk3 - dy - dl3 + dz + dm3);
        dm4 = delt*(-beta*(dy + dl3) - gamma*(dz + dm3));

        pt6[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt6[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt6[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
        pt6[3] = dx + (dk1 + 2*dk2 + 2*dk3 + dk4)/6;
        pt6[4] = dy + (dl1 + 2*dl2 + 2*dl3 + dl4)/6;
        pt6[5] = dz + (dm1 + 2*dm2 + 2*dm3 + dm4)/6;
    }

    protected static Point2D.Double project_2D(double x, double y, double z)
    {
        double phi_rad = project_phi*Math.PI/180;
        double theta_rad = project_theta*Math.PI/180;
        double psi_rad = project_psi*Math.PI/180;
        double xp =  x*Math.cos(phi_rad)
                  +  y*Math.sin(phi_rad);
        double yp = -x*Math.sin(phi_rad)*Math.cos(theta_rad)
                  +  y*Math.cos(phi_rad)*Math.cos(theta_rad)
                  +  z*Math.sin(theta_rad);
        double xpp =  Math.cos(psi_rad)*xp + Math.sin(psi_rad)*yp;
        double ypp = -Math.sin(psi_rad)*xp + Math.cos(psi_rad)*yp;
        if (skew_transform)              // transform first-order response into 'normal' form
        {
            // apply skew transform V_inverse from Book Chaos III, p. 66
            //double[][] xy2D = new double[][] {{ -1.7841497755881142 , 0.6713612244348759 }, { -6.973231134200858 , 2.0764691595691893 }};
            //double alpha = (xy2D[1][1] - xy2D[0][0])/2/xy2D[0][1];
            //double beta = Math.sqrt(-(xy2D[1][1] - xy2D[0][0])*(xy2D[1][1] - xy2D[0][0]) - 4*xy2D[0][1]*xy2D[1][0])/2/xy2D[0][1];
            //System.out.println("xy2D, " + xy2D[0][0] + ", " + xy2D[0][1] + ", " + xy2D[1][0] + ", " + xy2D[1][1] + ", " + alpha + ", " + beta);
            xp = (final_Re_V21 - final_Im_V21)*xpp - ypp;
            yp = (-final_Re_V21 - final_Im_V21)*xpp + ypp;
            xpp = -xp/2.0/final_Im_V21;
            ypp = -yp/2.0/final_Im_V21;
        }
        return new Point2D.Double(xpp, ypp);
    }

    protected static double project_zp(double x, double y, double z)      // transformed z'
    {
        double phi_rad = project_phi*Math.PI/180;
        double theta_rad = project_theta*Math.PI/180;
        return (x*Math.sin(phi_rad) - y*Math.cos(phi_rad))*Math.sin(theta_rad) + z*Math.cos(theta_rad);
    }

    protected static double invert_from_xp_yp(double xp, double yp, double zp, String to)     // back transform from (x', y', z')
    {
        double phi_rad = project_phi*Math.PI/180;
        double theta_rad = project_theta*Math.PI/180;
        double psi_rad = project_psi*Math.PI/180;
        double [][] T = new double[][] {{ Math.cos(phi_rad),
                                          Math.sin(phi_rad),
                                          0},
                                        {-Math.sin(phi_rad)*Math.cos(theta_rad),
                                          Math.cos(phi_rad)*Math.cos(theta_rad),
                                          Math.sin(theta_rad)},
                                        { Math.sin(phi_rad)*Math.sin(theta_rad),
                                         -Math.cos(phi_rad)*Math.sin(theta_rad),
                                          Math.cos(theta_rad)}};
        double xpp, ypp;
        if (skew_transform)
        {
            // apply skew transform V from Book Chaos III, p. 66
            //double[][] xy2D = new double[][] {{ -1.7841497755881142 , 0.6713612244348759 }, { -6.973231134200858 , 2.0764691595691893 }};
            //double alpha = (xy2D[1][1] - xy2D[0][0])/2/xy2D[0][1];
            //double beta = Math.sqrt(-(xy2D[1][1] - xy2D[0][0])*(xy2D[1][1] - xy2D[0][0]) - 4*xy2D[0][1]*xy2D[1][0])/2/xy2D[0][1];
            //System.out.println("xy2D, " + xy2D[0][0] + ", " + xy2D[0][1] + ", " + xy2D[1][0] + ", " + xy2D[1][1] + ", " + alpha + ", " + beta);
            xpp = xp + yp;
            ypp = (final_Re_V21 + final_Im_V21)*xp + (final_Re_V21 - final_Im_V21)*yp;
            xp = xpp;
            yp = ypp;
        }
        xpp = Math.cos(psi_rad)*xp - Math.sin(psi_rad)*yp;
        ypp = Math.sin(psi_rad)*xp + Math.cos(psi_rad)*yp;
        if (to.equals("x"))
            return T[0][0]*xpp + T[1][0]*ypp + T[2][0]*zp;
        else if (to.equals("y"))
            return T[0][1]*xpp + T[1][1]*ypp + T[2][1]*zp;
        else if (to.equals("z"))
            return T[0][2]*xpp + T[1][2]*ypp + T[2][2]*zp;
        else
            return Double.NaN;
    }

    protected static Point2D.Double project_stationary()            // stationary points (x', y')
    {
        double xstat = Math.sqrt((1/(1 + beta/gamma) - c)/a);
        double ystat = xstat/(1 + beta/gamma);
        double zstat = -beta*ystat/gamma;
        return project_2D(xstat, ystat, zstat);
    }

    protected static double calc_xdot(double x, double y, double z)      // Chua xdot
    {
        return alpha*(y - a*x*x*x - c*x);
    }

    protected static double calc_ydot(double x, double y, double z)      // Chua ydot
    {
        return x - y + z;
    }

    protected static double calc_zdot(double x, double y, double z)      // Chua zdot
    {
        return -beta*y - gamma*z;
    }

    protected static double calc_x2dot(double x, double y, double z)      // Chua xdot
    {
        return alpha*(calc_ydot(x, y, z) - (3*a*x*x + c)*calc_xdot(x, y, z));
    }

    protected static double calc_y2dot(double x, double y, double z)      // Chua ydot
    {
        return calc_xdot(x, y, z) - calc_ydot(x, y, z) + calc_zdot(x, y, z);
    }

    protected static double calc_z2dot(double x, double y, double z)      // Chua zdot
    {
        return -beta*calc_ydot(x, y, z) - gamma*calc_zdot(x, y, z);
    }

    protected static double calc_phi(double x, double y, double z)       // Euler phi
    {
        return Math.atan2(calc_xdot(x, y, z), -calc_ydot(x, y, z))*180/Math.PI;
    }

    protected static double calc_theta(double x, double y, double z)     // Euler theta
    {
        double v = Math.sqrt(calc_xdot(x, y, z)*calc_xdot(x, y, z) + calc_ydot(x, y, z)*calc_ydot(x, y, z) + calc_zdot(x, y, z)*calc_zdot(x, y, z));
        return Math.acos(calc_zdot(x, y, z)/v)*180/Math.PI;
    }
/*
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
*/
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

    protected static double solve_cubic(double p, double q, double r)
    {
        // see Math CRC book, page 392

        double cua = (3*q - p*p)/3;
        double cub = (2*p*p*p - 9*p*q + 27*r)/27;
        double cud = cub*cub/4 + cua*cua*cua/27;
        double ret;

//        System.out.println("\ncubic p,q,r = " + p + ", " + q + ", " + r);
        //System.out.println("\ncubic a,b,d = " + cua + ", " + cub + ", " + cud);
        if (cud < 0)
        {
            double myphi = Math.acos(-cub/2/Math.sqrt(-cua*cua*cua/27));
            //System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3));
            ret = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3) - p/3;
            if (Math.abs(ret) < 1)
                return ret;
            else
            {
                ret = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3) - p/3;
                if (Math.abs(ret) < 1)
                    return ret;
                else
                {
                    ret = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3) - p/3;
                    if (Math.abs(ret) < 1)
                        return ret;
                    else
                        return Double.NaN;
                }
            }
        }
        else
        {
            //System.out.println("1 cubic d > 0 : " + (Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3));
            return Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3;
        }
    }

    protected static void load_prefs()
    {
        try                                         // recall program properties
        {
            if (new File(System.getProperty("user.home"), "ChuaPrefs.ini").exists())
            {
                pgmProp.load(new FileInputStream(new File(System.getProperty("user.home"), "ChuaPrefs.ini")));
                type = pgmProp.getProperty("type", "phase");
                alpha = Double.parseDouble(pgmProp.getProperty("alpha", "0.2"));
                beta = Double.parseDouble(pgmProp.getProperty("beta", "0.2"));
                gamma = Double.parseDouble(pgmProp.getProperty("gamma", "2.5"));
                c = Double.parseDouble(pgmProp.getProperty("c", "0.144"));
                alpha_s = Double.parseDouble(pgmProp.getProperty("alpha_s", "2"));
                alpha_e = Double.parseDouble(pgmProp.getProperty("alpha_e", "NaN"));
                beta_s = Double.parseDouble(pgmProp.getProperty("beta_s", "2"));
                beta_e = Double.parseDouble(pgmProp.getProperty("beta_e", "NaN"));
                gamma_s = Double.parseDouble(pgmProp.getProperty("gamma_s", "2"));
                gamma_e = Double.parseDouble(pgmProp.getProperty("gamma_e", "NaN"));
                c_s = Double.parseDouble(pgmProp.getProperty("c_s", "0.1"));
                c_e = Double.parseDouble(pgmProp.getProperty("ca_e", "NaN"));
                x0 = Double.parseDouble(pgmProp.getProperty("x0", "0"));
                y0 = Double.parseDouble(pgmProp.getProperty("y0", "0"));
                z0 = Double.parseDouble(pgmProp.getProperty("z0", "0"));
                xc = Double.parseDouble(pgmProp.getProperty("xc", "0"));
                yc = Double.parseDouble(pgmProp.getProperty("yc", "0"));
                zc = Double.parseDouble(pgmProp.getProperty("zc", "1"));
                xmin = Double.parseDouble(pgmProp.getProperty("xmin", "-2"));
                xmax = Double.parseDouble(pgmProp.getProperty("xmax", "2"));
                ymin = Double.parseDouble(pgmProp.getProperty("ymin", "-2"));
                ymax = Double.parseDouble(pgmProp.getProperty("ymax", "2"));
                dx0du = Double.parseDouble(pgmProp.getProperty("dx0du", "0"));
                dy0du = Double.parseDouble(pgmProp.getProperty("dy0du", "0"));
                dz0du = Double.parseDouble(pgmProp.getProperty("dz0du", "0"));
                project_phi = Double.parseDouble(pgmProp.getProperty("project_phi", "0"));
                project_theta = Double.parseDouble(pgmProp.getProperty("project_theta", "0"));
                project_psi = Double.parseDouble(pgmProp.getProperty("project_psi", "0"));
                iTinit = Integer.parseInt(pgmProp.getProperty("iTinit", "0"));
            }
            else                                    // factory default
            {
                type = "phase";
                alpha = 0.2;
                beta = 0.2;
                gamma = 2.5;
                c = 0.144;
                alpha_s = 2;
                alpha_e = Double.NaN;
                beta_s = 2;
                beta_e = Double.NaN;
                gamma_s = 2;
                gamma_e = Double.NaN;
                c_s = 0.1;
                c_e = Double.NaN;
                x0 = 0;
                y0 = 0;
                z0 = 0;
                xc = 0;
                yc = 0;
                zc = 1;
                xmin = -2;
                xmax = 2;
                ymin = -2;
                ymax = 2;
                dx0du = 0;
                dy0du = 0;
                dz0du = 0;
                project_phi = 0;
                project_theta = 0;
                project_psi = 0;
                iTinit = 0;
            }
        }
        catch (IOException e)
            {System.out.println("error reading ChuaPrefs.ini : " + e);}
    }

    protected static void save_prefs()
    {
        pgmProp.setProperty("type", type);
        pgmProp.setProperty("alpha", "" + alpha);
        pgmProp.setProperty("beta", "" + beta);
        pgmProp.setProperty("gamma", "" + gamma);
        pgmProp.setProperty("c", "" + c);
        pgmProp.setProperty("alpha_s", "" + alpha_s);
        pgmProp.setProperty("alpha_e", "" + alpha_e);
        pgmProp.setProperty("beta_s", "" + beta_s);
        pgmProp.setProperty("beta_e", "" + beta_e);
        pgmProp.setProperty("gamma_s", "" + gamma_s);
        pgmProp.setProperty("gamma_e", "" + gamma_e);
        pgmProp.setProperty("c_s", "" + c_s);
        pgmProp.setProperty("c_e", "" + c_e);
        pgmProp.setProperty("x0", "" + x0);
        pgmProp.setProperty("y0", "" + y0);
        pgmProp.setProperty("z0", "" + z0);
        pgmProp.setProperty("xc", "" + xc);
        pgmProp.setProperty("yc", "" + yc);
        pgmProp.setProperty("zc", "" + zc);
        pgmProp.setProperty("xmin", "" + xmin);
        pgmProp.setProperty("xmax", "" + xmax);
        pgmProp.setProperty("ymin", "" + ymin);
        pgmProp.setProperty("ymax", "" + ymax);
        pgmProp.setProperty("dx0du", "" + dx0du);
        pgmProp.setProperty("dy0du", "" + dy0du);
        pgmProp.setProperty("dz0du", "" + dz0du);
        pgmProp.setProperty("project_phi", "" + project_phi);
        pgmProp.setProperty("project_theta", "" + project_theta);
        pgmProp.setProperty("project_psi", "" + project_psi);
        pgmProp.setProperty("iTinit", "" + iTinit);
        try
            {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "ChuaPrefs.ini"), "Chua System v" + VERSION_NO + " Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }

    protected static void save_PNG(BufferedImage img, String label)
    {
        try
        {
            ImageIO.write(img, "png", new File("\\Windows\\Temp", label + ".png"));
            System.out.println("save PNG: " + label);
        }
        catch (IOException e)
            {System.out.println("save_PNG error = " + e);}
    }
/*
    private static double[] gaussj(double[][] m, double[] v)
    {
        // copied from BSpline5.java
        // solve m*x = v
        // based on the routine gaussj() in "Numerical Recipes in C", p.39, W.H.Press
        // Gauss-Jordan elimination with full pivoting
        if (m.length != m[0].length) return null;
        if (m.length != v.length) return null;
        int i, j, k, l, ll, icol = 0, irow = 0;
        int n = m.length;
        double big, dum, pivinv;
        double[][] gaussa = new double[n][n];
        double[] gaussb = new double[n];
        int[] indxc = new int[n];
        int[] indxr = new int[n];
        int[] ipiv = new int[n];
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                gaussa[i][j] = m[i][j];
        for (i = 0; i < n; i++)
            gaussb[i] = v[i];
        for (i = 0; i < n; i++)
            ipiv[i] = 0;
        for (i = 0; i < n; i++)
        {
            big = 0;
            for (j = 0; j < n; j++)
                if (ipiv[j] != 1)
                    for (k = 0; k < n; k++)
                        if (ipiv[k] == 0)
                        {
                            if (Math.abs(gaussa[j][k]) >= big)
                            {
                                big = Math.abs(gaussa[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                        else if (ipiv[k] > 1)
                        {
                            System.out.println("gaussj error : Singular Matrix 1");
                            return null;
                        }
            ++(ipiv[icol]);
            if (irow != icol)
            {
                for (l = 0; l < n; l++) swapa (gaussa, irow, l, icol, l);
                swapv(gaussb, irow, icol);
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if (gaussa[icol][icol] == 0)
            {
                System.out.println("gaussj error : Singular Matrix 2");
                return null;
            }
            pivinv = 1/gaussa[icol][icol];
            gaussa[icol][icol] = 1;
            for (l = 0; l < n; l++) gaussa[icol][l] *= pivinv;
            gaussb[icol] *= pivinv;
            for (ll = 0; ll < n; ll++)
                if (ll != icol)
                {
                    dum = gaussa[ll][icol];
                    gaussa[ll][icol] = 0;
                    for (l = 0; l < n; l++) gaussa[ll][l] -= gaussa[icol][l]*dum;
                    gaussb[ll] -= gaussb[icol]*dum;
                }
        }
        for (l = n - 1; l >= 0; l--)
            if (indxr[l] != indxc[l])
                for (k = 0; k < n; k++)
                    swapa(gaussa, k, indxr[l], k, indxc[l]);
        //System.out.println("a inverse");
        //for (i = 0; i < n; i++)
        //{
        //    for (j = 0; j < n; j++)
        //        System.out.print(a[i][j] + ", ");
        //    System.out.println();
        //}
        //return a;
        //System.out.println("b");
        //for (i = 0; i < n; i++)
        //    System.out.println(b[i]);
        return gaussb;
    }

    private static void swapa(double[][] a, int i1, int j1, int i2, int j2)
    {
        // copied from BSpline5.java
        double temp = a[i1][j1];
        a[i1][j1] = a[i2][j2];
        a[i2][j2] = temp;
    }

    private static void swapv(double[] v, int i1, int i2)
    {
        // copied from BSpline5.java
        double temp = v[i1];
        v[i1] = v[i2];
        v[i2] = temp;
    }
*/
}
