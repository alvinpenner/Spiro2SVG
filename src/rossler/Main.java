
// RosslerSystem - Dec 2020 - Alvin Penner - penner@vaxxine.com

// simulate Rossler system of chaos
// see Strogatz, p. 377
// OLSEN AND DEGN, Chaos in biological systems
// xdot = -y - z
// ydot = x + a*y
// zdot = b + z*(x - c)

// this is file : \Documents\NetBeansProjects\RosslerSystem\src\rossler\Main.java

package rossler;

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
    protected static double a, b, c;                    // parameters
    protected static double astart, aend;               // bifurcate range
    protected static double bstart, bend;               // bifurcate range
    protected static double cstart, cend;               // bifurcate range
    protected static double x0, y0, z0;                 // initial values
    protected static double dx0dc, dy0dc, dz0dc;        // initial values
    protected static double xc, yc, zc;                 // x', y', z' crossover value for x-y scatter (transformed)
    protected static double xmin, xmax, ymin, ymax;     // (x', y') plot range for x_y_scatter (transformed)
    protected static double project_phi, project_theta; // Euler angles
    protected static double project_psi;
    protected static double final_x, final_y, final_z, final_delt;
    protected static int final_Period;
    protected static double final_Re_V21, final_Im_V21; // skew transform first-order response to 'normal' form
    protected static String type;
    protected static final String VERSION_NO = "0.1";
    protected static Rossler_y_vs_x plot_y_vs_x;
    protected static Rossler_z_bifurcate z_bifurcate;
    protected static Rossler_x_y_scatter x_y_scatter;
    protected static Euler_Slider euler_slider;
    private static double[][] arrinv = new double[][] {{ 2,-16, 16,-2},
                                                       {-1, 16, 16,-1},
                                                       {-2,  4, -4, 2},
                                                       { 1, -4, -4, 1}};
    private static double[] gen_x, gen_y, gen_z;
    private static double xdot, ydot, zdot, x2dot, y2dot, z2dot;
    private static double x3dot, y3dot, phi2dot;
    private static double phi, theta, phidot, thetadot;
    private static double[][] dTdphi, dTdtheta, dTd1, dTda, dTdxc, dTdz;
    private static double[][] M;                        // solve Mx + v = 0
    private static double[][][] Mout;                   // diagonal elements of M
    private static double[][] M2D = new double[2][2];   // 2-D coupling for transformed response
    private static double[] Moutdot00, Moutdot01;       // for the second-order d.e.
    private static double[] v;

    public static void main(String[] args)
    {
        load_prefs();
        if (type.equals("phase"))
        {
            plot_y_vs_x = new Rossler_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
            euler_slider = new Euler_Slider(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        }
        else if (type.equals("bifurcate"))
            z_bifurcate = new Rossler_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        else
        {
            x_y_scatter = new Rossler_x_y_scatter(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
            euler_slider = new Euler_Slider(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        }
        //double qa = 2, qb = 3.7654, qc = 1.276;
        //double qx0 = 1.23, qx1 = 2.15, qx2 = 7.21;
        //System.out.println(parabola(qx0, qx1, qx2, qa + qb*qx0 + qc*qx0*qx0, qa + qb*qx1 + qc*qx1*qx1, qa + qb*qx2 + qc*qx2*qx2, true));
        //System.out.println(parabola(qx0, qx1, qx2, qa + qb*qx0 + qc*qx0*qx0, qa + qb*qx1 + qc*qx1*qx1, qa + qb*qx2 + qc*qx2*qx2, false));
    }

    protected static void runge_kutta_rossler3(double[] pt3, double delt, double tempa, double tempb, double tempc)
    {
        double x = pt3[0];
        double y = pt3[1];
        double z = pt3[2];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;

        k1 = delt*(-y - z);
        l1 = delt*( x + tempa*y);
        m1 = delt*(tempb + z*(x - tempc));

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + tempa*(y + l1/2));
        m2 = delt*(tempb + (z + m1/2)*(x + k1/2 - tempc));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + tempa*(y + l2/2));
        m3 = delt*(tempb + (z + m2/2)*(x + k2/2 - tempc));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + tempa*(y + l3));
        m4 = delt*(tempb + (z + m3)*(x + k3 - tempc));

        pt3[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt3[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt3[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
    }

    protected static void runge_kutta_rossler6_ddc(double[] pt6, double delt, double tempc)
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

    protected static void runge_kutta_rossler6_ddu3(double[] pt6, double delt, double tempc)
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
        dm1 = delt*(dz*(x - tempc) + z*dx);

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + a*(y + l1/2));
        m2 = delt*(b + (z + m1/2)*(x + k1/2 - tempc));
        dk2 = delt*(-dy - dl1/2 - dz - dm1/2);
        dl2 = delt*(dx + dk1/2 + a*(dy + dl1/2));
        dm2 = delt*((dz + dm1/2)*(x + k1/2 - tempc) + (z + m1/2)*(dx + dk1/2));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + a*(y + l2/2));
        m3 = delt*(b + (z + m2/2)*(x + k2/2 - tempc));
        dk3 = delt*(-dy - dl2/2 - dz - dm2/2);
        dl3 = delt*(dx + dk2/2 + a*(dy + dl2/2));
        dm3 = delt*((dz + dm2/2)*(x + k2/2 - tempc) + (z + m2/2)*(dx + dk2/2));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + a*(y + l3));
        m4 = delt*(b + (z + m3)*(x + k3 - tempc));
        dk4 = delt*(-dy - dl3 - dz - dm3);
        dl4 = delt*(dx + dk3 + a*(dy + dl3));
        dm4 = delt*((dz + dm3)*(x + k3 - tempc) + (z + m3)*(dx + dk3));

        pt6[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt6[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt6[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
        pt6[3] = dx + (dk1 + 2*dk2 + 2*dk3 + dk4)/6;
        pt6[4] = dy + (dl1 + 2*dl2 + 2*dl3 + dl4)/6;
        pt6[5] = dz + (dm1 + 2*dm2 + 2*dm3 + dm4)/6;
    }

    protected static void runge_kutta_rossler6_ddu3_full(double[] pt6, double delt, double tempc)
    {
        // (x, y, z) = limit cycle
        // (dx, dy, dz) = deviation from (x, y, z), based on full nonlinear response
        //                                          (using original limit cycle)
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
        dm1 = delt*(dz*(x - tempc + dx) + z*dx);

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + a*(y + l1/2));
        m2 = delt*(b + (z + m1/2)*(x + k1/2 - tempc));
        dk2 = delt*(-dy - dl1/2 - dz - dm1/2);
        dl2 = delt*(dx + dk1/2 + a*(dy + dl1/2));
        dm2 = delt*((dz + dm1/2)*(x + k1/2 - tempc + dx + dk1/2) + (z + m1/2)*(dx + dk1/2));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + a*(y + l2/2));
        m3 = delt*(b + (z + m2/2)*(x + k2/2 - tempc));
        dk3 = delt*(-dy - dl2/2 - dz - dm2/2);
        dl3 = delt*(dx + dk2/2 + a*(dy + dl2/2));
        dm3 = delt*((dz + dm2/2)*(x + k2/2 - tempc + dx + dk2/2) + (z + m2/2)*(dx + dk2/2));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + a*(y + l3));
        m4 = delt*(b + (z + m3)*(x + k3 - tempc));
        dk4 = delt*(-dy - dl3 - dz - dm3);
        dl4 = delt*(dx + dk3 + a*(dy + dl3));
        dm4 = delt*((dz + dm3)*(x + k3 - tempc + dx + dk3) + (z + m3)*(dx + dk3));

        pt6[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt6[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt6[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
        pt6[3] = dx + (dk1 + 2*dk2 + 2*dk3 + dk4)/6;
        pt6[4] = dy + (dl1 + 2*dl2 + 2*dl3 + dl4)/6;
        pt6[5] = dz + (dm1 + 2*dm2 + 2*dm3 + dm4)/6;
    }

    protected static void runge_kutta_rossler6_ddu2(double[] pt6, double delt, double tempc)
    {
        double x = pt6[0];
        double y = pt6[1];
        double z = pt6[2];
        double dx = pt6[3];
        double dy = pt6[4];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;
        double dk1, dk2, dk3, dk4;
        double dl1, dl2, dl3, dl4;

        k1 = delt*(-y - z);
        l1 = delt*( x + a*y);
        m1 = delt*(b + z*(x - tempc));
        gen_M2D(x, y, z);
        dk1 = -delt*(M2D[0][0]*dx + M2D[0][1]*dy);
        dl1 = -delt*(M2D[1][0]*dx + M2D[1][1]*dy);

        k2 = delt*(-y - l1/2 - z - m1/2);
        l2 = delt*( x + k1/2 + a*(y + l1/2));
        m2 = delt*(b + (z + m1/2)*(x + k1/2 - tempc));
        gen_M2D(x + k1/2, y + l1/2, z + m1/2);
        dk2 = -delt*(M2D[0][0]*(dx + dk1/2) + M2D[0][1]*(dy + dl1/2));
        dl2 = -delt*(M2D[1][0]*(dx + dk1/2) + M2D[1][1]*(dy + dl1/2));

        k3 = delt*(-y - l2/2 - z - m2/2);
        l3 = delt*( x + k2/2 + a*(y + l2/2));
        m3 = delt*(b + (z + m2/2)*(x + k2/2 - tempc));
        gen_M2D(x + k2/2, y + l2/2, z + m2/2);
        dk3 = -delt*(M2D[0][0]*(dx + dk2/2) + M2D[0][1]*(dy + dl2/2));
        dl3 = -delt*(M2D[1][0]*(dx + dk2/2) + M2D[1][1]*(dy + dl2/2));

        k4 = delt*(-y - l3 - z - m3);
        l4 = delt*( x + k3 + a*(y + l3));
        m4 = delt*(b + (z + m3)*(x + k3 - tempc));
        gen_M2D(x + k3, y + l3, z + m3);
        dk4 = -delt*(M2D[0][0]*(dx + dk3) + M2D[0][1]*(dy + dl3));
        dl4 = -delt*(M2D[1][0]*(dx + dk3) + M2D[1][1]*(dy + dl3));

        pt6[0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        pt6[1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        pt6[2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
        pt6[3] = dx + (dk1 + 2*dk2 + 2*dk3 + dk4)/6;
        pt6[4] = dy + (dl1 + 2*dl2 + 2*dl3 + dl4)/6;
    }

    protected static void runge_kutta_rossler6_test(double[] pt6, double delt, double tempc)
    {
        double dTdc = -0.079172/576.9798;             // add compensation for dTau/dc

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

        pt6[3] += (-pt6[1] - pt6[2])*dTdc;              // compensate for dTdc
        pt6[4] += (pt6[0] + a*pt6[1])*dTdc;
        pt6[5] += (b + pt6[2]*(pt6[0] - tempc))*dTdc;
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

    protected static Point2D.Double project_stationary(String sign)    // stationary points (x', y')
    {
        double zstat = Math.sqrt(c*c - 4*a*b);
        if (sign.equals("+"))
            zstat = (c + zstat)/2/a;
        else
            zstat = (c - zstat)/2/a;
        double xstat = a*zstat;
        double ystat = -zstat;
        //System.out.println("pstat = " + xstat + ", " + ystat + ", " + zstat);
        return project_2D(xstat, ystat, zstat);
    }

    private static void gen_M2D(double x, double y, double z)
    {
        xdot = -y - z;
        ydot = x + a*y;
        zdot = b + z*(x - c);
        phi = Math.atan2(xdot, -ydot);
        theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));

        M2D[0][0] = -a*Math.sin(phi)*Math.sin(phi);
        M2D[0][1] =    Math.cos(phi)/Math.sin(theta);
        M2D[1][0] = -2.0*a*Math.sin(phi)*Math.cos(phi)*Math.cos(theta) - Math.cos(phi)/Math.sin(theta)
                    - (z - 1)*Math.cos(phi)*Math.sin(theta);
        M2D[1][1] = - a*Math.cos(phi)*Math.cos(phi)*Math.cos(theta)*Math.cos(theta)
                    - (x - c)*Math.sin(theta)*Math.sin(theta)
                    + (z - 1)*Math.sin(phi)*Math.sin(theta)*Math.cos(theta);
        //System.out.println("  , " + M2D[0][0] + ", " + M2D[0][1] + ", " + M2D[1][0] + ", " + M2D[1][1]);
    }

    protected static void gen_array(int index, int Period, double delt, double x, double y, double z)
    {
        // see book Chaos II, p. 35
        dTdtheta = new double[][] {{ 0,  0,  0},
                                   { 0,  0,  1},
                                   { 0, -1,  0}};
        double av_x = 0, av_y = 0, av_z = 0;     // accumulate averages to test Genesio Eq. 17
        if (index == 0)
        {
            gen_x = new double[Period];
            gen_y = new double[Period];
            gen_z = new double[Period];
            System.out.println("\nRossler diagonal elements of M");
            System.out.println("a      = " + a);
            System.out.println("b      = " + b);
            System.out.println("c      = " + c);
            System.out.println("Period = " + Period);
            System.out.println("delt   = " + delt);
            //System.out.print("d2Fdth2 = np.array([" + M[0][0]);
            //for (int j = 1; j < N; j++)
            //    System.out.print(", " + M[0][j]);
            //System.out.println("])");
            //System.out.print("d2Fdthdy = np.array([" + v[0]);
            //for (int j = 1; j < N; j++)
            //    System.out.print(", " + v[j]);
            //System.out.println("])");
        }
        gen_x[index] = x;
        gen_y[index] = y;
        gen_z[index] = z;
        if (index == Period - 1)
        {
            v = new double[3*Period];
            for (int i = 0; i < v.length; i++)
                v[i] = 0;
            Mout = new double[3][3][Period];
            Moutdot00 = new double[Period];
            Moutdot01 = new double[Period];
            //double dTaudc = 0.11636;             // add compensation for dTau/dc
            //double dTaudc = 0.1163616;             // add compensation for dTau/dc
            double v0 = Double.NaN, v1 = Double.NaN, vN2 = Double.NaN, vN1 = Double.NaN;
            //System.out.println("i, x, y, z, xdot, ydot, zdot, phi, theta, x2dot, y2dot, z2dot, phidot, thetadot, x3dot, y3dot, phi2dot");
            System.out.println("i, x, y, z, xdot, ydot, zdot, phi, theta");
            //System.out.println("i, M00, M01, M10, M11, cos(phi), sin(theta)");
            double xmin_temp = 999, xmax_temp = -999;
            double ymin_temp = 999, ymax_temp = -999;
            double zmin_temp = 999, zmax_temp = -999;
            double thmin_temp = 999, thmax_temp = -999;
            for (int i = 0; i < Period; i++)                        // i = time index
            {
                av_x += gen_x[i];
                av_y += gen_y[i];
                av_z += gen_z[i];
                xdot = -gen_y[i] - gen_z[i];
                ydot = gen_x[i] + a*gen_y[i];
                zdot = b + gen_z[i]*(gen_x[i] - c);
                if (i == 0) v0 = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
                else if (i == 1) v1 = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
                else if (i == Period - 2) vN2 = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
                else if (i == Period - 1) vN1 = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
                phi = Math.atan2(xdot, -ydot);
                theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));
                if (gen_x[i] > xmax_temp) xmax_temp = gen_x[i];
                if (gen_x[i] < xmin_temp) xmin_temp = gen_x[i];
                if (gen_y[i] > ymax_temp) ymax_temp = gen_y[i];
                if (gen_y[i] < ymin_temp) ymin_temp = gen_y[i];
                if (gen_z[i] > zmax_temp) zmax_temp = gen_z[i];
                if (gen_z[i] < zmin_temp) zmin_temp = gen_z[i];
                if (theta > thmax_temp) thmax_temp = theta;
                if (theta < thmin_temp) thmin_temp = theta;
                dTdphi = new double[][] {{                0,  Math.cos(theta), -Math.sin(theta)},
                                         { -Math.cos(theta),                0,               0},
                                         {  Math.sin(theta),                0,               0}};
                dTd1   = new double[][] {{                0, -Math.cos(theta) - Math.cos(phi)*Math.sin(theta), Math.sin(theta) - Math.cos(phi)*Math.cos(theta)},
                                         {  Math.cos(theta),  Math.sin(phi)*Math.sin(theta)*Math.cos(theta),  Math.sin(phi)*Math.cos(theta)*Math.cos(theta)},
                                         { -Math.sin(theta), -Math.sin(phi)*Math.sin(theta)*Math.sin(theta), -Math.sin(phi)*Math.sin(theta)*Math.cos(theta)}};
                dTda   = new double[][] {{  Math.sin(phi)*Math.sin(phi)                ,  Math.sin(phi)*Math.cos(phi)*Math.cos(theta)                , -Math.sin(phi)*Math.cos(phi)*Math.sin(theta)},
                                         {  Math.sin(phi)*Math.cos(phi)*Math.cos(theta),  Math.cos(phi)*Math.cos(phi)*Math.cos(theta)*Math.cos(theta), -Math.cos(phi)*Math.cos(phi)*Math.sin(theta)*Math.cos(theta)},
                                         { -Math.sin(phi)*Math.cos(phi)*Math.sin(theta), -Math.cos(phi)*Math.cos(phi)*Math.sin(theta)*Math.cos(theta),  Math.cos(phi)*Math.cos(phi)*Math.sin(theta)*Math.sin(theta)}};
                dTdxc  = new double[][] {{                0,                                0,                               0},
                                         {                0,  Math.sin(theta)*Math.sin(theta), Math.sin(theta)*Math.cos(theta)},
                                         {                0,  Math.sin(theta)*Math.cos(theta), Math.cos(theta)*Math.cos(theta)}};
                dTdz   = new double[][] {{                              0,                                              0,                                             0},
                                         {  Math.cos(phi)*Math.sin(theta), -Math.sin(phi)*Math.sin(theta)*Math.cos(theta), Math.sin(phi)*Math.sin(theta)*Math.sin(theta)},
                                         {  Math.cos(phi)*Math.cos(theta), -Math.sin(phi)*Math.cos(theta)*Math.cos(theta), Math.sin(phi)*Math.sin(theta)*Math.cos(theta)}};
                x2dot = -ydot - zdot;
                y2dot = xdot + a*ydot;
                z2dot = zdot*(gen_x[i] - c) + gen_z[i]*xdot;
                x3dot = -y2dot - z2dot;
                y3dot = x2dot + a*y2dot;
                phidot = (xdot*y2dot - ydot*x2dot)/(xdot*xdot + ydot*ydot);
                phi2dot = ((xdot*xdot + ydot*ydot)*(xdot*y3dot - ydot*x3dot) - 2.0*(xdot*y2dot - ydot*x2dot)*(xdot*x2dot + ydot*y2dot))
                        / (xdot*xdot + ydot*ydot)/(xdot*xdot + ydot*ydot);
                thetadot = ((x2dot*xdot + y2dot*ydot)*zdot - z2dot*(xdot*xdot + ydot*ydot))/
                           Math.sqrt(xdot*xdot + ydot*ydot)/(xdot*xdot + ydot*ydot + zdot*zdot);
                //System.out.println(i + ", " + gen_x[i] + ", " + gen_y[i] + ", " + gen_z[i] + ", " + xdot + ", " + ydot + ", " + zdot + ", " + phi + ", " + theta
                //                                                                           + ", " + x2dot + ", " + y2dot + ", " + z2dot + ", " + phidot + ", " + thetadot
                //                                                                           + ", " + x3dot + ", " + y3dot + ", " + phi2dot);
                //System.out.println(i + ", " + gen_x[i] + ", " + gen_y[i] + ", " + gen_z[i] + ", " + xdot + ", " + ydot + ", " + zdot);
                System.out.println(i + ", " + gen_x[i] + ", " + gen_y[i] + ", " + gen_z[i] + ", " + xdot + ", " + ydot + ", " + zdot + ", " + phi + ", " + theta);

                // generate diagonals of M, and v

                //phidot = 0;                                 // temporary fudge FIX FIX
                //thetadot = 0;                               // temporary fudge FIX FIX
                for (int vari = 0; vari < 3; vari++)                       // var = (dx'/dc, dy'/dc, dz'/dc)
                    for (int varj = 0; varj < 3; varj++)
                        Mout[vari][varj][i] = - dTdphi[vari][varj]*phidot - dTdtheta[vari][varj]*thetadot
                                              - dTd1[vari][varj] - dTda[vari][varj]*a
                                              - dTdxc[vari][varj]*(gen_x[i] - c) - dTdz[vari][varj]*gen_z[i];
                Moutdot00[i] = -2.0*a*Math.sin(phi)*Math.cos(phi)*phidot;
                Moutdot01[i] = -(phi2dot*Math.cos(theta) - phidot*Math.sin(theta)*thetadot
                             + thetadot*(Math.sin(theta) - Math.cos(phi)*Math.cos(theta) - a*Math.sin(phi)*Math.cos(phi)*Math.sin(theta))
                             + phidot*(Math.sin(phi)*Math.sin(theta) + a*Math.cos(phi)*Math.cos(phi)*Math.cos(theta) - a*Math.sin(phi)*Math.sin(phi)*Math.cos(theta)));
                //v[i + Period] = -Math.sin(theta)*gen_z[i];              // original code for d/dc
                //v[i + Period] = Math.sin(theta);                        // new code for d/db
                v[i] = Math.sin(phi)*gen_y[i];                          // new code for d/da
                v[i + Period] = Math.cos(theta)*Math.cos(phi)*gen_y[i]; // new code for d/da
                v[i + 2*Period] = -Math.cos(theta)*gen_z[i];            // not used anymore (part of original 3D system)
                //System.out.println(i + ", " + Mout[0][0][i] + ", " + Mout[0][1][i] + ", " + Mout[1][0][i] + ", " + Mout[1][1][i]);
                //gen_M2D(gen_x[i], gen_y[i], gen_z[i]);
            }
            //v[2*Period] += -dTaudc*(vN2 - 8.0*vN1)/delt/12.0;   // apply time drift in b.c.
            //v[2*Period + 1] += -dTaudc*vN1/delt/12.0;
            //v[3*Period - 2] += -dTaudc*v0/delt/12.0;
            //v[3*Period - 1] += dTaudc*(8.0*v0 - v1)/delt/12.0;
            System.out.println("\nrange, " + a + ", " + b + ", " + c + ", " + Period + ", " + delt + ", " + xmin_temp + ", " + xmax_temp + ", " + ymin_temp + ", " + ymax_temp + ", " + zmin_temp + ", " + zmax_temp + ", " + thmin_temp + ", " + thmax_temp + ", " + (av_x/Period + a - c));

            if (false)                                          // solve Mx = b in Java
            {
                M = new double[3*Period][3*Period];
                for (int i = 0; i < M.length; i++)
                    for (int j = 0; j < M.length; j++)
                        M[i][j] = 0;
                for (int i = 0; i < Period; i++)                // i = time index
                {
                    for (int vari = 0; vari < 3; vari++)        // var = (dx'/dc, dy'/dc, dz'/dc)
                        for (int varj = 0; varj < 3; varj++)
                            M[i + vari*Period][i + varj*Period] = Mout[vari][varj][i];
                    for (int vari = 0; vari < 3; vari++)        // d/dt of (dx'/dc, dy'/dc, dz'/dc)
                    {
                        M[i + vari*Period][(i + 1) % Period + vari*Period] += 8.0/delt/12.0;
                        M[i + vari*Period][(i - 1 + Period) % Period + vari*Period] += -8.0/delt/12.0;
                        M[i + vari*Period][(i + 2) % Period + vari*Period] += -1.0/delt/12.0;
                        M[i + vari*Period][(i - 2 + Period) % Period + vari*Period] += 1.0/delt/12.0;
                    }
                }
                double[] soln = gaussj(M, v);
                System.out.println("\nRossler soln: Ax = b");
                System.out.println("c      = " + c);
                System.out.println("Period = " + Period);
                System.out.println("delt   = " + delt);
                System.out.println("i, dx'/dc, dy'/dc, dz'/dc");
                for (int i = 0; i < Period; i++)                // i = time index
                    System.out.println(i + ", " + soln[i] + ", " + soln[i + Period] + ", " + soln[i + 2*Period]);
            }
            //System.out.println("M = ");
            //for (int i = 0; i < 3*Period; i++)
            //{
            //    System.out.print(i);
            //    for (int j = 0; j < 3*Period; j++)
            //        System.out.print(", " + M[i][j]);
            //    System.out.println();
            //}
            //System.out.println("v = ");
            //for (int i = 0; i < 3*Period; i++)
            //    System.out.println(i + ", " + v[i]);

            //if (true)                                         // test reduced 2-D version
            //{
            //    double[][] M2 = new double[2*Period][2*Period];
            //    double[] v2 = new double[2*Period];
            //    for (int i = 0; i < M2.length; i++)
            //        for (int j = 0; j < M2.length; j++)
            //            M2[i][j] = M[i][j];
            //    for (int i = 0; i < v2.length; i++)
            //        v2[i] = v[i];
            //    System.out.println("\nReduced 2-D solution of Rossler Mx = b");
            //    double[] soln2 = gaussj(M2, v2);
            //    for (int i = 0; i < Period; i++)                        // i = time index
            //        System.out.println(i + ", " + soln2[i] + ", " + soln2[i + Period]);
            //}

            if (false)                                           // solve first-order d.e. in Python
            {
                System.out.println("\nfinal (x y z) =," + a + ", " + b + ", " + c + ", " + gen_x[Period - 1] + ", " + gen_y[Period - 1] + ", " + gen_z[Period - 1]);
                System.out.println("\n# Rossler first-order d.e. elements of M (Python)");
                System.out.println("a_in   = " + a);
                System.out.println("b_in   = " + b);
                System.out.println("c_in   = " + c);
                System.out.println("delt   = " + delt);
                for (int vari = 0; vari < 2; vari++)            // var = (dx'/dc, dy'/dc, dz'/dc)
                    for (int varj = 0; varj < 2; varj++)
                    {
                        System.out.print("M" + vari + varj + " = np.array([" + Mout[vari][varj][0]);
                        for (int i = 1; i < Period; i++)            // i = time index
                            System.out.print("," + Mout[vari][varj][i]);
                        System.out.println("])");
                    }
                System.out.print("b0 = np.array([" + v[0]);
                for (int i = 1; i < Period; i++)
                    System.out.print(", " + v[i]);
                System.out.println("])");
                System.out.print("b1 = np.array([" + v[Period]);
                for (int i = 1; i < Period; i++)
                    System.out.print(", " + v[i + Period]);
                System.out.println("])");

                //System.out.println("\nindex,M00,M01,M10,M11,b");
                //for (int i = 0; i < Period; i++)
                //    System.out.println(i + ", " + Mout[0][0][i] + ", " + Mout[0][1][i] + ", " + Mout[1][0][i] + ", " + Mout[1][1][i] + ", " + v[i + Period]);
            }

            if (true)                                           // generate transformed Mout vs time
            {
                System.out.println("Rossler - Limit cycle and 2-D (x', y') (linear) coupling matrix");
                System.out.println("a_b_c , " + a + ", " + b + ", " + c);
                System.out.println("Period_delt , 0, " + Period + ", " + delt);
                System.out.println("i,x,y,z,phi,theta,phidot,thetadot,M[0][0],M[0][1],M[0][2],M[1][0],M[1][1],M[1][2],M[2][0],M[2][1],M[2][2]");
                for (int i = 0; i < Period; i++)            // i = time index
                    System.out.println(i + ", " + gen_x[i] + ", " + gen_y[i] + ", " + gen_z[i] + ",0.0,0.0,0.0,0.0, " + Mout[0][0][i] + ", " + Mout[0][1][i] + ", " + Mout[0][2][i] + ", " + Mout[1][0][i] + ", " + Mout[1][1][i] + ", " + Mout[1][2][i] + ", " + Mout[2][0][i] + ", " + Mout[2][1][i] + ", " + Mout[2][2][i]);
            }

            if (false)                                           // solve second-order d.e. in Python
            {
                System.out.println("\n# Rossler second-order d.e. elements of M");
                System.out.println("c      = " + c);
                System.out.println("delt   = " + delt);
                System.out.print("x2dot = np.array([");
                for (int i = 0; i < Period; i++)                // i = time index
                {
                    if (i > 0) System.out.print(",");
                    System.out.print(-Mout[0][1][i]);
                }
                System.out.println("])");
                System.out.print("xdot  = np.array([");
                for (int i = 0; i < Period; i++)                // i = time index
                {
                    if (i > 0) System.out.print(",");
                    System.out.print(Moutdot01[i] + Mout[0][0][i]*Mout[0][1][i] + Mout[0][1][i]*Mout[1][1][i]);
                }
                System.out.println("])");
                System.out.print("x     = np.array([");
                for (int i = 0; i < Period; i++)                // i = time index
                {
                    if (i > 0) System.out.print(",");
                    System.out.print(Moutdot00[i]*Mout[0][1][i] - Mout[0][0][i]*Moutdot01[i]
                                   + Mout[0][1][i]*(Mout[0][1][i]*Mout[1][0][i] - Mout[0][0][i]*Mout[1][1][i]));
                }
                System.out.println("])");
                System.out.print("b     = np.array([");
                for (int i = 0; i < Period; i++)                // i = time index
                {
                    if (i > 0) System.out.print(",");
                    System.out.print(-Mout[0][1][i]*Mout[0][1][i]*v[i + Period]);
                }
                System.out.println("])");
            }
        }
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
            if (new File(System.getProperty("user.home"), "RosslerPrefs.ini").exists())
            {
                pgmProp.load(new FileInputStream(new File(System.getProperty("user.home"), "RosslerPrefs.ini")));
                type = pgmProp.getProperty("type", "phase");
                a = Double.parseDouble(pgmProp.getProperty("a", "0.2"));
                b = Double.parseDouble(pgmProp.getProperty("b", "0.2"));
                c = Double.parseDouble(pgmProp.getProperty("c", "2.5"));
                astart = Double.parseDouble(pgmProp.getProperty("astart", "2"));
                aend = Double.parseDouble(pgmProp.getProperty("aend", "NaN"));
                bstart = Double.parseDouble(pgmProp.getProperty("bstart", "2"));
                bend = Double.parseDouble(pgmProp.getProperty("bend", "NaN"));
                cstart = Double.parseDouble(pgmProp.getProperty("cstart", "2"));
                cend = Double.parseDouble(pgmProp.getProperty("cend", "NaN"));
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
                dx0dc = Double.parseDouble(pgmProp.getProperty("dx0dc", "0"));
                dy0dc = Double.parseDouble(pgmProp.getProperty("dy0dc", "0"));
                dz0dc = Double.parseDouble(pgmProp.getProperty("dz0dc", "0"));
                project_phi = Double.parseDouble(pgmProp.getProperty("project_phi", "0"));
                project_theta = Double.parseDouble(pgmProp.getProperty("project_theta", "0"));
                project_psi = Double.parseDouble(pgmProp.getProperty("project_psi", "0"));
            }
            else                                    // factory default
            {
                type = "phase";
                a = 0.2;
                b = 0.2;
                c = 2.5;
                astart = 2;
                aend = Double.NaN;
                bstart = 2;
                bend = Double.NaN;
                cstart = 2;
                cend = Double.NaN;
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
                dx0dc = 0;
                dy0dc = 0;
                dz0dc = 0;
                project_phi = 0;
                project_theta = 0;
                project_psi = 0;
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
        pgmProp.setProperty("astart", "" + astart);
        pgmProp.setProperty("aend", "" + aend);
        pgmProp.setProperty("bstart", "" + bstart);
        pgmProp.setProperty("bend", "" + bend);
        pgmProp.setProperty("cstart", "" + cstart);
        pgmProp.setProperty("cend", "" + cend);
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
        pgmProp.setProperty("dx0dc", "" + dx0dc);
        pgmProp.setProperty("dy0dc", "" + dy0dc);
        pgmProp.setProperty("dz0dc", "" + dz0dc);
        pgmProp.setProperty("project_phi", "" + project_phi);
        pgmProp.setProperty("project_theta", "" + project_theta);
        pgmProp.setProperty("project_psi", "" + project_psi);
        try
            {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "RosslerPrefs.ini"), "Rossler System v" + VERSION_NO + " Prefs");}
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

    protected static double[][] multmm(double[][] m1, double[][] m2)
    {
        // matrix multiplication (for square matrices only)
        double[][] retVal = new double[m1.length][m1[0].length];

        for (int i = 0; i < m1.length; i++)
            for (int j = 0; j < m1.length; j++)
            {
                retVal[i][j] = 0;
                for (int k = 0; k < m1.length; k++)
                    retVal[i][j] += m1[i][k]*m2[k][j];
            }
        return retVal;
    }
}
