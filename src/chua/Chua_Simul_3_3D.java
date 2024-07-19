
package chua;

/*
 * this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Simul_3_3D.java
 * Feb 20, 2024, June 10, 2024
 * simulate a Poincare map of degree 3 in the x'-y'-z' space
 * based on model from Book "Averaging", p. 39
 * see also             java  : fit_cubic_response_full_3D()
 *                      Python: Chua_cylindrical_response.py
*/

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;

public class Chua_Simul_3_3D extends JDialog
{
    private static final BufferedImage image = new BufferedImage(480, 480, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static JButton btnCalc = new JButton("Calc");
    private static JButton btnReset = new JButton("Reset");
    private static JButton btnClear = new JButton("Clear");
    private static JCheckBox printChk = new JCheckBox("  print  ");
    private static JButton btnZerothOrder = new JButton("Generate Zeroth Order");

    private static JTextField[] txtCarr = {new JTextField("0.96775"),       // a
                                           new JTextField("0.25748"),       // b
                                           new JTextField("0.1"),           // c
                                           new JTextField("0.2"),           // d
                                           new JTextField("-0.05"),         // del r
                                           new JTextField("-0.8"),          // beta r
                                           new JTextField(""),              // del z
                                           new JTextField("-0.23"),         // beta z
                                           new JTextField(""),              // C r2z
                                           new JTextField("")};             // C zzz
    private static JTextField[] txtstart = {new JTextField("-0.001"), new JTextField("0"), new JTextField("0")};
    private static JTextField[] txtrange = {new JTextField("0.8"), new JTextField("0.8"), new JTextField("0.3")};
    private static double a, b, c, d, del_r, beta_r, del_z, beta_z, Cr2z, Czzz, alpha, angle;
    private static double[] pos = new double[3];                // position x, y, z
    private static int ix = 0, iy = 1;                          // axes to plot
    private static int iT = 0;

    public Chua_Simul_3_3D()
    {
        setTitle("Chua System - Simulate x'-y'-z' Scatter (cubic model)");
        setIconImage(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        setSize(900, 532);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());

        lblImage.setBorder(BorderFactory.createEtchedBorder());

        final JLabel[] lblarr = {new JLabel("a"),
                                 new JLabel("b"),
                                 new JLabel("c"),
                                 new JLabel("d"),
                                 new JLabel("del r"),
                                 new JLabel("beta r"),
                                 new JLabel("del z"),
                                 new JLabel("beta z"),
                                 new JLabel("C r2z"),
                                 new JLabel("C zzz")};
        JPanel spacerPanel = new JPanel();
        spacerPanel.setPreferredSize(new Dimension(300, 6));
        spacerPanel.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel.setOpaque(false);
        JPanel[] dataPanel = new JPanel[lblarr.length];
        for (int i = 0; i < dataPanel.length; i++)
        {
            dataPanel[i] = new JPanel();
            dataPanel[i].setPreferredSize(new Dimension(350, 22));
            dataPanel[i].setOpaque(false);
            lblarr[i].setPreferredSize(new Dimension(40, 20));
            dataPanel[i].add(lblarr[i]);
            txtCarr[i].setPreferredSize(new Dimension(85, 20));
            dataPanel[i].add(txtCarr[i]);
        }

        JPanel startPanel = new JPanel();
        startPanel.setOpaque(false);
        startPanel.setPreferredSize(new Dimension(350, 24));
        JLabel lblstart = new JLabel("start");
        lblstart.setPreferredSize(new Dimension(40, 20));
        startPanel.add(lblstart);
        for (int i = 0; i < txtstart.length; i++)
        {
            txtstart[i].setPreferredSize(new Dimension(85, 20));
            startPanel.add(txtstart[i]);
            pos[i] = Double.parseDouble(txtstart[i].getText());
        }
        JPanel rangePanel = new JPanel();
        rangePanel.setOpaque(false);
        rangePanel.setPreferredSize(new Dimension(350, 24));
        JLabel lblrange = new JLabel("range");
        lblrange.setPreferredSize(new Dimension(40, 20));
        rangePanel.add(lblrange);
        for (int i = 0; i < txtrange.length; i++)
        {
            txtrange[i].setPreferredSize(new Dimension(85, 20));
            rangePanel.add(txtrange[i]);
        }
        JPanel printPanel = new JPanel();
        printPanel.setOpaque(false);
        printPanel.setPreferredSize(new Dimension(300, 24));
        printChk.setOpaque(false);
        printPanel.add(printChk);

        final JPanel parmsPanel = new JPanel();
        parmsPanel.setBackground(new Color(200, 221, 242));
        for (int i = 0; i < dataPanel.length; i++)
            parmsPanel.add(dataPanel[i]);
        parmsPanel.add(spacerPanel);
        parmsPanel.add(btnCalc);
        parmsPanel.add(startPanel);
        parmsPanel.add(rangePanel);
        parmsPanel.add(btnReset);
        parmsPanel.add(btnClear);
        parmsPanel.add(printPanel);
        parmsPanel.add(btnZerothOrder);
        parmsPanel.setMaximumSize(new Dimension(350, 3000));
        parmsPanel.setPreferredSize(new Dimension(350, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        final JPanel scatterPanel = new JPanel();
        scatterPanel.setOpaque(false);
        scatterPanel.add(lblImage);

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().setBackground(new Color(200, 221, 242));
        getContentPane().add(parmsPanel);
        getContentPane().add(scatterPanel);

        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        btnCalc.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                btnCalc.setEnabled(false);
                for (int i = 0; i < txtCarr.length; i++)
                    if (txtCarr[i].getText().isEmpty())
                        txtCarr[i].setText("0");
                for (int j = 0; j < txtrange.length; j++)
                    if (txtrange[j].getText().isEmpty())
                        txtrange[j].setText("0");
                setTitle("Chua System - Simulate x'-y'-z' Scatter (cubic model) (" + ix + ", " + iy + ")");
                refresh_graph();
                btnCalc.setEnabled(true);
            }
        });

        btnReset.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                double ret = 0.0;                     // TEMPORARY CODE !!!!!!!!!!!!!!
                System.out.println("...................................");
                System.out.println("init ," + ret);
                for (int i = 0; i < 200; i++)
                    ret = calc_r_cubic(ret, "theta");
                calc_max_min_det(ret);

                for (int i = 0; i < txtstart.length; i++)
                    pos[i] = Double.parseDouble(txtstart[i].getText());
            }
        });

        btnClear.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                DC.clearRect(0, 0, image.getWidth(), image.getHeight());
                lblImage.repaint();
            }
        });

        btnZerothOrder.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                dump_zeroth_order();
            }
        });
    }

    private void refresh_graph()
    {
        PrintWriter fout = null;
        double[] range = new double[] {Double.parseDouble(txtrange[0].getText()), Double.parseDouble(txtrange[1].getText()), Double.parseDouble(txtrange[2].getText())};
        double tempx, tempy, r2;
        double oldx = 0, oldy = 0;
        int N = 1000;

        a = Double.parseDouble(txtCarr[0].getText());                   // first-order x-x
        b = Double.parseDouble(txtCarr[1].getText());                   // first-order y-x
        c = Double.parseDouble(txtCarr[2].getText());                   // first-order z-x
        d = Double.parseDouble(txtCarr[3].getText());                   // first-order z-y
        del_r =  Double.parseDouble(txtCarr[4].getText());              // third-order radial diagonal
        beta_r = Double.parseDouble(txtCarr[5].getText());              // third-order radial off-diagonal
        del_z  = Double.parseDouble(txtCarr[6].getText());              // third-order z diagonal
        beta_z = Double.parseDouble(txtCarr[7].getText());              // third-order z off-diagonal
        Cr2z   = Double.parseDouble(txtCarr[8].getText());              // z-map : r*r*z
        Czzz   = Double.parseDouble(txtCarr[9].getText());              // z-map : z*z*z

        alpha = Math.sqrt(a*a + b*b) - 1;                               // first-order alpha
        angle = Math.atan2(b, a);                                       // first-order angular response
        double rcubic = calc_r_cubic(0, "r");
        //System.out.println("r = " + calc_r_cubic("r"));
        //System.out.println("theta = " + angle*180/Math.PI + ", " + calc_r_cubic("theta"));
        //System.out.println("r1 = " + calc_r_cubic("r"));

        //rcubic = range[ix];                    // override cubic theory
        //double major  = Math.atan2(d, c);
        //System.out.println("major axis = " + rcubic + ", " + major*180/Math.PI);
        //System.out.println(alpha + ", " + theta*180/Math.PI);
        //System.out.println("start " + iT + ", " + pos[0] + ", " + pos[1] + ", " + pos[2]);

        String hdr = "Simul_3D_" + a + "_" + b + "_" + c + "_" + d + "_" + del_r + "_" + beta_r + "_" + del_z + "_" + beta_z + "_" + Cr2z + "_" + Czzz;

        if (printChk.isSelected())
            try
            {
                String fname = "C:\\Windows\\Temp\\" + hdr + ".csv";
                boolean fexists = new File(fname).exists();
                FileWriter fw = new FileWriter(fname, true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    //fout.println(hdr + ", " + ix + ", " + iy);
                    fout.println("     , a, b, c, d, del_r, beta_r, del_z, beta_z, Cr2z, Czzz");
                    fout.println("x_y_scatter, " + a + ", " + b + ", " + c + ", " + d + ", " + del_r + ", " + beta_r + ", " + del_z + ", " + beta_z + ", " + Cr2z + ", " + Czzz);
                    fout.println("scatter hdr, 0, " + Math.sqrt(a*a + b*b) + ", " + Math.atan2(b, a)*180/Math.PI);
                    fout.println("init x0 y0 , 0, 0, NaN, NaN");
                    fout.println("iT, x', y', z'");
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}
        image.setRGB(image.getWidth()/2, image.getHeight()/2, Color.BLACK.getRGB());

        for (int i = 0; i < N; i++)
        {
            oldx = pos[0];
            oldy = pos[1];
            r2 = pos[0]*pos[0] + pos[1]*pos[1];
            tempx = (1 + alpha)*pos[0] + r2*(del_r*pos[0] - beta_r*pos[1]) + pos[2]*pos[2]*(del_z*pos[0] - beta_z*pos[1]);
            tempy = (1 + alpha)*pos[1] + r2*(beta_r*pos[0] + del_r*pos[1]) + pos[2]*pos[2]*(beta_z*pos[0] + del_z*pos[1]);
            //tempy += 0.04*pos[2]*pos[2]*pos[0];
            //tempy += 0.04*pos[0]*pos[0]*pos[0];
            pos[2] = c*pos[0] + d*pos[1] + pos[2] + Cr2z*r2*pos[2] + Czzz*pos[2]*pos[2]*pos[2];
            pos[0] = Math.cos(angle)*tempx - Math.sin(angle)*tempy;
            pos[1] = Math.sin(angle)*tempx + Math.cos(angle)*tempy;

            if (image.getWidth()/2 + pos[ix]/range[ix]*image.getWidth()/2 < 0
            ||  image.getWidth()/2 + pos[ix]/range[ix]*image.getWidth()/2 > image.getWidth()
            ||  image.getHeight()/2 - pos[iy]/range[iy]*image.getHeight()/2 < 0
            ||  image.getHeight()/2 - pos[iy]/range[iy]*image.getHeight()/2 > image.getHeight())
                System.out.println("-------, " + iT + ", " + pos[0] + ", " + pos[1] + ", " + pos[2]);
            else
                image.setRGB(image.getWidth()/2 + (int) (pos[ix]/range[ix]*image.getWidth()/2), image.getHeight()/2  - (int) (pos[iy]/range[iy]*image.getHeight()/2), Color.BLACK.getRGB());
            if (fout != null)
                fout.println(iT + ", " + pos[0] + ", " + pos[1] + ", " + pos[2]);
            iT++;
        }
        //System.out.println(iT + ", " + Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1]) + ", " + oldx + ", " + oldy + ", " + pos[0] + ", " + pos[1]);
        System.out.println(iT + ", " + Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1]) + ", " + (Math.atan2(pos[1], pos[0])*180/Math.PI - Math.atan2(oldy, oldx)*180/Math.PI));
        DC.setColor(Color.blue);
        DC.drawLine(0, image.getHeight()/2, image.getWidth(), image.getHeight()/2);
        DC.drawLine(image.getWidth()/2, 0, image.getWidth()/2, image.getHeight());
        //DC.drawLine(image.getWidth()/2 - (int) (rcubic/range[ix]*image.getWidth()*Math.cos(major)/2), image.getWidth()/2 + (int) (rcubic/range[ix]*image.getWidth()*Math.sin(major)/2),
        //            image.getWidth()/2 + (int) (rcubic/range[ix]*image.getWidth()*Math.cos(major)/2), image.getWidth()/2 - (int) (rcubic/range[ix]*image.getWidth()*Math.sin(major))/2);    // draw major axis
        DC.setColor(Color.red);
        DC.drawOval(image.getWidth()/2 - (int) (rcubic/range[ix]*image.getWidth()/2), image.getHeight()/2 - (int) (rcubic/range[iy]*image.getWidth()/2),
                   (int) (rcubic/range[ix]*image.getWidth()), (int) (rcubic/range[iy]*image.getHeight()));
        lblImage.repaint();
        if (fout != null)
            fout.close();
    }

    private void dump_zeroth_order()
    {
        double R = 0.1693186731537729; // 0.75;
        double theta = 15.22750677088466; // 21.5;

        double sumcos, sumsin;
        double z, M21, M22;
        System.out.println("zeroth order, " + R + ", " + theta);
        theta = theta*Math.PI/180;
        double R_ampl = R*Math.sqrt((c*c + d*d)/2/(1 - Math.cos(theta)));
        try
        {
            String fname = "C:\\Windows\\Temp\\dump_zeroth.csv";
            FileWriter fw = new FileWriter(fname, true);
            PrintWriter fdump = new PrintWriter(fw);
            fdump.println("dump_zeroth , " + a + ", " + b + ", " + c + ", " + d + ", " + del_r + ", " + beta_r + ", " + del_z + ", " + beta_z);
            fdump.println("alpha_angle , " + alpha + ", " + angle*180/Math.PI + ", " + calc_r_cubic(0, "r"));
            fdump.println("R_theta_|R| , " + R + ", " + theta*180/Math.PI + ", " + R_ampl);
            fdump.println("N, x, y, z");
            for (int i = 0; i <= 200; i++)
            {
                sumcos = (Math.cos(i*theta) - Math.cos((i + 1)*theta))/2/(1 - Math.cos(theta));
                //sumcos = (1 - Math.cos(theta) + Math.cos(i*theta) - Math.cos((i + 1)*theta))/2/(1 - Math.cos(theta));
                sumsin = (Math.sin(i*theta) - Math.sin((i + 1)*theta))/2/(1 - Math.cos(theta));
                //sumsin = (Math.sin(theta) + Math.sin(i*theta) - Math.sin((i + 1)*theta))/2/(1 - Math.cos(theta));
                z = R*c*sumcos + R*d*sumsin;
                M21 = R*R*beta_r + z*z*beta_z;
                M22 = 1 + alpha + R*R*del_r + z*z*del_z;
                fdump.println(i + ", " + R*Math.cos(i*theta) + ", " + R*Math.sin(i*theta) + ", " + z + ", " + Math.atan2(M21, M22)*180/Math.PI);
            }
            fdump.close();
        }
        catch (java.io.IOException e)
            {System.out.println("dump_zeroth_order save error = " + e);}
    }

    private double calc_r_cubic(double del_theta, String type)
    {
        // original calc of cubic R^2 (see "Simplified Harmonic Balance" June 22)
        double rho2 = Math.sqrt(1 + beta_r*beta_r/del_r/del_r - beta_r*beta_r/del_r/del_r*(1 + alpha)*(1 + alpha));
        rho2 = (-(1 + alpha) + rho2)/(1 + beta_r*beta_r/del_r/del_r)/del_r;
        if (type.equals("r")) return Math.sqrt(rho2);   // return zeroth order 'r'

        // bootstrap calc of R^2 (see "Simplified Harmonic Balance" June 22)
        double z2_coeff = (c*c + d*d)/4/(1 - Math.cos(angle + del_theta));
        double del_new = del_r + z2_coeff*del_z;
        double beta_new = beta_r + z2_coeff*beta_z;
        rho2 = Math.sqrt(1 + beta_new*beta_new/del_new/del_new - beta_new*beta_new/del_new/del_new*(1 + alpha)*(1 + alpha));
        rho2 = (-(1 + alpha) + rho2)/(1 + beta_new*beta_new/del_new/del_new)/del_new;
        del_theta = Math.asin(beta_new*rho2);          // return updated del theta
        System.out.println(Math.sqrt(rho2) + ", " + (angle + del_theta)*180/Math.PI + ", " + del_theta);
        return del_theta;
    }

    private void calc_max_min_det(double del_theta)
    {
        System.out.println();
        double z2_coeff = (c*c + d*d)/4/(1 - Math.cos(angle + del_theta));
        //z2_coeff = z2_coeff;
        double del_new = del_r + z2_coeff*del_z;
        double beta_new = beta_r + z2_coeff*beta_z;
        double rho2 = Math.sqrt(1 + beta_new*beta_new/del_new/del_new - beta_new*beta_new/del_new/del_new*(1 + alpha)*(1 + alpha));
        rho2 = (-(1 + alpha) + rho2)/(1 + beta_new*beta_new/del_new/del_new)/del_new;
        System.out.println("average r," + Math.sqrt(rho2));
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                Chua_Simul_3_3D dlg = new Chua_Simul_3_3D();
            }
        });
    }
}
