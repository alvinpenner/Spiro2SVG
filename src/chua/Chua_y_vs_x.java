
package chua;

// plot y versus x
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import javax.swing.*;
import java.util.GregorianCalendar;
import java.io.*;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;

public class Chua_y_vs_x extends JDialog
{
    private boolean first = true;
    private final static double DEFAULT_DELT = -0.001; //-6.990205973471885E-5; // 0.0005; // -0.0008; // 0.0005; // -6.990205974363416E-5;    // -0.0005;
    private final static int N = 30000; // 480000; // 7000; // 480000;                // total # of iterations 480000
    private static double[] pt6_old;
    private static double eig, angle;                   // first-order response
    protected Path2D.Double path1 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);     // start band
    protected Path2D.Double path2 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);     // middle path
    protected Path2D.Double path3 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);     // end band
    protected Path2D.Double stataxis = new Path2D.Double(Path2D.WIND_NON_ZERO, 10);
    protected Line2D.Double xaxis = new Line2D.Double(0, 0, 0, 0);
    protected Line2D.Double yaxis = new Line2D.Double(0, 0, 0, 0);
    private static JCheckBox printChk = new JCheckBox("  print  ");
    private static JRadioButton phaseRadio = new JRadioButton("x-y phase");
    //private static JRadioButton ddcRadio = new JRadioButton("d/dc of phase");
    private static JRadioButton ddu3Radio = new JRadioButton("3-D d/du");
    //private static JRadioButton ddu2Radio = new JRadioButton("2-D d/du");
    private static JLabel lblxrange;
    private static JLabel lblyrange;
    private static JLabel lblzrange;
    private static JPanel phasePanel = new Plot_Phase_Panel();
    protected static String first_order_hdr;            // # iterations, eig, angle
    private static int iT = 0;                          // # iterations

    double[] zlist = new double[4];                     // previous z values
    double Told = 0;
    int Nfork = 1; // 23;                               // number of bifurcated branches
    double[] Tfork = new double[Nfork];                 // period per branch
    int Tindex = 0;                                     // number of peaks
    double delt = DEFAULT_DELT;

    public Chua_y_vs_x(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        Main.type = "phase";
        if (phaseRadio.isSelected())                            // normal phase space
            setTitle(" Chua System - phase y' vs. x'");
        else if (ddu3Radio.isSelected())
            setTitle(" Chua System - dy/du vs. dx/du (3-D)");  // Lyapunov response
        else
            setTitle(" Chua System - phase y' vs. x'");
        setIconImage(img);
        setSize(610 + 135, 493 + 135 + 65);
        setLocationByPlatform(true);

        final JLabel[] lbl = {new JLabel("alpha"),
                              new JLabel("beta"),
                              new JLabel("gamma"),
                              new JLabel("c"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0"),
                              new JLabel("dx0du"),
                              new JLabel("dy0du"),
                              new JLabel("dz0du"),
                              new JLabel("Period")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.alpha)),
                                  new JTextField(Double.toString(Main.beta)),
                                  new JTextField(Double.toString(Main.gamma)),
                                  new JTextField(Double.toString(Main.c)),
                                  new JTextField(Double.toString(Main.x0)),
                                  new JTextField(Double.toString(Main.y0)),
                                  new JTextField(Double.toString(Main.z0)),
                                  new JTextField(Double.toString(Main.dx0du)),
                                  new JTextField(Double.toString(Main.dy0du)),
                                  new JTextField(Double.toString(Main.dz0du)),
                                  new JTextField()};
        JPanel[] spacerPanel = new JPanel[4];
        JPanel[] dataPanel = new JPanel[lbl.length];
        JButton btnRun = new JButton("Run");
        final JCheckBox chkFast = new JCheckBox(" fast ");
        JButton btnFitLinear = new JButton("Fit Linear");
        JButton btnPerturbCoeff = new JButton("Perturb Coeffs.");
        JButton btnFitCubic = new JButton("Fit Cubic");
        chkFast.setOpaque(false);
        chkFast.setEnabled(false);

        for (int i = 0; i < spacerPanel.length; i++)
        {
            spacerPanel[i] = new JPanel();
            spacerPanel[i].setPreferredSize(new Dimension(110, 6));
            spacerPanel[i].setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
            spacerPanel[i].setOpaque(false);
        }
        for (int i = 0; i < dataPanel.length; i++)
        {
            dataPanel[i] = new JPanel();
            dataPanel[i].setPreferredSize(new Dimension(130, 22));
            dataPanel[i].setOpaque(false);
            lbl[i].setPreferredSize(new Dimension(45, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

        printChk.setOpaque(false);
        JPanel printPanel = new JPanel();
        printPanel.setPreferredSize(new Dimension(130, 24));
        printPanel.setOpaque(false);
        printPanel.add(printChk);
        printPanel.add(chkFast);

        phaseRadio.setOpaque(false);
        ddu3Radio.setOpaque(false);
        ButtonGroup group = new ButtonGroup();
        group.add(phaseRadio);
        group.add(ddu3Radio);
        phaseRadio.setSelected(true);
        JPanel radioPanel = new JPanel();
        radioPanel.setOpaque(false);
        radioPanel.setPreferredSize(new Dimension(120, 2*27));
        radioPanel.setBorder(BorderFactory.createEtchedBorder());
        radioPanel.setLayout(new BoxLayout(radioPanel, BoxLayout.Y_AXIS));
        radioPanel.add(phaseRadio);
        radioPanel.add(ddu3Radio);

        JPanel xrangePanel = new JPanel();
        xrangePanel.setPreferredSize(new Dimension(130, 20));
        xrangePanel.setOpaque(false);
        lblxrange = new JLabel("x = ");
        lblxrange.setPreferredSize(new Dimension(110, 18));
        xrangePanel.add(lblxrange);

        JPanel yrangePanel = new JPanel();
        yrangePanel.setPreferredSize(new Dimension(130, 20));
        yrangePanel.setOpaque(false);
        lblyrange = new JLabel("y = ");
        lblyrange.setPreferredSize(new Dimension(110, 18));
        yrangePanel.add(lblyrange);

        JPanel zrangePanel = new JPanel();
        zrangePanel.setPreferredSize(new Dimension(130, 20));
        zrangePanel.setOpaque(false);
        lblzrange = new JLabel("z = ");
        lblzrange.setPreferredSize(new Dimension(110, 18));
        zrangePanel.add(lblzrange);

        parmsPanel.removeAll();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(dataPanel[0]);
        parmsPanel.add(dataPanel[1]);
        parmsPanel.add(dataPanel[2]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(dataPanel[7]);
        parmsPanel.add(dataPanel[8]);
        parmsPanel.add(dataPanel[9]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(printPanel);
        parmsPanel.add(radioPanel);
        parmsPanel.add(xrangePanel);
        parmsPanel.add(yrangePanel);
        parmsPanel.add(zrangePanel);
        parmsPanel.add(spacerPanel[2]);
        parmsPanel.add(dataPanel[10]);
        parmsPanel.add(btnFitLinear);
        parmsPanel.add(btnPerturbCoeff);
        parmsPanel.add(btnFitCubic);

        parmsPanel.setMaximumSize(new Dimension(140, 3000));
        parmsPanel.setPreferredSize(new Dimension(140, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().add(parmsPanel);
        getContentPane().add(phasePanel);
        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        pt6_old = new double[] {Main.x0, Main.y0, Main.z0, Main.dx0du, Main.dy0du, Main.dz0du};

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                boolean changed = false;
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                if (phaseRadio.isSelected())                                // normal phase space
                    setTitle(" Chua System - phase y' vs. x' (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + Main.a + ", " + String.format(" %.6f", delt) + ")");
                else if (ddu3Radio.isSelected())
                    setTitle(" Chua System - dy/du vs. dx/du (3-D)");    // Lyapunov response
                else
                    setTitle(" Chua System - dy/du vs. dx/du (2-D)");    // Lyapunov response
                //if (Main.alpha != Double.parseDouble(txt[0].getText())) changed = true;
                Main.alpha = Double.parseDouble(txt[0].getText());
                //if (Main.beta != Double.parseDouble(txt[1].getText())) changed = true;
                Main.beta = Double.parseDouble(txt[1].getText());
                //if (Main.gamma != Double.parseDouble(txt[2].getText())) changed = true;
                Main.gamma = Double.parseDouble(txt[2].getText());
                //if (Main.c != Double.parseDouble(txt[3].getText())) changed = true;
                Main.c  = Double.parseDouble(txt[3].getText());
                if (Main.x0 != Double.parseDouble(txt[4].getText())) changed = true;
                Main.x0 = Double.parseDouble(txt[4].getText());
                if (Main.y0 != Double.parseDouble(txt[5].getText())) changed = true;
                Main.y0 = Double.parseDouble(txt[5].getText());
                if (Main.z0 != Double.parseDouble(txt[6].getText())) changed = true;
                Main.z0 = Double.parseDouble(txt[6].getText());
                //if (Main.dx0du != Double.parseDouble(txt[7].getText())) changed = true;
                Main.dx0du = Double.parseDouble(txt[7].getText());
                //if (Main.dy0du != Double.parseDouble(txt[8].getText())) changed = true;
                Main.dy0du = Double.parseDouble(txt[8].getText());
                //if (Main.dz0du != Double.parseDouble(txt[9].getText())) changed = true;
                Main.dz0du = Double.parseDouble(txt[9].getText());
                phase_space(changed, (int) Double.parseDouble(txt[10].getText()), chkFast.isSelected());
                //System.out.println("btnRun " + phasePanel.getSize());
                //System.out.println("btnRun " + parmsPanel.getSize());
            }
        });

        btnFitLinear.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                first_order_hdr = fit_linear_response();
                setTitle(" Chua System - phase y' vs. x' (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + String.format(" %.6f", Main.final_delt) + ")");
            }
        });

        btnPerturbCoeff.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                //Chua_Perturb_Torus.calc_coeff();     // calculate perturbation theory coeff (integral form)
                Chua_Perturb_d_e.calc_perturb_d_e();   // calculate perturbation theory coeff (d.e. form)
            }
        });

        btnFitCubic.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                //fit_cubic_response();               // original nonlinear response using a finite perturbation
                fit_cubic_response_full();          // new nonlinear response at fixed limit cycle Phi(t)
                //fit_cubic_response_projected_2D();  // (ABANDONED CODE) new 2D projected response at fixed limit cycle Phi(t)
            }
        });

        addWindowListener(new WindowAdapter() {
            @Override public void windowClosing(WindowEvent ev)
                {
                    Main.save_prefs();
                    Main.chua_euler_slider.dispose();
                }
            });

        for (JTextField tx: txt)
            tx.addKeyListener(new KeyAdapter()
            {
                @Override public void keyPressed(KeyEvent e)
                {
                    if (e.getKeyCode() == KeyEvent.VK_ESCAPE)
                    {
                        Main.save_prefs();
                        Main.chua_euler_slider.dispose();
                        Main.plot_y_vs_x.dispose();
                        Main.z_bifurcate = new Chua_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                }
            });
        //System.out.println("dataPanel = " + dataPanel[0].getSize());
    }

    private void phase_space(boolean ch, int Period, boolean fast)
    {
        double Tnew, Tsum;        // time of peak z
        GregorianCalendar now = new GregorianCalendar();
        PrintWriter fout = null;
        double[] pt6 = new double[] {Main.x0, Main.y0, Main.z0, Main.dx0du, Main.dy0du, Main.dz0du};
        double xmin, xmax, ymin, ymax, zmin, zmax;
        Point2D.Double pt2;                                 // projected (x', y') using Euler angles
        Point2D.Double pstat = Main.project_stationary();   // projected stationary point (x', y')
        //pstat = Main.project_2D(-0.37992427533176765 - 8.50122819e-01 , 0.000130106962529913 - 1.03983278e-02, 0.38005438 + 5.26481782e-01);
        //System.out.println("pstat = " + pstat.x + ", " + pstat.y + ", " + Main.project_phi + ", " + Main.project_theta);
        int Nloop = N;                                      // number of iterations
        String lblhdr;
        String fmt = "%.3f";
        if (fast)
            Nloop = 10*N;                                   // execute 10 times faster

        if (Period == 0)
            delt = DEFAULT_DELT;
        if (ch)
            iT = 0;
        else
            System.arraycopy(pt6_old, 0, pt6, 0, pt6.length);   // re-use previous run
        if (true)                                               // dump raw (x,y,z) to file
            try
            {
                String fname = "C:\\Windows\\Temp\\Chua_Output_" + Main.alpha + "_" + String.format("%.6f", delt) + ".csv";
                boolean fexists = new File(fname).exists();
                FileWriter fw = new FileWriter(fname, true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    fout.println("Chua x y z vs. i, " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c
                                                      + ", " + Period + ", " + delt + ", " + N
                                                      + ", " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi);
                    fout.println("iT, x', y', z'");
                    fout.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}

        if (ch || printChk.isSelected() || first)
        {
            if (phaseRadio.isSelected())
            {
                //System.out.println("Chua y' vs. x', " + now.getTime() + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Period + ", " + Main.project_phi + ", " + Main.project_theta + ", " + delt + ", " + N);
                System.out.print("Chua y' vs. x', " + now.getTime() + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                //System.out.println(", " + Main.calc_phi(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_theta(pt6[0], pt6[1], pt6[2]));
                System.out.println(", " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi);
            }
            else
                System.out.println("Chua dy/du vs. dx/du, " + now.getTime() + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Period + ", " + delt + ", " + N);
            //pt2 = Main.project_2D(0, 0.02454568025564129, 0.18455654708900665);      // temporary code debug only
            //System.out.println("cos coeff, " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(0, 0.02454568025564129, 0.18455654708900665));
            //System.out.println("a,b,c = " + );
            //double zt = (Main.c + Math.sqrt(Main.c*Main.c - 4*Main.a*Main.b))/2/Main.a;
            //System.out.println("P1 = (" + Main.a*zt + ", " + (-zt) + ", " + zt + ")");
            //zt = (Main.c - Math.sqrt(Main.c*Main.c - 4*Main.a*Main.b))/2/Main.a;
            //System.out.println("P2 = (" + Main.a*zt + ", " + (-zt) + ", " + zt + ")");
        }
        first = false;
        if (printChk.isSelected())
            if (phaseRadio.isSelected())
            {
                System.out.println("iT, x, y, z");
                System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                //System.out.println("start = ," + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
            }
            else
            {
                System.out.println("iT, x, y, z, dx/du, dy/du, dz/du");
                System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
            }
        if (phaseRadio.isSelected())                        // normal phase space
        {
            pt2 = Main.project_2D(pt6[0], pt6[1], pt6[2]);
            xmin = pt2.x;
            xmax = pt2.x;
            ymin = pt2.y;
            ymax = pt2.y;
            zmin = 0;
            zmax = 0;
        }
        else                                                    // ddc/ddy tangent phase space
        {
            xmin = pt6[3];
            xmax = pt6[3];
            ymin = pt6[4];
            ymax = pt6[4];
            zmin = pt6[5];
            zmax = pt6[5];
        }
        int Nstart = 0; //Nloop - 10000;
        int Nband = 2000; // 2000;
        path1.reset();
        path2.reset();
        path3.reset();
        path1.moveTo(xmin, ymin);
        for (int j = 0; j < Nloop; j++)
        {
            iT++;
            if (phaseRadio.isSelected())                        // normal phase space
            {
                Main.runge_kutta_chua3(pt6, delt);
                pt2 = Main.project_2D(pt6[0], pt6[1], pt6[2]);
                //if (j < 100)
                //    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                if (pt2.x > xmax) xmax = pt2.x;
                if (pt2.x < xmin) xmin = pt2.x;
                if (pt2.y > ymax) ymax = pt2.y;
                if (pt2.y < ymin) ymin = pt2.y;
                //if (pt2.x > xmax && pt2.x <  10) xmax = pt2.x;
                //if (pt2.x < xmin && pt2.x > -10) xmin = pt2.x;
                //if (pt2.y > ymax && pt2.y <  10) ymax = pt2.y;
                //if (pt2.y < ymin && pt2.y > -10) ymin = pt2.y;
                if (j == Nloop - Nband)                 // initiallize path3
                    path3.moveTo(pt2.x, pt2.y);
                else if (j > Nloop - Nband)             // draw path3 (blue)
                    path3.lineTo(pt2.x, pt2.y);
                else if (j == Nstart + Nband)           // initiallize path2
                    path2.moveTo(pt2.x, pt2.y);
                else if (j > Nstart + Nband)            // draw path2 (orange)
                    path2.lineTo(pt2.x, pt2.y);
                else if (j > Nstart)                    // draw path1 (green)
                    path1.lineTo(pt2.x, pt2.y);
                if (fout != null && printChk.isSelected())
                    //fout.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);                                 // org
                    fout.println(iT + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt6[0], pt6[1], pt6[2]));    // projected
                    //fout.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]
                    //                + ", " + Main.calc_xdot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_ydot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_zdot(pt6[0], pt6[1], pt6[2]));
                if (printChk.isSelected() && Period > 0 && j >= Nloop - Period - 1)     // transfer last cycle
                {
                    //System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);   // normal code, KEEP
                    System.out.println(iT + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt6[0], pt6[1], pt6[2]));   // temporary code (replace)
                    //System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]
                    //                      + ", " + Main.calc_xdot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_ydot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_zdot(pt6[0], pt6[1], pt6[2])
                    //                      + ", " + Main.calc_x2dot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_y2dot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_z2dot(pt6[0], pt6[1], pt6[2])
                    //                      + ", " + Main.calc_x3dot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_y3dot(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_z3dot(pt6[0], pt6[1], pt6[2]));
                }
                if (j == Nloop - 1)
                {
                    System.out.print("end = , " + iT + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                    System.out.println(", " + Main.calc_phi(pt6[0], pt6[1], pt6[2]) + ", " + Main.calc_theta(pt6[0], pt6[1], pt6[2]));
                    Main.final_Period = Period;
                    Main.final_delt = delt;
                    Main.final_x = pt6[0];
                    Main.final_y = pt6[1];
                    Main.final_z = pt6[2];
                }
            }
            else if (ddu3Radio.isSelected())               // Lyapunov response (3-D)
            {
                Main.runge_kutta_chua6_ddu3(pt6, delt);
                if (pt6[3] > xmax) xmax = pt6[3];
                if (pt6[3] < xmin) xmin = pt6[3];
                if (pt6[4] > ymax) ymax = pt6[4];
                if (pt6[4] < ymin) ymin = pt6[4];
                if (pt6[5] > zmax) zmax = pt6[5];
                if (pt6[5] < zmin) zmin = pt6[5];
                path1.lineTo(pt6[3], pt6[4]);
                if (printChk.isSelected() && j == Period - 1)                   // one-shot response after the first cycle
                {
                    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
                    double xdot = Main.calc_xdot(pt6[0], pt6[1], pt6[2]);
                    double ydot = Main.calc_ydot(pt6[0], pt6[1], pt6[2]);
                    double zdot = Main.calc_zdot(pt6[0], pt6[1], pt6[2]);
                    double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
                    System.out.println(iT + ", " + xdot/v + ", " + ydot/v + ", " + zdot/v
                                          + ", " + Main.calc_phi(pt6[0], pt6[1], pt6[2])
                                          + ", " + Main.calc_theta(pt6[0], pt6[1], pt6[2]));
                    System.out.println("ddu3 range, " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Period + ", " + delt + ", " + xmin + ", " + xmax + ", " + ymin + ", " + ymax + ", " + zmin + ", " + zmax);
                }
            }
            if (zlist[0] < zlist[1] && zlist[1] <= zlist[2] && zlist[2] > zlist[3] && zlist[3] > pt6[2])    // maximum of z coord
            {
                Tindex = (Tindex + 1) % Nfork;
                //Tnew = Main.parabola(iT - 3, iT - 2, iT - 1, zlist[1], zlist[2], zlist[3], false);
                Tnew = iT - 2 + Main.quarticT(zlist[0] - zlist[2], zlist[1] - zlist[2], zlist[3] - zlist[2], pt6[2] - zlist[2]); // use z peak
                //System.out.println("zlist, " + zlist[0] + ", " + zlist[1] + ", " + zlist[2] + ", " + zlist[3] + ", " + pt6[1]);
                Tfork[Tindex] = Tnew - Told;
                //System.out.println((iT - 2) + ", " + Main.c + ", " + zlist[2] + ", " + Tnew
                //                     + ", " + Main.parabola(iT - 3, iT - 2, iT - 1, zlist[1], zlist[2], zlist[3], true)
                //                     + ", " + Tfork[Tindex]);
                Told = Tnew;
                if (Tindex == Nfork - 1 && !printChk.isSelected())
                {
                    Tsum = 0;
                    System.out.print("sum = ," + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + (iT - 2));
                    for (double tf: Tfork)
                    {
                        System.out.print(", " + tf);
                        Tsum += tf;
                    }
                    System.out.println(", " + Tsum + ", " + Period + ", " + delt); // + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                    if (Period > 0 && true)
                        delt *= Tsum/Period;
                }
            }
            for (int k = 0; k < 3; k++)
                zlist[k] = zlist[k + 1];
            zlist[3] = pt6[2];
            if (Period > 0 && false)
                if (iT % Period == 0)
                    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
        }
        if (fout != null)
            fout.close();
        System.arraycopy(pt6, 0, pt6_old, 0, pt6.length);               // save
        AffineTransform at = new AffineTransform((phasePanel.getWidth() - 8)/(xmax - xmin),
                                          0, 0, -(phasePanel.getHeight() - 8)/(ymax - ymin),
                                                (xmax*4 + xmin*(4 - phasePanel.getWidth()))/(xmax - xmin),
                                               -(ymin*4 + ymax*(4 - phasePanel.getHeight()))/(ymax - ymin));
        path1.transform(at);
        path2.transform(at);
        path3.transform(at);
        stataxis.moveTo(-pstat.x, -pstat.y);    // this is a line joining 2 stationary points
        stataxis.lineTo( pstat.x,  pstat.y);    // this line will be slightly off-center because of the border used in 'at'
        stataxis.transform(at);

        //System.out.println("range, " + Main.c + ", " + xmin + ", " + xmax + ", " + ymin + ", " + ymax + ", " + zmin + ", " + zmax);
        lblhdr = phaseRadio.isSelected() ? "x" : "dx";
        lblxrange.setText(lblhdr + " = " + String.format(fmt, xmin) + ", " + String.format(fmt, xmax));
        lblhdr = phaseRadio.isSelected() ? "y" : "dy";
        lblyrange.setText(lblhdr + " = " + String.format(fmt, ymin) + ", " + String.format(fmt, ymax));
        lblhdr = phaseRadio.isSelected() ? "z" : "dz";
        lblzrange.setText(lblhdr + " = " + String.format(fmt, zmin) + ", " + String.format(fmt, zmax));
        xaxis = new Line2D.Double(0, ymax*phasePanel.getHeight()/(ymax - ymin), phasePanel.getWidth(), ymax*phasePanel.getHeight()/(ymax - ymin));
        yaxis = new Line2D.Double(-xmin*phasePanel.getWidth()/(xmax - xmin), 0, -xmin*phasePanel.getWidth()/(xmax - xmin), phasePanel.getHeight());
        phasePanel.repaint();
    }

    protected static String fit_linear_response()
    {
        // collect data to fit first-order response, after one cycle,
        // to a change in the (x, y, z) space, using runge_kutta_chua6_ddu3

        double[][] M = new double[2][2];
        //double[] xfer = new double[] {7.25, 16.0, 0.0, 1.0, -0.143, 500, 0.0038218641612062978, -0.37815340802372854, 5.423349576151025E-14, 0.37815340802378705};
        //double[] xfer = new double[] {101.0, 1499.25037, -0.51325, -1.0, 0.144, 200, -8.391111663733648E-4, -0.06280655416186308, -0.11453242046181626, 0.7194164907563163};
        //double[] xfer = new double[] {99.99, 1499.25037, -0.51325, -1.0, 0.144, 200, -8.388275344381413E-4, 0.21346406564909817, 0.07991823532663608, -3.3774899100816005};
        double[] xfer = new double[] {1609920000, 99.977, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990198943096418E-5, -0.18760431693110263, -0.1119262522179699, 1.48518071995642};
        //double[] xfer = new double[] {478080000, 99.978, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990201287968368E-5, -0.18760915811977436, -0.11192486733843215, 1.485180709880168};
        //double[] xfer = new double[] {347997600, 99.980, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990205973471885E-5, -0.18761359939889577, -0.11192288909528943, 1.4850862469573296};
        //double[] xfer = new double[] {347997700, 99.980, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990205973471885E-5, -0.11970240239014328, -0.117823496370916, 0.27107075945853604};
        //double[] xfer = new double[] {347997700, 99.980, 1499.25037, -0.51325, -1.0, 0.144, 24000, -6.990205973987019E-6, -0.11970240836566291, -0.11782349676998444, 0.2710706984438698};
        //double[] xfer = new double[] {347998900, 99.980, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990205973471885E-5, 0.033261747274911733, 0.11757542743865339, -0.17406255822626843};
        //double[] xfer = new double[] {347998900, 99.980, 1499.25037, -0.51325, -1.0, 0.144, 24000, -6.990205973360493E-6, 0.03326174478635965, 0.11757542749244691, -0.1740624919192874};
        //double[] xfer = new double[] {162720000, 99.980, 1499.25037, -0.51325, -1.0, 0.144, 9600, -1.747551493349202E-5, -0.11970240112953173, -0.11782349632166472, 0.27107074206366416};
        //double[] xfer = new double[] {1396320000, 99.983, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990212988270638E-5, -0.18762753711114005, -0.11191875877106697, 1.4850786625558816};
        //double[] xfer = new double[] {997920000, 99.986, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.99022002895545E-5, -0.18764188181492653, -0.11191464627359901, 1.4850788900013745};
        //double[] xfer = new double[] {757920000, 99.988, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990224713151071E-5, -0.1876514481504697, -0.11191187423625534, 1.4850786636760196};
        //double[] xfer = new double[] {530880000, 99.990, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.99022939734161E-5, -0.1876610226554186, -0.11190909945602992, 1.485078415215612};
        //double[] xfer = new double[] {320640000, 99.992, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990234080492677E-5, -0.18767053882372814, -0.11190633063345796, 1.4850781520490064};
        //double[] xfer = new double[] {1187520000, 99.994, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990238764870291E-5, -0.07442781388230205, -0.11741216153578302, -0.4929786723975686};
        //double[] xfer = new double[] {121440000, 99.9940, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990238765861372E-5, -0.18768009415178913, -0.11190356267981405, 1.4850779591762608};
        //double[] xfer = new double[] {237600000, 99.9944, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990239702576763E-5, -0.1876819989455841, -0.11190300881127653, 1.485077807549828};
        //double[] xfer = new double[] {69597600, 99.9948, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990240639707452E-5, -0.24848425467715657, -0.09786218088127782, 2.62850755558779};
        //double[] xfer = new double[] {180000000, 99.994810, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990240663172005E-5, -0.18768396485973995, -0.11190244052471447, 1.485078049480874};
        //double[] xfer = new double[] {423840000, 99.99482, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990240686565814E-5, -0.1876840026873076, -0.11190242816369231, 1.485077770460996};
        //double[] xfer = new double[] {290880000, 99.994821, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990240688952322E-5, -0.18768401920010822, -0.11190242521016512, 1.4850780500284413};
        //double[] xfer = new double[] {56640000, 99.994822, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990240691231722E-5, -0.18768400101962426, -0.11190242704156188, 1.485077745297063};
        //double[] xfer = new double[] {484320000, 99.99484, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990240733438612E-5, -0.18768409146260298, -0.11190240164480135, 1.4850777728383322};
        //double[] xfer = new double[] {14400000, 8.6956522, 14.2857, 0.0, 1.0, -0.42857, 2400, 7.991038476907487E-4, -0.6492164661497843, 0.21873831729699572, 3.17920947326913};
        //double[] xfer = new double[] {22080000, 8.6956522, 14.2857, 0.0, 1.0, -0.42857, 2400, 7.99103847690894E-4, -0.9833648715844555, -0.4448294298769549, 2.732230644119052};
        //double[] xfer = new double[] {19200000, 15.5279503, 28.3125708, 0.0, 1.0, -0.205, 2400, 5.201437828342687E-4, -0.8680072875863561, -0.3493391869076231, 5.300476743138388};
        //double[] xfer = new double[] {101280000, 100.0, 1499.25037, -0.51325, -1.0, 0.144, 2400, -6.990252819409822E-5, 0.12261973007393918, 0.11771544495390424, -0.32001502739612614};
        //double[] xfer = new double[] {7200000, -75.0187547, 31.8421053, -2.8947368, -1.0, 0.681, 2400, 3.561946267089565E-4, 0.8458471046904906, -0.012695584681198458, 1.0200109622035856};
        //double[] xfer = new double[] {12960000, 8.3333333, 16.0, 0.0, 1.0, -0.143, 2400, 7.354632501800597E-4, -0.5170548995280515, -0.6562639991620063, -1.0980044820904584};

        if (true)                                           // initiallize directly from data file, row Main.iT
        {
            String fDir = "\\APP\\Java\\ChuaOscillator\\skew_vs_time\\";
            //String fName = "xyz_input_2400_99.98";
            //String fName = "xyz_input_2400_99.9765";
            String fName = "xyz_input_24_99.9948";
            String str;
            try
            {
                BufferedReader istr = new BufferedReader(new FileReader(fDir + fName + ".csv"));
                try
                {
                    str = istr.readLine();                                          // read header
                    //System.out.println(str);
                    xfer = new double[] {Double.NaN,                                // iT
                                         Double.parseDouble(str.split(",")[2]),
                                         Double.parseDouble(str.split(",")[3]),
                                         Double.parseDouble(str.split(",")[4]),
                                         Double.parseDouble(str.split(",")[5]),
                                         Double.parseDouble(str.split(",")[6]),
                                         Double.parseDouble(str.split(",")[7]),
                                         Double.parseDouble(str.split(",")[8]),
                                         Double.NaN, Double.NaN, Double.NaN};       // x, y, z
                    istr.readLine();                                                // skip 1 line
                    for (int i = 0; i < Main.iTinit; i++)                           // skip lines
                        istr.readLine();
                    str = istr.readLine();                                          // read data
                    //System.out.println(str);
                    xfer[0]  = Double.parseDouble(str.split(",")[0]);
                    xfer[8]  = Double.parseDouble(str.split(",")[1]);
                    xfer[9]  = Double.parseDouble(str.split(",")[2]);
                    xfer[10] = Double.parseDouble(str.split(",")[3]);
                    istr.close();
                    Main.iTinit += 0;
                }
                catch (IOException e)
                {
                    System.out.println("read error : " + e.getMessage());
                    return "read error";
                }
            }
            catch (FileNotFoundException e)
            {
                System.out.println("file not found : " + e.getMessage());
                return "file not found";
            }
        }

        iT = (int) xfer[0];
        Main.alpha = xfer[1];
        Main.beta  = xfer[2];
        Main.gamma = xfer[3];
        Main.c     = xfer[5];
        Main.final_x = xfer[8];                 // emergency initiallization only
        Main.final_y = xfer[9];
        Main.final_z = xfer[10];
        Main.final_Period = (int) xfer[6];
        Main.final_delt = xfer[7];
        Main.final_delt *= -1;                  // temporarily reverse the direction of time

        //Main.final_Period *= 10;        // TEMPORARY CODE
        //Main.final_delt /= 10;          // TEMPORARY CODE

        Main.skew_transform = false;            // default = false, try to detect "non-uniformity"
        double incr = 1;
        double xdot = Main.calc_xdot(Main.final_x, Main.final_y, Main.final_z);
        double ydot = Main.calc_ydot(Main.final_x, Main.final_y, Main.final_z);
        double zdot = Main.calc_zdot(Main.final_x, Main.final_y, Main.final_z);
        double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
        Main.project_phi = Main.calc_phi(Main.final_x, Main.final_y, Main.final_z);
        Main.project_theta = Main.calc_theta(Main.final_x, Main.final_y, Main.final_z);
        double[] pt6 = new double[6];

        System.out.println("\nPython Chua output (fit_linear_response):");
        System.out.println("linear_hdr = \"\\n\\");
        System.out.println("ddu3 - measure linear response after one cycle\\n\\");
        System.out.println("incr            , " + incr + ", \\n\\");
        System.out.println("alpha_beta_gamma, " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", \\n\\");
        System.out.println("Period_delt     , " + Main.final_Period + ", " + Main.final_delt + ", \\n\\");
        System.out.println("x_y_z           , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ", \\n\\");
        System.out.println("vx_vy_vz        , " + xdot/v + ", " + ydot/v + ", " + zdot/v + ", \\n\\");
        System.out.println("phi_theta_psi   , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ", \\n\"");

        System.out.print("M = np.array([");
        for (int i = 0; i < 3; i++)         // initiallize dx/du, dy/du, or dz/du
        {
            if (i == 0)                                                     // initial (dxdu', dydu', dzdu')
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, incr*Main.invert_from_xp_yp(1, 0, 0, "x"), incr*Main.invert_from_xp_yp(1, 0, 0, "y"), incr*Main.invert_from_xp_yp(1, 0, 0, "z")};
            else if (i == 1)
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, incr*Main.invert_from_xp_yp(0, 1, 0, "x"), incr*Main.invert_from_xp_yp(0, 1, 0, "y"), incr*Main.invert_from_xp_yp(0, 1, 0, "z")};
            else
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, incr*xdot/v, incr*ydot/v, incr*zdot/v};
            //Point2D.Double pt2temp = Main.project_2D(pt6[3], pt6[4], pt6[5]);   // final (dxdu', dydu', dzdu')
            //double zptemp = Main.project_zp(pt6[3], pt6[4], pt6[5]);
            //System.out.println("init ," + i + ", " + pt2temp.x + ", " + pt2temp.y + ", " + zptemp);
            //System.out.println("0" + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);                // temporary code
            for (int k = 0; k < Main.final_Period; k++)                     // loop through one cycle
            {
                Main.runge_kutta_chua6_ddu3(pt6, Main.final_delt);
                //System.out.println((k + 1) + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);        // temporary code
            }
            Point2D.Double pt2 = Main.project_2D(pt6[3], pt6[4], pt6[5]);   // final (dxdu', dydu', dzdu')
            double zp = Main.project_zp(pt6[3], pt6[4], pt6[5]);
            System.out.print("[" + pt2.x + ", " + pt2.y + ", " + zp + "]");
            if (i < 2)
            {
                System.out.println(",");
                M[0][i] = pt2.x;
                M[1][i] = pt2.y;
            }
        }
        System.out.println("])");
        System.out.println("\nfinal x_y_z , " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
        eig = Math.sqrt(M[0][0]*M[1][1] - M[0][1]*M[1][0]);
        if (Math.abs(M[0][0] + M[1][1]) < 2*eig)
        {
            angle = Math.acos((M[0][0] + M[1][1])/2/eig);
            System.out.println("|eig|_angle , " + eig + ", " + angle*180/Math.PI + ", " + eig*Math.cos(angle) + ", " + -eig*Math.sin(angle));
            Main.final_Re_V21 = (M[1][1] - M[0][0])/2/M[0][1];          // define the skew-transform matrix
            Main.final_Im_V21 = Math.sqrt(-(M[1][1] - M[0][0])*(M[1][1] - M[0][0]) - 4*M[0][1]*M[1][0])/2/M[0][1];
            //System.out.println("test transform, " + iT + ", " + eig + ", " + angle*180/Math.PI + ", " + M[0][0] + ", " + M[0][1] + ", " + Main.final_Re_V21);
            //System.out.println("              , " + iT + ", " + eig + ", " + angle*180/Math.PI + ", " + M[1][0] + ", " + M[1][1] + ", " + Main.final_Im_V21);
            //System.out.println(" ");
            double sqr = Math.sqrt((M[0][1] + M[1][0])*(M[0][1] + M[1][0]) + (M[1][1] - M[0][0])*(M[1][1] - M[0][0]));
            System.out.println("ellipse phi_a/b = ," + Math.atan2(M[1][1] - M[0][0], M[0][1] + M[1][0])*180/Math.PI/2 + ", "
                                                     + Math.sqrt((M[0][1] - M[1][0] + sqr)/(M[0][1] - M[1][0] - sqr)));
            System.out.println("Re_V21 Im_V21 = , " + iT + ", " + Main.final_Re_V21 + ", " + Main.final_Im_V21);
//            System.out.println("summary, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + eig + ", " + angle*180/Math.PI + ", " + Main.final_Re_V21 + ", " + Main.final_Im_V21 + ", " + Math.atan2(M[1][1] - M[0][0], M[0][1] + M[1][0])*180/Math.PI/2 + ", " + Math.sqrt((M[0][1] - M[1][0] + sqr)/(M[0][1] - M[1][0] - sqr)));
            System.out.println("fit linear, " + iT + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + + Main.final_Period + ", " + Main.final_delt + ", " + eig + ", " + angle*180/Math.PI + ", " + eig*Math.cos(angle) + ", " + -eig*Math.sin(angle) + ", " + Main.final_Im_V21);
        }
        else
        {
            System.out.println("real eig, " + iT + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + + Main.final_Period + ", " + Main.final_delt
                                            + ", " + (M[0][0] + M[1][1] - Math.sqrt((M[1][1] - M[0][0])*(M[1][1] - M[0][0]) + 4*M[0][1]*M[1][0]))/2
                                            + ", " + (M[0][0] + M[1][1] + Math.sqrt((M[1][1] - M[0][0])*(M[1][1] - M[0][0]) + 4*M[0][1]*M[1][0]))/2);
        }
        //System.out.println("summary, " + M[0][0] + ", " + M[0][1] + ", " + Main.final_Re_V21);
        //System.out.println("summary, " + M[1][0] + ", " + M[1][1] + ", " + Main.final_Im_V21);
        return "" + iT + ", " + eig + ", " + angle*180/Math.PI;
    }

    private static void fit_cubic_response()
    {
        // collect Chua oscillator data to fit third-order response, after one cycle,
        // to a change in the (x', y') plane, perpendicular to velocity.
        // assume that the output has already been made uniform, by running 'fit_linear_response()'

        if (first_order_hdr == null)
        {
            System.out.println("'first_order_hdr' is not initialized");
            return;
        }
        //System.out.println("first-order hdr = '" + first_order_hdr + "'");
        Main.skew_transform = true;             // make the linear response "uniform"
        double incr = 0.001;     // 0.0001/2
        double[] pt3;
        int Nangl = 24;                         // # of angular positions (for radial grid)
        int Nrad = 3;                           // # of radii > 0 (for radial grid)
        //int Nrect = 3;                            // # of coordinates > 0 (for rectangular grid)

        Point2D.Double pt2 = Main.project_2D(Main.final_x, Main.final_y, Main.final_z); // initial (x', y')
        double zp = Main.project_zp(Main.final_x, Main.final_y, Main.final_z);          // setpoint for z'
        double xold = 0, yold = 0, zold = zp;
        Point2D.Double pt2old;          // projected, old, (x', y') using Euler angles
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        double zpold, zpnew;            // projected z'
        double xc, yc;                  // interpolated, projected (x', y')

        //System.out.println("projected zpdot = " + Main.project_zp(Main.calc_xdot(Main.final_x, Main.final_y, Main.final_z), Main.calc_ydot(Main.final_x, Main.final_y, Main.final_z), Main.calc_zdot(Main.final_x, Main.final_y, Main.final_z)));
        System.out.println("\nPython output (fit_cubic_response):");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Chua - Neimark-Sacker - measure cubic model response after one cycle\\n\\");
        System.out.println("incr_iT_eig_angle, " + incr + ", " + first_order_hdr + ", \\n\\");
        System.out.println("alpha_beta_gamma , " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt      , " + Main.final_Period + ", " + Main.final_delt + ", " + Nangl + ", " + Nrad + ",\\n\\");
        //System.out.println("Period_delt      , " + Main.final_Period + ", " + Main.final_delt + ", " + Nrect + ",\\n\\");
        System.out.println("x_y_z            , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi    , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("x'_y'_z'         , " + pt2.x + ", " + pt2.y + ", " + zp + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',tc\"");
        System.out.print("data = np.array([");
        //for (int i = -Nrect; i < Nrect + 1; i++)          // increment x' by i*incr
        //    for (int j = -Nrect; j < Nrect + 1; j++)      // increment y' by j*incr
        for (int i = 0; i < Nrad + 1; i++)      // increment r' radius by i*incr
            for (int j = 0; j < Nangl; j++)     // increment angle by 360/Nposn degrees
                if (i > 0 || j == 0)
            {
                pt3 = new double[] {Main.final_x + Main.invert_from_xp_yp(incr*i*Math.cos(j*2*Math.PI/Nangl), incr*i*Math.sin(j*2*Math.PI/Nangl), 0, "x"),
                                    Main.final_y + Main.invert_from_xp_yp(incr*i*Math.cos(j*2*Math.PI/Nangl), incr*i*Math.sin(j*2*Math.PI/Nangl), 0, "y"),
                                    Main.final_z + Main.invert_from_xp_yp(incr*i*Math.cos(j*2*Math.PI/Nangl), incr*i*Math.sin(j*2*Math.PI/Nangl), 0, "z")};
                //pt3 = new double[] {Main.final_x + Main.invert_from_xp_yp(incr*i, incr*j, 0, "x"),
                //                    Main.final_y + Main.invert_from_xp_yp(incr*i, incr*j, 0, "y"),
                //                    Main.final_z + Main.invert_from_xp_yp(incr*i, incr*j, 0, "z")};
                //System.out.println("org   , " + i + ", " + j + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                //if (i == 0 && j == 0 || i == 1 && j == 0)
                //    System.out.println("start -   , " + i + ", " + j + ", " + 0 + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                // loop through one cycle
                for (int k = 0; k < Main.final_Period + 0*40; k++) // add 10 iterations just to be sure it crosses (KEEP)
                {
                    //if (k == 0 || k == 10 || k == 20)     // temporary debugging code
                    //    System.out.println("debug, " + k + ", " + i + ", " + j + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                    Main.runge_kutta_chua3(pt3, Main.final_delt);
                    //if (k + 1 == 24*1000)                       // temporary test code for in-flight response at k
                    //{
                    //    pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                    //    zpnew = Main.project_zp(pt3[0], pt3[1], pt3[2]);
                    //    System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");
                    //    if (i < Nrad || j < Nangl - 1)
                    //        System.out.println(",");
                    //}
                    //if ((i == 0 && j == 0 || i == 1 && j == 0) && ((k + 1) % 10 == 0))
                    //    System.out.println("---    , " + i + ", " + j + ", " + (k + 1) + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                        //System.out.println((k+1) + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + Main.calc_phi(pt3[0], pt3[1], pt3[2]) + ", " + Main.calc_theta(pt3[0], pt3[1], pt3[2]));
                    if (false && k > 1 && (Main.project_zp(xold, yold, zold) - zp)*Main.final_delt <= 0 && (Main.project_zp(pt3[0], pt3[1], pt3[2]) - zp)*Main.final_delt > 0)
                    {
                        pt2old = Main.project_2D(xold, yold, zold);
                        pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                        zpold = Main.project_zp(xold, yold, zold);
                        zpnew = Main.project_zp(pt3[0], pt3[1], pt3[2]);
                        //System.out.println("zp = ," + zpold + ", " + zpnew);
                        //System.out.println("old = ," + xold + ", " + yold + ", " + zold);
                        //System.out.println("new = ," + i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                        xc = ((zpnew - zp)*pt2old.x + (zp - zpold)*pt2new.x)/(zpnew - zpold);
                        yc = ((zpnew - zp)*pt2old.y + (zp - zpold)*pt2new.y)/(zpnew - zpold);
                        System.out.print("[" + i + ", " + j + ", " + xc + ", " + yc + ", " + (k + (zp - zpold)/(zpnew - zpold)) + "]");
                        //if (i < Nrect || j < Nrect)
                        //    System.out.println(",");
                        if (i < Nrad || j < Nangl - 1)
                            System.out.println(",");
                    }
                    xold = pt3[0];
                    yold = pt3[1];
                    zold = pt3[2];
                }
                if (!false)                        // use last point (synchronize in time)
                {
                    pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                    zpnew = Main.project_zp(pt3[0], pt3[1], pt3[2]);
                    System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");
                    //if (i < Nrect || j < Nrect)
                    //    System.out.println(",");
                    if (i < Nrad || j < Nangl - 1)
                        System.out.println(",");
                }
                if (false)                        // fit unprojected data (org coordinates)
                {
                    //System.out.print("[" + i + ", " + j + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + "]");
                    //System.out.print("[" + i + ", " + j + ", " + Main.calc_xdot(pt3[0], pt3[1], pt3[2]) + ", " + Main.calc_ydot(pt3[0], pt3[1], pt3[2]) + ", " + Main.calc_zdot(pt3[0], pt3[1], pt3[2]) + "]");
                    System.out.print("[" + i + ", " + j + ", " + Main.calc_x2dot(pt3[0], pt3[1], pt3[2]) + ", " + Main.calc_y2dot(pt3[0], pt3[1], pt3[2]) + ", " + Main.calc_z2dot(pt3[0], pt3[1], pt3[2]) + "]");
                    //if (i < Nrect || j < Nrect)
                    //    System.out.println(",");
                    if (i < Nrad || j < Nangl - 1)
                        System.out.println(",");
                }
            }
        System.out.println("])");
    }

    private static void fit_cubic_response_full()
    {
        // collect Chua oscillator data to fit full nonlinear response, after one cycle,
        // to a change in the (x', y') plane, perpendicular to velocity.
        // assume that the output has already been made uniform, by running 'fit_linear_response()'
        // use 'Main.runge_kutta_chua6_ddu3_full' which leaves the original limit cycle unchanged

        if (first_order_hdr == null)
        {
            System.out.println("'first_order_hdr' is not initialized");
            return;
        }
        //System.out.println("first-order hdr = '" + first_order_hdr + "'");
        Main.skew_transform = true;             // make the linear response "uniform"
        double incr = 0.001;     // 0.0001/2
        double[] pt6 = new double[6];
        int Nangl = 24;                         // # of angular positions (for radial grid)
        int Nrad = 3;                           // # of radii > 0 (for radial grid)

        Point2D.Double pt2 = Main.project_2D(Main.final_x, Main.final_y, Main.final_z); // initial (x', y')
        double zp = Main.project_zp(Main.final_x, Main.final_y, Main.final_z);          // setpoint for z'
        double xold = 0, yold = 0, zold = 0;
        Point2D.Double pt2old;          // projected, old, (x', y') using Euler angles
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        double zpold, zpnew;            // projected z'
        double xc, yc;                  // interpolated, projected (x', y')

        //System.out.println("projected zpdot = " + Main.project_zp(Main.calc_xdot(Main.final_x, Main.final_y, Main.final_z), Main.calc_ydot(Main.final_x, Main.final_y, Main.final_z), Main.calc_zdot(Main.final_x, Main.final_y, Main.final_z)));
        System.out.println("\nPython output (fit_cubic_response_full):");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Chua - Neimark-Sacker - measure cubic model response after one cycle\\n\\");
        System.out.println("incr_iT_eig_angle, " + incr + ", " + first_order_hdr + ", \\n\\");
        System.out.println("alpha_beta_gamma , " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt      , " + Main.final_Period + ", " + Main.final_delt + ", " + Nangl + ", " + Nrad + ",\\n\\");
        System.out.println("x_y_z            , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi    , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("x'_y'_z'         , " + pt2.x + ", " + pt2.y + ", " + zp + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',tc\"");
        System.out.print("data = np.array([");

        // an aside to test 'Main.runge_kutta_chua6_ddu3_full'

        //pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z,       // displace in x' direction
        //                    incr*Main.invert_from_xp_yp(1, 0, 0, "x"),
        //                    incr*Main.invert_from_xp_yp(1, 0, 0, "y"),
        //                    incr*Main.invert_from_xp_yp(1, 0, 0, "z")};
        //System.out.println("start  , " + 0 + ", " + 0 + ", " + 0 + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
        //for (int k = 0; k < 10*Main.final_Period; k++)                      // run 10 loops of limit cycle (TEMPORARY)
        //{
        //    Main.runge_kutta_chua6_ddu3_full(pt6, Main.final_delt);
        //    System.out.println("---    , " + 0 + ", " + 0 + ", " + (k + 1) + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
        //}                                                                   // END of aside
        for (int i = 0; i < Nrad + 1; i++)      // increment r' radius by i*incr
            for (int j = 0; j < Nangl; j++)     // increment angle by 360/Nposn degrees
                if (i > 0 || j == 0)
            {
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z,
                                    incr*i*Main.invert_from_xp_yp(Math.cos(j*2*Math.PI/Nangl), Math.sin(j*2*Math.PI/Nangl), 0, "x"),
                                    incr*i*Main.invert_from_xp_yp(Math.cos(j*2*Math.PI/Nangl), Math.sin(j*2*Math.PI/Nangl), 0, "y"),
                                    incr*i*Main.invert_from_xp_yp(Math.cos(j*2*Math.PI/Nangl), Math.sin(j*2*Math.PI/Nangl), 0, "z")};
                //System.out.println("org   , " + i + ", " + j + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                //if (i == 1 && j == 0)
                //    System.out.println("start -   , " + i + ", " + j + ", " + 0 + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
                // loop through one cycle
                for (int k = 0; k < Main.final_Period + 0*40; k++) // add 10 iterations just to be sure it crosses (KEEP)
                {
                    Main.runge_kutta_chua6_ddu3_full(pt6, Main.final_delt);
                    //if (k + 1 == 24*1000)                       // temporary test code for in-flight response at k
                    //{
                    //    pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                    //    zpnew = Main.project_zp(pt3[0], pt3[1], pt3[2]);
                    //    System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");
                    //    if (i < Nrad || j < Nangl - 1)
                    //        System.out.println(",");
                    //}
                    //if ((i == 1 && j == 0) && ((k + 1) % 10 == 0))
                    //    System.out.println("---    , " + i + ", " + j + ", " + (k + 1) + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
                    if (false && k > 1 && Main.project_zp(xold, yold, zold)*Main.final_delt <= 0 && Main.project_zp(pt6[3], pt6[4], pt6[5])*Main.final_delt > 0)
                    {
                        pt2old = Main.project_2D(xold, yold, zold);
                        pt2new = Main.project_2D(pt6[3], pt6[4], pt6[5]);
                        zpold = Main.project_zp(xold, yold, zold);
                        zpnew = Main.project_zp(pt6[3], pt6[4], pt6[5]);
                        //System.out.println("zp = ," + zpold + ", " + zpnew);
                        //System.out.println("old = ," + xold + ", " + yold + ", " + zold);
                        //System.out.println("new = ," + i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                        xc = (zpnew*pt2old.x - zpold*pt2new.x)/(zpnew - zpold);
                        yc = (zpnew*pt2old.y - zpold*pt2new.y)/(zpnew - zpold);
                        System.out.print("[" + i + ", " + j + ", " + xc + ", " + yc + ", " + (k - zpold/(zpnew - zpold)) + "]");
                        if (i < Nrad || j < Nangl - 1)
                            System.out.println(",");
                    }
                    xold = pt6[3];
                    yold = pt6[4];
                    zold = pt6[5];
                }
                if (!false)                        // use last point (synchronize in time)
                {
                    pt2new = Main.project_2D(pt6[3], pt6[4], pt6[5]);
                    zpnew = Main.project_zp(pt6[3], pt6[4], pt6[5]);
                    System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");
                    if (i < Nrad || j < Nangl - 1)
                        System.out.println(",");
                }
            }
        System.out.println("])");
    }
/*
    private static void fit_cubic_response_projected_2D()
    {
        // collect Chua oscillator data to fit full nonlinear response, after one cycle,
        // to a change in the (x', y') plane, perpendicular to velocity.
        // assume that the output has already been made uniform, by running 'fit_linear_response()' at the start pt.
        // use 'Main.runge_kutta_chua6_ddu3_projected_2D' which leaves the original limit cycle unchanged
        // define plane perpendicular to the velocity vector at every point to produce a 2D decoupled response

        if (first_order_hdr == null)
        {
            System.out.println("'first_order_hdr' is not initialized");
            return;
        }
        Main.skew_transform = true;             // make the linear response "uniform"
        double incr = 0.001;
        double[] pt6 = new double[6];
        int Nangl = 24;                         // # of angular positions (for radial grid)
        int Nrad = 3;                           // # of radii > 0 (for radial grid)

        Point2D.Double pt2 = Main.project_2D(Main.final_x, Main.final_y, Main.final_z); // initial (x', y')
        double zp = Main.project_zp(Main.final_x, Main.final_y, Main.final_z);          // setpoint for z'

        System.out.println("\nPython output (fit_cubic_response_projected_2D):");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Chua - Neimark-Sacker - measure cubic model response after one cycle\\n\\");
        System.out.println("incr_iT_eig_angle, " + incr + ", " + first_order_hdr + ", \\n\\");
        System.out.println("alpha_beta_gamma , " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt      , " + Main.final_Period + ", " + Main.final_delt + ", " + Nangl + ", " + Nrad + ",\\n\\");
        System.out.println("x_y_z            , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi    , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("x'_y'_z'         , " + pt2.x + ", " + pt2.y + ", " + zp + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',tc\"");
        System.out.print("data = np.array([");

        for (int i = 0; i < Nrad + 1; i++)      // increment r' radius by i*incr
            for (int j = 0; j < Nangl; j++)     // increment angle by 360/Nposn degrees
                if (i > 0 || j == 0)
            {
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z,   // original coordinates (x, y, z)
                                    incr*i*Math.cos(j*2*Math.PI/Nangl),         // projected coordinates (dx', dy', dz')
                                    incr*i*Math.sin(j*2*Math.PI/Nangl),
                                    incr*i*0};
                // loop through one cycle
                for (int k = 0; k < Main.final_Period; k++)             // perform t-sync ONLY
                {
                    Main.runge_kutta_chua6_ddu3_projected_2D(pt6, Main.final_delt);
                }
                if (true)                                               // use last point (synchronize in time)
                {
                    System.out.print("[" + i + ", " + j + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5] + "]");
                    if (i < Nrad || j < Nangl - 1)
                        System.out.println(",");
                }
            }
        System.out.println("])");
    }
*/
 }

class Plot_Phase_Panel extends JPanel
{
    @Override public void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;

        if (Main.plot_y_vs_x == null)
        {
            System.out.println("plot_y_vs_x null in paintComponent");
            return;
        }
        g2.setPaint(new Color(0, 128, 0));      // green
        g2.draw(Main.plot_y_vs_x.path1);
        g2.setPaint(new Color(255, 127, 39));   // orange
        g2.draw(Main.plot_y_vs_x.path2);
        g2.setPaint(new Color(0, 0, 192));      // blue
        g2.draw(Main.plot_y_vs_x.path3);
        g2.setPaint(Color.BLUE);
        g2.draw(Main.plot_y_vs_x.xaxis);
        g2.draw(Main.plot_y_vs_x.yaxis);
        g2.setPaint(Color.black);
        g2.draw(Main.plot_y_vs_x.stataxis);
    }
}
