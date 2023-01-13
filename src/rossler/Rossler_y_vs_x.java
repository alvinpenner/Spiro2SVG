
package rossler;

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
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;

public class Rossler_y_vs_x extends JDialog
{
    private boolean first = true;
    private final static double DEFAULT_DELT = 0.02;    // 0.005; // .001223734471716991;
    private final static int N = 480000;                // total # of iterations 160000
    //private final static int N = 6624*75;             // total # of iterations 160000
    private static double[] pt6_old;
    private static double eig, angle;                   // first-order response
    protected Path2D.Double path1 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    //protected Path2D.Double path2 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Line2D.Double xaxis = new Line2D.Double(0, 0, 0, 0);
    protected Line2D.Double yaxis = new Line2D.Double(0, 0, 0, 0);
    private static JCheckBox printChk = new JCheckBox("  print  ");
    private static JRadioButton phaseRadio = new JRadioButton("x-y phase");
    private static JRadioButton ddcRadio = new JRadioButton("d/dc of phase");
    private static JRadioButton ddu3Radio = new JRadioButton("3-D d/du");
    private static JRadioButton ddu2Radio = new JRadioButton("2-D d/du");
    private static JLabel lblxrange;
    private static JLabel lblyrange;
    private static JLabel lblzrange;
    private static JPanel phasePanel = new Plot_Phase_Panel();

    int iT = 0;                                     // # iterations
    double[] zlist = new double[4];                 // previous z values
    double Told = 0;
    int Nfork = 1; // 23;                                  // number of bifurcated branches
    double[] Tfork = new double[Nfork];             // period per branch
    int Tindex = 0;                                 // number of peaks
    double delt = DEFAULT_DELT;

    public Rossler_y_vs_x(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        Main.type = "phase";
        if (phaseRadio.isSelected())                            // normal phase space
            setTitle(" Rossler System - phase y' vs. x'");
        else if (ddcRadio.isSelected())                           // ddy phase space
            setTitle(" Rossler System - dy/dc vs. dx/dc");
        else if (ddu3Radio.isSelected())
            setTitle(" Rossler System - dy/du vs. dx/du (3-D)");  // Lyapunov response
        else if (ddu2Radio.isSelected())
            setTitle(" Rossler System - dy/du vs. dx/du (2-D)");  // Lyapunov response
        else
            setTitle(" Rossler System - phase y' vs. x'");
        setIconImage(img);
        setSize(610 + 135, 493 + 135 + 65);
        setLocationByPlatform(true);

        final JLabel[] lbl = {new JLabel("a"),
                              new JLabel("b"),
                              new JLabel("c"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0"),
                              new JLabel("dx0dc"),
                              new JLabel("dy0dc"),
                              new JLabel("dz0dc"),
                              new JLabel("Period")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.a)),
                                  new JTextField(Double.toString(Main.b)),
                                  new JTextField(Double.toString(Main.c)),
                                  new JTextField(Double.toString(Main.x0)),
                                  new JTextField(Double.toString(Main.y0)),
                                  new JTextField(Double.toString(Main.z0)),
                                  new JTextField(Double.toString(Main.dx0dc)),
                                  new JTextField(Double.toString(Main.dy0dc)),
                                  new JTextField(Double.toString(Main.dz0dc)),
                                  new JTextField()};
        JPanel[] spacerPanel = new JPanel[4];
        JPanel[] dataPanel = new JPanel[lbl.length];
        JButton btnRun = new JButton("Run");
        JButton btnFitLinear = new JButton("Fit Linear");
        JButton btnPerturbCoeff = new JButton("Perturb Coeffs.");
        JButton btnFitCubic = new JButton("Fit Cubic");

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
            lbl[i].setPreferredSize(new Dimension(40, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

        printChk.setOpaque(false);
        JPanel printPanel = new JPanel();
        printPanel.setPreferredSize(new Dimension(130, 24));
        printPanel.setOpaque(false);
        printPanel.add(printChk);

        phaseRadio.setOpaque(false);
        ddcRadio.setOpaque(false);
        ddu3Radio.setOpaque(false);
        ddu2Radio.setOpaque(false);
        ButtonGroup group = new ButtonGroup();
        group.add(phaseRadio);
        group.add(ddcRadio);
        group.add(ddu3Radio);
        group.add(ddu2Radio);
        phaseRadio.setSelected(true);
        JPanel radioPanel = new JPanel();
        radioPanel.setOpaque(false);
        radioPanel.setPreferredSize(new Dimension(120, 4*26));
        radioPanel.setBorder(BorderFactory.createEtchedBorder());
        radioPanel.setLayout(new BoxLayout(radioPanel, BoxLayout.Y_AXIS));
        radioPanel.add(phaseRadio);
        radioPanel.add(ddcRadio);
        radioPanel.add(ddu3Radio);
        radioPanel.add(ddu2Radio);

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
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(dataPanel[7]);
        parmsPanel.add(dataPanel[8]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(printPanel);
        parmsPanel.add(radioPanel);
        parmsPanel.add(xrangePanel);
        parmsPanel.add(yrangePanel);
        parmsPanel.add(zrangePanel);
        parmsPanel.add(spacerPanel[2]);
        parmsPanel.add(dataPanel[9]);
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
        pt6_old = new double[] {Main.x0, Main.y0, Main.z0, Main.dx0dc, Main.dy0dc, Main.dz0dc};

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                boolean changed = false;
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                if (phaseRadio.isSelected())                                // normal phase space
                    setTitle(" Rossler System - phase y' vs. x' (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + String.format(" %.6f", delt) + ")");
                else if (ddcRadio.isSelected())                             // ddc phase space
                    setTitle(" Rossler System - dy/dc vs. dx/dc");
                else if (ddu3Radio.isSelected())
                    setTitle(" Rossler System - dy/du vs. dx/du (3-D)");    // Lyapunov response
                else
                    setTitle(" Rossler System - dy/du vs. dx/du (2-D)");    // Lyapunov response
                //if (Main.a != Double.parseDouble(txt[0].getText())) changed = true;
                Main.a = Double.parseDouble(txt[0].getText());
                //if (Main.b != Double.parseDouble(txt[1].getText())) changed = true;
                Main.b = Double.parseDouble(txt[1].getText());
                //if (Main.c != Double.parseDouble(txt[2].getText())) changed = true;
                Main.c = Double.parseDouble(txt[2].getText());
                if (Main.x0 != Double.parseDouble(txt[3].getText())) changed = true;
                Main.x0 = Double.parseDouble(txt[3].getText());
                if (Main.y0 != Double.parseDouble(txt[4].getText())) changed = true;
                Main.y0 = Double.parseDouble(txt[4].getText());
                if (Main.z0 != Double.parseDouble(txt[5].getText())) changed = true;
                Main.z0 = Double.parseDouble(txt[5].getText());
                //if (Main.dx0dc != Double.parseDouble(txt[6].getText())) changed = true;
                Main.dx0dc = Double.parseDouble(txt[6].getText());
                //if (Main.dy0dc != Double.parseDouble(txt[7].getText())) changed = true;
                Main.dy0dc = Double.parseDouble(txt[7].getText());
                //if (Main.dz0dc != Double.parseDouble(txt[8].getText())) changed = true;
                Main.dz0dc = Double.parseDouble(txt[8].getText());
                phase_space(changed, (int) Double.parseDouble(txt[9].getText()));
                //System.out.println("btnRun " + phasePanel.getSize());
                //System.out.println("btnRun " + parmsPanel.getSize());
            }
        });

        btnFitLinear.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                fit_linear_response();
                setTitle(" Rossler System - phase y' vs. x' (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + String.format(" %.6f", Main.final_delt) + ")");
            }
        });

        btnPerturbCoeff.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                //fit_linear_uniform();
                Perturb_Torus.calc_coeff();     // calculate perturbation theory coeff
            }
        });

        btnFitCubic.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                fit_cubic_response();
            }
        });

        addWindowListener(new WindowAdapter() {
            @Override public void windowClosing(WindowEvent ev)
                {
                    Main.save_prefs();
                    Main.euler_slider.dispose();
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
                        Main.euler_slider.dispose();
                        Main.plot_y_vs_x.dispose();
                        Main.z_bifurcate = new Rossler_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                }
            });
        //System.out.println("dataPanel = " + dataPanel[0].getSize());
    }

    private void phase_space(boolean ch, int Period)
    {
        double Tnew, Tsum;        // time of peak z
        GregorianCalendar now = new GregorianCalendar();
        PrintWriter fout = null;
        double[] pt6 = new double[] {Main.x0, Main.y0, Main.z0, Main.dx0dc, Main.dy0dc, Main.dz0dc};
        double xmin, xmax, ymin, ymax, zmin, zmax;
        Point2D.Double pt2;                             // projected (x', y') using Euler angles
        int Nloop = N;                                  // number of iterations
        String lblhdr;
        String fmt = "%.3f";

        if (Period == 0)
            delt = DEFAULT_DELT;
        if (ch)
            iT = 0;
        else
            System.arraycopy(pt6_old, 0, pt6, 0, pt6.length);   // re-use previous run
        if (false)
            try
            {
                String fname = "C:\\Windows\\Temp\\Rossler_Output_" + Main.a + "_" + Main.b + "_" + Main.c + "_" + String.format("%.6f", delt) + ".csv";
                boolean fexists = new File(fname).exists();
                FileWriter fw = new FileWriter(fname, true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    fout.println("Rossler x y z vs. i, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + N);
                    fout.println("iT,x,y,z,xdot,ydot,zdot,v_phi,v_theta,x2dot,y2dot,z2dot,kx,ky,kz,v_psi");
                    fout.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Rossler_Output.csv save error = " + e);}

        if (ch || printChk.isSelected() || first)
        {
            if (phaseRadio.isSelected())
                System.out.println("Rossler y' vs. x', " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + Main.project_phi + ", " + Main.project_theta + ", " + delt + ", " + N);
            else if (ddcRadio.isSelected())
                System.out.println("Rossler dy/dc vs. dx/dc, " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + N);
            else
                System.out.println("Rossler dy/du vs. dx/du, " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + N);
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
                //Point2D.Double pt2test = Main.project_2D(pt6[0], pt6[1], pt6[2]); // test code
                //System.out.println("primed 0, " + ", " + pt2test.x + ", " + pt2test.y + ", " + Main.project_zp(pt6[0], pt6[1], pt6[2]));
                //pt6[0] += Main.invert_from_xp_yp(0.0, 0.001, "x");            // test code only
                //pt6[1] += Main.invert_from_xp_yp(0.0, 0.001, "y");            // test code only
                //pt6[2] += Main.invert_from_xp_yp(0.0, 0.001, "z");            // test code only
                //pt6[0] += 0.001*(-0.858330234);
                //pt6[1] += 0.001*0.509240728;
                //pt6[2] += 0.001*0.062794029;
                //pt2test = Main.project_2D(pt6[0], pt6[1], pt6[2]);              // test code
                //System.out.println("primed 1, " + ", " + pt2test.x + ", " + pt2test.y + ", " + Main.project_zp(pt6[0], pt6[1], pt6[2]));
                System.out.println("iT, x, y, z");
                System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                //System.out.println("start = ," + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
            }
            else
            {
                System.out.println("iT, x, y, z, dx/dc, dy/dc, dz/dc");
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
        path1.reset();
        path1.moveTo(xmin, ymin);
        for (int j = 0; j < Nloop; j++)
        {
            iT++;
            if (phaseRadio.isSelected())                        // normal phase space
            {
                Main.runge_kutta_rossler3(pt6, delt, Main.a, Main.b, Main.c);
                pt2 = Main.project_2D(pt6[0], pt6[1], pt6[2]);
                if (pt2.x > xmax) xmax = pt2.x;
                if (pt2.x < xmin) xmin = pt2.x;
                if (pt2.y > ymax) ymax = pt2.y;
                if (pt2.y < ymin) ymin = pt2.y;
                path1.lineTo(pt2.x, pt2.y);
                //if (j <= 1000)
                //    System.out.println(j + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                if (fout != null && true)
                    gen_rdot_r2dot(fout, pt6[0], pt6[1], pt6[2]);
                if (printChk.isSelected() && Period > 0 && j >= Nloop - Period) // transfer last cycle
                    Main.gen_array(j - Nloop + Period, Period, delt, pt6[0], pt6[1], pt6[2]);
                //if (printChk.isSelected() && Period == 0 && j >= Nloop - 288*30)
                //    System.out.println(iT + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                if (j == Nloop - 1)
                {
                    //PrintWriter pout = new PrintWriter(System.out);
                    System.out.print("end = ," + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                    double xdot = -pt6[1] - pt6[2];
                    double ydot =  pt6[0] + Main.a*pt6[1];
                    double zdot =  Main.b + pt6[2]*(pt6[0] - Main.c);
                    System.out.println(", " + Math.atan2(xdot, -ydot)*180/Math.PI
                                     + ", " + Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot))*180/Math.PI);
                    //Point2D.Double pt2test = Main.project_2D(pt6[0], pt6[1], pt6[2]);   // test code
                    //System.out.println("primed 2, " + ", " + pt2test.x + ", " + pt2test.y + ", " + Main.project_zp(pt6[0], pt6[1], pt6[2]));
                    Main.final_Period = Period;
                    Main.final_delt = delt;
                    Main.final_x = pt6[0];
                    Main.final_y = pt6[1];
                    Main.final_z = pt6[2];
                    //gen_rdot_r2dot(fout, pt6[0], pt6[1], pt6[2]);
                    //pout.close();
                }
            }
            else if (ddcRadio.isSelected())                        // ddc tangent phase space
            {
                Main.runge_kutta_rossler6_ddc(pt6, delt, Main.c);
                if (pt6[3] > xmax) xmax = pt6[3];
                if (pt6[3] < xmin) xmin = pt6[3];
                if (pt6[4] > ymax) ymax = pt6[4];
                if (pt6[4] < ymin) ymin = pt6[4];
                if (pt6[5] > zmax) zmax = pt6[5];
                if (pt6[5] < zmin) zmin = pt6[5];
                path1.lineTo(pt6[3], pt6[4]);
                if (printChk.isSelected() && j >= Nloop - 2001)
                    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
            }
            else if (ddu3Radio.isSelected())               // Lyapunov response (3-D)
            {
                Main.runge_kutta_rossler6_ddu3(pt6, delt, Main.c);
                if (pt6[3] > xmax) xmax = pt6[3];
                if (pt6[3] < xmin) xmin = pt6[3];
                if (pt6[4] > ymax) ymax = pt6[4];
                if (pt6[4] < ymin) ymin = pt6[4];
                if (pt6[5] > zmax) zmax = pt6[5];
                if (pt6[5] < zmin) zmin = pt6[5];
                path1.lineTo(pt6[3], pt6[4]);
                if (printChk.isSelected() && j == Period - 1)                   // one-shot response after the first cycle
                //if (printChk.isSelected() && j == N - 1)      // fix fix bug bug test code             // one-shot response after the first cycle
                {
                    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
                    double xdot = -pt6[1] - pt6[2];
                    double ydot =  pt6[0] + Main.a*pt6[1];
                    double zdot =  Main.b + pt6[2]*(pt6[0] - Main.c);
                    double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
                    System.out.println(iT + ", " + xdot/v + ", " + ydot/v + ", " + zdot/v
                                          + ", " + Math.atan2(xdot, -ydot)*180/Math.PI
                                          + ", " + Math.acos(zdot/v)*180/Math.PI);
                    System.out.println("ddu3 range, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + xmin + ", " + xmax + ", " + ymin + ", " + ymax + ", " + zmin + ", " + zmax);
                }
                if (fout != null && true && j <= 30000)
                    fout.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
                //if (printChk.isSelected() && Period > 0 && j >= Nloop - 2*Period) // continuous response over last cycle
                //    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
            }
            else if (ddu2Radio.isSelected())               // Lyapunov response (2-D)
            {
                Main.runge_kutta_rossler6_ddu2(pt6, delt, Main.c);
                if (pt6[3] > xmax) xmax = pt6[3];
                if (pt6[3] < xmin) xmin = pt6[3];
                if (pt6[4] > ymax) ymax = pt6[4];
                if (pt6[4] < ymin) ymin = pt6[4];
                zmax = 0;
                zmin = 0;
                path1.lineTo(pt6[3], pt6[4]);
                //if (printChk.isSelected() && Period > 0 && j >= Nloop - Period) // continuous response over last cycle
                //    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
                if (fout != null && true && j >= (30000 - 1) && j < 60000)
                    fout.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4]);
            }
            if (zlist[0] < zlist[1] && zlist[1] <= zlist[2] && zlist[2] > zlist[3] && zlist[3] > pt6[2])
            {
                Tindex = (Tindex + 1) % Nfork;
                //Tnew = Main.parabola(iT - 3, iT - 2, iT - 1, zlist[1], zlist[2], zlist[3], false);
                Tnew = iT - 2 + Main.quarticT(zlist[0] - zlist[2], zlist[1] - zlist[2], zlist[3] - zlist[2], pt6[2] - zlist[2]);
                //System.out.println("zlist, " + zlist[0] + ", " + zlist[1] + ", " + zlist[2] + ", " + zlist[3] + ", " + pt6[2]);
                Tfork[Tindex] = Tnew - Told;
                //System.out.println((iT - 2) + ", " + Main.c + ", " + zlist[2] + ", " + Tnew
                //                     + ", " + Main.parabola(iT - 3, iT - 2, iT - 1, zlist[1], zlist[2], zlist[3], true)
                //                     + ", " + Tfork[Tindex]);
                Told = Tnew;
                if (Tindex == Nfork - 1 && !printChk.isSelected())
                {
                    Tsum = 0;
                    System.out.print("sum = ," + Main.a + ", " + Main.b + ", " + Main.c + ", " + (iT - 2));
                    for (double tf: Tfork)
                    {
                        System.out.print(", " + tf);
                        Tsum += tf;
                    }
                    System.out.println(", " + Tsum + ", " + Period + ", " + delt); // + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                    if (Period > 0)
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

    protected static void fit_linear_response()
    {
        // collect data to fit first-order response, after one cycle,
        // to a change in the (x, y, z) space, using runge_kutta_rossler6_ddu3

        double[][] M = new double[2][2];
        //double[] xfer = new double[] {0.8475, 0.6, 2.0, 2400, 0.002045685610183913, 2.7792546555929794, -1.2277175829196465, 2.9101772988508485};
        //double[] xfer = new double[] {0.84190, 0.59, 2.0, 2400, 0.0020321339745981207, 2.32059658980927, -0.8485003749147986, 3.515665813849232};
        //double[] xfer = new double[] {0.84192, 0.59, 2.0, 2400, 0.0020320893106168156, 2.84772250612556, -1.2630065015865268, 2.9076520146213745};
        //double[] xfer = new double[] {0.848, 0.6018, 2.0, 2400, 0.0020491121506988103, 3.074569698559594, -1.8633039841465269, 1.8026220073366606};
        //double[] xfer = new double[] {0.8484, 0.6018, 2.0, 2500, 0.0019663861357734525, 2.354907117802819, -0.892886695902259, 3.3797137892723295};
        //double[] xfer = new double[] {0.849, 0.6018, 2.0, 2400, 0.0020471164547495166, 2.741538523820649, -1.197415664183056, 2.951294913869168};
        //double[] xfer = new double[] {0.8493, 0.6018, 2.0, 2400, 0.0020465094608093123, 3.013772875361961, -1.5742013537827195, 2.30657089928223};
        //double[] xfer = new double[] {0.8483, 0.6018, 2.0, 2400, 0.0020485178260750037, 3.0173316830218306, -1.5789207124022282, 2.2953108811103995};
        //double[] xfer = new double[] {0.8492, 0.6018, 2.0, 2500, 0.0019648437369485885, 2.34511116788451, -0.8889898118445804, 3.38661715735611};
        //double[] xfer = new double[] {0.8476, 0.6018, 2.0, 2500, 0.0019679029217163522, 0.09536624520817999, -1.579300031673738, 0.32683653610549074};
        //double[] xfer = new double[] {0.8486, 0.6018, 2.0, 2500, 0.0019660029788507332, 2.352444315329206, -0.8918955133940338, 3.3814619126442818};
        //double[] xfer = new double[] {0.85, 0.6018, 2.0, 2400, 0.0020450777377757036, 2.2702741689085024, -2.4815181841064016, 0.6853905939587694};
        //double[] xfer = new double[] {0.6154, 0.6, 1.25, 2500, 0.001958255038616749, 0.943436061311421, -0.8801316636610237, 1.3939058197338807};
        //double[] xfer = new double[] {0.6152, 0.6, 1.25, 2400, 0.0020404301233040785, 1.1021570263574632, -1.0697358812222264, 1.2582782496952685};
        //double[] xfer = new double[] {0.61535, 0.6, 1.25, 2400, 0.0020399945248164965, 1.0408280966384564, -0.9731256685650206, 1.3390957887692634};
        //double[] xfer = new double[] {0.61455, 0.6, 1.25, 2400, 0.002042301013948591, 1.045805952303036, -0.9724592171925505, 1.33788244951892};
        //double[] xfer = new double[] {0.6142, 0.6, 1.25, 2400, 0.002043297429904907, 1.0479769853341947, -0.9722469347640545, 1.3372799587873525};
        //double[] xfer = new double[] {0.61436, 0.6, 1.25, 2400, 0.002042842863915326, 0.9585504551259202, -0.8846192577991049, 1.3909892432494375};
        //double[] xfer = new double[] {0.615, 0.6, 1.25, 2400, 0.002041008652432798, 1.1053031693677278, -1.0747243349326034, 1.2530723349641284};
        //double[] xfer = new double[] {0.613, 0.6, 1.25, 1000.0, -0.00491197867733379, 1.0803402502099697, -1.0065396774688258, 1.3078939309021518};
        //double[] xfer = new double[] {0.614, 0.6, 1.25, 2500, 0.0019621089010658493, 0.1526369391508889, -0.8648954643212351, 0.8659408980845922};
        //double[] xfer = new double[] {0.6156, 0.6, 1.25, 2400, 0.0020392652462394504, 0.3224917019058715, -1.150200373750082, 0.7301173767526574};
//        //double[] xfer = new double[] {0.613615, 0.6, 1.25, 2500, 0.0019631483769110302, 0.14895150728061218, -0.8615515776624472, 0.8647183170080577};
        //double[] xfer = new double[] {0.613615, 0.6, 1.25, 2400, 0.0020449462259348994, 0.14904317050291274, -0.8569829974422363, 0.8689988111468419};
        //double[] xfer = new double[] {0.61362, 0.6, 1.25, 2400, 0.002044932220316744, 0.15212507446654686, -0.8291541692640159, 0.8971886609340222};
        double[] xfer = new double[] {0.61362, 0.6, 1.25, 2400, 0.0020449322203147774, 0.8253714931854378, -0.7906394525837482, 1.415558702706121};
        //double[] xfer = new double[] {0.613612788, 0.6, 1.25, 2500, 0.0019631483769110302, 0.14895150728061218, -0.8615515776624472, 0.8647183170080577};
        //double[] xfer = new double[] {0.6155, 0.6, 1.25, 2400, 0.00203955745286074, 1.0727371201489548, -1.0179184419750524, 1.3046631553824795};
        //double[] xfer = new double[] {0.6155, 0.6, 1.25, 2400, 0.00203955745286074, 0.6062959213425169, -1.320068706969935, 0.7425192346784949};
        //double[] xfer = new double[] {0.613615, 0.6, 1.25, 2400, 0.0020449462259342285, 1.034450798924672, -0.9516501511131881, 1.3506453448473772};
        //double[] xfer = new double[] {0.6155, 0.6, 1.25, 2400, 0.00203955745286074, 1.0727371201489548, -1.0179184419750524, 1.3046631553824795};
        //double[] xfer = new double[] {0.6155, 0.6, 1.25, 2400*2, 0.00203955745286074/2, 0.9980865296029361, -0.9283243327845676, 1.3690760769937402};
        //double[] xfer = new double[] {0.61465, 0.6, 1.25, 2400, 0.002042014926546927, 1.073065448368517, -1.010296410713409, 1.3088652634252895};
        //double[] xfer = new double[] {0.61445, 0.6, 1.25, 2400, 0.0020425864773399565, 1.0742955273651251, -1.0103094893478959, 1.3083703787600285};
        //double[] xfer = new double[] {0.61415, 0.6, 1.25, 2400, 0.0020434391601832, 0.30280341144849077, -1.140230754847078, 0.7220977493765032};

        Main.a = xfer[0];
        Main.b = xfer[1];
        Main.c = xfer[2];
        Main.final_x = xfer[5];                 // emergency initiallization only
        Main.final_y = xfer[6];
        Main.final_z = xfer[7];
        Main.final_Period = (int) xfer[3];
        Main.final_delt = xfer[4];
        Main.final_delt *= -1;                  // temporarily reverse the direction of time

        Main.skew_transform = false;            // default = false, try to detect "non-uniformity"
        double incr = 1;
        double xdot = -Main.final_y - Main.final_z;
        double ydot =  Main.final_x + Main.a*Main.final_y;
        double zdot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
        Main.project_phi = Math.atan2(xdot, -ydot)*180/Math.PI;
        Main.project_theta = Math.acos(zdot/v)*180/Math.PI;
        double[] pt6 = new double[6];

        System.out.println("\nPython output (fit_linear_response):");
        System.out.println("linear_hdr = \"\\n\\");
        System.out.println("ddu3 - measure linear response after one cycle\\n\\");
        System.out.println("incr         , " + incr + ", \\n\\");
        System.out.println("a_b_c        , " + Main.a + ", " + Main.b + ", " + Main.c + ", \\n\\");
        System.out.println("Period_delt  , " + Main.final_Period + ", " + Main.final_delt + ", \\n\\");
        System.out.println("x_y_z        , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ", \\n\\");
        System.out.println("vx_vy_vz     , " + xdot/v + ", " + ydot/v + ", " + zdot/v + ", \\n\\");
        System.out.println("phi_theta_psi, " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ", \\n\"");

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
            for (int k = 0; k < Main.final_Period; k++)                     // loop through one cycle
                Main.runge_kutta_rossler6_ddu3(pt6, Main.final_delt, Main.c);
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
            double sqr = Math.sqrt((M[0][1] + M[1][0])*(M[0][1] + M[1][0]) + (M[1][1] - M[0][0])*(M[1][1] - M[0][0]));
            System.out.println("ellipse phi_a/b = ," + Math.atan2(M[1][1] - M[0][0], M[0][1] + M[1][0])*180/Math.PI/2 + ", "
                                                     + Math.sqrt((M[0][1] - M[1][0] + sqr)/(M[0][1] - M[1][0] - sqr)));
            System.out.println("Re_V21, Im_V21 = , " + Main.project_psi + ", " + Main.final_Re_V21 + ", " + Main.final_Im_V21);
//            System.out.println("summary, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + eig + ", " + angle*180/Math.PI + ", " + Main.final_Re_V21 + ", " + Main.final_Im_V21 + ", " + Math.atan2(M[1][1] - M[0][0], M[0][1] + M[1][0])*180/Math.PI/2 + ", " + Math.sqrt((M[0][1] - M[1][0] + sqr)/(M[0][1] - M[1][0] - sqr)));
            System.out.println("summary, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + + Main.final_Period + ", " + Main.final_delt + ", " + eig + ", " + angle*180/Math.PI + ", " + eig*Math.cos(angle) + ", " + -eig*Math.sin(angle));
        }
        else
        {
            System.out.println("real eigenvalues, " + (M[0][0] + M[1][1] - Math.sqrt((M[1][1] - M[0][0])*(M[1][1] - M[0][0]) + 4*M[0][1]*M[1][0]))/2
                                             + ", " + (M[0][0] + M[1][1] + Math.sqrt((M[1][1] - M[0][0])*(M[1][1] - M[0][0]) + 4*M[0][1]*M[1][0]))/2);
        }
        //System.out.println("summary, " + M[0][0] + ", " + M[0][1] + ", " + Main.final_Re_V21);
        //System.out.println("summary, " + M[1][0] + ", " + M[1][1] + ", " + Main.final_Im_V21);

        // test transforms (test code, to be deleted)

//        Main.project_phi = 72.1724028032559;
//        Main.project_theta = 77.37490911440288;
//        double testx = 1.73, testy = -2.41, testz = 7.28;
//        Point2D.Double testpt2 = Main.project_2D(testx, testy, testz);
//        double testzp = Main.project_zp(testx, testy, testz);
//        System.out.println("Euler = " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi);
//        System.out.println("orig  = " + testx + ", " + testy + ", " + testz);
//        System.out.println("prime = " + testpt2.x + ", " + testpt2.y + ", " + testzp);
//        System.out.println("back  = " + Main.invert_from_xp_yp(testpt2.x, testpt2.y, testzp, "x") + ", " + Main.invert_from_xp_yp(testpt2.x, testpt2.y, testzp, "y") + ", " + Main.invert_from_xp_yp(testpt2.x, testpt2.y, testzp, "z"));
    }

    protected static void fit_linear_uniform_obsolete()          // OBSOLETE, May 25, 2022 (replaced with class Perturb_Torus)
    {
        // calculate first-order response, during one cycle,
        // assume that the output has already been made uniform, by running 'fit_linear_response()'
        // generate file output of a first-order S matrix

        if (Main.final_Re_V21 == 0 || Main.final_Im_V21 == 0 || Main.project_phi == 0 || Main.project_theta == 0
        ||  Main.final_Period == 0 || eig == 0 || angle == 0 || Main.skew_transform)
        {
            System.out.println("Bad data in 'fit_linear_uniform()'");
            return;
        }
        double[] pt6 = new double[6];
        Point2D.Double pt2;
        Main.skew_transform = true;                 // make the linear response "uniform"
        double xdot = -Main.final_y - Main.final_z;
        double ydot =  Main.final_x + Main.a*Main.final_y;
        double zdot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
        //double deleig = -Math.log(eig)/Main.final_Period;
        //double delangle = (angle + 2*Math.PI)/Main.final_Period;
        //double delx, dely;
        //System.out.println(Math.exp(deleig) + ", " + delangle);
        try
        {
            String fname = "C:\\Windows\\Temp\\Rossler_S_Matrix_" + Main.a + "_" + Main.b + "_" + Main.c + "_" + String.format("%.6f", Main.final_delt) + ".csv";
            FileWriter fw = new FileWriter(fname, false);
            PrintWriter fout = new PrintWriter(fw);
            fout.println("Rossler S matrix, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Main.final_Period + ", " + Main.final_delt);
            fout.println("phi_theta_psi, " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi);
            fout.println("|eig|_angle, " + eig + ", " + angle*180/Math.PI + ", " + eig*Math.cos(angle) + ", " + -eig*Math.sin(angle));
            fout.println("Re_V21_Im_V21, " + Main.final_Re_V21 + ", " + Main.final_Im_V21);
            fout.println("");
            fout.println("i,pt6[0],pt6[1],pt6[2],dxp,dyp,dzp");
            for (int i = 0; i < 3; i++)         // initiallize dx/du, dy/du, or dz/du
            {
                if (i == 0)                                                     // initial (dxdu', dydu', dzdu')
                    pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, Main.invert_from_xp_yp(1, 0, 0, "x"), Main.invert_from_xp_yp(1, 0, 0, "y"), Main.invert_from_xp_yp(1, 0, 0, "z")};
                else if (i == 1)
                    pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, Main.invert_from_xp_yp(0, 1, 0, "x"), Main.invert_from_xp_yp(0, 1, 0, "y"), Main.invert_from_xp_yp(0, 1, 0, "z")};
                else
                    pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, xdot/v, ydot/v, zdot/v};
                fout.println("axis_" + i);
                pt2 = Main.project_2D(pt6[3], pt6[4], pt6[5]);   // final (dxdu', dydu', dzdu')
                fout.println("0" + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt6[3], pt6[4], pt6[5]));
                for (int k = 1; k <= Main.final_Period; k++)                     // loop through one cycle
                {
                    Main.runge_kutta_rossler6_ddu3(pt6, Main.final_delt, Main.c);
                    pt2 = Main.project_2D(pt6[3], pt6[4], pt6[5]);
                    //delx = Math.exp(k*deleig)*( pt2.x*Math.cos(k*delangle) + pt2.y*Math.sin(k*delangle));
                    //dely = Math.exp(k*deleig)*(-pt2.x*Math.sin(k*delangle) + pt2.y*Math.cos(k*delangle));
                    //fout.println(k + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + delx + ", " + dely + ", " + Main.project_zp(pt6[3], pt6[4], pt6[5]));
                    fout.println(k + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt6[3], pt6[4], pt6[5]));
                }
            }
            fout.close();
        }
        catch (java.io.IOException e)
            {System.out.println("Rossler_Output.csv save error = " + e);}
    }

    private static void fit_cubic_response()
    {
        // collect data to fit third-order response, after one cycle,
        // to a change in the (x', y') plane, perpendicular to velocity

        //double[] xfer = new double[] {0.849, 0.6018, 2.0, 1000, 0.004913079491639286, 2.0079106745938855, -0.6989335762557917, 3.5403687512600683};
        //Main.a = xfer[0];
        //Main.b = xfer[1];
        //Main.c = xfer[2];
        //Main.final_x = xfer[5];                 // emergency initiallization only
        //Main.final_y = xfer[6];
        //Main.final_z = xfer[7];
        //Main.final_Period = xfer[3];
        //Main.final_delt = xfer[4];

        Main.skew_transform = true;             // make the linear response "uniform"
        double incr = 0.001;
        double x0dot = -Main.final_y - Main.final_z;
        double y0dot =  Main.final_x + Main.a*Main.final_y;
        double z0dot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        double[] pt3;

        //Main.project_phi = Math.atan2(xdot, -ydot)*180/Math.PI;
        //Main.project_theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot))*180/Math.PI;
        Point2D.Double pt2 = Main.project_2D(Main.final_x, Main.final_y, Main.final_z); // initial (x', y')
        double zp = Main.project_zp(Main.final_x, Main.final_y, Main.final_z);          // setpoint for z'
        double xold = 0, yold = 0, zold = zp;
        Point2D.Double pt2old;          // projected, old, (x', y') using Euler angles
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        double zpold, zpnew;            // projected z'
        double xc, yc;                  // interpolated, projected (x', y')
        //System.out.println("pt2 = " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ", " + pt2);

        System.out.println("\nPython output (fit_cubic_response):");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Neimark-Sacker - measure cubic model response after one cycle\\n\\");
        System.out.println("incr        , " + incr + ", " + "\\n\\");
        System.out.println("a_b_c       , " + Main.a + ", " + Main.b + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt , " + Main.final_Period + ", " + Main.final_delt + ",\\n\\");
        System.out.println("x_y_z       , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta   , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("x'_y'_z'    , " + pt2.x + ", " + pt2.y + ", " + zp + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',tc\"");
        System.out.print("data = np.array([");
//        for (int i = -2; i < 3; i++)        // increment x' by i*incr
//            for (int j = -2; j < 3; j++)    // increment y' by j*incr
        for (int i = 0; i < 3; i++)         // increment r' radius by i*incr
            for (int j = 0; j < 12; j++)    // increment angle by 30 degrees
                if (i > 0 || j == 0)
            {
                pt3 = new double[] {Main.final_x + Main.invert_from_xp_yp(incr*i*Math.cos(j*Math.PI/6), incr*i*Math.sin(j*Math.PI/6), 0, "x"),
                                    Main.final_y + Main.invert_from_xp_yp(incr*i*Math.cos(j*Math.PI/6), incr*i*Math.sin(j*Math.PI/6), 0, "y"),
                                    Main.final_z + Main.invert_from_xp_yp(incr*i*Math.cos(j*Math.PI/6), incr*i*Math.sin(j*Math.PI/6), 0, "z")};
                //System.out.println("org   , " + i + ", " + j + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                // loop through one cycle
                for (int k = 0; k < Main.final_Period + 1*40; k++) // add 10 iterations just to be sure it crosses
                {
                    Main.runge_kutta_rossler3(pt3, Main.final_delt, Main.a, Main.b, Main.c);
                    //if (k < 2 || k > Main.final_Period - 3)
                    //    System.out.println("---   , " + (k + 1) + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                    if (true && k > 1 && (Main.project_zp(xold, yold, zold) - zp)*Main.final_delt <= 0 && (Main.project_zp(pt3[0], pt3[1], pt3[2]) - zp)*Main.final_delt > 0)
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
                        if (i < 2 || j < 11)
                            System.out.println(",");
                    }
                    xold = pt3[0];
                    yold = pt3[1];
                    zold = pt3[2];
                }
                if (false)                          // use last point (synchronize in time)
                {
                    pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                    zpnew = Main.project_zp(pt3[0], pt3[1], pt3[2]);
                    System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");
                    if (i < 2 || j < 11)
                        System.out.println(",");
                }
                if (false)                           // calculate theoretical z_sync (Chaos IV, p. 31)
                {
                    double xdot = -pt3[1] - pt3[2];
                    double ydot =  pt3[0] + Main.a*pt3[1];
                    double zdot =  Main.b + pt3[2]*(pt3[0] - Main.c);
                    double x2dot = -ydot - zdot;
                    double y2dot =  xdot + Main.a*ydot;
                    double z2dot =  zdot*(pt3[0] - Main.c) + pt3[2]*xdot;
                    double A = (x0dot*x2dot + y0dot*y2dot + z0dot*z2dot)/2;
                    double B = x0dot*xdot + y0dot*ydot + z0dot*zdot;
                    double C = x0dot*(pt3[0] - Main.final_x) + y0dot*(pt3[1] - Main.final_y) + z0dot*(pt3[2] - Main.final_z);
                    double alpha = (-B + Math.sqrt(B*B - 4*A*C))/2/A;
                    pt2new = Main.project_2D(pt3[0] + alpha*xdot + alpha*alpha*x2dot/2, pt3[1] + alpha*ydot + alpha*alpha*y2dot/2, pt3[2] + alpha*zdot + alpha*alpha*z2dot/2);
                    zpnew = Main.project_zp(pt3[0] + alpha*xdot + alpha*alpha*x2dot/2, pt3[1] + alpha*ydot + alpha*alpha*y2dot/2, pt3[2] + alpha*zdot + alpha*alpha*z2dot/2);
                    //System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + (Main.final_Period + alpha/Main.final_delt) + "]");
                    System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");
                    if (i < 2 || j < 11)
                        System.out.println(",");
                }
            }
        System.out.println("])");
    }

    private void gen_rdot_r2dot(PrintWriter fout, double x, double y, double z)
    {
        // generate rdot and r2dot, velocity and acceleration
        // calculate curvature vector as a cross-product
        // generate Euler angles of velocity and curvature (phi, theta)

        fout.print(iT + ", " + x + ", " + y + ", " + z);
        double xdot = -y - z;
        double ydot =  x + Main.a*y;
        double zdot =  Main.b + z*(x - Main.c);
        fout.print(", " + xdot + ", " + ydot + ", " + zdot);
        double phi = Math.atan2(xdot, -ydot);
        double theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));
        fout.print(", " + phi*180/Math.PI + ", " + theta*180/Math.PI);

        double x2dot = -ydot - zdot;
        double y2dot =  xdot + Main.a*ydot;
        double z2dot =  zdot*(x - Main.c) + z*xdot;
        fout.print(", " + x2dot + ", " + y2dot + ", " + z2dot);
        double kx = ydot*z2dot - zdot*y2dot;
        double ky = zdot*x2dot - xdot*z2dot;
        double kz = xdot*y2dot - ydot*x2dot;
        fout.print(", " + kx + ", " + ky + ", " + kz);
        //fout.print(", " + Math.asin(kz/Math.sqrt(kx*kx + ky*ky + kz*kz)/Math.sin(theta))*180/Math.PI);
        //fout.print(", " + Math.asin((-kx*Math.sin(phi) + ky*Math.cos(phi))/Math.sqrt(kx*kx + ky*ky + kz*kz)/Math.cos(theta))*180/Math.PI);
        double psi = Math.acos((kx*Math.cos(phi) + ky*Math.sin(phi))/Math.sqrt(kx*kx + ky*ky + kz*kz));
        //Main.project_phi = phi*180/Math.PI;
        //Main.project_theta = theta*180/Math.PI;
        //Main.project_psi = psi*180/Math.PI;
        //Point2D.Double pt2 = Main.project_2D(x, y, z);
        //fout.print(", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(x, y, z));
        fout.println(", " + psi*180/Math.PI);
        //fout.println(", " + Math.atan2(kx, -ky)*180/Math.PI + ", " + Math.acos(kz/Math.sqrt(kx*kx + ky*ky + kz*kz))*180/Math.PI);

        // test code
/*
        Main.project_phi = phi*180/Math.PI;
        Main.project_theta = theta*180/Math.PI;
        Main.project_psi = psi*180/Math.PI;
        Main.zc = Main.project_zp(x, y, z);
        Point2D.Double pt2 = Main.project_2D(x, y, z);
        Main.xc = pt2.x;
        Main.yc = pt2.y;
        fout.println("x'_y'_z', " + Main.xc + "," + Main.yc + "," + Main.zc);
        //pt2 = Main.project_2D(0, 1, 0, psi);
        //fout.println(pt2.x + "," + pt2.y + "," + Main.project_zp(0, 1, 0));
        //pt2 = Main.project_2D(0, 0, 1, psi);
        //fout.println(pt2.x + "," + pt2.y + "," + Main.project_zp(0, 0, 1));
*/
}

    private double dtdc(double[] pt6)          // two-dimensional (dx, dy)
    {
        double xdot = -pt6[1] - pt6[2];
        double ydot = pt6[0] + Main.a*pt6[1];
        return -(pt6[3]*xdot + pt6[4]*ydot)/(xdot*xdot + ydot*ydot);
    }

    private double dtdc3(double[] pt6)          // three-dimensional (dx, dy, dz)
    {
        double xdot = -pt6[1] - pt6[2];
        double ydot = pt6[0] + Main.a*pt6[1];
        double zdot = Main.b + pt6[2]*(pt6[0] - Main.c);
        return -(pt6[3]*xdot + pt6[4]*ydot + pt6[5]*zdot)/(xdot*xdot + ydot*ydot + zdot*zdot);
    }
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
        g2.setPaint(new Color(255, 127, 39));
        g2.draw(Main.plot_y_vs_x.path1);
        //g2.setPaint(new Color(0, 0, 0));
        //g2.setStroke(new BasicStroke(2));
        //g2.draw(Main.plot_y_vs_x.path2);
        g2.setPaint(Color.BLUE);
        g2.draw(Main.plot_y_vs_x.xaxis);
        g2.draw(Main.plot_y_vs_x.yaxis);
    }
}
