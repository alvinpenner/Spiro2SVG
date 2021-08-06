
package rossler;

// plot y versus x
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import javax.swing.*;
import java.util.GregorianCalendar;
import java.io.PrintWriter;
import java.io.FileWriter;

public class Rossler_y_vs_x extends JDialog
{
    private boolean first = true;
    private final static double DEFAULT_DELT = 0.02; // 0.02;
    private final static int N = 500000; //50000;            // total # of iterations 160000
    private static double[] pt6_old;
    protected Path2D.Double path1 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    //protected Path2D.Double path2 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Line2D.Double xaxis = new Line2D.Double(0, 0, 0, 0);
    protected Line2D.Double yaxis = new Line2D.Double(0, 0, 0, 0);
    private static JCheckBox printChk = new JCheckBox("  print  ");
    private static JCheckBox ddyChk = new JCheckBox("d/dc of phase");
    private static JCheckBox delddyChk = new JCheckBox("del d/dc");
    private static JLabel lblxrange;
    private static JLabel lblyrange;
    private static JLabel lblzrange;
    private static JPanel phasePanel = new Plot_Phase_Panel();
    private static PrintWriter out = null;

    int iT = 0;                                     // # iterations
    double[] zlist = new double[4];                 // previous z values
    double Told = 0;
    int Nfork = 1;                                  // number of bifurcated branches
    double[] Tfork = new double[Nfork];             // period per branch
    int Tindex = 0;                                 // number of peaks
    double delt = DEFAULT_DELT;

    public Rossler_y_vs_x(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        Main.type = "phase";
        if (delddyChk.isSelected())
            ddyChk.setSelected(false);
        if (!ddyChk.isSelected() && !delddyChk.isSelected())    // normal phase space
            setTitle(" Rossler System - phase y vs. x");
        else if (ddyChk.isSelected())                           // ddy phase space
            setTitle(" Rossler System - dy/dc vs. dx/dc");
        else
            setTitle(" Rossler System - del dy/dc vs. del dx/dc");  // compensated ddy
        setIconImage(img);
        setSize(610 + 100, 493 + 100);
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
            dataPanel[i].setPreferredSize(new Dimension(130, 24));
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

        ddyChk.setOpaque(false);
        delddyChk.setOpaque(false);

        JPanel xrangePanel = new JPanel();
        xrangePanel.setPreferredSize(new Dimension(130, 24));
        xrangePanel.setOpaque(false);
        lblxrange = new JLabel("x = ");
        lblxrange.setPreferredSize(new Dimension(110, 18));
        xrangePanel.add(lblxrange);

        JPanel yrangePanel = new JPanel();
        yrangePanel.setPreferredSize(new Dimension(130, 24));
        yrangePanel.setOpaque(false);
        lblyrange = new JLabel("y = ");
        lblyrange.setPreferredSize(new Dimension(110, 18));
        yrangePanel.add(lblyrange);

        JPanel zrangePanel = new JPanel();
        zrangePanel.setPreferredSize(new Dimension(130, 24));
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
        parmsPanel.add(ddyChk);
        parmsPanel.add(delddyChk);
        parmsPanel.add(spacerPanel[2]);
        parmsPanel.add(xrangePanel);
        parmsPanel.add(yrangePanel);
        parmsPanel.add(zrangePanel);
        parmsPanel.add(spacerPanel[3]);
        parmsPanel.add(dataPanel[9]);

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
                if (delddyChk.isSelected())
                    ddyChk.setSelected(false);
                if (!ddyChk.isSelected() && !delddyChk.isSelected())    // normal phase space
                    setTitle(" Rossler System - phase y vs. x");
                else if (ddyChk.isSelected())                           // ddy phase space
                    setTitle(" Rossler System - dy/dc vs. dx/dc");
                else
                    setTitle(" Rossler System - del dy/dc vs. del dx/dc");  // compensated ddy
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
                if (Main.dx0dc != Double.parseDouble(txt[6].getText())) changed = true;
                Main.dx0dc = Double.parseDouble(txt[6].getText());
                if (Main.dy0dc != Double.parseDouble(txt[7].getText())) changed = true;
                Main.dy0dc = Double.parseDouble(txt[7].getText());
                if (Main.dz0dc != Double.parseDouble(txt[8].getText())) changed = true;
                Main.dz0dc = Double.parseDouble(txt[8].getText());
                phase_space(changed, (int) Double.parseDouble(txt[9].getText()));
                //System.out.println("btnRun " + phasePanel.getSize());
                //System.out.println("btnRun " + parmsPanel.getSize());
            }
        });

        addWindowListener(new WindowAdapter() {
            @Override public void windowClosing(WindowEvent ev)
                {
                    Main.save_prefs();
                    if (out != null)
                        out.close();
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
                        Main.plot_y_vs_x.dispose();
                        Main.z_bifurcate = new Rossler_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                }
            });

        if (false)
            try
            {
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\Rossler_Output_" + Main.c + ".csv", false);
                out = new PrintWriter(fw);
                out.println("iT, x, y, z");
            }
            catch (java.io.IOException e)
                {System.out.println("Rossler_Output.csv save error = " + e);}
    }

    private void phase_space(boolean ch, int Period)
    {
        double Tnew, Tsum;        // time of peak z
        GregorianCalendar now = new GregorianCalendar();
        double[] pt6 = new double[] {Main.x0, Main.y0, Main.z0, Main.dx0dc, Main.dy0dc, Main.dz0dc};
        double xmin, xmax, ymin, ymax, zmin, zmax;
        double delx, dely;                              // only for compensated ddy runs
        int Nloop = N;                                  // number of iterations
        String lblhdr;
        String fmt = "%.4f";

        if (out != null)
            out.println("Rossler y vs. x, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + delt + ", " + N);
        if (Period == 0)
            delt = DEFAULT_DELT;
        //else
        //    Nloop = 1000*Period;
        if (ch)
            iT = 0;
        else
            System.arraycopy(pt6_old, 0, pt6, 0, pt6.length);   // re-use previous run
        if (ch || printChk.isSelected() || first)
        {
            if (!ddyChk.isSelected() && !delddyChk.isSelected())
                System.out.println("Rossler y vs. x, " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + N);
            else
                System.out.println("Rossler dy/dc vs. dx/dc, " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Period + ", " + delt + ", " + N);
            //System.out.println("a,b,c = " + );
            double zt = (Main.c + Math.sqrt(Main.c*Main.c - 4*Main.a*Main.b))/2/Main.a;
            System.out.println("P1 = (" + Main.a*zt + ", " + (-zt) + ", " + zt + ")");
            zt = (Main.c - Math.sqrt(Main.c*Main.c - 4*Main.a*Main.b))/2/Main.a;
            System.out.println("P2 = (" + Main.a*zt + ", " + (-zt) + ", " + zt + ")");
        }
        first = false;
        if (printChk.isSelected())
            if (!ddyChk.isSelected() && !delddyChk.isSelected())
            {
                System.out.println("iT, x, y, z");
                System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
            }
            else
            {
                System.out.println("iT, x, y, z, dx/dc, dy/dc, dz/dc, dt/dc");
                //System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);
            }
        if (!ddyChk.isSelected() && !delddyChk.isSelected())    // normal phase space
        {
            xmin = pt6[0];
            xmax = pt6[0];
            ymin = pt6[1];
            ymax = pt6[1];
            zmin = pt6[2];
            zmax = pt6[2];
        }
        else if (ddyChk.isSelected())                           // ddy phase space
        {
            xmin = pt6[3];
            xmax = pt6[3];
            ymin = pt6[4];
            ymax = pt6[4];
            zmin = pt6[5];
            zmax = pt6[5];
        }
        else                                                    // compensated ddy
        {
            xmin = pt6[3] + (-pt6[1] - pt6[2])*dtdc(pt6);
            xmax = xmin;
            ymin = pt6[4] + (pt6[0] + Main.a*pt6[1])*dtdc(pt6);
            ymax = ymin;
            zmin = 0;
            zmax = 0;
        }
        path1.reset();
        path1.moveTo(xmin, ymin);
        for (int j = 0; j < Nloop; j++)
        {
            iT++;
            if (!ddyChk.isSelected() && !delddyChk.isSelected())    // normal phase space
            {
                Main.runge_kutta_rossler3(pt6, delt, Main.a, Main.c);
                if (pt6[0] > xmax) xmax = pt6[0];
                if (pt6[0] < xmin) xmin = pt6[0];
                if (pt6[1] > ymax) ymax = pt6[1];
                if (pt6[1] < ymin) ymin = pt6[1];
                if (pt6[2] > zmax) zmax = pt6[2];
                if (pt6[2] < zmin) zmin = pt6[2];
                path1.lineTo(pt6[0], pt6[1]);
                //if (j <= 5760)
                //    System.out.println(j + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                if (out != null)
                    out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
                if (printChk.isSelected() && Period > 0 && j >= Nloop - Period)
                    Main.gen_array(j - Nloop + Period, Period, delt, pt6[0], pt6[1], pt6[2]);
                //if (printChk.isSelected() && Period == 0 && j >= Nloop - 288*30)
                //    System.out.println(iT + ", " + Period + ", " + delt + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2]);
            }
            else if (ddyChk.isSelected())                           // ddy phase space
            {
                Main.runge_kutta_rossler6(pt6, delt, Main.c);
                if (pt6[3] > xmax) xmax = pt6[3];
                if (pt6[3] < xmin) xmin = pt6[3];
                if (pt6[4] > ymax) ymax = pt6[4];
                if (pt6[4] < ymin) ymin = pt6[4];
                if (pt6[5] > zmax) zmax = pt6[5];
                if (pt6[5] < zmin) zmin = pt6[5];
                path1.lineTo(pt6[3], pt6[4]);
                if (printChk.isSelected() && j >= Nloop - 2001)
                    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5] + ", " + dtdc(pt6));
            }
            else
            {
                Main.runge_kutta_rossler6(pt6, delt, Main.c);       // compensated ddy
                delx = pt6[3] + (-pt6[1] - pt6[2])*dtdc(pt6);
                dely = pt6[4] + (pt6[0] + Main.a*pt6[1])*dtdc(pt6);
                if (delx > xmax) xmax = delx;
                if (delx < xmin) xmin = delx;
                if (dely > ymax) ymax = dely;
                if (dely < ymin) ymin = dely;
                path1.lineTo(delx, dely);
                if (printChk.isSelected() && j >= Nloop - 2001)
                    System.out.println(iT + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5] + ", " + dtdc(pt6));
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
                    System.out.print("sum = ," + Main.c + ", "+ (iT - 2));
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
        System.arraycopy(pt6, 0, pt6_old, 0, pt6.length);               // save
        AffineTransform at = new AffineTransform((phasePanel.getWidth() - 8)/(xmax - xmin),
                                          0, 0, -(phasePanel.getHeight() - 8)/(ymax - ymin),
                                                (xmax*4 + xmin*(4 - phasePanel.getWidth()))/(xmax - xmin),
                                               -(ymin*4 + ymax*(4 - phasePanel.getHeight()))/(ymax - ymin));
        path1.transform(at);
        //System.out.println("range, " + Main.c + ", " + xmin + ", " + xmax + ", " + ymin + ", " + ymax + ", " + zmin + ", " + zmax);
        lblhdr = !ddyChk.isSelected() ? "x" : "dx";
        lblxrange.setText(lblhdr + " = " + String.format(fmt, xmin) + ", " + String.format(fmt, xmax));
        lblhdr = !ddyChk.isSelected() ? "y" : "dy";
        lblyrange.setText(lblhdr + " = " + String.format(fmt, ymin) + ", " + String.format(fmt, ymax));
        lblhdr = !ddyChk.isSelected() ? "z" : "dz";
        lblzrange.setText(lblhdr + " = " + String.format(fmt, zmin) + ", " + String.format(fmt, zmax));
        xaxis = new Line2D.Double(0, ymax*phasePanel.getHeight()/(ymax - ymin), phasePanel.getWidth(), ymax*phasePanel.getHeight()/(ymax - ymin));
        yaxis = new Line2D.Double(-xmin*phasePanel.getWidth()/(xmax - xmin), 0, -xmin*phasePanel.getWidth()/(xmax - xmin), phasePanel.getHeight());
        phasePanel.repaint();
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
