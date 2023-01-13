
package chua;

// bifurcation plot of Tau versus c
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.List;
import javax.swing.*;

public class Chua_z_bifurcate extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(600, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static double save_alpha, save_beta, save_gamma, save_c;        // backup and restore from other dialogs
    protected static JButton btnRun = new JButton("Run");
    protected static JButton btnClear = new JButton("Clr");
    protected static double scan_start, scan_end;
    //protected static double Tmin = 1.8, Tmax = 2.8;
    protected static double Tmin = 0.1675, Tmax = 0.1681;

    public Chua_z_bifurcate(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        final JPanel bifurcatePanel = new JPanel();
        Main.type = "bifurcate";
        setTitle(" Chua Oscillator - Tau bifurcate");
        setIconImage(img);
        setSize(770, 540);
        setLocationByPlatform(true);
        save_alpha = Main.alpha;
        save_beta = Main.beta;
        save_gamma = Main.gamma;
        save_c = Main.c;

        final JLabel[] lbl = {new JLabel("alpha s"),
                              new JLabel("alpha e"),
                              new JLabel("beta s"),
                              new JLabel("beta e"),
                              new JLabel("gamma s"),
                              new JLabel("gamma e"),
                              new JLabel("c s"),
                              new JLabel("c e"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.alpha_s)),
                                  new JTextField(Double.toString(Main.alpha_e)),
                                  new JTextField(Double.toString(Main.beta_s)),
                                  new JTextField(Double.toString(Main.beta_e)),
                                  new JTextField(Double.toString(Main.gamma_s)),
                                  new JTextField(Double.toString(Main.gamma_e)),
                                  new JTextField(Double.toString(Main.c_s)),
                                  new JTextField(Double.toString(Main.c_e)),
                                  new JTextField(Double.toString(Main.x0)),
                                  new JTextField(Double.toString(Main.y0)),
                                  new JTextField(Double.toString(Main.z0))};
        JPanel[] spacerPanel = new JPanel[2];
        JPanel[] dataPanel = new JPanel[lbl.length];

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
            dataPanel[i].setOpaque(false);
            lbl[i].setPreferredSize(new Dimension(55, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

        JPanel rangePanel = new JPanel();
        rangePanel.setOpaque(false);
        JLabel lblrange = new JLabel("T = " + Tmin + " - " + Tmax);
        lblrange.setPreferredSize(new Dimension(110, 20));
        rangePanel.add(lblrange);

        JPanel posnPanel = new JPanel();
        posnPanel.setOpaque(false);
        final JLabel lblposn = new JLabel();
        lblposn.setBorder(BorderFactory.createEtchedBorder());
        lblposn.setPreferredSize(new Dimension(110, 20));
        posnPanel.add(lblposn);

        parmsPanel.removeAll();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(dataPanel[0]);
        parmsPanel.add(dataPanel[1]);
        parmsPanel.add(dataPanel[2]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(dataPanel[7]);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[8]);
        parmsPanel.add(dataPanel[9]);
        parmsPanel.add(dataPanel[10]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(btnClear);
        parmsPanel.add(rangePanel);
        parmsPanel.add(posnPanel);
        parmsPanel.setMaximumSize(new Dimension(140, 3000));
        parmsPanel.setPreferredSize(new Dimension(140, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        int ulx = 20 + 450;
        int uly = 350;
        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        DC.setFont(new Font( "SansSerif", Font.BOLD, 12 ));
        DC.setColor(new Color(40, 40, 136));
        DC.drawLine(ulx, uly - 20, ulx + 10, uly - 20);
        DC.drawString("increasing", ulx + 20, uly - 15);
        DC.setColor(new Color(255, 128, 0));
        DC.drawLine(ulx, uly + 0, ulx + 10, uly + 0);
        DC.drawString("decreasing", ulx + 20, uly + 5);
//        DC.setColor(new Color(0, 0, 0));
//        DC.drawLine(ulx, uly + 20, ulx + 10, uly + 20);
//        DC.drawString("cubic model", ulx + 20, uly + 25);
        DC.setColor(Color.white);

        lblImage.setBorder(BorderFactory.createEtchedBorder());
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                if (scan_end > scan_start)
                    lblposn.setText(String.format(" %.3f", scan_start + e.getX()*(scan_end - scan_start)/image.getWidth()) + ", " + String.format("%.4f", Tmax + e.getY()*(Tmin - Tmax)/image.getHeight()));
                else
                    lblposn.setText(String.format(" %.3f", scan_end + e.getX()*(scan_start - scan_end)/image.getWidth()) + ", " + String.format("%.4f", Tmax + e.getY()*(Tmin - Tmax)/image.getHeight()));
            }
        });
        bifurcatePanel.add(lblImage);

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(parmsPanel, BorderLayout.WEST);
        getContentPane().add(bifurcatePanel, BorderLayout.EAST);
        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                Main.alpha_s = Double.parseDouble(txt[0].getText());
                Main.alpha_e = Double.parseDouble(txt[1].getText());
                Main.beta_s  = Double.parseDouble(txt[2].getText());
                Main.beta_e  = Double.parseDouble(txt[3].getText());
                Main.gamma_s = Double.parseDouble(txt[4].getText());
                Main.gamma_e = Double.parseDouble(txt[5].getText());
                Main.c_s = Double.parseDouble(txt[6].getText());
                Main.c_e = Double.parseDouble(txt[7].getText());
                Main.x0 = Double.parseDouble(txt[8].getText());
                Main.y0 = Double.parseDouble(txt[9].getText());
                Main.z0 = Double.parseDouble(txt[10].getText());
                Main.alpha = Main.alpha_s;                          // temporary for scanning only
                Main.beta = Main.beta_s;                            // temporary for scanning only
                Main.gamma = Main.gamma_s;                          // temporary for scanning only
                Main.c = Main.c_s;                                  // temporary for scanning only
                if (!Double.isNaN(Main.alpha_e))
                {
                    setTitle(" Chua System - Tau bifurcate vs. alpha (" + Main.a + ", " + BifurcateActivity.delt + ")");
                    scan_start = Main.alpha_s;
                    scan_end = Main.alpha_e;
                }
                else if (!Double.isNaN(Main.beta_e))
                {
                    setTitle(" Chua System - Tau bifurcate vs. beta (" + Main.a + ", " + BifurcateActivity.delt + ")");
                    scan_start = Main.beta_s;
                    scan_end = Main.beta_e;
                }
                else if (!Double.isNaN(Main.gamma_e))
                {
                    setTitle(" Chua System - Tau bifurcate vs. gamma (" + Main.a + ", " + BifurcateActivity.delt + ")");
                    scan_start = Main.gamma_s;
                    scan_end = Main.gamma_e;
                }
                else if (!Double.isNaN(Main.c_e))
                {
                    setTitle(" Chua System - Tau bifurcate vs. c (" + Main.a + ", " + BifurcateActivity.delt + ")");
                    scan_start = Main.c_s;
                    scan_end = Main.c_e;
                }
                else
                {
                    setTitle(" Chua System - Bad Data, Do Not Use!");
                    System.out.println("Chua System - Bad Data, Do Not Use!");
                }
                btnRun.setEnabled(false);
                BifurcateActivity activity = new BifurcateActivity();
                activity.execute();
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

        addWindowListener(new WindowAdapter() {
            @Override public void windowClosing(WindowEvent ev)
                {
                    Main.alpha = save_alpha;                            // restore original
                    Main.beta = save_beta;
                    Main.gamma = save_gamma;
                    Main.c = save_c;
                    Main.save_prefs();
                }
            });

        for (JTextField tx: txt)
            tx.addKeyListener(new KeyAdapter()
            {
                @Override public void keyPressed(KeyEvent e)
                {
                    if (e.getKeyCode() == KeyEvent.VK_ESCAPE)
                    {
                        Main.alpha = save_alpha;                        // restore original
                        Main.beta = save_beta;
                        Main.gamma = save_gamma;
                        Main.c = save_c;
                        Main.save_prefs();
                        Main.z_bifurcate.dispose();
                        Main.x_y_scatter = new Chua_x_y_scatter(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                        Main.chua_euler_slider = new Chua_Euler_Slider(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                    if (e.getKeyCode() == KeyEvent.VK_F1)       // save a PNG file using the F1 key
                    {
                        if (Main.z_bifurcate != null && Main.z_bifurcate.isShowing())
                            Main.save_PNG(image, String.format("Chua_bifurcate_%.2f_%.2f_%.1f_%.3f", scan_start, scan_end, Main.a, Main.c));
                        return;
                    }
                }
            });
    }
}

class BifurcateActivity extends SwingWorker<Void, Point>
{
    protected static double delt = 0.0005;                        // normally negative
    private static Color clr;
    private static double[] pt3;

    public BifurcateActivity()
    {
        pt3 = new double[] {Main.x0, Main.y0, Main.z0};
        if (Chua_z_bifurcate.scan_end > Chua_z_bifurcate.scan_start)
            clr = new Color(40, 40, 136);
        else
            clr = new Color(255, 128, 0);
    }

    protected Void doInBackground() throws Exception
    {
        // use peaks in y, pt3[1], to detect period
        int N = 400000;                             // # of iterations per c
        int Ninit = 0*1152;                      // # approx 10 Tx cycles to initiallize limit cycle
        double scan;                                // horizontal axis
        double[] zlist = new double[4];             // previous z values
        double Told = 0, Tnew = 0;                  // time of peak z
        double Tdiff = Tnew - Told;
        boolean first;                              // first peak at c

        for (int i = 0; i < 1000; i++)              // pre-initialize
            Main.runge_kutta_chua3(pt3, delt);
        //System.out.println("init 1 = " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        for (int i = 0; i < Chua_z_bifurcate.image.getWidth(); i++)
        {
            first = true;
            scan = Chua_z_bifurcate.scan_start + i*(Chua_z_bifurcate.scan_end - Chua_z_bifurcate.scan_start)/Chua_z_bifurcate.image.getWidth();
            for (int j = 0; j < N; j++)
            {
                if (!Double.isNaN(Main.alpha_e))
                    Main.alpha = scan;
                else if(!Double.isNaN(Main.beta_e))
                    Main.beta = scan;
                else if(!Double.isNaN(Main.gamma_e))
                    Main.gamma = scan;
                else
                    Main.c = scan;
                Main.runge_kutta_chua3(pt3, delt);
                if (i == 300 && j == 10)
                    System.out.println(i + ", " + j + ", " + scan + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                if (zlist[0] < zlist[1] && zlist[1] <= zlist[2] && zlist[2] > zlist[3] && zlist[3] > pt3[1])
                {
                    Tnew = j - 2 + Main.quarticT(zlist[0] - zlist[2], zlist[1] - zlist[2], zlist[3] - zlist[2], pt3[1] - zlist[2]);
                    //System.out.println("zlist, " + i + ", " + scan + ", " + (j - 2) + ", " + Tnew + ", " + zlist[0] + ", " + zlist[1] + ", " + zlist[2] + ", " + zlist[3] + ", " + pt3[1]);
                    //System.out.println(i + ", " + j + ", " + Tdiff + ", " + Told + ", " + Tnew + ", " + zlist[2]);
                    //if (Tnew - Told > 0.0 || i < 10)              // temporary code fix fix
                    //if (Tnew - Told < 0.0 || Tnew - Told > 0.6*Tdiff)
                    //if (Tnew - Told < 0.0 || zlist[2] > 50)
                    if (first)
                    {
                        Told = Tnew;
                        Tdiff = 0;
                        first = false;
                    }
                    else //if (zlist[2] > zold/3 || c < 40) //if (Tnew - Told > 0.7*Tdiff)
                    {
                        Tdiff = Tnew - Told;
                        if (delt*Tdiff > Chua_z_bifurcate.Tmin && delt*Tdiff < Chua_z_bifurcate.Tmax && j > Ninit)
                            //if ((i > 20 && Rossler_z_bifurcate.scan_end > Rossler_z_bifurcate.scan_start)
                            //||  (i < 640 && Rossler_z_bifurcate.scan_end < Rossler_z_bifurcate.scan_start))
                            publish(new Point(i, (int) (Chua_z_bifurcate.image.getHeight()*(Chua_z_bifurcate.Tmax - delt*(Tnew - Told))/(Chua_z_bifurcate.Tmax - Chua_z_bifurcate.Tmin))));
                        //if (c == 3)
                        //    System.out.println(c + ", " + (j - 2) + ", " + Tnew + ", " + (Tnew - Told));
                        Told = Tnew;
                    }
                }
                for (int k = 0; k < 3; k++)
                    zlist[k] = zlist[k + 1];
                zlist[3] = pt3[1];
            }
            Chua_z_bifurcate.lblImage.repaint();
        }

        // paint a theoretical envelope (normally disabled)
/*
        clr = new Color(0, 0, 0);
        double av = -4.895;
        double start = (0.61563 - Chua_z_bifurcate.scan_start)/(Chua_z_bifurcate.scan_end - Chua_z_bifurcate.scan_start)*Chua_z_bifurcate.image.getWidth();
        double end = (0.613613 - Chua_z_bifurcate.scan_start)/(Chua_z_bifurcate.scan_end - Chua_z_bifurcate.scan_start)*Chua_z_bifurcate.image.getWidth();
        double b = 4.92;
        for (int i = 0; i < Chua_z_bifurcate.image.getWidth(); i++)
            if (i >= start && i <= end)
            {
                double alpha = 0.02043*(i - start)/(end - start);
                double theo = 0.4*Math.sqrt(1 + alpha - Math.sqrt(1 + b*b - b*b*(1 + alpha)*(1 + alpha)));
                //System.out.println(i + ", " + alpha + ", " + (- av + theo));
                if (- av + theo < -Chua_z_bifurcate.Tmin)
                    publish(new Point(i, (int) (Chua_z_bifurcate.image.getHeight()*(Chua_z_bifurcate.Tmax - av + theo)/(Chua_z_bifurcate.Tmax - Chua_z_bifurcate.Tmin))));
                if (- av - theo > -Chua_z_bifurcate.Tmax)
                    publish(new Point(i, (int) (Chua_z_bifurcate.image.getHeight()*(Chua_z_bifurcate.Tmax - av - theo)/(Chua_z_bifurcate.Tmax - Chua_z_bifurcate.Tmin))));
            }
*/
        return null;
    }

    @Override protected void process(List<Point> points)
    {
        for (Point listpt : points)
        {
            //System.out.println(listpt.x + ", " + listpt.y);
            if (Chua_z_bifurcate.scan_end > Chua_z_bifurcate.scan_start)
                Chua_z_bifurcate.image.setRGB(listpt.x, listpt.y, clr.getRGB());
            else                                // plot in reverse
                Chua_z_bifurcate.image.setRGB(Chua_z_bifurcate.image.getWidth() - 1 - listpt.x, listpt.y, clr.getRGB());
        }
    }

    @Override protected void done()
    {
        Chua_z_bifurcate.lblImage.repaint();
        Chua_z_bifurcate.btnRun.setEnabled(true);
    }
}
