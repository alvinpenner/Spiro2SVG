
package chua;

// x-y scatter plot at a crossover value of z

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.awt.geom.Point2D;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.util.List;
import javax.swing.*;

public class Chua_x_y_scatter extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(400, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    protected static JButton btnRun = new JButton("Run");
    private static JButton btnClear = new JButton("Clr");
    protected static JCheckBox printChk = new JCheckBox("  print  ");
    private static JButton btnUniform = new JButton("Set-Uniform");
    protected static double del_xp = 0, del_yp = 0;
    protected static String first_order_hdr;                // # iterations, eig, angle
    private boolean first = true;

    public Chua_x_y_scatter(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        final JPanel scatterPanel = new JPanel();
        Main.type = "scatter";
        setTitle(" Chua Oscillator - x'-y' Scatter (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + ScatterActivity.delt + ")");
        setIconImage(img);
        setSize(570, 450 + 24 + 28);
        setLocationByPlatform(true);
        String fmt = "%.3f";

        final JLabel[] lbl = {new JLabel("alpha"),
                              new JLabel("beta"),
                              new JLabel("gamma"),
                              new JLabel("c"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0"),
                              new JLabel("zc'")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.alpha)),
                                  new JTextField(Double.toString(Main.beta)),
                                  new JTextField(Double.toString(Main.gamma)),
                                  new JTextField(Double.toString(Main.c)),
                                  new JTextField(Double.toString(Main.x0)),
                                  new JTextField(Double.toString(Main.y0)),
                                  new JTextField(Double.toString(Main.z0)),
                                  new JTextField(Double.toString(Main.zc))};
        JPanel[] spacerPanel = new JPanel[3];
        JPanel[] dataPanel = new JPanel[lbl.length];

        for (int i = 0; i < spacerPanel.length; i++)
        {
            spacerPanel[i] = new JPanel();
            spacerPanel[i].setPreferredSize(new Dimension(110, 1));
            spacerPanel[i].setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
            spacerPanel[i].setOpaque(false);
        }
        for (int i = 0; i < dataPanel.length; i++)
        {
            dataPanel[i] = new JPanel();
            dataPanel[i].setOpaque(false);
            dataPanel[i].setPreferredSize(new Dimension(125, 24));
            lbl[i].setPreferredSize(new Dimension(45, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

        printChk.setOpaque(false);
        JPanel xrange = new JPanel();
        xrange.setOpaque(false);
        xrange.setPreferredSize(new Dimension(120, 20));
        xrange.add(new JLabel("dx'"));
        final JTextField txtxmin = new JTextField(String.format(fmt, Main.xmin));
        final JTextField txtxmax = new JTextField(String.format(fmt, Main.xmax));
        txtxmin.setPreferredSize(new Dimension(40, 18));
        txtxmax.setPreferredSize(new Dimension(40, 18));
        xrange.add(txtxmin);
        xrange.add(txtxmax);

        JPanel yrange = new JPanel();
        yrange.setOpaque(false);
        yrange.setPreferredSize(new Dimension(120, 20));
        yrange.add(new JLabel("dy'"));
        final JTextField txtymin = new JTextField(String.format(fmt, Main.ymin));
        final JTextField txtymax = new JTextField(String.format(fmt, Main.ymax));
        txtymin.setPreferredSize(new Dimension(40, 18));
        txtymax.setPreferredSize(new Dimension(40, 18));
        yrange.add(txtymin);
        yrange.add(txtymax);

        JPanel posnPanel = new JPanel();                // live display of (x',y')
        posnPanel.setOpaque(false);
        final JLabel lblposn = new JLabel();
        lblposn.setBorder(BorderFactory.createEtchedBorder());
        lblposn.setPreferredSize(new Dimension(105, 20));
        posnPanel.add(lblposn);

        JPanel pnl_del_xp_yp = new JPanel();                // increment initial (x',y')
        pnl_del_xp_yp.setOpaque(false);
        pnl_del_xp_yp.setPreferredSize(new Dimension(120, 20));
        pnl_del_xp_yp.add(new JLabel("del '"));
        final JTextField txt_del_xp = new JTextField("" + del_xp);
        final JTextField txt_del_yp = new JTextField("" + del_yp);
        txt_del_xp.setPreferredSize(new Dimension(40, 18));
        txt_del_yp.setPreferredSize(new Dimension(40, 18));
        pnl_del_xp_yp.add(txt_del_xp);
        pnl_del_xp_yp.add(txt_del_yp);

        parmsPanel.removeAll();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(dataPanel[0]);
        parmsPanel.add(dataPanel[1]);
        parmsPanel.add(dataPanel[2]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[7]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(btnClear);
        parmsPanel.add(printChk);
        parmsPanel.add(xrange);
        parmsPanel.add(yrange);
        parmsPanel.add(posnPanel);
        parmsPanel.add(pnl_del_xp_yp);
        parmsPanel.add(spacerPanel[2]);
        parmsPanel.add(btnUniform);            // make first-order response uniform
        parmsPanel.setMaximumSize(new Dimension(140, 3000));
        parmsPanel.setPreferredSize(new Dimension(140, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());

        lblImage.setBorder(BorderFactory.createEtchedBorder());
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                lblposn.setText(String.format(" %.4f", Main.xmin + e.getX()*(Main.xmax - Main.xmin)/image.getWidth())
                        + "," + String.format(" %.4f", Main.ymax + e.getY()*(Main.ymin - Main.ymax)/image.getHeight()));
            }
        });
        scatterPanel.add(lblImage);

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(parmsPanel, BorderLayout.WEST);
        getContentPane().add(scatterPanel, BorderLayout.EAST);
        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (first_order_hdr == null)
                {
                    System.out.println("'first_order_hdr' is not initialized");
                    return;
                }
                boolean changed = false;
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                //if (Main.alpha != Double.parseDouble(txt[0].getText())) changed = true;
                Main.alpha = Double.parseDouble(txt[0].getText());
                //if (Main.beta != Double.parseDouble(txt[1].getText())) changed = true;
                Main.beta = Double.parseDouble(txt[1].getText());
                //if (Main.gamma != Double.parseDouble(txt[2].getText())) changed = true;
                Main.gamma = Double.parseDouble(txt[2].getText());
                //if (Main.c != Double.parseDouble(txt[3].getText())) changed = true;
                Main.c = Double.parseDouble(txt[3].getText());
                if (Main.x0 != Double.parseDouble(txt[4].getText())) changed = true;
                Main.x0 = Double.parseDouble(txt[4].getText());
                if (Main.y0 != Double.parseDouble(txt[5].getText())) changed = true;
                Main.y0 = Double.parseDouble(txt[5].getText());
                if (Main.z0 != Double.parseDouble(txt[6].getText())) changed = true;
                Main.z0 = Double.parseDouble(txt[6].getText());
                //if (Main.zc != Double.parseDouble(txt[7].getText())) changed = true;
                Main.zc = Double.parseDouble(txt[7].getText());         // transformed coordinates
                Main.xmin = Double.parseDouble(txtxmin.getText());      // transformed coordinates
                Main.xmax = Double.parseDouble(txtxmax.getText());      // transformed coordinates
                Main.ymin = Double.parseDouble(txtymin.getText());      // transformed coordinates
                Main.ymax = Double.parseDouble(txtymax.getText());      // transformed coordinates
                if (del_xp != Double.parseDouble(txt_del_xp.getText())) changed = true;
                del_xp = Double.parseDouble(txt_del_xp.getText());
                if (del_yp != Double.parseDouble(txt_del_yp.getText())) changed = true;
                del_yp = Double.parseDouble(txt_del_yp.getText());
                setTitle(" Chua Oscillator - x'-y' Scatter (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + Main.a + ", " + ScatterActivity.delt + ")");
                btnRun.setEnabled(false);
                ScatterActivity activity = new ScatterActivity(changed || first);
                activity.execute();
                first = false;
                DC.setColor(Color.BLUE);
                //Point2D.Double pstat = Main.project_stationary();   // projected stationary point (x', y')
                //System.out.println("pstat = " + pstat);
                //DC.drawLine((int) ((-pstat.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth()),
                //            (int) ((-pstat.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()),
                //            (int) ( (pstat.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth()),
                //            (int) ( (pstat.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()));
                //double kx = Main.calc_ydot(Main.x0, Main.y0, Main.z0)*Main.calc_z2dot(Main.x0, Main.y0, Main.z0) - Main.calc_zdot(Main.x0, Main.y0, Main.z0)*Main.calc_y2dot(Main.x0, Main.y0, Main.z0);
                //double ky = Main.calc_zdot(Main.x0, Main.y0, Main.z0)*Main.calc_x2dot(Main.x0, Main.y0, Main.z0) - Main.calc_xdot(Main.x0, Main.y0, Main.z0)*Main.calc_z2dot(Main.x0, Main.y0, Main.z0);
                //double kz = Main.calc_xdot(Main.x0, Main.y0, Main.z0)*Main.calc_y2dot(Main.x0, Main.y0, Main.z0) - Main.calc_ydot(Main.x0, Main.y0, Main.z0)*Main.calc_x2dot(Main.x0, Main.y0, Main.z0);
                //Point2D.Double pt2 = Main.project_2D(kx, ky, kz);
                //double kangle = Math.atan2(pt2.y, pt2.x);
                //kangle = 58*Math.PI/180;
                //System.out.println("project k = " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(kx, ky, kz) + ", " + kangle*180/Math.PI);
                Point2D.Double pt2 = Main.project_2D(Main.x0, Main.y0, Main.z0);
                //DC.drawLine((int) ((pt2.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth()),
                //            (int) ((pt2.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()),
                //            (int) ((pt2.x + (Main.xmax - Main.xmin)/10*Math.cos(kangle) - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth()),
                //            (int) ((pt2.y + (Main.xmax - Main.xmin)/10*Math.sin(kangle) - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()));
                DC.drawLine((int) ((pt2.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth() - 2),
                            (int) ((pt2.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight() - 2),
                            (int) ((pt2.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth() + 2),
                            (int) ((pt2.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()) + 2);
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
        btnUniform.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                first_order_hdr = Chua_y_vs_x.fit_linear_response();
                txt[0].setText(Double.toString(Main.alpha));
                txt[1].setText(Double.toString(Main.beta));
                txt[2].setText(Double.toString(Main.gamma));
                txt[3].setText(Double.toString(Main.c));
                txt[4].setText(String.format("%.8f", Main.final_x));
                txt[5].setText(String.format("%.8f", Main.final_y));
                txt[6].setText(String.format("%.8f", Main.final_z));
                txt[7].setText(String.format("%.8f", Main.project_zp(Main.final_x, Main.final_y, Main.final_z)));
                Main.skew_transform = true;             // make the linear response "uniform"
                setTitle(" Chua Oscillator - x'-y' Scatter (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + Main.a + ", " + ScatterActivity.delt + ")");
                //System.out.println("ScatterActivity = " + ScatterActivity.pt3_cross[0] + ", " + ScatterActivity.pt3_cross[1] + ", " + ScatterActivity.pt3_cross[2]);
                //Point2D.Double pt2_cross = Main.project_2D(ScatterActivity.pt3_cross[0], ScatterActivity.pt3_cross[1], ScatterActivity.pt3_cross[2]);
                //System.out.println("ScatterActivity = " + pt2_cross.x + ", " + pt2_cross.y + ", " + Main.project_zp(ScatterActivity.pt3_cross[0], ScatterActivity.pt3_cross[1], ScatterActivity.pt3_cross[2]));
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
                        Main.x_y_scatter.dispose();
                        Main.plot_y_vs_x = new Chua_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                    if (e.getKeyCode() == KeyEvent.VK_F1)       // save a PNG file using the F1 key
                    {
                        if (Main.x_y_scatter != null && Main.x_y_scatter.isShowing())
                            Main.save_PNG(image, String.format("Chua_scatter_%f_%.2f_%.2f", Main.alpha, Main.beta, Main.gamma));
                        return;
                    }
                }
            });
        //System.out.println("dataPanel = " + dataPanel[0].getSize());
    }
}

class ScatterActivity extends SwingWorker<Void, Point>
{
    protected static double delt = 0.00005;
    //protected static double[] Cx = new double[] {-1.7144498681885263 , 2.887590649924726 , 0.41865398546281907};
    //protected static double[] Cy = new double[] {-2.7302926626819404 , -0.9454671597617976 , -1.5553130222501819};
    //protected static double[] Cx = new double[] {-1.9854121867626764 , 5.278816440613675 , 1.1677876751189247};
    //protected static double[] Cy = new double[] {-4.781366025646555 , -2.8490240722558156 , -1.573431311366903};
    //protected static double[] Cx = new double[] {-1.107155942014133 , 8.7312104139989 , -34.17683230319818};
    //protected static double[] Cy = new double[] {-0.005320900593378551 , 43.475293577200446 , 9.063606245173178};
    protected static double[] Cx = new double[] {0, 0, 0};
    protected static double[] Cy = new double[] {0, 0, 0};
    protected static Point2D.Double pt2_proj;           // projected (x0, y0, z0)
    protected static double[] pt3_cross = new double[3];
    private static double[] pt3 = new double[3];
    private static double xold = 0, yold = 0, zold = Main.zc;
    private static Color clr = Color.BLACK;            // new Color(40, 40, 136);
    private static PrintWriter fout = null;
    private static int iter;                             // # iterations

    public ScatterActivity(boolean chg)
    {
        if (true)
            try
            {
                //FileWriter fw = new FileWriter("C:\\Windows\\Temp\\Chua_scatter_" + Chua_x_y_scatter.first_order_hdr.split(",")[0] + "_" + Main.alpha + "_" + Main.beta + "_" + Main.gamma + "_" + Main.a + "_" + Main.c + ".csv", true);
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\Chua_scatter_" + Chua_x_y_scatter.first_order_hdr.split(",")[0] + "_" + Chua_x_y_scatter.del_xp + "_" + Chua_x_y_scatter.del_yp + "_" + Main.alpha + ".csv", true);
                fout = new PrintWriter(fw);
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}

        if (chg)
        {
            iter = 0;
            pt2_proj = Main.project_2D(Main.x0, Main.y0, Main.z0);          // projected (x0, y0, z0)
            //System.out.println("ScatterActivity 0, " + Main.project_phi + ", " + Main.project_theta + ", " + Main.x0 + ", " + Main.y0 + ", " +  Main.z0 + ", " + pt2_proj);
            System.out.println("reset, " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + delt + ", " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ", " + Main.zc);
            System.out.println(" ,iter,x,y,z,x',y',z',tc,phi,theta");
            System.arraycopy(new double[] {Main.x0, Main.y0, Main.z0}, 0, pt3, 0, pt3.length);
            Point2D.Double pt2 = Main.project_2D(pt3[0], pt3[1], pt3[2]);
            System.out.println("init  , " + Chua_x_y_scatter.first_order_hdr + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]));
            //System.out.println("test data org ," + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2.x + ", " + pt2.y);
            // back-transform del (x',y') into del (x,y,z)
            pt3[0] += Main.invert_from_xp_yp(Chua_x_y_scatter.del_xp, Chua_x_y_scatter.del_yp, 0, "x");
            pt3[1] += Main.invert_from_xp_yp(Chua_x_y_scatter.del_xp, Chua_x_y_scatter.del_yp, 0, "y");
            pt3[2] += Main.invert_from_xp_yp(Chua_x_y_scatter.del_xp, Chua_x_y_scatter.del_yp, 0, "z");
            //System.out.println("ScatterActivity 1, " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + Main.project_2D(pt3[0], pt3[1], pt3[2]));

            if (fout != null)
            {
                Point2D.Double pstat = Main.project_stationary();           // projected stationary point (x', y')
                fout.println("           , alpha, beta, gamma, a, c, delt, phi, theta, psi");
                fout.println("x_y_scatter, " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + delt + ", " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi);
                fout.println("scatter hdr, " + Chua_x_y_scatter.first_order_hdr);
                fout.println("init x0 y0 ," + pt2.x + ", " + pt2.y + ", " + pstat.x + ", " + pstat.y);
                fout.println("iter       , x', y'");
            }
        }
        //System.out.println("ScatterActivity 1, " + chg + ", " + Main.x0 + ", " + Main.y0 + ", " +  Main.z0 + ", " + pt2_proj);
        //System.out.println("ScatterActivity 2, " + chg + ", " + Main.x0 + ", " + Main.y0 + ", " +  Main.z0 + ", " + Main.project_2D(pt3[0], pt3[1], pt3[2]));
    }

    protected Void doInBackground() throws Exception
    {
        int Nloop = 5*10000*2000;
        Point2D.Double pt2old;          // projected, old, (x', y') using Euler angles
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        Point2D.Double ztrans;
        double zpold, zpnew;            // projected z'
        double xc, yc;                  // interpolated, projected (x', y')
        double xplt, yplt;              // transformed coord

        System.out.println("print = " + Chua_x_y_scatter.printChk.isSelected() + " at " + iter);
        for (int i = 0; i < Nloop; i++)
        {
            iter++;
            Main.runge_kutta_chua3(pt3, delt);
            //System.out.println("old = " + zold + ", " + zp(xold, yold, zold));
            //System.out.println("new = " + pt3[2] + ", " + zp(pt3[0], pt3[1], pt3[2]));
            //if ((Main.project_zp(xold, yold, zold) - Main.zc)*delt < 0 && (Main.project_zp(pt3[0], pt3[1], pt3[2]) - Main.zc)*delt >= 0)
            //System.out.println(i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
            //System.out.println("test z , " + iter + ", " + xold + ", " + yold + ", " + zold + ", " + Main.project_zp(xold, yold, zold) + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]) + ", " + Main.zc);
            if (true && (Main.project_zp(xold, yold, zold) - Main.zc)*delt <= 0 && (Main.project_zp(pt3[0], pt3[1], pt3[2]) - Main.zc)*delt > 0)
            {
                pt2old = Main.project_2D(xold, yold, zold);
                pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                zpold = Main.project_zp(xold, yold, zold);
                zpnew = Main.project_zp(pt3[0], pt3[1], pt3[2]);
                //System.out.println("zp = ," + zpold + ", " + zpnew);
                //System.out.println("old = ," + xold + ", " + yold + ", " + zold);
                //System.out.println("new = ," + i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
                xc = ((zpnew - Main.zc)*pt2old.x + (Main.zc - zpold)*pt2new.x)/(zpnew - zpold);
                yc = ((zpnew - Main.zc)*pt2old.y + (Main.zc - zpold)*pt2new.y)/(zpnew - zpold);
                // transform to remove quadratic component (xc, yc)
                //System.out.println("ScatterActivity in , " + i + ", " + xc + ", " + yc + ", " + pt2_proj.x + ", " + pt2_proj.y);
                //xplt = transform_quadratic("x", xc, yc);
                //yplt = transform_quadratic("y", xc, yc);
                ztrans = Chua_Simul_3.transform(xc - pt2_proj.x, yc - pt2_proj.y);
                xplt = ztrans.x + pt2_proj.x;
                yplt = ztrans.y + pt2_proj.y;
                //xplt = xc;
                //yplt = yc;

                //System.out.println("ScatterActivity out, " + i + ", " + xplt + ", " + yplt + ", " + pt2_proj.x + ", " + pt2_proj.y);
                if (true)
                {
                    System.out.println("z inter , " + iter + ", " + ((zpnew - Main.zc)*xold + (Main.zc - zpold)*pt3[0])/(zpnew - zpold)
                                                           + ", " + ((zpnew - Main.zc)*yold + (Main.zc - zpold)*pt3[1])/(zpnew - zpold)
                                                           + ", " + ((zpnew - Main.zc)*zold + (Main.zc - zpold)*pt3[2])/(zpnew - zpold)
                                                           + ", " + xplt + ", " + yplt + ", " + (iter + (Main.zc - zpold)/(zpnew - zpold)));
                }
                pt3_cross[0] = pt3[0];
                pt3_cross[1] = pt3[1];
                pt3_cross[2] = pt3[2];
                if (fout != null && Chua_x_y_scatter.printChk.isSelected())
                    fout.println(iter + ", " + xc + ", " + yc);
                if (xplt > Main.xmin && xplt < Main.xmax && yplt > Main.ymin && yplt < Main.ymax)
                    publish(new Point((int) ((xplt - Main.xmin)/(Main.xmax - Main.xmin)*Chua_x_y_scatter.image.getWidth()),
                                      (int) ((yplt - Main.ymax)/(Main.ymin - Main.ymax)*Chua_x_y_scatter.image.getHeight())));
            }
            if (false && iter % 500 == 0)
            {
                pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                System.out.println("t cross , " + iter + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2new.x + ", " + pt2new.y + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]));
                if (fout != null)
                    fout.println(iter + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2new.x + ", " + pt2new.y + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]));
                if (pt2new.x > Main.xmin && pt2new.x < Main.xmax && pt2new.y > Main.ymin && pt2new.y < Main.ymax)
                    publish(new Point((int) ((pt2new.x - Main.xmin)/(Main.xmax - Main.xmin)*Chua_x_y_scatter.image.getWidth()),
                                      (int) ((pt2new.y - Main.ymax)/(Main.ymin - Main.ymax)*Chua_x_y_scatter.image.getHeight())));
            }
            xold = pt3[0];
            yold = pt3[1];
            zold = pt3[2];
            Chua_x_y_scatter.lblImage.repaint();
            //System.out.println(i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        }
        //System.out.println("final , " + iter + ", " + pt3_cross[0] + ", " + pt3_cross[1] + ", " + pt3_cross[2]);
        return null;
    }

//    protected Point2D.Double transform_z(double xt, double yt)
//    {
        //System.out.println("in  , " + pt2_proj + ", " + x + ", " + y);
//        xt = xt - pt2_proj.x;
//        yt = yt - pt2_proj.y;
        //ztrans = Chua_Simul_3.transform(xt, yt);

        //System.out.println("disp, " + pt2_proj + ", " + x + ", " + y);

//            return xt - Cx[0]*xt*xt - Cx[1]*xt*yt - Cx[2]*yt*yt + pt2_proj.x;
//        else
//            return yt - Cy[0]*xt*xt - Cy[1]*xt*yt - Cy[2]*yt*yt + pt2_proj.y;
//    }

    protected double transform_quadratic(String type, double xt, double yt)
    {
        //System.out.println("in  , " + pt2_proj + ", " + x + ", " + y);
        xt = xt - pt2_proj.x;
        yt = yt - pt2_proj.y;
        //System.out.println("disp, " + pt2_proj + ", " + x + ", " + y);
        if (type.equals("x"))
            return xt - Cx[0]*xt*xt - Cx[1]*xt*yt - Cx[2]*yt*yt + pt2_proj.x;
        else
            return yt - Cy[0]*xt*xt - Cy[1]*xt*yt - Cy[2]*yt*yt + pt2_proj.y;
    }

    @Override protected void process(List<Point> points)
    {
        for (Point listpt : points)
        {
            //System.out.println(listpt.x + ", " + listpt.y + ", " + clr.getRGB());
            Chua_x_y_scatter.image.setRGB(listpt.x, listpt.y, clr.getRGB());
        }
    }

    @Override protected void done()
    {
        Chua_x_y_scatter.btnRun.setEnabled(true);
        if (fout != null)
            fout.close();
    }
}
