
package rossler;

// x-y scatter plot at a crossover value of z

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.awt.geom.Point2D;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.util.List;
import javax.swing.*;

public class Rossler_x_y_scatter extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(400, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    protected static JButton btnRun = new JButton("Run");
    private static JButton btnClear = new JButton("Clr");
    private static JButton btnUniform = new JButton("Set-Uniform");
    protected static double del_xp = 0, del_yp = 0;
    private boolean first = true;

    public Rossler_x_y_scatter(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        final JPanel scatterPanel = new JPanel();
        Main.type = "scatter";
        setTitle(" Rossler System - x'-y' Scatter (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + ScatterActivity.delt + ")");
        setIconImage(img);
        setSize(570, 450);
        setLocationByPlatform(true);
        String fmt = "%.3f";

        final JLabel[] lbl = {new JLabel("a"),
                              new JLabel("b"),
                              new JLabel("c"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0"),
                              new JLabel("zc'")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.a)),
                                  new JTextField(Double.toString(Main.b)),
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
            lbl[i].setPreferredSize(new Dimension(40, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

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
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(btnClear);
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
                boolean changed = false;
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
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
                //if (Main.zc != Double.parseDouble(txt[6].getText())) changed = true;
                Main.zc = Double.parseDouble(txt[6].getText());         // transformed coordinates
                Main.xmin = Double.parseDouble(txtxmin.getText());      // transformed coordinates
                Main.xmax = Double.parseDouble(txtxmax.getText());      // transformed coordinates
                Main.ymin = Double.parseDouble(txtymin.getText());      // transformed coordinates
                Main.ymax = Double.parseDouble(txtymax.getText());      // transformed coordinates
                if (del_xp != Double.parseDouble(txt_del_xp.getText())) changed = true;
                del_xp = Double.parseDouble(txt_del_xp.getText());
                if (del_yp != Double.parseDouble(txt_del_yp.getText())) changed = true;
                del_yp = Double.parseDouble(txt_del_yp.getText());
                setTitle(" Rossler System - x'-y' Scatter (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + ScatterActivity.delt + ")");
                btnRun.setEnabled(false);
                ScatterActivity activity = new ScatterActivity(changed || first);
                activity.execute();
                first = false;
                DC.setColor(Color.BLUE);
                Point2D.Double pstat1 = Main.project_stationary("+");   // projected stationary point (x', y')
                Point2D.Double pstat2 = Main.project_stationary("-");   // projected stationary point (x', y')
                System.out.println("pstat1 = " + pstat1);
                System.out.println("pstat2 = " + pstat2);
                DC.drawLine((int) ((pstat1.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth()),
                            (int) ((pstat1.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()),
                            (int) ((pstat2.x - Main.xmin)/(Main.xmax - Main.xmin)*image.getWidth()),
                            (int) ((pstat2.y - Main.ymax)/(Main.ymin - Main.ymax)*image.getHeight()));
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
                Rossler_y_vs_x.fit_linear_response();
                txt[0].setText(Double.toString(Main.a));
                txt[1].setText(Double.toString(Main.b));
                txt[2].setText(Double.toString(Main.c));
                txt[3].setText(String.format("%.8f", Main.final_x));
                txt[4].setText(String.format("%.8f", Main.final_y));
                txt[5].setText(String.format("%.8f", Main.final_z));
                txt[6].setText(String.format("%.8f", Main.project_zp(Main.final_x, Main.final_y, Main.final_z)));
                Main.skew_transform = true;             // make the linear response "uniform"
                setTitle(" Rossler System - x'-y' Scatter (" + String.format("%.4f", Main.project_phi) + ", " + String.format("%.4f", Main.project_theta) + ", " + Double.toString(Main.project_psi) + ") (" + ScatterActivity.delt + ")");
                //System.out.println("ScatterActivity = " + ScatterActivity.pt3_cross[0] + ", " + ScatterActivity.pt3_cross[1] + ", " + ScatterActivity.pt3_cross[2]);
                //Point2D.Double pt2_cross = Main.project_2D(ScatterActivity.pt3_cross[0], ScatterActivity.pt3_cross[1], ScatterActivity.pt3_cross[2]);
                //System.out.println("ScatterActivity = " + pt2_cross.x + ", " + pt2_cross.y + ", " + Main.project_zp(ScatterActivity.pt3_cross[0], ScatterActivity.pt3_cross[1], ScatterActivity.pt3_cross[2]));
                //txt[6].setText(String.format(" %.4f", Main.project_zp(ScatterActivity.pt3_cross[0], ScatterActivity.pt3_cross[1], ScatterActivity.pt3_cross[2])));      // refresh zc'
                //txtxmin.setText(String.format(" %.3f", pt2_cross.x - (Main.xmax - Main.xmin)/2));
                //txtxmax.setText(String.format(" %.3f", pt2_cross.x + (Main.xmax - Main.xmin)/2));
                //txtymin.setText(String.format(" %.3f", pt2_cross.y - (Main.ymax - Main.ymin)/2));
                //txtymax.setText(String.format(" %.3f", pt2_cross.y + (Main.ymax - Main.ymin)/2));
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
                        Main.x_y_scatter.dispose();
                        Main.plot_y_vs_x = new Rossler_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                    if (e.getKeyCode() == KeyEvent.VK_F1)       // save a PNG file using the F1 key
                    {
                        if (Main.x_y_scatter != null && Main.x_y_scatter.isShowing())
                            Main.save_PNG(image, String.format("Rossler_scatter_%f_%.2f_%.2f", Main.a, Main.b, Main.c));
                        return;
                    }
                }
            });
        //System.out.println("dataPanel = " + dataPanel[0].getSize());
    }
}

class ScatterActivity extends SwingWorker<Void, Point>
{
    protected static double delt = -0.002; // -0.02;
    protected static double[] pt3_cross = new double[3];
    private static double[] pt3 = new double[3];
    private static double xold = 0, yold = 0, zold = Main.zc;
    private static Color clr = Color.BLACK;            // new Color(40, 40, 136);
    private static PrintWriter fout = null;
    private static int iT;                             // # iterations

    public ScatterActivity(boolean chg)
    {
        if (chg)
        {
            iT = 0;
            System.out.println("reset, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ", " + Main.zc + ", " + delt);
            System.out.println(" ,iT,x,y,z,x',y',z',tc,phi,theta");
            System.arraycopy(new double[] {Main.x0, Main.y0, Main.z0}, 0, pt3, 0, pt3.length);
            Point2D.Double pt2 = Main.project_2D(pt3[0], pt3[1], pt3[2]);
            System.out.println("init  , " + iT + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]));
            // back-transform del (x',y') into del (x,y,z)
            pt3[0] += Main.invert_from_xp_yp(Rossler_x_y_scatter.del_xp, Rossler_x_y_scatter.del_yp, 0, "x");
            pt3[1] += Main.invert_from_xp_yp(Rossler_x_y_scatter.del_xp, Rossler_x_y_scatter.del_yp, 0, "y");
            pt3[2] += Main.invert_from_xp_yp(Rossler_x_y_scatter.del_xp, Rossler_x_y_scatter.del_yp, 0, "z");
            // test velocity direction test code only, to be deleted
/*
            System.out.print("test1 , " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
            double xdot = -pt3[1] - pt3[2];
            double ydot =  pt3[0] + Main.a*pt3[1];
            double zdot =  Main.b + pt3[2]*(pt3[0] - Main.c);
            System.out.println(", " + xdot + ", " + ydot + ", " + zdot + ", " + Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));
            Point2D.Double pt2 = Main.project_2D(xdot, ydot, zdot);
            System.out.println("test2 , " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(xdot, ydot, zdot));
            System.out.println(Math.atan2(xdot, -ydot)*180/Math.PI
                      + ", " + Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot))*180/Math.PI);
*/
        }

        if (false)
            try
            {
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\x_y_scatter_Output_" + Main.a + "_" + Main.b + "_" + Main.c + "_" + Main.project_phi + "_" + Main.project_theta + "_" + delt + ".csv", true);
                fout = new PrintWriter(fw);
                //fout.println("x_y_scatter, " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Main.project_phi + ", " + Main.project_theta + ", " + delt);
                //fout.println("iT,x,y,z,xdot,ydot,zdot,v_phi,v_theta,x2dot,y2dot,z2dot,kx,ky,kz,k_phi,k_theta");
            }
            catch (java.io.IOException e)
                {System.out.println("Rossler_Output.csv save error = " + e);}
    }

    protected Void doInBackground() throws Exception
    {
        int Nloop = 10000*100;
        Point2D.Double pt2old;          // projected, old, (x', y') using Euler angles
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        double zpold, zpnew;            // projected z'
        double xc, yc;                  // interpolated, projected (x', y')

        for (int i = 0; i < Nloop; i++)
        {
            iT++;
            Main.runge_kutta_rossler3(pt3, delt, Main.a, Main.b, Main.c);
            //System.out.println("old = " + zold + ", " + zp(xold, yold, zold));
            //System.out.println("new = " + pt3[2] + ", " + zp(pt3[0], pt3[1], pt3[2]));
            //if ((Main.project_zp(xold, yold, zold) - Main.zc)*delt < 0 && (Main.project_zp(pt3[0], pt3[1], pt3[2]) - Main.zc)*delt >= 0)
            //System.out.println(i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
            //System.out.println("test z , " + iT + ", " + xold + ", " + yold + ", " + zold + ", " + Main.project_zp(xold, yold, zold) + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]) + ", " + Main.zc);
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
                //xdot = -yc - Main.zc;
                //ydot =  xc + Main.a*yc;
                //zdot =  Main.b + Main.zc*(xc - Main.c);
                //System.out.println("cross , " + iT + ", " + xc + ", " + yc + ", " + (iT + (Main.zc - zpold)/(zpnew - zpold))
                //                       + ", " + xdot + ", " + ydot + ", " + zdot
                //                       + ", " + Math.atan2(xdot, -ydot)*180/Math.PI
                //                       + ", " + Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot))*180/Math.PI);
                //System.out.println("z cross , " + iT + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + xc + ", " + yc + ", " + (iT + (Main.zc - zpold)/(zpnew - zpold))); // original code
                System.out.println("z inter , " + iT + ", " + ((zpnew - Main.zc)*xold + (Main.zc - zpold)*pt3[0])/(zpnew - zpold)
                                                     + ", " + ((zpnew - Main.zc)*yold + (Main.zc - zpold)*pt3[1])/(zpnew - zpold)
                                                     + ", " + ((zpnew - Main.zc)*zold + (Main.zc - zpold)*pt3[2])/(zpnew - zpold)
                                                     + ", " + xc + ", " + yc + ", " + (iT + (Main.zc - zpold)/(zpnew - zpold)));
                pt3_cross[0] = pt3[0];
                pt3_cross[1] = pt3[1];
                pt3_cross[2] = pt3[2];
                if (fout != null)
                    fout.println(iT + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + xc + ", " + yc);
                if (xc > Main.xmin && xc < Main.xmax && yc > Main.ymin && yc < Main.ymax)
                    publish(new Point((int) ((xc - Main.xmin)/(Main.xmax - Main.xmin)*Rossler_x_y_scatter.image.getWidth()),
                                      (int) ((yc - Main.ymax)/(Main.ymin - Main.ymax)*Rossler_x_y_scatter.image.getHeight())));
                //Rossler_x_y_scatter.lblImage.repaint();
            }
            if (false && iT % 500 == 0)
            {
                pt2new = Main.project_2D(pt3[0], pt3[1], pt3[2]);
                System.out.println("t cross , " + iT + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2new.x + ", " + pt2new.y + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]));
                if (fout != null)
                    fout.println(iT + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2] + ", " + pt2new.x + ", " + pt2new.y + ", " + Main.project_zp(pt3[0], pt3[1], pt3[2]));
                if (pt2new.x > Main.xmin && pt2new.x < Main.xmax && pt2new.y > Main.ymin && pt2new.y < Main.ymax)
                    publish(new Point((int) ((pt2new.x - Main.xmin)/(Main.xmax - Main.xmin)*Rossler_x_y_scatter.image.getWidth()),
                                      (int) ((pt2new.y - Main.ymax)/(Main.ymin - Main.ymax)*Rossler_x_y_scatter.image.getHeight())));
            }
            xold = pt3[0];
            yold = pt3[1];
            zold = pt3[2];
            Rossler_x_y_scatter.lblImage.repaint();
            //System.out.println(i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        }
        //System.out.println("final , " + pt3_cross[0] + ", " + pt3_cross[1] + ", " + pt3_cross[2]);
        return null;
    }

    @Override protected void process(List<Point> points)
    {
        for (Point listpt : points)
        {
            //System.out.println(listpt.x + ", " + listpt.y + ", " + clr.getRGB());
            Rossler_x_y_scatter.image.setRGB(listpt.x, listpt.y, clr.getRGB());
        }
    }

    @Override protected void done()
    {
        Rossler_x_y_scatter.btnRun.setEnabled(true);
        if (fout != null)
            fout.close();
    }
}
