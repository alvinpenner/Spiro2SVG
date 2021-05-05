
package rossler;

// bifurcation plot of Tau versus c
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.List;
import javax.swing.*;

public class Rossler_z_bifurcate extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(600, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    protected static JButton btnRun = new JButton("Run");
    protected static JButton btnClear = new JButton("Clr");
    protected static double Tmin = 2, Tmax = 8;

    public Rossler_z_bifurcate(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        final JPanel bifurcatePanel = new JPanel();
        Main.type = "bifurcate";
        setTitle(" Rossler System - Tau bifurcate vs. c");
        setIconImage(img);
        setSize(770, 450);
        setLocationByPlatform(true);

        final JLabel[] lbl = {new JLabel("a"),
                              new JLabel("b"),
                              new JLabel("c start"),
                              new JLabel("c end"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.a)),
                                  new JTextField(Double.toString(Main.b)),
                                  new JTextField(Double.toString(Main.cstart)),
                                  new JTextField(Double.toString(Main.cend)),
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
            lbl[i].setPreferredSize(new Dimension(40, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

        JPanel rangePanel = new JPanel();
        rangePanel.setOpaque(false);
        JLabel lblrange = new JLabel("Tau = " + Tmin + " - " + Tmax);
        lblrange.setPreferredSize(new Dimension(95, 20));
        rangePanel.add(lblrange);

        JPanel posnPanel = new JPanel();
        posnPanel.setOpaque(false);
        final JLabel lblposn = new JLabel();
        lblposn.setBorder(BorderFactory.createEtchedBorder());
        lblposn.setPreferredSize(new Dimension(95, 20));
        posnPanel.add(lblposn);

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
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(btnClear);
        parmsPanel.add(rangePanel);
        parmsPanel.add(posnPanel);
        parmsPanel.setMaximumSize(new Dimension(140, 3000));
        parmsPanel.setPreferredSize(new Dimension(140, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        int ulx = 20;
        int uly = 350;
        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        DC.setFont(new Font( "SansSerif", Font.BOLD, 12 ));
        DC.setColor(new Color(40, 40, 136));
        DC.drawLine(ulx, uly, ulx + 10, uly);
        DC.drawString("increasing", ulx + 20, uly + 5);
        DC.setColor(new Color(255, 128, 0));
        DC.drawLine(ulx, uly + 20, ulx + 10, uly + 20);
        DC.drawString("decreasing", ulx + 20, uly + 25);
        DC.setColor(Color.white);

        lblImage.setBorder(BorderFactory.createEtchedBorder());
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                if (Main.cend > Main.cstart)
                    lblposn.setText(String.format(" %.4f", Main.cstart + e.getX()*(Main.cend - Main.cstart)/image.getWidth()) + ", " + String.format("%.3f", Tmax + e.getY()*(Tmin - Tmax)/image.getHeight()));
                else
                    lblposn.setText(String.format(" %.4f", Main.cend + e.getX()*(Main.cstart - Main.cend)/image.getWidth()) + ", " + String.format("%.3f", Tmax + e.getY()*(Tmin - Tmax)/image.getHeight()));
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
                Main.a = Double.parseDouble(txt[0].getText());
                Main.b = Double.parseDouble(txt[1].getText());
                Main.cstart = Double.parseDouble(txt[2].getText());
                Main.cend = Double.parseDouble(txt[3].getText());
                Main.x0 = Double.parseDouble(txt[4].getText());
                Main.y0 = Double.parseDouble(txt[5].getText());
                Main.z0 = Double.parseDouble(txt[6].getText());
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
                {Main.save_prefs();}
            });

        for (JTextField tx: txt)
            tx.addKeyListener(new KeyAdapter()
            {
                @Override public void keyPressed(KeyEvent e)
                {
                    if (e.getKeyCode() == KeyEvent.VK_ESCAPE)
                    {
                        Main.save_prefs();
                        Main.z_bifurcate.dispose();
                        Main.plot_y_vs_x = new Rossler_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                    if (e.getKeyCode() == KeyEvent.VK_F1)       // save a PNG file using the F1 key
                    {
                        if (Main.z_bifurcate != null && Main.z_bifurcate.isShowing())
                            Main.save_PNG(image, String.format("Rossler_bifurcate_%.2f_%.2f", Main.cstart, Main.cend));
                        return;
                    }
                }
            });
    }
}

class BifurcateActivity extends SwingWorker<Void, Point>
{
    Color clr;
    private static double[] pt3;

    public BifurcateActivity()
    {
        pt3 = new double[] {Main.x0, Main.y0, Main.z0};
        if (Main.cend > Main.cstart)
            clr = new Color(40, 40, 136);
        else
            clr = new Color(255, 128, 0);
    }

    protected Void doInBackground() throws Exception
    {
        int N = 100000;                             // # of iterations per c
        double delt = 0.02;
        double c;                                   // horizontal axis
        double[] zlist = new double[4];             // previous z values
        double Told = 0, Tnew = 0;                  // time of peak z

        for (int i = 0; i < 1000; i++)              // pre-initialize
            Main.runge_kutta_rossler3(pt3, delt, Main.cstart);
        //System.out.println("init 1 = " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        for (int i = 0; i < Rossler_z_bifurcate.image.getWidth(); i++)
        {
            c = Main.cstart + i*(Main.cend - Main.cstart)/Rossler_z_bifurcate.image.getWidth();
            for (int j = 0; j < N; j++)
            {
                Main.runge_kutta_rossler3(pt3, delt, c);
                if (zlist[0] < zlist[1] && zlist[1] <= zlist[2] && zlist[2] > zlist[3] && zlist[3] > pt3[2])
                {
                    Tnew = j - 2 + Main.quarticT(zlist[0] - zlist[2], zlist[1] - zlist[2], zlist[3] - zlist[2], pt3[2] - zlist[2]);
                    //System.out.println("zlist, " + zlist[0] + ", " + zlist[1] + ", " + zlist[2] + ", " + zlist[3] + ", " + pt6[2]);
                    if (delt*(Tnew - Told) > Rossler_z_bifurcate.Tmin && delt*(Tnew - Told) < Rossler_z_bifurcate.Tmax)
                        publish(new Point(i, (int) (Rossler_z_bifurcate.image.getHeight()*(Rossler_z_bifurcate.Tmax - delt*(Tnew - Told))/(Rossler_z_bifurcate.Tmax - Rossler_z_bifurcate.Tmin))));
                    //if (c == 3)
                    //    System.out.println(c + ", " + (j - 2) + ", " + Tnew + ", " + (Tnew - Told));
                    Told = Tnew;
                }
                for (int k = 0; k < 3; k++)
                    zlist[k] = zlist[k + 1];
                zlist[3] = pt3[2];
            }
            Rossler_z_bifurcate.lblImage.repaint();
        }
        return null;
    }

    @Override protected void process(List<Point> points)
    {
        for (Point listpt : points)
        {
            //System.out.println(listpt.x + ", " + listpt.y);
            if (Main.cend > Main.cstart)
                Rossler_z_bifurcate.image.setRGB(listpt.x, listpt.y, clr.getRGB());
            else                                // plot in reverse
                Rossler_z_bifurcate.image.setRGB(Rossler_z_bifurcate.image.getWidth() - 1 - listpt.x, listpt.y, clr.getRGB());
        }
    }

    @Override protected void done()
    {
        Rossler_z_bifurcate.btnRun.setEnabled(true);
    }
}
