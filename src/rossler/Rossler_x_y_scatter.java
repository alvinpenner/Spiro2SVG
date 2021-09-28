
package rossler;

// x-y scatter plot at a crossover value of z

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.List;
import javax.swing.*;

public class Rossler_x_y_scatter extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(400, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
//    protected static JCheckBox reverseChk = new JCheckBox("reverse time");
    protected static JButton btnRun = new JButton("Run");
    protected static JButton btnClear = new JButton("Clr");

    public Rossler_x_y_scatter(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        final JPanel scatterPanel = new JPanel();
        Main.type = "scatter";
        setTitle(" Rossler System - x-y Scatter (" + ScatterActivity.delt + ")");
        setIconImage(img);
        setSize(570, 450);
        setLocationByPlatform(true);

        final JLabel[] lbl = {new JLabel("a"),
                              new JLabel("b"),
                              new JLabel("c"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0"),
                              new JLabel("zc")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.a)),
                                  new JTextField(Double.toString(Main.b)),
                                  new JTextField(Double.toString(Main.c)),
                                  new JTextField(Double.toString(Main.x0)),
                                  new JTextField(Double.toString(Main.y0)),
                                  new JTextField(Double.toString(Main.z0)),
                                  new JTextField(Double.toString(Main.zc))};
        JPanel[] spacerPanel = new JPanel[2];
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
            lbl[i].setPreferredSize(new Dimension(40, 18));
            dataPanel[i].add(lbl[i]);
            txt[i].setPreferredSize(new Dimension(70, 18));
            dataPanel[i].add(txt[i]);
        }

//        reverseChk.setOpaque(false);

        JPanel xrange = new JPanel();
        xrange.setOpaque(false);
        xrange.setPreferredSize(new Dimension(120, 20));
        xrange.add(new JLabel("dx"));
        final JTextField txtxmin = new JTextField(Double.toString(Main.xmin));
        final JTextField txtxmax = new JTextField(Double.toString(Main.xmax));
        txtxmin.setPreferredSize(new Dimension(40, 18));
        txtxmax.setPreferredSize(new Dimension(40, 18));
        xrange.add(txtxmin);
        xrange.add(txtxmax);

        JPanel yrange = new JPanel();
        yrange.setOpaque(false);
        yrange.setPreferredSize(new Dimension(120, 20));
        yrange.add(new JLabel("dy"));
        final JTextField txtymin = new JTextField(Double.toString(Main.ymin));
        final JTextField txtymax = new JTextField(Double.toString(Main.ymax));
        txtymin.setPreferredSize(new Dimension(40, 18));
        txtymax.setPreferredSize(new Dimension(40, 18));
        yrange.add(txtymin);
        yrange.add(txtymax);

        JPanel posnPanel = new JPanel();
        posnPanel.setOpaque(false);
        final JLabel lblposn = new JLabel();
        lblposn.setBorder(BorderFactory.createEtchedBorder());
        lblposn.setPreferredSize(new Dimension(105, 20));
        posnPanel.add(lblposn);

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
//        parmsPanel.add(reverseChk);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(btnClear);
        parmsPanel.add(xrange);
        parmsPanel.add(yrange);
        parmsPanel.add(posnPanel);
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
                lblposn.setText(String.format(" %.3f", Main.xmin + e.getX()*(Main.xmax - Main.xmin)/image.getWidth())
                        + "," + String.format(" %.3f", Main.ymax + e.getY()*(Main.ymin - Main.ymax)/image.getHeight()));
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
                if (Main.a != Double.parseDouble(txt[0].getText())) changed = true;
                Main.a = Double.parseDouble(txt[0].getText());
                if (Main.b != Double.parseDouble(txt[1].getText())) changed = true;
                Main.b = Double.parseDouble(txt[1].getText());
                if (Main.c != Double.parseDouble(txt[2].getText())) changed = true;
                Main.c = Double.parseDouble(txt[2].getText());
                if (Main.x0 != Double.parseDouble(txt[3].getText())) changed = true;
                Main.x0 = Double.parseDouble(txt[3].getText());
                if (Main.y0 != Double.parseDouble(txt[4].getText())) changed = true;
                Main.y0 = Double.parseDouble(txt[4].getText());
                if (Main.z0 != Double.parseDouble(txt[5].getText())) changed = true;
                Main.z0 = Double.parseDouble(txt[5].getText());
                if (Main.zc != Double.parseDouble(txt[6].getText())) changed = true;
                Main.zc = Double.parseDouble(txt[6].getText());
                Main.xmin = Double.parseDouble(txtxmin.getText());
                Main.xmax = Double.parseDouble(txtxmax.getText());
                Main.ymin = Double.parseDouble(txtymin.getText());
                Main.ymax = Double.parseDouble(txtymax.getText());
                btnRun.setEnabled(false);
                ScatterActivity activity = new ScatterActivity(changed);
                activity.execute();
            }
        });
        btnClear.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                DC.clearRect(0, 0, image.getWidth(), image.getHeight());
                lblImage.repaint();
                //ScatterActivity.pt3 = new double[] {Main.x0, Main.y0, Main.z0};
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
    }
}

class ScatterActivity extends SwingWorker<Void, Point>
{
    protected static double[] pt3 = new double[] {Main.x0, Main.y0, Main.z0};
    protected static double delt = -0.02;
    private static Color clr = Color.BLACK;             // new Color(40, 40, 136);

    public ScatterActivity(boolean chg)
    {
        if (chg)
            System.arraycopy(new double[] {Main.x0, Main.y0, Main.z0}, 0, pt3, 0, pt3.length);
//        if (Rossler_x_y_scatter.reverseChk.isSelected() && delt > 0)
//            delt = -delt;
//        System.out.println(Rossler_x_y_scatter.reverseChk.isSelected() + ", " + delt);
    }

    protected Void doInBackground() throws Exception
    {
        int Nloop = 10000;
        double xc, yc;                      // interpolated (x,y) at z = zc
        double xold = 0, yold = 0, zold = Main.zc;

        System.out.println("init, " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        for (int i = 0; i < Nloop; i++)
        {
            Main.runge_kutta_rossler3(pt3, delt, Main.a, Main.b, Main.c);
            if ((zold - Main.zc)*delt < 0 && (pt3[2] - Main.zc)*delt >= 0)
            {
                xc = ((pt3[2] - Main.zc)*xold + (Main.zc - zold)*pt3[0])/(pt3[2] - zold);
                yc = ((pt3[2] - Main.zc)*yold + (Main.zc - zold)*pt3[1])/(pt3[2] - zold);
                System.out.println("cross , " + xc + ", " + yc + ", " + Main.zc);
                if (xc > Main.xmin && xc < Main.xmax && yc > Main.ymin && yc < Main.ymax)
                    publish(new Point((int) ((xc - Main.xmin)/(Main.xmax - Main.xmin)*Rossler_x_y_scatter.image.getWidth()),
                                      (int) ((yc - Main.ymax)/(Main.ymin - Main.ymax)*Rossler_x_y_scatter.image.getHeight())));
                //Rossler_x_y_scatter.lblImage.repaint();
            }
            xold = pt3[0];
            yold = pt3[1];
            zold = pt3[2];
            Rossler_x_y_scatter.lblImage.repaint();
            //System.out.println(i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        }
        //Rossler_x_y_scatter.lblImage.repaint();
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
    }
}
