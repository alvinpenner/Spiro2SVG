
package zeeman;

// plot a bifurcation diagram of theta as fxn of y0
// implement Runge-Kutta: see Froberg page 269
// for animation, see CoreJava v2ch06, ProgressBarTest.java

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.List;
import javax.swing.*;

public class Bifurcate extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(600, 400, BufferedImage.TYPE_3BYTE_BGR);
    protected static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    protected static JPanel bifurcatePanel = new JPanel();
    protected static JCheckBox initChk = new JCheckBox("constant init");
    private static JPanel parmsPanel = new JPanel();
    protected static JButton btnRun = new JButton("Run");
    protected static JButton btnClear = new JButton("Clr");
    protected static double thmin = 1.5, thmax = 4.8;
    protected static boolean activitycancel = false;
    //protected static double thmin = 0, thmax = 6.28;

    public Bifurcate(Image img)
    {
        JLabel[] lbl = {new JLabel("A"),
                        new JLabel("θ0"),
                        new JLabel("ω0"),
                        new JLabel("x0"),
                        new JLabel("xa"),
                        new JLabel("ystart"),
                        new JLabel("yend"),
                        new JLabel("c"),
                        new JLabel("Tx"),
                        new JLabel("φ0")};
        final JTextField[] txt = {new JTextField(Double.toString(main.A)),
                                  new JTextField(String.format("%.2f", main.theta0*180/Math.PI)),
                                  new JTextField(String.format("%.3f", main.w0)),
                                  new JTextField(Double.toString(main.x0)),
                                  new JTextField(String.format("%.2f", main.xa)),
                                  new JTextField(String.format("%.3f", main.ystart)),
                                  new JTextField(String.format("%.3f", main.yend)),
                                  new JTextField(Double.toString(main.c)),
                                  new JTextField(Double.toString(main.Tx)),
                                  new JTextField(String.format("%.2f", main.phi0*180/Math.PI))};
        JPanel[] spacerPanel = new JPanel[2];
        JPanel[] dataPanel = new JPanel[10];

        setTitle(" Dynamic Zeeman - Bifurcate diagram (θ vs. y)");
        setIconImage(img);
        setSize(755, 530);
        setLocationByPlatform(true);

        for (int i = 0; i < spacerPanel.length; i++)
        {
            spacerPanel[i] = new JPanel();
            spacerPanel[i].setPreferredSize(new Dimension(85, 6));
            spacerPanel[i].setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
            spacerPanel[i].setOpaque(false);
        }
        for (int i = 0; i < dataPanel.length; i++)
        {
            dataPanel[i] = new JPanel();
            dataPanel[i].setOpaque(false);
            if (i < 7)
                lbl[i].setPreferredSize(new Dimension(35, 18));
            else
                lbl[i].setPreferredSize(new Dimension(25, 18));
            dataPanel[i].add(lbl[i]);
            if (i < 7)
                txt[i].setPreferredSize(new Dimension(50, 18));
            else
                txt[i].setPreferredSize(new Dimension(60, 18));
            dataPanel[i].add(txt[i]);
        }
        txt[0].setEditable(false);                              // A distance
        initChk.setOpaque(false);

        JPanel thrangePanel = new JPanel();
        thrangePanel.setOpaque(false);
        JLabel lblthrange = new JLabel("θ = " + thmin + " - " + thmax);
        lblthrange.setPreferredSize(new Dimension(85, 20));
        thrangePanel.add(lblthrange);

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
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(initChk);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[7]);
        parmsPanel.add(dataPanel[8]);
        parmsPanel.add(dataPanel[9]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(btnClear);
        parmsPanel.add(thrangePanel);
        parmsPanel.add(posnPanel);
        parmsPanel.setMaximumSize(new Dimension(125, 3000));
        parmsPanel.setPreferredSize(new Dimension(125, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        DC.setFont(new Font( "SansSerif", Font.BOLD, 12 ));
        DC.setColor(new Color(0, 160, 0));
        DC.drawLine(480, 55, 490, 55);
        DC.drawString("constant init", 500, 60);
        DC.setColor(new Color(40, 40, 136));
        DC.drawLine(480, 75, 490, 75);
        DC.drawString("increasing", 500, 80);
        DC.setColor(new Color(255, 128, 0));
        DC.drawLine(480, 95, 490, 95);
        DC.drawString("decreasing", 500, 100);
        DC.setColor(Color.white);

        lblImage.setBorder(BorderFactory.createEtchedBorder());
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                if (main.yend > main.ystart)
                    lblposn.setText(String.format(" %.4f", main.ystart + e.getX()*(main.yend - main.ystart)/image.getWidth()) + ", " + String.format("%.3f", thmin + e.getY()*(thmax - thmin)/image.getHeight()));
                else
                    lblposn.setText(String.format(" %.4f", main.yend + e.getX()*(main.ystart - main.yend)/image.getWidth()) + ", " + String.format("%.3f", thmin + e.getY()*(thmax - thmin)/image.getHeight()));
            }
        });
        bifurcatePanel.add(lblImage);

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(parmsPanel, BorderLayout.WEST);
        getContentPane().add(bifurcatePanel, BorderLayout.EAST);
        setVisible(true);

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                main.theta0 = Math.PI*Double.parseDouble(txt[1].getText())/180;   // radians
                main.w0 = Double.parseDouble(txt[2].getText());                   // radians/sec
                main.x0 = Double.parseDouble(txt[3].getText());                   // units of R
                main.xa = Double.parseDouble(txt[4].getText());                   // units of R
                main.ystart = Double.parseDouble(txt[5].getText());               // units of R
                main.yend = Double.parseDouble(txt[6].getText());                 // units of R
                main.c = Double.parseDouble(txt[7].getText());                    // moment of inertia
                main.Tx = Double.parseDouble(txt[8].getText());                   // period of x motion
                main.phi0 = Math.PI*Double.parseDouble(txt[9].getText())/180;     // radians
                btnRun.setEnabled(false);
                Bifurcate.activitycancel = false;
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
    }
}

class BifurcateActivity extends SwingWorker<Void, Point>
{
    final int Nper = 100;                       // # of iterations per Tx
    final int NCycle = 256;                     // # of Tx cycles to execute/record per y
    final double delt = main.Tx/Nper;
    int Ninit = 100;                            // # of Tx cycles to initiallize limit cycle (if 0, then use previous run)
    Color clr;
    double tempmin = 7, tempmax = -1;           // range of theta
    Point2D.Double pt = new Point2D.Double(main.theta0, main.w0);   // phase-space point (theta, w)
    double y;                                   // distance to forcing function
    private int itime = 0;

    public BifurcateActivity()
    {
        if (!Bifurcate.initChk.isSelected())    // override the factory default, use last values
            Ninit = 0;
        if (Ninit > 0)                          // use factory default for initiallizing theta, w
            clr = new Color(0, 160, 0);
        else if (main.yend > main.ystart)       // use previous run (forward)
            clr = new Color(40, 40, 136);
        else                                    // use previous run (reverse)
            clr = new Color(255, 128, 0);
    }

    protected Void doInBackground() throws Exception
    {
        for (int i = 0; i < Bifurcate.image.getWidth(); i++)
        {
            y = main.ystart + i*(main.yend - main.ystart)/Bifurcate.image.getWidth();
            if (Ninit > 0)                                  // use factory default
            {
                pt.x = main.theta0;
                pt.y = main.w0;
            }
            //System.out.println(i + ", " + y + ", " + pt.x + ", " + pt.y);
            for (int j = 0; j < NCycle + Ninit; j++)
            {
                //System.out.println(i + ", " + (x0 + xa*Math.cos(phi0)) + ", " + y + ", " + pt.x + ", " + pt.y);
                for (int k = 0; k < Nper; k++)
                {
                    pt = main.runge_kutta(itime, pt.x, pt.y, delt, y);
                    //if (j == NCycle - 1) System.out.print(", " + (int) (100*pt.y));
                    itime++;
                }
                if (pt.x > Bifurcate.thmin && pt.x < Bifurcate.thmax && j > Ninit)
                    publish(new Point(i, (int) (Bifurcate.image.getHeight()*(pt.x - Bifurcate.thmin)/(Bifurcate.thmax - Bifurcate.thmin))));
                if (pt.x > tempmax && j > Ninit) tempmax = pt.x;
                if (pt.x < tempmin && j > Ninit) tempmin = pt.x;
                //System.out.print(", " + (int) (100*pt.y));
            }
            publish(new Point(i, (int) (Bifurcate.image.getHeight() - 1))); // cursor at the bottom
            Bifurcate.lblImage.repaint();
            if (Bifurcate.activitycancel)
                break;
        }
        System.out.println("min->max = " + tempmin + ", " + tempmax);
        return null;
    }

    @Override protected void process(List<Point> points)
    {
        for (Point listpt : points)
        {
            if (main.yend > main.ystart)
                Bifurcate.image.setRGB(listpt.x, listpt.y, clr.getRGB());
            else                                // plot in reverse
                Bifurcate.image.setRGB(Bifurcate.image.getWidth() - 1 - listpt.x, listpt.y, clr.getRGB());
        }
    }

    @Override protected void done()
    {
        Bifurcate.btnRun.setEnabled(true);
        Bifurcate.activitycancel = false;
        Bifurcate.DC.drawLine(0, Bifurcate.image.getHeight() - 1, Bifurcate.image.getWidth() - 1, Bifurcate.image.getHeight() - 1); // clear cursor at the bottom
        Bifurcate.lblImage.repaint();
    }
}
