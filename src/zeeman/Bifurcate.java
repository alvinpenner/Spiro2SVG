
package zeeman;

// plot a bifurcation diagram of theta as fxn of y0
// implement Runge-Kutta: see Froberg page 269
// for animation, see CoreJava v2ch06, ProgressBarTest.java

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import javax.swing.*;

public class Bifurcate extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(600, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static JPanel bifurcatePanel = new JPanel();
    private static JPanel parmsPanel = new JPanel();
    private static JTextField txtA = new JTextField();
    private static JTextField txtTheta0 = new JTextField();
    private static JTextField txtw0 = new JTextField();
    private static JTextField txtx0 = new JTextField();
    private static JTextField txtxa = new JTextField();
    private static JTextField txtymin = new JTextField();
    private static JTextField txtymax = new JTextField();
    private static double theta0;                           // static variables
    protected static double xa, ymin, ymax;                 // static variables
    protected static double c, x0, w0, Tx, phi0;            // dynamic variables
    private double thmin = 1.5, thmax = 4.8;

    public Bifurcate(Image img, double xorg)
    {
        JPanel spacerPanel1 = new JPanel();
        JPanel spacerPanel2 = new JPanel();
        final JTextField txtc = new JTextField();
        final JTextField txtTx = new JTextField();
        final JTextField txtphi0 = new JTextField();
        JButton btnRun = new JButton("Run ");

        setTitle(" Dynamic Zeeman - Bifurcate diagram");
        setIconImage(img);
        setSize(770, 500);
        setLocationByPlatform(true);

        JPanel APanel = new JPanel();
        APanel.setOpaque(false);
        JLabel lblA = new JLabel("A");
        lblA.setPreferredSize(new Dimension(40, 18));
        APanel.add(lblA);
        txtA.setPreferredSize(new Dimension(45, 18));
        txtA.setText(Double.toString(staticComponent.A));
        txtA.setEditable(false);
        APanel.add(txtA);

        JPanel theta0Panel = new JPanel();
        theta0Panel.setOpaque(false);
        JLabel lblTheta0 = new JLabel("θ0");
        lblTheta0.setPreferredSize(new Dimension(40, 18));
        theta0Panel.add(lblTheta0);
        txtTheta0.setPreferredSize(new Dimension(45, 18));
        //txtTheta0.setText(String.format("%.2f", staticComponent.theta*180/Math.PI));
        txtTheta0.setText("0.0");
        theta0Panel.add(txtTheta0);

        JPanel w0Panel = new JPanel();
        w0Panel.setOpaque(false);
        JLabel lblw0 = new JLabel("ω0");
        lblw0.setPreferredSize(new Dimension(40, 18));
        w0Panel.add(lblw0);
        txtw0.setPreferredSize(new Dimension(45, 18));
        txtw0.setText(Double.toString(w0));
        w0Panel.add(txtw0);

        JPanel x0Panel = new JPanel();
        x0Panel.setOpaque(false);
        JLabel lblx0 = new JLabel("x0");
        lblx0.setPreferredSize(new Dimension(40, 18));
        x0Panel.add(lblx0);
        txtx0.setPreferredSize(new Dimension(45, 18));
        txtx0.setText(Double.toString(x0));
        x0Panel.add(txtx0);

        JPanel xaPanel = new JPanel();
        xaPanel.setOpaque(false);
        JLabel lblxa = new JLabel("xa");
        lblxa.setPreferredSize(new Dimension(40, 18));
        xaPanel.add(lblxa);
        txtxa.setPreferredSize(new Dimension(45, 18));
        txtxa.setText(String.format("%.2f", xorg));
        xaPanel.add(txtxa);

        JPanel yminPanel = new JPanel();
        yminPanel.setOpaque(false);
        JLabel lblymin = new JLabel("ymin");
        lblymin.setPreferredSize(new Dimension(40, 18));
        yminPanel.add(lblymin);
        txtymin.setPreferredSize(new Dimension(45, 18));
        txtymin.setText(String.format("%.3f", ymin));
        yminPanel.add(txtymin);

        JPanel ymaxPanel = new JPanel();
        ymaxPanel.setOpaque(false);
        JLabel lblymax = new JLabel("ymax");
        lblymax.setPreferredSize(new Dimension(40, 18));
        ymaxPanel.add(lblymax);
        txtymax.setPreferredSize(new Dimension(45, 18));
        txtymax.setText(String.format("%.3f", ymax));
        ymaxPanel.add(txtymax);

        spacerPanel1.setPreferredSize(new Dimension(85, 6));
        spacerPanel1.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel1.setOpaque(false);

        JPanel cPanel = new JPanel();
        cPanel.setOpaque(false);
        JLabel lblc = new JLabel("c");
        lblc.setPreferredSize(new Dimension(40, 18));
        cPanel.add(lblc);
        txtc.setPreferredSize(new Dimension(45, 18));
        txtc.setText(Double.toString(c));
        cPanel.add(txtc);

        JPanel TxPanel = new JPanel();
        TxPanel.setOpaque(false);
        JLabel lblTx = new JLabel("Tx");
        lblTx.setPreferredSize(new Dimension(40, 18));
        TxPanel.add(lblTx);
        txtTx.setPreferredSize(new Dimension(45, 18));
        txtTx.setText(Double.toString(Tx));
        TxPanel.add(txtTx);

        JPanel phi0Panel = new JPanel();
        phi0Panel.setOpaque(false);
        JLabel lblphi0 = new JLabel("φ0");
        lblphi0.setPreferredSize(new Dimension(40, 18));
        phi0Panel.add(lblphi0);
        txtphi0.setPreferredSize(new Dimension(45, 18));
        txtphi0.setText(String.format("%.2f", phi0*180/Math.PI));
        phi0Panel.add(txtphi0);

        spacerPanel2.setPreferredSize(new Dimension(85, 6));
        spacerPanel2.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel2.setOpaque(false);

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
        parmsPanel.add(APanel);
        parmsPanel.add(theta0Panel);
        parmsPanel.add(w0Panel);
        parmsPanel.add(x0Panel);
        parmsPanel.add(xaPanel);
        parmsPanel.add(yminPanel);
        parmsPanel.add(ymaxPanel);
        parmsPanel.add(spacerPanel1);
        parmsPanel.add(cPanel);
        parmsPanel.add(TxPanel);
        parmsPanel.add(phi0Panel);
        parmsPanel.add(spacerPanel2);
        parmsPanel.add(btnRun);
        parmsPanel.add(thrangePanel);
        parmsPanel.add(posnPanel);
        parmsPanel.setMaximumSize(new Dimension(120, 3000));
        parmsPanel.setPreferredSize(new Dimension(120, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());

        final JLabel lblImage = new JLabel(new ImageIcon(image));
        lblImage.setBorder(BorderFactory.createEtchedBorder());
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                lblposn.setText(String.format(" %.4f", ymin + e.getX()*(ymax - ymin)/image.getWidth()) + ", " + String.format("%.3f", thmin + e.getY()*(thmax - thmin)/image.getHeight()));
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
                theta0 = Math.PI*Double.parseDouble(txtTheta0.getText())/180;   // radians
                w0 = Double.parseDouble(txtw0.getText());                       // radians/sec
                x0 = Double.parseDouble(txtx0.getText());                       // units of R
                xa = Double.parseDouble(txtxa.getText());                       // units of R
                ymin = Double.parseDouble(txtymin.getText());                   // units of R
                ymax = Double.parseDouble(txtymax.getText());                   // units of R
                c = Double.parseDouble(txtc.getText());                         // moment of inertia
                Tx = Double.parseDouble(txtTx.getText());                       // period of x motion
                phi0 = Math.PI*Double.parseDouble(txtphi0.getText())/180;       // radians
                plot_bifurcate();
            }
        });
    }

    private void plot_bifurcate()
    {
        final int Nper = 100;                        // # of iterations per Tx
        final int NCycle = 128;                       // # of Tx cycles to execute
        final double delt = Tx/Nper;
        Point2D.Double pt = new Point2D.Double(theta0, w0);
        Color clr = new Color(40, 40, 136);
        double tempmin = 7, tempmax = -1;
        double y;
        double t = 0;

        //DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        for (int i = 0; i < image.getWidth(); i++)
        {
            y = ymin + i*(ymax - ymin)/image.getWidth();
            //System.out.println(i + ", " + y);
            for (int j = 0; j < NCycle; j++)
            {
                //System.out.println(i + ", " + (x0 + xa*Math.cos(phi0)) + ", " + y + ", " + pt.x + ", " + pt.y);
                for (int k = 0; k < Nper; k++)
                {
                    pt = runge_kutta(t, pt.x, pt.y, delt, y);
                    t += delt;
                }
                if (pt.x > thmin && pt.x < thmax)
                    if (ymax > ymin)
                        image.setRGB(i, (int) (image.getHeight()*(pt.x - thmin)/(thmax - thmin)), clr.getRGB());
                    else                                // plot in reverse
                        image.setRGB(image.getWidth() - 1 - i, (int) (image.getHeight()*(pt.x - thmin)/(thmax - thmin)), clr.getRGB());
                if (pt.x > tempmax) tempmax = pt.x;
                if (pt.x < tempmin) tempmin = pt.x;
            }
            bifurcatePanel.repaint();
        }
        System.out.println("min-max = " + tempmin + ", " + tempmax);
    }

    private static Point2D.Double runge_kutta(double t, double th, double w, double delt, double y)
    {
        double x;
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;                      // fix fix bug bug replace ymin with a variable y

        x = x0 + xa*Math.cos(2*Math.PI*t/Tx + phi0);
        k1 = delt*w;
        l1 = delt*(-c*main.calc_dFdth(th, staticComponent.A, x, y) -w);

        x = x0 + xa*Math.cos(2*Math.PI*(t + delt/2)/Tx + phi0);
        k2 = delt*(w + l1/2);
        l2 = delt*(-c*main.calc_dFdth(th + k1/2, staticComponent.A, x, y) -(w + l1/2));

        k3 = delt*(w + l2/2);
        l3 = delt*(-c*main.calc_dFdth(th + k2/2, staticComponent.A, x, y) -(w + l2/2));

        x = x0 + xa*Math.cos(2*Math.PI*(t + delt)/Tx + phi0);
        k4 = delt*(w + l3);
        l4 = delt*(-c*main.calc_dFdth(th + k3, staticComponent.A, x, y) -(w + l3));

        return new Point2D.Double(th + (k1 + 2*k2 + 2*k3 + k4)/6, w + (l1 + 2*l2 + 2*l3 + l4)/6);
    }
}
