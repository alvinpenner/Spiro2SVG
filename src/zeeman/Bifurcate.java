
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
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    protected static JPanel bifurcatePanel = new JPanel();
    private static JPanel parmsPanel = new JPanel();
    protected static JButton btnRun = new JButton("Run");
    protected static JButton btnClear = new JButton("Clr");
    protected static double thmin = 1.5, thmax = 4.8;

    public Bifurcate(Image img)
    {
        JPanel spacerPanel1 = new JPanel();
        JPanel spacerPanel2 = new JPanel();
        final JTextField txtA = new JTextField();
        final JTextField txtTheta0 = new JTextField();
        final JTextField txtw0 = new JTextField();
        final JTextField txtx0 = new JTextField();
        final JTextField txtxa = new JTextField();
        final JTextField txtc = new JTextField();
        final JTextField txtTx = new JTextField();
        final JTextField txtphi0 = new JTextField();
        final JTextField txtystart = new JTextField();
        final JTextField txtyend = new JTextField();

        setTitle(" Dynamic Zeeman - Bifurcate diagram (θ vs. y)");
        setIconImage(img);
        setSize(755, 500);
        setLocationByPlatform(true);

        JPanel APanel = new JPanel();
        APanel.setOpaque(false);
        JLabel lblA = new JLabel("A");
        lblA.setPreferredSize(new Dimension(35, 18));
        APanel.add(lblA);
        txtA.setPreferredSize(new Dimension(50, 18));
        txtA.setText(Double.toString(main.A));
        txtA.setEditable(false);
        APanel.add(txtA);

        JPanel theta0Panel = new JPanel();
        theta0Panel.setOpaque(false);
        JLabel lblTheta0 = new JLabel("θ0");
        lblTheta0.setPreferredSize(new Dimension(35, 18));
        theta0Panel.add(lblTheta0);
        txtTheta0.setPreferredSize(new Dimension(50, 18));
        txtTheta0.setText(String.format("%.2f", main.theta0*180/Math.PI));
        theta0Panel.add(txtTheta0);

        JPanel w0Panel = new JPanel();
        w0Panel.setOpaque(false);
        JLabel lblw0 = new JLabel("ω0");
        lblw0.setPreferredSize(new Dimension(35, 18));
        w0Panel.add(lblw0);
        txtw0.setPreferredSize(new Dimension(50, 18));
        txtw0.setText(Double.toString(main.w0));
        w0Panel.add(txtw0);

        JPanel x0Panel = new JPanel();
        x0Panel.setOpaque(false);
        JLabel lblx0 = new JLabel("x0");
        lblx0.setPreferredSize(new Dimension(35, 18));
        x0Panel.add(lblx0);
        txtx0.setPreferredSize(new Dimension(50, 18));
        txtx0.setText(Double.toString(main.x0));
        x0Panel.add(txtx0);

        JPanel xaPanel = new JPanel();
        xaPanel.setOpaque(false);
        JLabel lblxa = new JLabel("xa");
        lblxa.setPreferredSize(new Dimension(35, 18));
        xaPanel.add(lblxa);
        txtxa.setPreferredSize(new Dimension(50, 18));
        txtxa.setText(String.format("%.2f", main.xa));
        xaPanel.add(txtxa);

        JPanel ystartPanel = new JPanel();
        ystartPanel.setOpaque(false);
        JLabel lblystart = new JLabel("ystart");
        lblystart.setPreferredSize(new Dimension(35, 18));
        ystartPanel.add(lblystart);
        txtystart.setPreferredSize(new Dimension(50, 18));
        txtystart.setText(String.format("%.3f", main.ystart));
        ystartPanel.add(txtystart);

        JPanel yendPanel = new JPanel();
        yendPanel.setOpaque(false);
        JLabel lblyend = new JLabel("yend");
        lblyend.setPreferredSize(new Dimension(35, 18));
        yendPanel.add(lblyend);
        txtyend.setPreferredSize(new Dimension(50, 18));
        txtyend.setText(String.format("%.3f", main.yend));
        yendPanel.add(txtyend);

        spacerPanel1.setPreferredSize(new Dimension(85, 6));
        spacerPanel1.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel1.setOpaque(false);

        JPanel cPanel = new JPanel();
        cPanel.setOpaque(false);
        JLabel lblc = new JLabel("c");
        lblc.setPreferredSize(new Dimension(25, 18));
        cPanel.add(lblc);
        txtc.setPreferredSize(new Dimension(60, 18));
        txtc.setText(Double.toString(main.c));
        cPanel.add(txtc);

        JPanel TxPanel = new JPanel();
        TxPanel.setOpaque(false);
        JLabel lblTx = new JLabel("Tx");
        lblTx.setPreferredSize(new Dimension(25, 18));
        TxPanel.add(lblTx);
        txtTx.setPreferredSize(new Dimension(60, 18));
        txtTx.setText(Double.toString(main.Tx));
        TxPanel.add(txtTx);

        JPanel phi0Panel = new JPanel();
        phi0Panel.setOpaque(false);
        JLabel lblphi0 = new JLabel("φ0");
        lblphi0.setPreferredSize(new Dimension(25, 18));
        phi0Panel.add(lblphi0);
        txtphi0.setPreferredSize(new Dimension(60, 18));
        txtphi0.setText(String.format("%.2f", main.phi0*180/Math.PI));
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
        parmsPanel.add(ystartPanel);
        parmsPanel.add(yendPanel);
        parmsPanel.add(spacerPanel1);
        parmsPanel.add(cPanel);
        parmsPanel.add(TxPanel);
        parmsPanel.add(phi0Panel);
        parmsPanel.add(spacerPanel2);
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
        DC.setColor(new Color(40, 40, 136));
        DC.drawLine(480, 55, 490, 55);
        DC.drawString("increasing", 500, 60);
        DC.setColor(new Color(255, 128, 0));
        DC.drawLine(480, 75, 490, 75);
        DC.drawString("decreasing", 500, 80);

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
                main.theta0 = Math.PI*Double.parseDouble(txtTheta0.getText())/180;   // radians
                main.w0 = Double.parseDouble(txtw0.getText());                       // radians/sec
                main.x0 = Double.parseDouble(txtx0.getText());                       // units of R
                main.xa = Double.parseDouble(txtxa.getText());                       // units of R
                main.c = Double.parseDouble(txtc.getText());                         // moment of inertia
                main.Tx = Double.parseDouble(txtTx.getText());                       // period of x motion
                main.phi0 = Math.PI*Double.parseDouble(txtphi0.getText())/180;       // radians
                main.ystart = Double.parseDouble(txtystart.getText());               // units of R
                main.yend = Double.parseDouble(txtyend.getText());                   // units of R
                btnRun.setEnabled(false);
                BifurcateActivity activity = new BifurcateActivity(main.theta0, main.w0);
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
        //addKeyListener(new KeyAdapter()
        //{
        //    @Override public void keyPressed(KeyEvent e)
        //    {
        //        System.out.println("bifurcate key event");
        //    }
        //});
    }
}

class BifurcateActivity extends SwingWorker<Void, Point>
{
    final int Nper = 100;                       // # of iterations per Tx
    final int NCycle = 256;                     // # of Tx cycles to execute per y
    final double delt = main.Tx/Nper;
    Color clr;
    double tempmin = 7, tempmax = -1;           // range of theta
    Point2D.Double pt;                          // phase-space point (theta, w)
    double y;                                   // distance to forcing function
    private int itime = 0;

    public BifurcateActivity(double passtheta0, double passw0)
    {
        pt = new Point2D.Double(passtheta0, passw0);
        if (main.yend > main.ystart)
            clr = new Color(40, 40, 136);
        else
            clr = new Color(255, 128, 0);
    }

    protected Void doInBackground() throws Exception
    {
        for (int i = 0; i < Bifurcate.image.getWidth(); i++)
        {
            y = main.ystart + i*(main.yend - main.ystart)/Bifurcate.image.getWidth();
            for (int j = 0; j < NCycle; j++)
            {
                //System.out.println(i + ", " + (x0 + xa*Math.cos(phi0)) + ", " + y + ", " + pt.x + ", " + pt.y);
                for (int k = 0; k < Nper; k++)
                {
                    pt = main.runge_kutta(itime, pt.x, pt.y, delt, y);
                    itime++;
                }
                if (pt.x > Bifurcate.thmin && pt.x < Bifurcate.thmax)
                    publish(new Point(i, (int) (Bifurcate.image.getHeight()*(pt.x - Bifurcate.thmin)/(Bifurcate.thmax - Bifurcate.thmin))));
                if (pt.x > tempmax) tempmax = pt.x;
                if (pt.x < tempmin) tempmin = pt.x;
            }
            Bifurcate.lblImage.repaint();
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
    }
}
