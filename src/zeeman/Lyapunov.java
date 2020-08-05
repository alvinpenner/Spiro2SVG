
package zeeman;

// calculate the largest Lyapunov exponent as fxn of y0
// at each y0, run NCycle cycles of length Tx to initiallize
// then perturb both theta and w
// then run NPerturb cycles of length Tx to determine response

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.List;
import javax.swing.*;

public class Lyapunov extends JDialog
{
    protected static final BufferedImage image = new BufferedImage(600, 400, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    protected static final JLabel lblImage = new JLabel(new ImageIcon(image));
    protected static JPanel lyapunovPanel = new JPanel();
    private static JPanel parmsPanel = new JPanel();
    protected static JButton btnRun = new JButton("Run");
    protected static double lymin = -.4, lymax = .2;

    public Lyapunov(Image img)
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

        setTitle(" Dynamic Zeeman - Lyapunov Exponent vs. y");
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

        JPanel lyrangePanel = new JPanel();
        lyrangePanel.setOpaque(false);
        JLabel lbllyrange = new JLabel("ly = " + lymin + " - " + lymax);
        lbllyrange.setPreferredSize(new Dimension(85, 20));
        lyrangePanel.add(lbllyrange);

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
        parmsPanel.add(lyrangePanel);
        parmsPanel.add(posnPanel);
        parmsPanel.setMaximumSize(new Dimension(125, 3000));
        parmsPanel.setPreferredSize(new Dimension(125, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        DC.setBackground(Color.white);
        DC.setPaint(Color.black);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        DC.draw(new Line2D.Double(0, image.getHeight()*(0 - lymax)/(lymin - lymax),
                   image.getWidth(), image.getHeight()*(0 - lymax)/(lymin - lymax)));
        lblImage.setBorder(BorderFactory.createEtchedBorder());
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                if (main.yend > main.ystart)
                    lblposn.setText(String.format(" %.4f", main.ystart + e.getX()*(main.yend - main.ystart)/image.getWidth()) + ", " + String.format("%.3f", lymax + e.getY()*(lymin - lymax)/image.getHeight()));
                else
                    lblposn.setText(String.format(" %.4f", main.yend + e.getX()*(main.ystart - main.yend)/image.getWidth()) + ", " + String.format("%.3f", lymax + e.getY()*(lymin - lymax)/image.getHeight()));
            }
        });
        lyapunovPanel.add(lblImage);

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(parmsPanel, BorderLayout.WEST);
        getContentPane().add(lyapunovPanel, BorderLayout.EAST);
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
                LyapunovActivity activity = new LyapunovActivity(main.theta0, main.w0);
                activity.execute();
            }
        });
    }
}

class LyapunovActivity extends SwingWorker<Void, Point>
{
    final int Nper = 100;                       // # of iterations per Tx
    final int NCycle = 30;                      // # of Tx cycles before perturb
    final int NPerturb = 15;			// # of Tx cycles after perturb
    final double delt = main.Tx/Nper;
    final Color clr = new Color(40, 40, 136);
    double tempmin = 10, tempmax = -10;         // range of lyapunov exponent
    Point2D.Double pt;                          // phase-space point (theta, w)
    double y;                                   // distance to forcing function
    private int itime = 0;                      // time index

    public LyapunovActivity(double passtheta0, double passw0)
    {
        pt = new Point2D.Double(passtheta0, passw0);
    }

    protected Void doInBackground() throws Exception
    {
        Point2D.Double ptN;                     // phase-space point after NCycle
        Point2D.Double ptN0;                    // unperturbed     after extra NPerturb cycles
        Point2D.Double ptNth;                   // theta perturbed after extra NPerturb cycles
        Point2D.Double ptNw;                    // w perturbed     after extra NPerturb cycles
        double deltheta = 0.01;
        double delw = 0.01;
        double delin = Math.sqrt(deltheta*deltheta + delw*delw);
        double delout, exp;
        for (int i = 0; i < Lyapunov.image.getWidth(); i++)
        {
            y = main.ystart + i*(main.yend - main.ystart)/Lyapunov.image.getWidth();
            for (int j = 0; j < NCycle; j++)
            {
                for (int k = 0; k < Nper; k++)
                {
                    pt = main.runge_kutta(itime, pt.x, pt.y, delt, y);
                    itime++;
                }
            }
            ptN = new Point2D.Double(pt.x, pt.y);       // unperturbed pt after NCycle

            for (int j = 0; j < NPerturb; j++)
                for (int k = 0; k < Nper; k++)
                {
                    pt = main.runge_kutta(itime, pt.x, pt.y, delt, y);
                    itime++;
                }
            ptN0 = new Point2D.Double(pt.x, pt.y);      // unperturbed after NPerturb

            pt.x = ptN.x + deltheta;                    // perturb theta
            pt.y = ptN.y;
            for (int j = 0; j < NPerturb; j++)
                for (int k = 0; k < Nper; k++)
                {
                    pt = main.runge_kutta(itime, pt.x, pt.y, delt, y);
                    itime++;
                }
            ptNth = new Point2D.Double(pt.x, pt.y);     // theta perturbed after NPerturb

            pt.x = ptN.x;
            pt.y = ptN.y + delw;                        // perturb w
            for (int j = 0; j < NPerturb; j++)
                for (int k = 0; k < Nper; k++)
                {
                    pt = main.runge_kutta(itime, pt.x, pt.y, delt, y);
                    itime++;
                }
            ptNw = new Point2D.Double(pt.x, pt.y);      // w perturbed after NPerturb
            //System.out.println(i + ", " + y + ", " + ptN.x + ", " + ptN.y + ", " + ptN0.x + ", " + ptN0.y + ", " + ptNth.x + ", " + ptNth.y + ", " + ptNw.x + ", " + ptNw.y);
            delout = Math.sqrt((ptN0.x - ptNth.x)*(ptN0.x - ptNth.x) + (ptN0.y - ptNth.y)*(ptN0.y - ptNth.y)
                             + (ptN0.x - ptNw.x)*(ptN0.x - ptNw.x) + (ptN0.y - ptNw.y)*(ptN0.y - ptNw.y));
            exp = Math.log(delout/delin)/NPerturb/main.Tx;
            if (exp > Lyapunov.lymin && exp < Lyapunov.lymax)
                publish(new Point(i, (int) (Lyapunov.image.getHeight()*(exp - Lyapunov.lymax)/(Lyapunov.lymin - Lyapunov.lymax))));
            if (exp > tempmax) tempmax = exp;
            if (exp < tempmin) tempmin = exp;
            Lyapunov.lblImage.repaint();

            pt.x = ptN.x;                               // revert to common pt
            pt.y = ptN.y;
        }
        System.out.println("exp min->max = " + tempmin + ", " + tempmax);
        return null;
    }

    @Override protected void process(List<Point> points)
    {
        for (Point listpt : points)
        {
            if (main.yend > main.ystart)
                Lyapunov.image.setRGB(listpt.x, listpt.y, clr.getRGB());
            else                                // plot in reverse
                Lyapunov.image.setRGB(Lyapunov.image.getWidth() - 1 - listpt.x, listpt.y, clr.getRGB());
        }
    }

    @Override protected void done()
    {
        Lyapunov.btnRun.setEnabled(true);
    }
}
