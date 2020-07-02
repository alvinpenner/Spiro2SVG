
package zeeman;

// plot w versus theta for dynamic simulation
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Path2D;
import javax.swing.*;
import java.util.GregorianCalendar;

public class PhaseSpace extends JDialog
{
    private final static int N = 21000;                 // total # of iterations
    private final static int Nper = 100;                // # of iterations per Tx
    protected static int NLimit;                        // # of Tx per limit cycle (to be determined)
    protected Path2D.Double path1 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Path2D.Double path2 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Path2D.Double path3 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected static JPanel phasePanel = new Plot_Phase_Panel();
    private static JPanel parmsPanel = new JPanel();
    private static JTextField txtA = new JTextField();
    private static JTextField txtTheta0 = new JTextField();
    private static JTextField txtw0 = new JTextField();
    private static JTextField txtx0 = new JTextField();
    private static JTextField txtxa = new JTextField();
    private static JTextField txty0 = new JTextField();
    private static JCheckBox complementChk = new JCheckBox("complement", true);
    private static JCheckBox printChk = new JCheckBox("print");
    private static JLabel lblThetarange;
    private static JLabel lblwrange;
    private static double theta0;                           // static variables
    protected static double xa, y0;                         // static variables
    protected static double c, x0, w0, Tx, phi0;            // dynamic variables

    public PhaseSpace(Image img, double xorg, double yorg)
    {
        //phasePanel.setBackground(Color.white);
        JPanel spacerPanel1 = new JPanel();
        JPanel spacerPanel2 = new JPanel();
        JPanel spacerPanel3 = new JPanel();
        final JTextField txtc = new JTextField();
        final JTextField txtTx = new JTextField();
        final JTextField txtphi0 = new JTextField();
        JButton btnRun = new JButton("Run");
        String[] NLimitdata = {" 0", " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", " 10"};
        final JComboBox NLimitCombo = new JComboBox(NLimitdata);

        setTitle(" Dynamic Zeeman - w vs. theta");
        setIconImage(img);
        setSize(680, 580);
        setLocationByPlatform(true);

        JPanel APanel = new JPanel();
        APanel.setOpaque(false);
        JLabel lblA = new JLabel("A");
        lblA.setPreferredSize(new Dimension(25, 18));
        APanel.add(lblA);
        txtA.setPreferredSize(new Dimension(55, 18));
        txtA.setText(Double.toString(staticComponent.A));
        txtA.setEditable(false);
        APanel.add(txtA);

        JPanel theta0Panel = new JPanel();
        theta0Panel.setOpaque(false);
        JLabel lblTheta0 = new JLabel("θ0");
        lblTheta0.setPreferredSize(new Dimension(25, 18));
        theta0Panel.add(lblTheta0);
        txtTheta0.setPreferredSize(new Dimension(55, 18));
        //txtTheta0.setText(String.format("%.2f", staticComponent.theta*180/Math.PI));
        txtTheta0.setText("0.0");
        theta0Panel.add(txtTheta0);

        JPanel w0Panel = new JPanel();
        w0Panel.setOpaque(false);
        JLabel lblw0 = new JLabel("ω0");
        lblw0.setPreferredSize(new Dimension(25, 18));
        w0Panel.add(lblw0);
        txtw0.setPreferredSize(new Dimension(55, 18));
        txtw0.setText(Double.toString(w0));
        w0Panel.add(txtw0);

        JPanel x0Panel = new JPanel();
        x0Panel.setOpaque(false);
        JLabel lblx0 = new JLabel("x0");
        lblx0.setPreferredSize(new Dimension(25, 18));
        x0Panel.add(lblx0);
        txtx0.setPreferredSize(new Dimension(55, 18));
        txtx0.setText(Double.toString(x0));
        x0Panel.add(txtx0);

        JPanel xaPanel = new JPanel();
        xaPanel.setOpaque(false);
        JLabel lblxa = new JLabel("xa");
        lblxa.setPreferredSize(new Dimension(25, 18));
        xaPanel.add(lblxa);
        txtxa.setPreferredSize(new Dimension(55, 18));
        txtxa.setText(String.format("%.2f", xorg));
        xaPanel.add(txtxa);

        JPanel y0Panel = new JPanel();
        y0Panel.setOpaque(false);
        JLabel lbly0 = new JLabel("y0");
        lbly0.setPreferredSize(new Dimension(25, 18));
        y0Panel.add(lbly0);
        txty0.setPreferredSize(new Dimension(55, 18));
        txty0.setText(String.format("%.5f", yorg));
        y0Panel.add(txty0);

        spacerPanel1.setPreferredSize(new Dimension(80, 6));
        spacerPanel1.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel1.setOpaque(false);

        JPanel cPanel = new JPanel();
        cPanel.setOpaque(false);
        JLabel lblc = new JLabel("c");
        lblc.setPreferredSize(new Dimension(25, 18));
        cPanel.add(lblc);
        txtc.setPreferredSize(new Dimension(55, 18));
        txtc.setText(Double.toString(c));
        cPanel.add(txtc);

        JPanel TxPanel = new JPanel();
        TxPanel.setOpaque(false);
        JLabel lblTx = new JLabel("Tx");
        lblTx.setPreferredSize(new Dimension(25, 18));
        TxPanel.add(lblTx);
        txtTx.setPreferredSize(new Dimension(55, 18));
        txtTx.setText(Double.toString(Tx));
        TxPanel.add(txtTx);

        JPanel phi0Panel = new JPanel();
        phi0Panel.setOpaque(false);
        JLabel lblphi0 = new JLabel("φ0");
        lblphi0.setPreferredSize(new Dimension(25, 18));
        phi0Panel.add(lblphi0);
        txtphi0.setPreferredSize(new Dimension(55, 18));
        txtphi0.setText(String.format("%.2f", phi0*180/Math.PI));
        phi0Panel.add(txtphi0);

        spacerPanel2.setPreferredSize(new Dimension(80, 6));
        spacerPanel2.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel2.setOpaque(false);

        spacerPanel3.setPreferredSize(new Dimension(80, 6));
        spacerPanel3.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel3.setOpaque(false);

        JPanel thetarangePanel = new JPanel();
        thetarangePanel.setOpaque(false);
        lblThetarange = new JLabel("θ = " + String.format("%.1f", 0.0) + " - " + String.format("%.1f", 360.0));
        lblThetarange.setPreferredSize(new Dimension(95, 18));
        thetarangePanel.add(lblThetarange);

        JPanel wrangePanel = new JPanel();
        wrangePanel.setOpaque(false);
        lblwrange = new JLabel("ω = " + String.format("%.2f", -1.0) + " - " + String.format("%.2f", 1.0));
        lblwrange.setPreferredSize(new Dimension(95, 18));
        wrangePanel.add(lblwrange);

        JPanel NLimitPanel = new JPanel();
        NLimitPanel.setOpaque(false);
        JLabel lblLimit = new JLabel("# limit ");
        NLimitPanel.add (lblLimit);
        NLimitCombo.setSelectedIndex(NLimit);
        NLimitPanel.add(NLimitCombo);

        complementChk.setOpaque(false);
        printChk.setOpaque(false);

        parmsPanel.removeAll();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(APanel);
        parmsPanel.add(theta0Panel);
        parmsPanel.add(w0Panel);
        parmsPanel.add(x0Panel);
        parmsPanel.add(xaPanel);
        parmsPanel.add(y0Panel);
        parmsPanel.add(spacerPanel1);
        parmsPanel.add(cPanel);
        parmsPanel.add(TxPanel);
        parmsPanel.add(phi0Panel);
        parmsPanel.add(spacerPanel2);
        parmsPanel.add(btnRun);
        parmsPanel.add(NLimitPanel);
        parmsPanel.add(complementChk);
        parmsPanel.add(printChk);
        parmsPanel.add(spacerPanel3);
        parmsPanel.add(thetarangePanel);
        parmsPanel.add(wrangePanel);
        parmsPanel.setMaximumSize(new Dimension(120, 3000));
        parmsPanel.setPreferredSize(new Dimension(120, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().add(parmsPanel);
        getContentPane().add(phasePanel);
        setVisible(true);

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                theta0 = Math.PI*Double.parseDouble(txtTheta0.getText())/180;   // radians
                w0 = Double.parseDouble(txtw0.getText());                       // radians/sec
                x0 = Double.parseDouble(txtx0.getText());                       // units of R
                xa = Double.parseDouble(txtxa.getText());                       // units of R
                y0 = Double.parseDouble(txty0.getText());                       // units of R
                c = Double.parseDouble(txtc.getText());                         // moment of inertia
                Tx = Double.parseDouble(txtTx.getText());                       // period of x motion
                phi0 = Math.PI*Double.parseDouble(txtphi0.getText())/180;       // radians
                NLimit = NLimitCombo.getSelectedIndex();                        // number of Tx in limit cycle
                phase_space();
            }
        });
    }

    private void phase_space()
    {
        final double delt = Tx/Nper;
        Point2D.Double pt = new Point2D.Double(theta0, w0);
        double thmin = theta0, thmax = theta0;
        double wmin = w0, wmax = w0;
        double tempx;

        if (printChk.isSelected())
        {
            GregorianCalendar now = new GregorianCalendar();
            System.out.println("\nDynamic Zeeman: " + now.getTime());
            System.out.printf("A, %.1f\ntheta0, %f\nw0, %.2f\nx0, %.2f\nxa, %.2f\ny0, %.5f\nc, %f\nTx, %.2f\nphi0, %f\n\n", staticComponent.A, theta0, w0, x0, xa, y0, c, Tx, phi0);
            System.out.println("t, x, theta, w, rhs");
            tempx = x0 + xa*Math.cos(phi0);
            System.out.println("0.0, " + tempx + ", " + pt.x + ", " + pt.y + ", " + (-c*main.calc_dFdth(pt.x, staticComponent.A, tempx, y0) - pt.y));
        }
        path1.reset();              // transient path
        path2.reset();              // limit cycle
        path1.moveTo(pt.x, pt.y);
        System.out.println(" , , " + y0);
        for (int i = 0; i < N; i++)
        {
            if ((i/Nper)*Nper == i)
                System.out.println(i + ", " + (x0 + xa*Math.cos(phi0)) + ", " + pt.x + ", " + pt.y);
            pt = runge_kutta(i*delt, pt.x, pt.y, delt);
            if (i == N - NLimit*Nper - 1)
            {
                path2.moveTo(pt.x, pt.y);
                if (NLimit > 0)
                {
                    thmax = pt.x;
                    thmin = pt.x;
                    wmax = pt.y;
                    wmin = pt.y;
                }
            }
            if (i < N - NLimit*Nper)
                path1.lineTo(pt.x, pt.y);
            else
                path2.lineTo(pt.x, pt.y);
            if (pt.x > thmax) thmax = pt.x;
            if (pt.x < thmin) thmin = pt.x;
            if (pt.y > wmax) wmax = pt.y;
            if (pt.y < wmin) wmin = pt.y;
            tempx = x0 + xa*Math.cos(2*Math.PI*(i + 1)*delt/Tx + phi0);
            if (printChk.isSelected()) System.out.println((i + 1)*delt + ", " + tempx + ", " + pt.x + ", " + pt.y + ", " + (-c*main.calc_dFdth(pt.x, staticComponent.A, tempx, y0) - pt.y));
        }
        if (2*Math.PI - thmax < thmin) thmin = 2*Math.PI - thmax;
        if (2*Math.PI - thmin > thmax) thmax = 2*Math.PI - thmin;
        if (-wmax < wmin) wmin = -wmax;
        if (-wmin > wmax) wmax = -wmin;
        AffineTransform at = new AffineTransform(phasePanel.getWidth()/(thmax - thmin),
                                          0, 0, -phasePanel.getHeight()/(wmax - wmin),
                                                -thmin*phasePanel.getWidth()/(thmax - thmin),
                                                 wmax*phasePanel.getHeight()/(wmax - wmin));
        path1.transform(at);
        path2.transform(at);
        path3.reset();                              // complementary limit cycle
        if (complementChk.isSelected() && NLimit > 0)
        {
            pt.x = 2*Math.PI - pt.x;
            pt.y = -pt.y;
            path3.moveTo(pt.x, pt.y);               // use the endpoint of path2 as a start point
            if (printChk.isSelected())
            {
                System.out.println("complement");
                System.out.println("0.0, " + ", " + pt.x + ", " + pt.y);
            }
            for (int i = 0; i < NLimit*Nper; i++)
            {
                pt = runge_kutta(Tx/2 + i*delt, pt.x, pt.y, delt);
                path3.lineTo(pt.x, pt.y);
                if (printChk.isSelected()) System.out.println((i + 1)*delt + ", " + ", " + pt.x + ", " + pt.y);
            }
            path3.transform(at);
        }
        phasePanel.repaint();
        lblThetarange.setText("θ = " + String.format("%.1f", 180*thmin/Math.PI) + " - " + String.format("%.1f", 180*thmax/Math.PI));
        lblwrange.setText("ω = " + String.format("%.2f", wmin) + " - " + String.format("%.2f", wmax));
    }

    private static Point2D.Double runge_kutta(double t, double th, double w, double delt)
    {
        double x;
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;

        x = x0 + xa*Math.cos(2*Math.PI*t/Tx + phi0);
        k1 = delt*w;
        l1 = delt*(-c*main.calc_dFdth(th, staticComponent.A, x, y0) -w);

        x = x0 + xa*Math.cos(2*Math.PI*(t + delt/2)/Tx + phi0);
        k2 = delt*(w + l1/2);
        l2 = delt*(-c*main.calc_dFdth(th + k1/2, staticComponent.A, x, y0) -(w + l1/2));

        k3 = delt*(w + l2/2);
        l3 = delt*(-c*main.calc_dFdth(th + k2/2, staticComponent.A, x, y0) -(w + l2/2));

        x = x0 + xa*Math.cos(2*Math.PI*(t + delt)/Tx + phi0);
        k4 = delt*(w + l3);
        l4 = delt*(-c*main.calc_dFdth(th + k3, staticComponent.A, x, y0) -(w + l3));

        return new Point2D.Double(th + (k1 + 2*k2 + 2*k3 + k4)/6, w + (l1 + 2*l2 + 2*l3 + l4)/6);
    }
}

class Plot_Phase_Panel extends JPanel
{
    @Override public void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;

        g2.setPaint(new Color(255, 127, 39));
        g2.draw(staticComponent.phase_dlg.path1);
        g2.setPaint(new Color(0, 0, 0));
        g2.setStroke(new BasicStroke(2));
        g2.draw(staticComponent.phase_dlg.path2);
        g2.setPaint(new Color(0, 224, 0));
        g2.setStroke(new BasicStroke(1));
        g2.draw(staticComponent.phase_dlg.path3);
    }
}
