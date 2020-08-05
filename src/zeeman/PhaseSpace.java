
package zeeman;

// plot w versus theta for dynamic simulation
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Path2D;
import javax.swing.*;
import java.util.GregorianCalendar;

public class PhaseSpace extends JDialog
{
    private final static int N = 159600;             // total # of iterations
    private final static int Nper = 100;            // # of iterations per Tx (assume even)
    protected Path2D.Double path1 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Path2D.Double path2 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Path2D.Double path3 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected static JPanel phasePanel = new Plot_Phase_Panel();
    protected static double wmin, wmax;
    private static double[] finalphase;
    private static JPanel parmsPanel = new JPanel();
    private static JCheckBox complementChk = new JCheckBox("complement", true);
    private static JCheckBox printChk = new JCheckBox("print");
    private static JCheckBox ddyChk = new JCheckBox("d/dy of phase");
    private static JLabel lblThetarange;
    private static JLabel lblwrange;

    public PhaseSpace(Image img)
    {
        //phasePanel.setBackground(Color.white);
        JPanel spacerPanel1 = new JPanel();
        JPanel spacerPanel2 = new JPanel();
        JPanel spacerPanel3 = new JPanel();
        final JTextField txtA = new JTextField();
        final JTextField txtTheta0 = new JTextField();
        final JTextField txtw0 = new JTextField();
        final JTextField txtx0 = new JTextField();
        final JTextField txtxa = new JTextField();
        final JTextField txty0 = new JTextField();
        final JTextField txtc = new JTextField();
        final JTextField txtTx = new JTextField();
        final JTextField txtphi0 = new JTextField();
        final JTextField txtdTheta0 = new JTextField();
        final JTextField txtdw0 = new JTextField();
        JButton btnRun = new JButton("Run");
        String[] NLimitdata = {" 0", " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", " 10"};
        final JComboBox NLimitCombo = new JComboBox(NLimitdata);

        if (!ddyChk.isSelected())
            setTitle(" Dynamic Zeeman - ω vs. θ");
        else
            setTitle(" Dynamic Zeeman - dω/dy vs. dθ/dy");
        setIconImage(img);
        setSize(786, 639);
        setLocationByPlatform(true);

        JPanel APanel = new JPanel();
        APanel.setOpaque(false);
        JLabel lblA = new JLabel("A");
        lblA.setPreferredSize(new Dimension(40, 18));
        APanel.add(lblA);
        txtA.setPreferredSize(new Dimension(70, 18));
        txtA.setText(Double.toString(main.A));
        txtA.setEditable(false);
        APanel.add(txtA);

        JPanel theta0Panel = new JPanel();
        theta0Panel.setOpaque(false);
        JLabel lblTheta0 = new JLabel("θ0");
        lblTheta0.setPreferredSize(new Dimension(40, 18));
        theta0Panel.add(lblTheta0);
        txtTheta0.setPreferredSize(new Dimension(70, 18));
        txtTheta0.setText(String.format("%.4f", main.theta0*180/Math.PI));
        theta0Panel.add(txtTheta0);

        JPanel w0Panel = new JPanel();
        w0Panel.setOpaque(false);
        JLabel lblw0 = new JLabel("ω0");
        lblw0.setPreferredSize(new Dimension(40, 18));
        w0Panel.add(lblw0);
        txtw0.setPreferredSize(new Dimension(70, 18));
        txtw0.setText(Double.toString(main.w0));
        w0Panel.add(txtw0);

        JPanel x0Panel = new JPanel();
        x0Panel.setOpaque(false);
        JLabel lblx0 = new JLabel("x0");
        lblx0.setPreferredSize(new Dimension(40, 18));
        x0Panel.add(lblx0);
        txtx0.setPreferredSize(new Dimension(70, 18));
        txtx0.setText(Double.toString(main.x0));
        x0Panel.add(txtx0);

        JPanel xaPanel = new JPanel();
        xaPanel.setOpaque(false);
        JLabel lblxa = new JLabel("xa");
        lblxa.setPreferredSize(new Dimension(40, 18));
        xaPanel.add(lblxa);
        txtxa.setPreferredSize(new Dimension(70, 18));
        txtxa.setText(String.format("%.2f", main.xa));
        xaPanel.add(txtxa);

        JPanel y0Panel = new JPanel();
        y0Panel.setOpaque(false);
        JLabel lbly0 = new JLabel("y0");
        lbly0.setPreferredSize(new Dimension(40, 18));
        y0Panel.add(lbly0);
        txty0.setPreferredSize(new Dimension(70, 18));
        txty0.setText(String.format("%.6f", main.y0));
        y0Panel.add(txty0);

        spacerPanel1.setPreferredSize(new Dimension(110, 6));
        spacerPanel1.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel1.setOpaque(false);

        JPanel cPanel = new JPanel();
        cPanel.setOpaque(false);
        JLabel lblc = new JLabel("c");
        lblc.setPreferredSize(new Dimension(40, 18));
        cPanel.add(lblc);
        txtc.setPreferredSize(new Dimension(70, 18));
        txtc.setText(Double.toString(main.c));
        cPanel.add(txtc);

        JPanel TxPanel = new JPanel();
        TxPanel.setOpaque(false);
        JLabel lblTx = new JLabel("Tx");
        lblTx.setPreferredSize(new Dimension(40, 18));
        TxPanel.add(lblTx);
        txtTx.setPreferredSize(new Dimension(70, 18));
        txtTx.setText(Double.toString(main.Tx));
        TxPanel.add(txtTx);

        JPanel phi0Panel = new JPanel();
        phi0Panel.setOpaque(false);
        JLabel lblphi0 = new JLabel("φ0");
        lblphi0.setPreferredSize(new Dimension(40, 18));
        phi0Panel.add(lblphi0);
        txtphi0.setPreferredSize(new Dimension(70, 18));
        txtphi0.setText(String.format("%.1f", main.phi0*180/Math.PI));
        phi0Panel.add(txtphi0);

        JPanel dtheta0Panel = new JPanel();
        dtheta0Panel.setOpaque(false);
        JLabel lbldTheta0 = new JLabel("dθ0dy");
        lbldTheta0.setPreferredSize(new Dimension(40, 18));
        dtheta0Panel.add(lbldTheta0);
        txtdTheta0.setPreferredSize(new Dimension(70, 18));
        txtdTheta0.setText(Double.toString(main.dtheta0dy));
        dtheta0Panel.add(txtdTheta0);

        JPanel dw0Panel = new JPanel();
        dw0Panel.setOpaque(false);
        JLabel lbldw0 = new JLabel("dω0dy");
        lbldw0.setPreferredSize(new Dimension(40, 18));
        dw0Panel.add(lbldw0);
        txtdw0.setPreferredSize(new Dimension(70, 18));
        txtdw0.setText(Double.toString(main.dw0dy));
        dw0Panel.add(txtdw0);

        spacerPanel2.setPreferredSize(new Dimension(110, 6));
        spacerPanel2.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel2.setOpaque(false);

        spacerPanel3.setPreferredSize(new Dimension(110, 6));
        spacerPanel3.setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
        spacerPanel3.setOpaque(false);

        JPanel thetarangePanel = new JPanel();
        thetarangePanel.setOpaque(false);
        lblThetarange = new JLabel("θ = " + String.format("%.1f", 0.0) + ", " + String.format("%.1f", 360.0));
        lblThetarange.setPreferredSize(new Dimension(150, 18));
        thetarangePanel.add(lblThetarange);

        JPanel wrangePanel = new JPanel();
        wrangePanel.setOpaque(false);
        lblwrange = new JLabel("ω = " + String.format("%.2f", -1.0) + ", " + String.format("%.2f", 1.0));
        lblwrange.setPreferredSize(new Dimension(150, 18));
        wrangePanel.add(lblwrange);

        JPanel NLimitPanel = new JPanel();
        NLimitPanel.setOpaque(false);
        JLabel lblLimit = new JLabel("# limit ");
        lblLimit.setPreferredSize(new Dimension(60, 18));
        NLimitPanel.add (lblLimit);
        NLimitCombo.setSelectedIndex(main.NLimit);
        NLimitPanel.add(NLimitCombo);

        complementChk.setOpaque(false);
        printChk.setOpaque(false);
        ddyChk.setOpaque(false);

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
        parmsPanel.add(dtheta0Panel);
        parmsPanel.add(dw0Panel);
        parmsPanel.add(spacerPanel2);
        parmsPanel.add(btnRun);
        parmsPanel.add(NLimitPanel);
        parmsPanel.add(complementChk);
        parmsPanel.add(printChk);
        parmsPanel.add(ddyChk);
        parmsPanel.add(spacerPanel3);
        parmsPanel.add(thetarangePanel);
        parmsPanel.add(wrangePanel);
        parmsPanel.setMaximumSize(new Dimension(170, 3000));
        parmsPanel.setPreferredSize(new Dimension(170, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().add(parmsPanel);
        getContentPane().add(phasePanel);
        setVisible(true);

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (!ddyChk.isSelected())
                    setTitle(" Dynamic Zeeman - ω vs. θ");
                else
                    setTitle(" Dynamic Zeeman - dω/dy vs. dθ/dy");
                main.theta0 = Math.PI*Double.parseDouble(txtTheta0.getText())/180;   // radians
                main.w0 = Double.parseDouble(txtw0.getText());                       // radians/sec
                main.x0 = Double.parseDouble(txtx0.getText());                       // units of R
                main.xa = Double.parseDouble(txtxa.getText());                       // units of R
                main.y0 = Double.parseDouble(txty0.getText());                       // units of R
                main.c = Double.parseDouble(txtc.getText());                         // moment of inertia
                main.Tx = Double.parseDouble(txtTx.getText());                       // period of x motion
                main.phi0 = Math.PI*Double.parseDouble(txtphi0.getText())/180;       // radians
                main.dtheta0dy = Double.parseDouble(txtdTheta0.getText());           // radians
                main.dw0dy = Double.parseDouble(txtdw0.getText());                   // radians/sec
                main.NLimit = NLimitCombo.getSelectedIndex();                        // number of Tx in limit cycle
                phase_space();
            }
        });
    }

    private void phase_space()
    {
        final double delt = main.Tx/Nper;
        Point2D.Double pt = new Point2D.Double(main.theta0, main.w0);
        double[] pt4 = new double[] {main.theta0, main.w0, main.dtheta0dy, main.dw0dy};   // theta, w, dthetady, dwdy
        double thmin, thmax;
        double tempx;

        if (printChk.isSelected())
        {
            GregorianCalendar now = new GregorianCalendar();
            System.out.println("\nDynamic Zeeman: " + now.getTime());
            System.out.printf("A, %.1f\ntheta0, %f\nw0, %.2f\nx0, %.2f\nxa, %.2f\ny0, %.5f\nc, %f\nTx, %.2f\nphi0, %f\n\n", main.A, main.theta0, main.w0, main.x0, main.xa, main.y0, main.c, main.Tx, main.phi0);
            System.out.println("t, x, theta, w, rhs");
            tempx = main.x0 + main.xa*Math.cos(main.phi0);
            System.out.println("0.0, " + tempx + ", " + pt.x + ", " + pt.y + ", " + (-main.c*main.calc_dFdth(pt.x, main.A, tempx, main.y0) - pt.y));
        }
        path1.reset();              // transient path
        path2.reset();              // limit cycle
        if (!ddyChk.isSelected())
        {
            path1.moveTo(pt.x, pt.y);
            thmin = pt.x;
            thmax = pt.x;
            wmin = pt.y;
            wmax = pt.y;
        }
        else
        {
            if (finalphase != null)
                System.arraycopy(finalphase, 0, pt4, 0, pt4.length);
            path1.moveTo(pt4[2], pt4[3]);
            thmin = pt4[2];
            thmax = pt4[2];
            wmin = pt4[3];
            wmax = pt4[3];
        }
        for (int i = 0; i < N; i++)
        {
            //if ((i/Nper)*Nper == i)
            //    System.out.println(i + ", " + (x0 + xa*Math.cos(phi0)) + ", " + pt.x + ", " + pt.y);
            if (!ddyChk.isSelected())
                pt = main.runge_kutta(i, pt.x, pt.y, delt, main.y0);
            else
                pt = main.runge_kutta_4(i, pt4, delt, main.y0);
            if (i == N - main.NLimit*Nper - 1)
            {
                path2.moveTo(pt.x, pt.y);
                if (main.NLimit > 0)
                {
                    thmax = pt.x;
                    thmin = pt.x;
                    wmax = pt.y;
                    wmin = pt.y;
                }
            }
            if (i < N - main.NLimit*Nper)
                path1.lineTo(pt.x, pt.y);
            else
                path2.lineTo(pt.x, pt.y);
            if (pt.x > thmax) thmax = pt.x;
            if (pt.x < thmin) thmin = pt.x;
            if (pt.y > wmax) wmax = pt.y;
            if (pt.y < wmin) wmin = pt.y;
            if (printChk.isSelected() && i >= N - main.NLimit*Nper - 1)
            {
                tempx = main.x0 + main.xa*Math.cos(2*Math.PI*(i + 1)*delt/main.Tx + main.phi0);
                if (!ddyChk.isSelected())
                    System.out.println((i + 1)*delt + ", " + tempx + ", " + pt.x + ", " + pt.y + ", " + (-main.c*main.calc_dFdth(pt.x, main.A, tempx, main.y0) - pt.y));
                else
                    System.out.println((i + 1)*delt + ", " + tempx + ", " + pt4[0] + ", " + pt4[1] + ", " + pt4[2] + ", " + pt4[3]);
            }
        }
        if (!ddyChk.isSelected())
        {
            lblThetarange.setText("θ = " + String.format("%.1f", 180*thmin/Math.PI) + ", " + String.format("%.1f", 180*thmax/Math.PI));
            lblwrange.setText("ω = " + String.format("%.2f", wmin) + ", " + String.format("%.2f", wmax));
            System.out.println("final phase, " + main.y0 + ", " + main.NLimit + ", " + (thmax - thmin) + ", " + (wmax - wmin) + ", " + 180*pt.x/Math.PI + ", " + pt.y);
        }
        else
        {
            lblThetarange.setText("dθdy = " + String.format("%.2f", thmin) + ", " + String.format("%.2f", thmax));
            lblwrange.setText("dωdy = " + String.format("%.2f", wmin) + ", " + String.format("%.2f", wmax));
            finalphase = new double[] {pt4[0], pt4[1], pt4[2], pt4[3]};
            System.out.println("final ddy, " + main.y0 + ", " + main.NLimit + ", " + (thmax - thmin) + ", " + (wmax - wmin) + ", " + 180*pt4[0]/Math.PI + ", " + pt4[1] + ", " + pt4[2] + ", " + pt4[3]);
        }
        if (complementChk.isSelected())
        {
            if (2*Math.PI - thmax < thmin) thmin = 2*Math.PI - thmax;
            if (2*Math.PI - thmin > thmax) thmax = 2*Math.PI - thmin;
            if (-wmax < wmin) wmin = -wmax;
            if (-wmin > wmax) wmax = -wmin;
        }
        AffineTransform at = new AffineTransform(phasePanel.getWidth()/(thmax - thmin),
                                          0, 0, -phasePanel.getHeight()/(wmax - wmin),
                                                -thmin*phasePanel.getWidth()/(thmax - thmin),
                                                 wmax*phasePanel.getHeight()/(wmax - wmin));
        path1.transform(at);
        path2.transform(at);
        path3.reset();                              // complementary limit cycle
        if (complementChk.isSelected() && main.NLimit > 0 && !ddyChk.isSelected())  // only for 'normal' phase space, not d/dy
        {
            pt.x = 2*Math.PI - pt.x;
            pt.y = -pt.y;
            path3.moveTo(pt.x, pt.y);               // use the endpoint of path2 as a start point
            if (printChk.isSelected())
            {
                System.out.println("complement");
                System.out.println("0.0, " + ", " + pt.x + ", " + pt.y);
            }
            for (int i = 0; i < main.NLimit*Nper; i++)
            {
                pt = main.runge_kutta(Nper/2 + i, pt.x, pt.y, delt, main.y0);
                path3.lineTo(pt.x, pt.y);
                if (printChk.isSelected()) System.out.println((i + 1)*delt + ", " + ", " + pt.x + ", " + pt.y);
            }
            path3.transform(at);
        }
        phasePanel.repaint();
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
        g2.setPaint(new Color(0, 128, 128));
        g2.setStroke(new BasicStroke(1));
        g2.draw(staticComponent.phase_dlg.path3);
        g2.setPaint(Color.BLUE);
        g2.draw(new Line2D.Double(0, getHeight()*(0 - PhaseSpace.wmax)/(PhaseSpace.wmin - PhaseSpace.wmax),
                         getWidth(), getHeight()*(0 - PhaseSpace.wmax)/(PhaseSpace.wmin - PhaseSpace.wmax)));
    }
}
