
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
    private final static int N = 160000; // 161000;    // total # of iterations 160000
    //private final static int N = 300;
    private final static int Nper = 100; //100;        // # of iterations per Tx (assume even)
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
    private static double[][] M;       // for use by gaussj
    private static double[] v;         // for use by gaussj
    private static double[] sol;       // for use by gaussj

    public PhaseSpace(Image img)
    {
        //phasePanel.setBackground(Color.white);
        JLabel[] lbl = {new JLabel("A"),
                        new JLabel("θ0"),
                        new JLabel("ω0"),
                        new JLabel("x0"),
                        new JLabel("xa"),
                        new JLabel("y0"),
                        new JLabel("c"),
                        new JLabel("Tx"),
                        new JLabel("φ0"),
                        new JLabel("dθ0dy"),
                        new JLabel("dω0dy")};
        final JTextField[] txt = {new JTextField(Double.toString(main.A)),
                                  new JTextField(String.format("%.4f", main.theta0*180/Math.PI)),
                                  new JTextField(String.format("%.4f", main.w0)),
                                  new JTextField(Double.toString(main.x0)),
                                  new JTextField(String.format("%.2f", main.xa)),
                                  new JTextField(String.format("%.6f", main.y0)),
                                  new JTextField(Double.toString(main.c)),
                                  new JTextField(Double.toString(main.Tx)),
                                  new JTextField(String.format("%.1f", main.phi0*180/Math.PI)),
                                  new JTextField(String.format("%.4f", main.dtheta0dy)),
                                  new JTextField(String.format("%.4f", main.dw0dy))};
        JPanel[] spacerPanel = new JPanel[3];
        JPanel[] dataPanel = new JPanel[11];
        JButton btnRun = new JButton("Run");
        String[] NLimitdata = {" 0", " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", " 10", " 11", " 12", " 13", " 14", " 15", " 16"};
        final JComboBox NLimitCombo = new JComboBox(NLimitdata);

        if (!ddyChk.isSelected())
            setTitle(" Dynamic Zeeman - ω vs. θ");
        else
            setTitle(" Dynamic Zeeman - dω/dy vs. dθ/dy");
        setIconImage(img);
        setSize(786, 639);
        setLocationByPlatform(true);

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
        txt[0].setEditable(false);                              // A distance

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
        parmsPanel.add(dataPanel[0]);
        parmsPanel.add(dataPanel[1]);
        parmsPanel.add(dataPanel[2]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(dataPanel[7]);
        parmsPanel.add(dataPanel[8]);
        parmsPanel.add(dataPanel[9]);
        parmsPanel.add(dataPanel[10]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(NLimitPanel);
        parmsPanel.add(complementChk);
        parmsPanel.add(printChk);
        parmsPanel.add(ddyChk);
        parmsPanel.add(spacerPanel[2]);
        parmsPanel.add(thetarangePanel);
        parmsPanel.add(wrangePanel);
        parmsPanel.setMaximumSize(new Dimension(170, 3000));
        parmsPanel.setPreferredSize(new Dimension(170, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().add(parmsPanel);
        getContentPane().add(phasePanel);
        setVisible(true);
        System.out.println("PhaseSpace, " + N + ", " + Nper + ", " + (float) N/Nper);

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                if (!ddyChk.isSelected())
                    setTitle(" Dynamic Zeeman - ω vs. θ");
                else
                    setTitle(" Dynamic Zeeman - dω/dy vs. dθ/dy");
                main.theta0 = Math.PI*Double.parseDouble(txt[1].getText())/180;     // radians
                main.w0 = Double.parseDouble(txt[2].getText());                     // radians/sec
                main.x0 = Double.parseDouble(txt[3].getText());                     // units of R
                main.xa = Double.parseDouble(txt[4].getText());                     // units of R
                main.y0 = Double.parseDouble(txt[5].getText());                     // units of R
                main.c = Double.parseDouble(txt[6].getText());                      // moment of inertia
                main.Tx = Double.parseDouble(txt[7].getText());                     // period of x motion
                main.phi0 = Math.PI*Double.parseDouble(txt[8].getText())/180;       // radians
                main.dtheta0dy = Double.parseDouble(txt[9].getText());              // radians
                main.dw0dy = Double.parseDouble(txt[10].getText());                 // radians/sec
                main.NLimit = NLimitCombo.getSelectedIndex();                       // number of Tx in limit cycle
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

        M = new double[Nper*main.NLimit][Nper*main.NLimit];     // for use by gaussj
        v = new double[Nper*main.NLimit];                       // for use by gaussj

        if (printChk.isSelected())
        {
            GregorianCalendar now = new GregorianCalendar();
            System.out.println("\nDynamic Zeeman: " + now.getTime());
            System.out.printf("A, %.1f\ntheta0, %f\nw0, %.3f\nx0, %.3f\nxa, %.3f\ny0, %.5f\nc, %f\nTx, %.3f\nphi0, %f\n\n", main.A, main.theta0, main.w0, main.x0, main.xa, main.y0, main.c, main.Tx, main.phi0);
            System.out.println("t, x, theta, w, rhs");
            tempx = main.x0 + main.xa*Math.cos(main.phi0);
            System.out.println("0.0, " + tempx + ", " + pt.x + ", " + pt.y + ", " + (-main.c*main.calc_dFdth(pt.x, main.A, tempx, main.y0) - pt.y));
        }
        path1.reset();              // transient path
        path2.reset();              // limit cycle
        if (!ddyChk.isSelected())
        {
            if (finalphase != null)
            {
                pt.x = finalphase[0];
                pt.y = finalphase[1];
            }
            path1.moveTo(pt.x, pt.y);
            if (N == main.NLimit*Nper)
                path2.moveTo(pt.x, pt.y);               // in case N = main.NLimit*Nper
            thmin = pt.x;
            thmax = pt.x;
            wmin = pt.y;
            wmax = pt.y;
        }
        else
        {
            if (finalphase != null)
                if (finalphase.length == pt4.length)
                    System.arraycopy(finalphase, 0, pt4, 0, pt4.length);
                else if (finalphase.length == 2)
                    System.arraycopy(finalphase, 0, pt4, 0, 2);
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
                {
                    //System.out.println((i + 1)*delt + ", " + tempx + ", " + pt.x + ", " + pt.y + ", " + (-main.c*main.calc_dFdth(pt.x, main.A, tempx, main.y0) - pt.y));
                    //System.out.println((i + 1)*delt + ", " + tempx + ", " + pt.x + ", " + pt.y + ", " + main.c*main.calc_d2Fdth2(pt.x, main.A, tempx, main.y0));
                    gen_array(i - N + main.NLimit*Nper + 1, main.NLimit*Nper, main.c*main.calc_d2Fdth2(pt.x, main.A, tempx, main.y0), main.c*main.calc_d2Fdthdy(pt.x, main.A, tempx, main.y0), delt);
                    //gen_vector(i - N + main.NLimit*Nper + 1, main.NLimit*Nper, main.c*main.calc_d2Fdth2(pt.x, main.A, tempx, main.y0), main.c*main.calc_d2Fdthdy(pt.x, main.A, tempx, main.y0), delt);
                }
                else
                    System.out.println((i + 1)*delt + ", " + tempx + ", " + pt4[0] + ", " + pt4[1] + ", " + pt4[2] + ", " + pt4[3]);
            }
        }
        if (!ddyChk.isSelected())
        {
            lblThetarange.setText("θ = " + String.format("%.1f", 180*thmin/Math.PI) + ", " + String.format("%.1f", 180*thmax/Math.PI));
            lblwrange.setText("ω = " + String.format("%.2f", wmin) + ", " + String.format("%.2f", wmax));
            finalphase = new double[] {pt.x, pt.y};
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

    private void gen_array(int i, int N, double d2Fdth2, double d2Fdthdy, double delt)
    {
        // this is for Python, to generate eigenvalues of M
        // M[0][j] = d2Fdth2, v[j] = d2Fdthdy
        if (i < N)                  // fill the array M
        {
            M[0][i] = d2Fdth2;
            v[i] = d2Fdthdy;
        }
        else                        // print the array
        {
            System.out.println("# diagonal elements of M");
            System.out.println("x0 = " + main.x0);
            System.out.println("y0 = " + main.y0 + " # N = " + N);
            System.out.println("delt = " + delt);
            System.out.print("d2Fdth2 = np.array([" + M[0][0]);
            for (int j = 1; j < N; j++)
                System.out.print(", " + M[0][j]);
            System.out.println("])");
            System.out.print("d2Fdthdy = np.array([" + v[0]);
            for (int j = 1; j < N; j++)
                System.out.print(", " + v[j]);
            System.out.println("])");
        }
    }

    private void gen_array_Excel(int i, int N, double d2Fdth2, double d2Fdthdy, double delt)
    {
        // this is for Excel
        // generate an NxN array and Nx1 vector suitable for calculating dtheta/dy
        // OBSOLETE - quadratic formula
        if (i >= N) return;
        for (int j = 0; j < N; j++)
            if (j == i)
                System.out.print((d2Fdth2 - 2.0/delt/delt) + ",");
            else if (j == (N + i - 1) % N)
                System.out.print((-0.5/delt + 1.0/delt/delt) + ",");
            else if (j == (i + 1) % N)
                System.out.print((0.5/delt + 1.0/delt/delt) + ",");
            else
                System.out.print("0,");
        System.out.println("," + (-d2Fdthdy));
    }

    private void gen_vector(int i, int N, double d2Fdth2, double d2Fdthdy, double delt)
    {
        // this is to solve directly for dthetady using gaussj
        // and using a previously pre-converged solution for theta, w versus t
        if (i < N)                  // fill the matrix M
        {
            for (int j = 0; j < N; j++)
                if (j == i)
                    M[i][j] = d2Fdth2 - 30.0/delt/delt/12;
                else if (j == (N + i - 1) % N)
                    M[i][j] = -8.0/delt/12 + 16.0/delt/delt/12;
                else if (j == (i + 1) % N)
                    M[i][j] =  8.0/delt/12 + 16.0/delt/delt/12;
                else if (j == (N + i - 2) % N)
                    M[i][j] =  1.0/delt/12 - 1.0/delt/delt/12;
                else if (j == (i + 2) % N)
                    M[i][j] = -1.0/delt/12 - 1.0/delt/delt/12;
                else
                    M[i][j] = 0.0;
            v[i] = -d2Fdthdy;
        }
        else                        // solve M*sol = v
        {
            System.out.println("gaussj solution for theta");
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                    System.out.print(M[j][k] + ", ");
                System.out.println();
            }
            sol = main.gaussj(M, v);
            System.out.println("i, dtheta/dy, dw/dy");
            for (int j = 0; j < N; j++)
                System.out.println(j + ", " + sol[j] + ", " + (-sol[(j + 2) % N] + 8.0*sol[(j + 1) % N] - 8.0*sol[(j + N - 1) % N] + sol[(j + N - 2) % N])/12/delt);
        }
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
