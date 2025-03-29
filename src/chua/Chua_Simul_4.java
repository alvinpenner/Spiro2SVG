
package chua;

/*
 * this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Simul_4.java.
 * simulate a Poincare map of degree 4 in the x'-y' plane
 * based on standard cubic model, Book IV, page 62 and page 69
*/

import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;
import java.util.Properties;

public class Chua_Simul_4 extends JDialog
{
    private static Properties simulProp = new Properties();
    private static final BufferedImage image = new BufferedImage(480, 480, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static JButton btnCalc = new JButton("Calc");
    private static JButton btnReset = new JButton("Reset");
    private static JButton btnClear = new JButton("Clear");
    private static JCheckBox invertChk = new JCheckBox(" invert ");
    private static JCheckBox printChk = new JCheckBox("  print  ");
    private static final String C_or_g = "C";       // use Cx + iCy parms or use gij

    private static JTextField txtalpha = new JTextField("0.00004");
    private static JTextField txttheta = new JTextField("10.75");
    private static JTextField txta = new JTextField("-5");
    private static JTextField txtb = new JTextField("5");
    private static JTextField[][] txtCarr = {{new JTextField(), new JTextField("3.0E10")},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()}};
    private static double[][] Carr = new double[5][2];
    private static Point2D.Double[] g_4 = new Point2D.Double[5];           // quartic model
    private static Point2D.Double[] h_4 = new Point2D.Double[5];           // quartic transform
    private static JTextField txtstart = new JTextField("0");
    private static JTextField txtrange = new JTextField("0");
    private static double alpha;
    private static double costheta;
    private static double sintheta;
    private static double x0, y0;
    private static int iT = 0;

    public Chua_Simul_4()
    {
        setTitle("Chua System - Simulate x'-y' Scatter (4th degree)");
        setIconImage(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        setSize(780, 532);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        lblImage.setBorder(BorderFactory.createEtchedBorder());

        final JLabel[] lblarr = {new JLabel(C_or_g + " 40"),
                                 new JLabel(C_or_g + " 31"),
                                 new JLabel(C_or_g + " 22"),
                                 new JLabel(C_or_g + " 13"),
                                 new JLabel(C_or_g + " 04")};
        JPanel[] spacerPanel = new JPanel[2];
        for (int i = 0; i < spacerPanel.length; i++)
        {
            spacerPanel[i] = new JPanel();
            spacerPanel[i].setPreferredSize(new Dimension(190, 6));
            spacerPanel[i].setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
            spacerPanel[i].setOpaque(false);
        }
        JPanel[] dataPanel = new JPanel[lblarr.length];
        for (int i = 0; i < dataPanel.length; i++)
        {
            dataPanel[i] = new JPanel();
            dataPanel[i].setPreferredSize(new Dimension(230, 22));
            dataPanel[i].setOpaque(false);
            lblarr[i].setPreferredSize(new Dimension(30, 18));
            dataPanel[i].add(lblarr[i]);
            txtCarr[i][0].setPreferredSize(new Dimension(90, 18));
            txtCarr[i][1].setPreferredSize(new Dimension(90, 18));
            dataPanel[i].add(txtCarr[i][0]);
            dataPanel[i].add(txtCarr[i][1]);
        }

        JPanel alphaPanel = new JPanel();
        alphaPanel.setOpaque(false);
        alphaPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lblalpha = new JLabel("alpha");
        lblalpha.setPreferredSize(new Dimension(50, 18));
        alphaPanel.add(lblalpha);
        txtalpha.setPreferredSize(new Dimension(70, 18));
        alphaPanel.add(txtalpha);

        JPanel thetaPanel = new JPanel();
        thetaPanel.setOpaque(false);
        thetaPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lbltheta = new JLabel("theta");
        lbltheta.setPreferredSize(new Dimension(50, 18));
        thetaPanel.add(lbltheta);
        txttheta.setPreferredSize(new Dimension(70, 18));
        thetaPanel.add(txttheta);

        JPanel aPanel = new JPanel();
        aPanel.setOpaque(false);
        aPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lbla = new JLabel("a");
        lbla.setPreferredSize(new Dimension(50, 18));
        aPanel.add(lbla);
        txta.setPreferredSize(new Dimension(70, 18));
        aPanel.add(txta);

        JPanel bPanel = new JPanel();
        bPanel.setOpaque(false);
        bPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lblb = new JLabel("b");
        lblb.setPreferredSize(new Dimension(50, 18));
        bPanel.add(lblb);
        txtb.setPreferredSize(new Dimension(70, 18));
        bPanel.add(txtb);

        JPanel startPanel = new JPanel();
        startPanel.setOpaque(false);
        startPanel.setPreferredSize(new Dimension(185, 24));
        JLabel lblstart = new JLabel("start");
        lblstart.setPreferredSize(new Dimension(60, 18));
        startPanel.add(lblstart);
        txtstart.setPreferredSize(new Dimension(70, 18));
        //txtstart.setEnabled(false);
        startPanel.add(txtstart);

        JPanel rangePanel = new JPanel();
        rangePanel.setOpaque(false);
        rangePanel.setPreferredSize(new Dimension(185, 24));
        JLabel lblrange = new JLabel("range");
        lblrange.setPreferredSize(new Dimension(60, 18));
        rangePanel.add(lblrange);
        txtrange.setPreferredSize(new Dimension(70, 18));
        rangePanel.add(txtrange);

        JPanel printPanel = new JPanel();
        printPanel.setOpaque(false);
        printPanel.setPreferredSize(new Dimension(185, 24));
        invertChk.setOpaque(false);
        printPanel.add(invertChk);
        printChk.setOpaque(false);
        printPanel.add(printChk);

        final JPanel parmsPanel = new JPanel();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(alphaPanel);
        parmsPanel.add(thetaPanel);
        parmsPanel.add(aPanel);
        parmsPanel.add(bPanel);
        parmsPanel.add(spacerPanel[0]);
        for (int i = 0; i < dataPanel.length; i++)
            parmsPanel.add(dataPanel[i]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnCalc);
        parmsPanel.add(startPanel);
        parmsPanel.add(rangePanel);
        parmsPanel.add(btnReset);
        parmsPanel.add(btnClear);
        parmsPanel.add(printPanel);
        parmsPanel.setMaximumSize(new Dimension(250, 3000));
        parmsPanel.setPreferredSize(new Dimension(250, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        final JPanel scatterPanel = new JPanel();
        scatterPanel.setOpaque(false);
        scatterPanel.add(lblImage);

        load_prefs();
        x0 = Double.parseDouble(txtstart.getText());
        y0 = Double.parseDouble(txtstart.getText());
        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().setBackground(new Color(200, 221, 242));
        getContentPane().add(parmsPanel);
        getContentPane().add(scatterPanel);

        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        btnCalc.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                btnCalc.setEnabled(false);
                for (int i = 0; i < txtCarr.length; i++)
                    for (int j = 0; j < 2; j++)
                        if (txtCarr[i][j].getText().isEmpty())
                            txtCarr[i][j].setText("0");
                refresh_graph();
                btnCalc.setEnabled(true);
            }
        });

        btnReset.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                x0 = Double.parseDouble(txtstart.getText());
                y0 = Double.parseDouble(txtstart.getText());
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
                {
                    save_prefs();
                }
            });
    }

    protected static void init()
    {
        // convert from Cij (without factorials)
        // to gij (with factorials)
        // see Book "Averaging" p.57
        final double[][] C_g_4_R = new double[][] {{24,  0,-24,  0, 24},
                                                   {24,  0,  0,  0,-24},
                                                   {24,  0,  8,  0, 24},
                                                   {24,  0,  0,  0,-24},
                                                   {24,  0,-24,  0, 24}};
        final double[][] C_g_4_I = new double[][] {{ 0,-24,  0, 24,  0},
                                                   { 0,-12,  0,-12,  0},
                                                   { 0,  0,  0,  0,  0},
                                                   { 0, 12,  0, 12,  0},
                                                   { 0, 24,  0,-24,  0}};
        for (int i = 0; i < C_g_4_R.length; i++)
            for (int j = 0; j < C_g_4_R.length; j++)
            {
                C_g_4_R[i][j] /= 16.0;
                C_g_4_I[i][j] /= 16.0;
            }

        Point2D.Double mu10 = new Point2D.Double((1 + alpha)*costheta, (1 + alpha)*sintheta); // first-order response
        //System.out.println("linear a b, " + mu10);
        Point2D.Double mu01 = Main.conjugate(mu10);                  // conjugate of mu
        Point2D.Double mu20 = Main.multiply(mu10, mu10);
        Point2D.Double mu11 = Main.multiply(mu10, mu01);
        Point2D.Double mu02 = Main.multiply(mu01, mu01);
        Point2D.Double mu40 = Main.multiply(mu20, mu20);
        Point2D.Double mu31 = Main.multiply(mu20, mu11);
        Point2D.Double mu22 = Main.multiply(mu11, mu11);
        Point2D.Double mu13 = Main.multiply(mu02, mu11);
        Point2D.Double mu04 = Main.multiply(mu02, mu02);

        // define Taylor coeff 'gij' using Kuznetsov notation and a 2D point for complex numbers
        for (int i = 0; i < g_4.length; i++)
        {
            g_4[i] = new Point2D.Double(0.0, 0.0);
            for (int j = 0; j < g_4.length; j++)
            {
                g_4[i].x += C_g_4_R[i][j]*Carr[j][0] - C_g_4_I[i][j]*Carr[j][1];
                g_4[i].y += C_g_4_R[i][j]*Carr[j][1] + C_g_4_I[i][j]*Carr[j][0];
            }
        }
        System.out.print("g_4");
        for (int i = 0; i < g_4.length; i++)
            System.out.print(", " + g_4[i]);
        System.out.println();

        // coeff 'h' of nonlinear transform

        Point2D.Double htemp = new Point2D.Double();            // (Kuznetsov p. 151)
        htemp = Main.subtract(mu40, mu10);
        h_4[0] = Main.divide(g_4[0], htemp);
        htemp = Main.subtract(mu31, mu10);
        h_4[1] = Main.divide(g_4[1], htemp);
        htemp = Main.subtract(mu22, mu10);
        h_4[2] = Main.divide(g_4[2], htemp);
        htemp = Main.subtract(mu13, mu10);
        h_4[3] = Main.divide(g_4[3], htemp);
        htemp = Main.subtract(mu04, mu10);
        h_4[4] = Main.divide(g_4[4], htemp);
        System.out.print("h_4");
        for (int i = 0; i < h_4.length; i++)
            System.out.print(", " + h_4[i]);
        System.out.println();
        //System.out.println("test transform : " + transform(-0.2374, 1.825));
    }

    protected static Point2D.Double transform(double x, double y)
    {
        // transform z to w
        Point2D.Double z10 = new Point2D.Double(x, y);
        Point2D.Double z01 = Main.conjugate(z10);
        Point2D.Double z20 = Main.multiply(z10, z10);
        Point2D.Double z11 = Main.multiply(z10, z01);
        Point2D.Double z02 = Main.multiply(z01, z01);
        Point2D.Double z40 = Main.multiply(z20, z20);
        Point2D.Double z31 = Main.multiply(z20, z11);
        Point2D.Double z22 = Main.multiply(z11, z11);
        Point2D.Double z13 = Main.multiply(z02, z11);
        Point2D.Double z04 = Main.multiply(z02, z02);

        z10 = Main.subtract(z10, Main.multiply(1.0/24.0, h_4[0], z40));
        z10 = Main.subtract(z10, Main.multiply( 1.0/6.0, h_4[1], z31));
        z10 = Main.subtract(z10, Main.multiply( 1.0/4.0, h_4[2], z22));
        z10 = Main.subtract(z10, Main.multiply( 1.0/6.0, h_4[3], z13));
        z10 = Main.subtract(z10, Main.multiply(1.0/24.0, h_4[4], z04));
        //z10 = Main.add(z10, Main.multiply(inv_3[0], z30));
        //z10 = Main.add(z10, Main.multiply(inv_3[1], z21));
        //z10 = Main.add(z10, Main.multiply(inv_3[2], z12));
        //z10 = Main.add(z10, Main.multiply(inv_3[3], z03));
        return z10;
    }

    private void refresh_graph()
    {
        PrintWriter fout = null;
        alpha = Double.parseDouble(txtalpha.getText());
        costheta = Math.cos(Double.parseDouble(txttheta.getText())*Math.PI/180);
        sintheta = Math.sin(Double.parseDouble(txttheta.getText())*Math.PI/180);
        double a = Double.parseDouble(txta.getText());
        double b = Double.parseDouble(txtb.getText());
        double range = Double.parseDouble(txtrange.getText());
        Point2D.Double ztrans;
        double xplt = 0, yplt = 0;              // transformed coord
        double x, y, xtemp, ytemp;
        int N = 100000;

        System.out.println("refresh_graph " + invertChk.isSelected());
        if (C_or_g.startsWith("C"))
        {
            for (int i = 0; i < txtCarr.length; i++)    // default to Cxy input
            {
                Carr[i][0] = Double.parseDouble(txtCarr[i][0].getText());   // Cx
                Carr[i][1] = Double.parseDouble(txtCarr[i][1].getText());   // Cy
            }
            init();                                     // use Cxy to define gij, hij
        }
        else                // override Cxy input with gij input (May 19/23 looseleaf)
        {
            // convert from gij (with factorials)
            // to Cij (without factorials)
            final double[][] coeffR = new double[][] {{1,  0, -6,  0,  1},   // for g40
                                                      {1,  0,  0,  0, -1},
                                                      {1,  0,  2,  0,  1},
                                                      {1,  0,  0,  0, -1},
                                                      {1,  0, -6,  0,  1}};  // for g04
            final double[][] coeffI = new double[][] {{0,  4,  0, -4,  0},
                                                      {0,  2,  0,  2,  0},
                                                      {0,  0,  0,  0,  0},
                                                      {0, -2,  0, -2,  0},
                                                      {0, -4,  0,  4,  0}};
            System.out.println("\nCarr from gi");
            System.out.println("i, Carr[i][0], Carr[i][1]");
            for (int i = 0; i < txtCarr.length; i++)            // default to Cxy input
            {
                Carr[i][0] = 0;   // Cx
                Carr[i][1] = 0;   // Cy
                for (int j = 0; j < txtCarr.length; j++)
                {
                    Carr[i][0] += coeffR[j][i]*Double.parseDouble(txtCarr[j][0].getText()) - coeffI[j][i]*Double.parseDouble(txtCarr[j][1].getText());
                    Carr[i][1] += coeffI[j][i]*Double.parseDouble(txtCarr[j][0].getText()) + coeffR[j][i]*Double.parseDouble(txtCarr[j][1].getText());
                }
                System.out.println(i + ", " + Carr[i][0] + ", " + Carr[i][1]);
            }
        }
        if (printChk.isSelected())
            try
            {
                boolean fexists = new File("C:\\Windows\\Temp\\Chua_Simul_4" + C_or_g + "_" + txtalpha.getText() + "_" + txttheta.getText() + "_" + txta.getText() + "_" + txtb.getText() + ".csv").exists();
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\Chua_Simul_4" + C_or_g + "_" + txtalpha.getText() + "_" + txttheta.getText() + "_" + txta.getText() + "_" + txtb.getText() + ".csv", true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    fout.println("Cxi, " + Carr[0][0] + ", " + Carr[1][0] + ", " + Carr[2][0] + ", " + Carr[3][0] + ", " + Carr[4][0] + ", NaN");
                    fout.println("Cyi, " + Carr[0][1] + ", " + Carr[1][1] + ", " + Carr[2][1] + ", " + Carr[3][1] + ", " + Carr[4][1] + ", NaN");
                    fout.println("scatter hdr, NaN, " + (1 + alpha) + ", " + txttheta.getText());
                    fout.println("init x0 y0, 0, 0, NaN, NaN");
                    fout.println("iter       , x', y'");
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}
        //image.setRGB(image.getWidth()/2, image.getHeight()/2, Color.BLACK.getRGB());
        for (int i = 0; i < N; i++)
        {
            xtemp = (1 + alpha)*x0 + (x0*x0 + y0*y0)*(a*x0 - b*y0);
            ytemp = (1 + alpha)*y0 + (x0*x0 + y0*y0)*(b*x0 + a*y0);
            //xtemp += -.5*x0*x0*x0 + 20*x0*x0*y0 - 0*x0*y0*y0 + 0*y0*y0*y0;     // fudge some quadratics
            //ytemp +=  20*x0*x0*x0 - .5*x0*x0*y0 + 20*x0*y0*y0 - 5*y0*y0*y0;     // fudge some quadratics
            xtemp += Carr[0][0]*x0*x0*x0*x0 + Carr[1][0]*x0*x0*x0*y0 + Carr[2][0]*x0*x0*y0*y0 + Carr[3][0]*x0*y0*y0*y0 + Carr[4][0]*y0*y0*y0*y0;
            ytemp += Carr[0][1]*x0*x0*x0*x0 + Carr[1][1]*x0*x0*x0*y0 + Carr[2][1]*x0*x0*y0*y0 + Carr[3][1]*x0*y0*y0*y0 + Carr[4][1]*y0*y0*y0*y0;
            x = costheta*xtemp - sintheta*ytemp;
            y = sintheta*xtemp + costheta*ytemp;
            //x += Carr[0][0]*x0*x0*x0*x0 + Carr[1][0]*x0*x0*x0*y0 + Carr[2][0]*x0*x0*y0*y0 + Carr[3][0]*x0*y0*y0*y0 + Carr[4][0]*y0*y0*y0*y0;
            //y += Carr[0][1]*x0*x0*x0*x0 + Carr[1][1]*x0*x0*x0*y0 + Carr[2][1]*x0*x0*y0*y0 + Carr[3][1]*x0*y0*y0*y0 + Carr[4][1]*y0*y0*y0*y0;
            if (invertChk.isSelected())
            {
                ztrans = transform(x, y);
                xplt = ztrans.x;
                yplt = ztrans.y;
            }
            else
            {
                xplt = x;
                yplt = y;
            }
            if (image.getWidth()/2 + xplt/range*image.getWidth()/2 < 0
            ||  image.getWidth()/2 + xplt/range*image.getWidth()/2 > image.getWidth()
            ||  image.getHeight()/2 - yplt/range*image.getHeight()/2 < 0
            ||  image.getHeight()/2 - yplt/range*image.getHeight()/2 > image.getHeight())
                System.out.println(xplt + ", " + yplt);
            else
                image.setRGB(image.getWidth()/2 + (int) (xplt/range*image.getWidth()/2), image.getHeight()/2  - (int) (yplt/range*image.getHeight()/2), Color.BLACK.getRGB());
            x0 = x;
            y0 = y;
            if (fout != null)
                fout.println(iT + ", " + xplt + ", " + yplt);
            iT++;
        }
        int r = (int) (Math.sqrt(xplt*xplt + yplt*yplt)/range*image.getWidth()/2);
        //DC.setColor(Color.blue);
        DC.setColor(new Color(192, 128, 96));
        DC.drawLine(0, image.getHeight()/2, image.getWidth(), image.getHeight()/2);
        DC.drawLine(image.getWidth()/2, 0, image.getWidth()/2, image.getHeight());
        DC.drawOval(image.getHeight()/2 - r, image.getWidth()/2 - r, 2*r, 2*r);
        DC.drawLine((int) (image.getWidth()/2*(1 - costheta)), (int) (image.getHeight()/2*(1 + sintheta)), (int) (image.getWidth()/2*(1 + costheta)), (int) (image.getHeight()/2*(1 - sintheta)));
        lblImage.repaint();
        if (fout != null)
            fout.close();
        //calc_response_anal(a, b, alpha, x0, y0);
    }

    private static void calc_response_incr(double a, double b, double alpha, double x_in, double y_in)
    {
        // calculate the response, numerically, to a slight perturbation, after one cycle
        double[] x_step = {0, -1, 1,  0, 0};
        double[] y_step = {0,  0, 0, -1, 1};
        double incr = 0.00001;
        double x, y, xtemp, ytemp;

        //x_in = 0.0023009946092825175;
        //y_in = 0.0012750025504145068;
        System.out.println();
        for (int i = 0; i < x_step.length; i++)
        {
            x = x_in + incr*x_step[i];
            y = y_in + incr*y_step[i];
            xtemp = (1 + alpha)*x + (x*x + y*y)*(a*x - b*y);
            ytemp = (1 + alpha)*y + (x*x + y*y)*(b*x + a*y);
            xtemp += Carr[0][0]*x*x*x*x + Carr[1][0]*x*x*x*y + Carr[2][0]*x*x*y*y + Carr[3][0]*x*y*y*y + Carr[4][0]*y*y*y*y;
            ytemp += Carr[0][1]*x*x*x*x + Carr[1][1]*x*x*x*y + Carr[2][1]*x*x*y*y + Carr[3][1]*x*y*y*y + Carr[4][1]*y*y*y*y;
            System.out.println("calc_numer , " + x + ", " + y + ", " + (costheta*xtemp - sintheta*ytemp) + ", " + (sintheta*xtemp + costheta*ytemp));
        }
    }

    private static void calc_response_anal(double a, double b, double alpha, double x_in, double y_in)
    {
        // calculate the response, analytically, to a perturbation, after one cycle
        double A11, A12, A21, A22;

        System.out.println();
        A11 = 1 + alpha + 2*(a*x_in - b*y_in)*x_in + a*(x_in*x_in + y_in*y_in);
        A12 =             2*(a*x_in - b*y_in)*y_in - b*(x_in*x_in + y_in*y_in);
        A21 =             2*(b*x_in + a*y_in)*x_in + b*(x_in*x_in + y_in*y_in);
        A22 = 1 + alpha + 2*(b*x_in + a*y_in)*y_in + a*(x_in*x_in + y_in*y_in);
        A11 += 4*Carr[0][0]*x_in*x_in*x_in + 3*Carr[1][0]*x_in*x_in*y_in + 2*Carr[2][0]*x_in*y_in*y_in + 1*Carr[3][0]*y_in*y_in*y_in;
        A12 += 1*Carr[1][0]*x_in*x_in*x_in + 2*Carr[2][0]*x_in*x_in*y_in + 3*Carr[3][0]*x_in*y_in*y_in + 4*Carr[4][0]*y_in*y_in*y_in;
        A21 += 4*Carr[0][1]*x_in*x_in*x_in + 3*Carr[1][1]*x_in*x_in*y_in + 2*Carr[2][1]*x_in*y_in*y_in + 1*Carr[3][1]*y_in*y_in*y_in;
        A22 += 1*Carr[1][1]*x_in*x_in*x_in + 2*Carr[2][1]*x_in*x_in*y_in + 3*Carr[3][1]*x_in*y_in*y_in + 4*Carr[4][1]*y_in*y_in*y_in;
        System.out.println("calc_anal  , " + x_in + ", " + y_in + ", " + (costheta*A11 - sintheta*A21) + ", " + (costheta*A12 - sintheta*A22));
        System.out.println("           ,   ,   , " +                     (sintheta*A11 + costheta*A21) + ", " + (sintheta*A12 + costheta*A22));
    }

    private static void load_prefs()
    {
        try                                         // recall simulation properties
        {
            if (new File(System.getProperty("user.home"), "ChuaSimul4.ini").exists())
            {
                simulProp.load(new FileInputStream(new File(System.getProperty("user.home"), "ChuaSimul4.ini")));
                txtalpha.setText(simulProp.getProperty("alpha", "0.00004"));
                txttheta.setText(simulProp.getProperty("theta", "10.75"));
                txta.setText(simulProp.getProperty("a", "-5"));
                txtb.setText(simulProp.getProperty("b", "5"));
                for (int i = 0; i < txtCarr.length; i++)
                    txtCarr[i][0].setText(simulProp.getProperty("Carr_" + i + "0", "0"));
                for (int i = 0; i < txtCarr.length; i++)
                    txtCarr[i][1].setText(simulProp.getProperty("Carr_" + i + "1", "0"));
                txtstart.setText(simulProp.getProperty("start", "0.00001"));
                txtrange.setText(simulProp.getProperty("range", "0.005"));
            }
        }
        catch (IOException e)
            {System.out.println("error reading ChuaSimul4.ini : " + e);}
    }

    private static void save_prefs()
    {
        simulProp.setProperty("alpha", txtalpha.getText());
        simulProp.setProperty("theta", txttheta.getText());
        simulProp.setProperty("a", txta.getText());
        simulProp.setProperty("b", txtb.getText());
        for (int i = 0; i < txtCarr.length; i++)
            simulProp.setProperty("Carr_" + i + "0", txtCarr[i][0].getText());
        for (int i = 0; i < txtCarr.length; i++)
            simulProp.setProperty("Carr_" + i + "1", txtCarr[i][1].getText());
        simulProp.setProperty("start", txtstart.getText());
        simulProp.setProperty("range", txtrange.getText());
        try
            {simulProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "ChuaSimul4.ini"), "Chua Simulate Degree 4 Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                Chua_Simul_4 dlg = new Chua_Simul_4();
            }
        });
    }
}
