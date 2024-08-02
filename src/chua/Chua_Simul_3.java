
package chua;

/*
 * this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Simul_3.java.
 * simulate a Poincare map of degree 3 in the x'-y' plane
 * based on standard cubic model, Book IV, page 62 and page 69.

 * separately: convert real, 2D, cubic, coeff Cxy to complex gij - July, 2024
 * multiply complex gij to produce 2D point
 * implement, and invert, a cubic transform in z-space, using gij
 * to emulate a normalized cubic model of a torus
 * 'hdr' must be the output of 'Chua_N_response_circ.py'
 * Note: for Cxy we use a series expansion in (x, y) without any factorials
 * Note: for gij, hij, we use standard Taylor expansion with factorials
 */

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;
import java.util.Properties;

public class Chua_Simul_3 extends JDialog
{
    //private static String hdr = " 69597600 ,  99.9948 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990240639707452E-5 , 1.2809579565997563e-06 , 0.9826226599215292 , -0.18659758937064014 , -1.103029973043194 , 8.462874528763841 , -34.24146413987307 , 8.728918049582369 , 272.23195311295456 , -1730.6613177460904 , -1879.1440596659322 , -1.0364954821513763e-06 , 0.18683689988920496 , 0.9821346041621322 , -0.004294046470095772 , 43.389922533105086 , 9.147510786513596 , 44.33005570148138 , 1047.6157267475264 , -440.64698522425294 , 3090.601610076812";
    //private static String hdr = " 69597600 ,  99.9948 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990240639707452E-5 , 0.11185872467139167 , 0.9826267056297189 , -0.18659642970817192 , -1.107155942014133 , 8.7312104139989 , -34.17683230319818 , 9.36956946945483 , 232.6286221124231 , -1698.9717975602443 , -1868.6516702469346 , 0.027948518804287947 , 0.18683698202314097 , 0.9821347500967078 , -0.005320900593378551 , 43.475293577200446 , 9.063606245173178 , 44.52993222294985 , 1034.1542723056987 , -415.05667095352396 , 3083.976801673444";
    //private static String hdr = " 347999900 ,  99.98 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990205973471885E-5 , 1.3917989234171478e-07 , 0.982507316281533 , -0.18646289595781568 , 2.0254511106548967 , 0.2462332337496824 , -18.866347271505557 , 6.09069155938875 , -256.79226952744136 , -601.7640282013914 , 91.04610381992342 , -2.499974395182974e-07 , 0.1864439547726416 , 0.9824552395821273 , 9.108590542550685 , 23.934163656379955 , 2.276766550194862 , 181.00960970931143 , 635.0614869130052 , 907.5156477005373 , 914.4589182330253";
    private static String hdr = " 347997600 ,  99.98 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990205973471885E-5 , 0.10988280519552579 , 0.9825569079471931 , -0.18616685775475939 , -1.258623242989341 , 9.215500941711188 , -34.71353701828044 , 4.932463163344982 , 278.0045531672077 , -1642.4571438608352 , -2005.693670826398 , 0.029299664078839736 , 0.18643930591184257 , 0.9820440804963968 , -0.7295355591727333 , 43.909843164825794 , 9.333720416981011 , 28.50130065994334 , 1083.498286323412 , -576.6878374054445 , 3137.806271392609";
    private static double[][] C_g_2_R = new double[][] {{ 1, 0,-1},     // convert C to g
                                                        { 1, 0, 1},     // quadratic, real
                                                        { 1, 0,-1}};
    private static double[][] C_g_2_I = new double[][] {{ 0,-1, 0},     // see 'Chua_2D_cubic_variable_c1.py'
                                                        { 0, 0, 0},
                                                        { 0, 1, 0}};
    private static double[][] C_g_3_R = new double[][] {{ 3, 0,-3, 0},
                                                        { 3, 0, 1, 0},
                                                        { 3, 0, 1, 0},
                                                        { 3, 0,-3, 0}};
    private static double[][] C_g_3_I = new double[][] {{ 0,-3, 0, 3},
                                                        { 0,-1, 0,-3},
                                                        { 0, 1, 0, 3},
                                                        { 0, 3, 0,-3}};
    private static Point2D.Double[] h_2 = new Point2D.Double[3];        // quadratic response
    private static Point2D.Double[] inv_3 = new Point2D.Double[4];      // inverse cubic response

    private static Properties simulProp = new Properties();
    private static final BufferedImage image = new BufferedImage(480, 480, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static JButton btnCalc = new JButton("Calc");
    private static JButton btnReset = new JButton("Reset");
    private static JButton btnClear = new JButton("Clear");
    private static JCheckBox printChk = new JCheckBox("  print  ");

    private static JTextField txtalpha = new JTextField();
    private static JTextField txttheta = new JTextField();
    private static JTextField txta = new JTextField();
    private static JTextField txtb = new JTextField();
    private static JTextField[][] txtCarr = {{new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()}};
    private static double[][] Carr = new double[7][2];
    private static JTextField txtstart = new JTextField("0");
    private static JTextField txtrange = new JTextField("0");
    private static double x0, y0;
    private static int iT = 0;

    public Chua_Simul_3()
    {
        setTitle("Chua System - Simulate x'-y' Scatter (cubic map)");
        setIconImage(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        setSize(780, 544);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        lblImage.setBorder(BorderFactory.createEtchedBorder());

        final JLabel[] lblarr = {new JLabel("C 20"),
                                 new JLabel("C 11"),
                                 new JLabel("C 02"),
                                 new JLabel("C 30"),
                                 new JLabel("C 21"),
                                 new JLabel("C 12"),
                                 new JLabel("C 03")};
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
        lblalpha.setPreferredSize(new Dimension(40, 18));
        alphaPanel.add(lblalpha);
        txtalpha.setPreferredSize(new Dimension(90, 18));
        txtalpha.setEditable(false);
        alphaPanel.add(txtalpha);

        JPanel thetaPanel = new JPanel();
        thetaPanel.setOpaque(false);
        thetaPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lbltheta = new JLabel("theta");
        lbltheta.setPreferredSize(new Dimension(40, 18));
        thetaPanel.add(lbltheta);
        txttheta.setPreferredSize(new Dimension(90, 18));
        txttheta.setEditable(false);
        thetaPanel.add(txttheta);

        JPanel aPanel = new JPanel();
        aPanel.setOpaque(false);
        aPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lbla = new JLabel("a");
        lbla.setPreferredSize(new Dimension(40, 18));
        aPanel.add(lbla);
        txta.setPreferredSize(new Dimension(90, 18));
        aPanel.add(txta);

        JPanel bPanel = new JPanel();
        bPanel.setOpaque(false);
        bPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lblb = new JLabel("b");
        lblb.setPreferredSize(new Dimension(40, 18));
        bPanel.add(lblb);
        txtb.setPreferredSize(new Dimension(90, 18));
        bPanel.add(txtb);

        JPanel startPanel = new JPanel();
        startPanel.setOpaque(false);
        startPanel.setPreferredSize(new Dimension(185, 24));
        JLabel lblstart = new JLabel("start");
        lblstart.setPreferredSize(new Dimension(60, 18));
        startPanel.add(lblstart);
        txtstart.setPreferredSize(new Dimension(70, 18));
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
        if (true) load_cubic_fit();
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
        System.out.println("call init from Chua_Simul_3");
        String iTstr  = hdr.split(",")[0].trim();
        double Chua_alpha  = Double.parseDouble(hdr.split(",")[1]);
        double Chua_beta   = Double.parseDouble(hdr.split(",")[2]);
        double Chua_gamma  = Double.parseDouble(hdr.split(",")[3]);
        double Chua_a      = Double.parseDouble(hdr.split(",")[4]);
        double Chua_c      = Double.parseDouble(hdr.split(",")[5]);
        double Period = Double.parseDouble(hdr.split(",")[6]);
        double delt   = Double.parseDouble(hdr.split(",")[7]);
        double[] Cx = new double[10];
        double[] Cy = new double[10];
        for (int i = 0; i < Cx.length; i++)
            Cx[i] = Double.parseDouble(hdr.split(",")[8 + i]);
        for (int i = 0; i < Cy.length; i++)
            Cy[i] = Double.parseDouble(hdr.split(",")[18 + i]);
        for (int i = 0; i < C_g_2_R.length; i++)
            for (int j = 0; j < C_g_2_R.length; j++)
            {
                C_g_2_R[i][j] /= 2.0;
                C_g_2_I[i][j] /= 2.0;
            }
        for (int i = 0; i < C_g_3_R.length; i++)
            for (int j = 0; j < C_g_3_R.length; j++)
            {
                C_g_3_R[i][j] /= 4.0;
                C_g_3_I[i][j] /= 4.0;
            }

        // define Taylor coeff 'g' using Kuznetsov notation and a 2D point for complex numbers

        Point2D.Double mu10 = new Point2D.Double(Cx[1], Cy[1]); // first-order response
        Point2D.Double mu01 = conjugate(mu10);                  // conjugate of mu
        Point2D.Double mu20 = multiply(mu10, mu10);
        Point2D.Double mu11 = multiply(mu10, mu01);
        Point2D.Double mu02 = multiply(mu01, mu01);
        Point2D.Double mu30 = multiply(mu20, mu10);
        Point2D.Double mu21 = multiply(mu20, mu01);
        Point2D.Double mu12 = multiply(mu11, mu01);
        Point2D.Double mu03 = multiply(mu02, mu01);
        Point2D.Double[] g_2 = new Point2D.Double[3];           // quadratic response
        Point2D.Double[] g_3 = new Point2D.Double[4];           // cubic response
        for (int i = 0; i < g_2.length; i++)
        {
            g_2[i] = new Point2D.Double(0.0, 0.0);
            for (int j = 0; j < g_2.length; j++)
            {
                g_2[i].x += C_g_2_R[i][j]*Cx[j + 3] - C_g_2_I[i][j]*Cy[j + 3];
                g_2[i].y += C_g_2_R[i][j]*Cy[j + 3] + C_g_2_I[i][j]*Cx[j + 3];
            }
        }
        for (int i = 0; i < g_3.length; i++)
        {
            g_3[i] = new Point2D.Double(0.0, 0.0);
            for (int j = 0; j < g_3.length; j++)
            {
                g_3[i].x += C_g_3_R[i][j]*Cx[j + 6] - C_g_3_I[i][j]*Cy[j + 6];
                g_3[i].y += C_g_3_R[i][j]*Cy[j + 6] + C_g_3_I[i][j]*Cx[j + 6];
            }
        }
        //System.out.println("convert_Cxy_to_gij:" + hdr);
        System.out.println("convert_Cxy_to_gij: " + iTstr + ", " + Chua_alpha + ", " + Chua_beta + ", " + Chua_gamma + ", " + Chua_a + ", " + Chua_c + ", " + Period + ", " + delt);
        System.out.print("Cx ");
        for (int i = 0; i < Cx.length; i++)
            System.out.print(", " + Cx[i]);
        System.out.println();
        System.out.print("Cy ");
        for (int i = 0; i < Cy.length; i++)
            System.out.print(", " + Cy[i]);
        System.out.println();

        System.out.print("g_2");
        for (int i = 0; i < g_2.length; i++)
            System.out.print(", " + g_2[i]);
        System.out.println();
        System.out.print("g_3");
        for (int i = 0; i < g_3.length; i++)
            System.out.print(", " + g_3[i]);
        System.out.println();

        // coeff 'h' of nonlinear transform

        Point2D.Double htemp = new Point2D.Double();            // (Kuznetsov p. 151)
        htemp = subtract(mu20, mu10);
        h_2[0] = divide(g_2[0], htemp);
        htemp = subtract(mu11, mu10);
        h_2[1] = divide(g_2[1], htemp);
        htemp = subtract(mu02, mu10);
        h_2[2] = divide(g_2[2], htemp);
        System.out.print("h_2");
        for (int i = 0; i < h_2.length; i++)
            System.out.print(", " + h_2[i]);
        System.out.println();

        Point2D.Double[] h_3 = new Point2D.Double[4];           // cubic response
        htemp = subtract(mu30, mu10);                           // (Kuznetsov p. 152)
        h_3[0] = divide(g_3[0], htemp);
        htemp = subtract(mu21, mu10);
        h_3[1] = divide(g_3[1], htemp);
        htemp = subtract(mu12, mu10);
        h_3[2] = divide(g_3[2], htemp);
        htemp = subtract(mu03, mu10);
        h_3[3] = divide(g_3[3], htemp);
        System.out.print("h_3");
        for (int i = 0; i < h_3.length; i++)
            System.out.print(", " + h_3[i]);
        System.out.println();

        Point2D.Double htemp2 = new Point2D.Double();           // terms A, B, C, D
        Point2D.Double htemp3 = new Point2D.Double();           // N_S bifurc, Dec 21, 2021
        Point2D.Double htemp4 = new Point2D.Double();
        htemp  = multiply(0.5, multiply(h_2[0], h_2[0]));
        htemp2 = multiply(0.5, multiply(h_2[1], conjugate(h_2[2])));
        inv_3[0] = add(htemp, htemp2);                              // term A

        htemp  = multiply(1.5, multiply(h_2[0], h_2[1]));
        htemp2 = multiply(0.5, multiply(h_2[2], conjugate(h_2[2])));
        htemp3 = multiply(h_2[1], conjugate(h_2[1]));
        inv_3[1] = add(htemp, add(htemp2, htemp3));                 // term B

        htemp  = multiply(0.5, multiply(h_2[0], h_2[2]));
        htemp2 = multiply(h_2[1], h_2[1]);
        htemp3 = multiply(0.5, multiply(h_2[1], conjugate(h_2[0])));
        htemp4 = multiply(h_2[2], conjugate(h_2[1]));
        inv_3[2] = add(htemp, add(htemp2, add(htemp3, htemp4)));    // term C

        htemp  = multiply(0.5, multiply(h_2[1], h_2[2]));
        htemp2 = multiply(0.5, multiply(h_2[2], conjugate(h_2[0])));
        inv_3[3] = add(htemp, htemp2);                              // term D

        System.out.print("inv_3");
        for (int i = 0; i < inv_3.length; i++)
            System.out.print(", " + inv_3[i]);
        System.out.println();
        //System.out.println("test transform : " + transform(new Point2D.Double(-0.2374, 1.825)));
    }

    protected static Point2D.Double transform(double x, double y)
    {
        // transform z to w
        Point2D.Double z10 = new Point2D.Double(x, y);
        Point2D.Double z01 = conjugate(z10);
        Point2D.Double z20 = multiply(z10, z10);
        Point2D.Double z11 = multiply(z10, z01);
        Point2D.Double z02 = multiply(z01, z01);
        Point2D.Double z30 = multiply(z20, z10);
        Point2D.Double z21 = multiply(z20, z01);
        Point2D.Double z12 = multiply(z11, z01);
        Point2D.Double z03 = multiply(z02, z01);

        z10 = subtract(z10, multiply(0.5, h_2[0], z20));
        z10 = subtract(z10, multiply(h_2[1], z11));
        z10 = subtract(z10, multiply(0.5, h_2[2], z02));
        z10 = add(z10, multiply(inv_3[0], z30));
        z10 = add(z10, multiply(inv_3[1], z21));
        z10 = add(z10, multiply(inv_3[2], z12));
        z10 = add(z10, multiply(inv_3[3], z03));
        return z10;
    }

    private static Point2D.Double conjugate(Point2D.Double p1)
    {
        // complex conjugate of p1
        return new Point2D.Double(p1.x, -p1.y);
    }

    private static Point2D.Double add(Point2D.Double p1, Point2D.Double p2)
    {
        // add (x1 + i*y1) + (x2 + i*y2)
        return new Point2D.Double(p1.x + p2.x, p1.y + p2.y);
    }

    private static Point2D.Double subtract(Point2D.Double p1, Point2D.Double p2)
    {
        // subtract (x1 + i*y1) - (x2 + i*y2)
        return new Point2D.Double(p1.x - p2.x, p1.y - p2.y);
    }

    private static Point2D.Double multiply(double d, Point2D.Double p1)
    {
        // multiply d*(x2 + i*y2)
        return new Point2D.Double(d*p1.x, d*p1.y);
    }

    private static Point2D.Double multiply(Point2D.Double p1, Point2D.Double p2)
    {
        // multiply (x1 + i*y1)*(x2 + i*y2)
        return new Point2D.Double(p1.x*p2.x - p1.y*p2.y, p1.x*p2.y + p2.x*p1.y);
    }

    private static Point2D.Double multiply(double d, Point2D.Double p1, Point2D.Double p2)
    {
        // multiply d*(x1 + i*y1)*(x2 + i*y2)
        return multiply(d, multiply(p1, p2));
    }

    private static Point2D.Double divide(Point2D.Double p1, Point2D.Double p2)
    {
        // divide (x1 + i*y1)/(x2 + i*y2)
        double abs = p2.x*p2.x + p2.y*p2.y;
        return new Point2D.Double((p1.x*p2.x + p1.y*p2.y)/abs, (-p1.x*p2.y + p2.x*p1.y)/abs);
    }

    private void refresh_graph()
    {
        PrintWriter fout = null;
        double a = Double.parseDouble(txta.getText());
        double b = Double.parseDouble(txtb.getText());
        double range = Double.parseDouble(txtrange.getText());
        double x, y;
        int N = 100000;

        //double alpha = Math.sqrt(a*a + b*b) - 1;            // first-order response (redundant)
        //double costheta = a/(alpha + 1);                    //                      (redundant)
        //double sintheta = b/(alpha + 1);                    //                      (redundant)
        for (int i = 0; i < txtCarr.length; i++)
        {
            Carr[i][0] = Double.parseDouble(txtCarr[i][0].getText());   // Cx
            Carr[i][1] = Double.parseDouble(txtCarr[i][1].getText());   // Cy
        }
        if (printChk.isSelected())
            try
            {
                boolean fexists = new File("C:\\Windows\\Temp\\Chua_Simul_3_" + txtalpha.getText() + "_" + txttheta.getText() + "_" + txta.getText() + "_" + txtb.getText() + ".csv").exists();
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\Chua_Simul_3_" + txtalpha.getText() + "_" + txttheta.getText() + "_" + txta.getText() + "_" + txtb.getText() + ".csv", true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    fout.println("Cxi, " + Carr[0][0] + ", " + Carr[1][0] + ", " + Carr[2][0] + ", " + Carr[3][0] + ", " + Carr[4][0] + ", " + Carr[5][0] + ", " + Carr[6][0]);
                    fout.println("Cyi, " + Carr[0][1] + ", " + Carr[1][1] + ", " + Carr[2][1] + ", " + Carr[3][1] + ", " + Carr[4][1] + ", " + Carr[5][1] + ", " + Carr[6][1]);
                    fout.println("scatter hdr, NaN, " + Math.sqrt(a*a + b*b) + ", " + txttheta.getText());
                    fout.println("init x0 y0, 0, 0, NaN, NaN");
                    fout.println("iter       , x', y'");
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}
        for (int i = 0; i < N; i++)
        {
            x = a*x0 - b*y0
              + Carr[0][0]*x0*x0 + Carr[1][0]*x0*y0 + Carr[2][0]*y0*y0
              + Carr[3][0]*x0*x0*x0 + Carr[4][0]*x0*x0*y0 + Carr[5][0]*x0*y0*y0 + Carr[6][0]*y0*y0*y0;
            y = b*x0 + a*y0
              + Carr[0][1]*x0*x0 + Carr[1][1]*x0*y0 + Carr[2][1]*y0*y0
              + Carr[3][1]*x0*x0*x0 + Carr[4][1]*x0*x0*y0 + Carr[5][1]*x0*y0*y0 + Carr[6][1]*y0*y0*y0;
            //x = costheta*xtemp - sintheta*ytemp;
            //y = sintheta*xtemp + costheta*ytemp;
            if (image.getWidth()/2 + x/range*image.getWidth()/2 < 0
            ||  image.getWidth()/2 + x/range*image.getWidth()/2 > image.getWidth()
            ||  image.getHeight()/2 - y/range*image.getHeight()/2 < 0
            ||  image.getHeight()/2 - y/range*image.getHeight()/2 > image.getHeight())
                System.out.println(x + ", " + y);
            else
                image.setRGB(image.getWidth()/2 + (int) (x/range*image.getWidth()/2), image.getHeight()/2  - (int) (y/range*image.getHeight()/2), Color.BLACK.getRGB());
            x0 = x;
            y0 = y;
            if (fout != null)
                fout.println(iT + ", " + x0 + ", " + y0);
            iT++;
        }
        //int r = (int) (Math.sqrt(x0*x0 + y0*y0)/range*image.getWidth()/2);
        //DC.setColor(Color.blue);
        DC.setColor(new Color(192, 128, 96));
        DC.drawLine(0, image.getHeight()/2, image.getWidth(), image.getHeight()/2);
        DC.drawLine(image.getWidth()/2, 0, image.getWidth()/2, image.getHeight());
        //DC.drawOval(image.getHeight()/2 - r, image.getWidth()/2 - r, 2*r, 2*r);
        //DC.drawLine((int) (image.getWidth()/2*(1 - costheta)), (int) (image.getHeight()/2*(1 + sintheta)), (int) (image.getWidth()/2*(1 + costheta)), (int) (image.getHeight()/2*(1 - sintheta)));
        lblImage.repaint();
        if (fout != null)
            fout.close();
    }

    private static void load_cubic_fit()
    {
        double a = Double.parseDouble(hdr.split(",")[9]);
        txta.setText(("" + a).substring(0, 10));
        double b = Double.parseDouble(hdr.split(",")[19]);
        txtb.setText(("" + b).substring(0, 10));
        txtalpha.setText(String.format("%.7f", Math.sqrt(a*a + b*b) - 1));      // first-order response
        txttheta.setText(String.format("%.5f", Math.atan2(b, a)*180/Math.PI));  // (redundant)
        for (int i = 0; i < txtCarr.length; i++)
            txtCarr[i][0].setText(hdr.split(",")[11 + i].trim().substring(0, 10));
        for (int i = 0; i < txtCarr.length; i++)
            txtCarr[i][1].setText(hdr.split(",")[21 + i].trim().substring(0, 10));
    }

    private static void load_prefs()
    {
        try                                         // recall simulation properties
        {
            if (new File(System.getProperty("user.home"), "ChuaSimul3.ini").exists())
            {
                simulProp.load(new FileInputStream(new File(System.getProperty("user.home"), "ChuaSimul3.ini")));
                double a = Double.parseDouble(simulProp.getProperty("a", "0.98"));
                txta.setText("" + a);
                double b = Double.parseDouble(simulProp.getProperty("b", "0.2"));
                txtb.setText("" + b);
                txtalpha.setText(String.format("%.7f", Math.sqrt(a*a + b*b) - 1));      // first-order response
                txttheta.setText(String.format("%.5f", Math.atan2(b, a)*180/Math.PI));  // (redundant)
                for (int i = 0; i < txtCarr.length; i++)
                    txtCarr[i][0].setText(simulProp.getProperty("Carr_" + i + "0", "0"));
                for (int i = 0; i < txtCarr.length; i++)
                    txtCarr[i][1].setText(simulProp.getProperty("Carr_" + i + "1", "0"));
                txtstart.setText(simulProp.getProperty("start", "0.00001"));
                txtrange.setText(simulProp.getProperty("range", "0.005"));
            }
            else
                System.out.println("ChuaSimul3.ini : file not found");
        }
        catch (IOException e)
            {System.out.println("error reading ChuaSimul3.ini : " + e);}
    }

    private static void save_prefs()
    {
        simulProp.setProperty("a", txta.getText());
        simulProp.setProperty("b", txtb.getText());
        for (int i = 0; i < txtCarr.length; i++)
            simulProp.setProperty("Carr_" + i + "0", txtCarr[i][0].getText());
        for (int i = 0; i < txtCarr.length; i++)
            simulProp.setProperty("Carr_" + i + "1", txtCarr[i][1].getText());
        simulProp.setProperty("start", txtstart.getText());
        simulProp.setProperty("range", txtrange.getText());
        try
            {simulProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "ChuaSimul3.ini"), "Chua Simulate Degree 3 Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                Chua_Simul_3.init();
                //Chua_Simul_3 dlg = new Chua_Simul_3();
            }
        });
    }
}
