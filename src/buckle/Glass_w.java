
// solve the Glass - 'Hopf bifurcation' problem
// as per Langford 1977: "Numerical Solution of Bifurcation Problems"
// use Runge-Kutta to solve a two-point boundary value problem
// with constrained cyclic behaviour at t = 2*pi
// see loose leaf: 'Implementation of Glass Model', June 9

package buckle;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import javax.swing.*;

// this is file : \Documents\NetBeansProjects\ChuaOscillator\src\buckle\Glass_w.java

public class Glass_w extends JDialog
{
    private static final BufferedImage image = new BufferedImage(480, 480, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static JTextField txtradius = new JTextField("11");
    private static JTextField txtx = new JTextField();
    private static JTextField txty = new JTextField();
    private static JButton btnCalc = new JButton(" Calc ");
    private static JButton btnClear = new JButton(" Clear ");

    private static final int N = 128;                            // number of steps in (0, 1)
    private static final double beta = Math.sqrt(3);            // eigenvalue at bifurc
    private static double eps = 0.1;
    private static double eta = 0;                              // mu = 4 + eps*eps*eta
    private static double tau = 0;                              // T = 2*pi*(1 + eps*eps*tau)/beta
    private static double[][] T_inv = new double[][] {{ -1.0/beta,  1.0/beta, -1.0/beta},
                                                      {  2.0/beta,  1.0/beta, -1.0/beta},
                                                      {  0       , -1.0     , -1.0}};
    private static double[][] T     = new double[][] {{ -1.0/beta,  1.0/beta,  0},
                                                      {  1.0/beta,  0.5/beta, -0.5},
                                                      { -1.0/beta, -0.5/beta, -0.5}};
    //private static double g1_sin, g1_cos, g2_sin, g2_cos;       // coeff for boundary solution
    private static double [][] phi = new double[N + 1][3];      // linear response
    private static double [][] w = new double[N + 1][3];        // incremental response (plus final b.c. test value)
    private static double [][] P = new double[N + 1][3];        // intermediate storage (forcing function)
    private static double wx0 = 0, wy0 = 0.1;                    // start point in w plane

    public Glass_w()
    {
        setTitle("Glass_Hopf-Bifurcation System - w limit cycle");
        setIconImage(Toolkit.getDefaultToolkit().getImage(chua.Main.class.getResource("images/icon.gif")));
        setSize(820, 535);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        lblImage.setBorder(BorderFactory.createEtchedBorder());

        JPanel[] spacerPanel = new JPanel[2];
        for (int i = 0; i < spacerPanel.length; i++)
        {
            spacerPanel[i] = new JPanel();
            spacerPanel[i].setPreferredSize(new Dimension(240, 10));
            spacerPanel[i].setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
            spacerPanel[i].setOpaque(false);
        }
        JPanel radiusPanel = new JPanel();
        radiusPanel.setOpaque(false);
        radiusPanel.setPreferredSize(new Dimension(160, 24));
        JLabel lblradius = new JLabel("radius");
        lblradius.setPreferredSize(new Dimension(40, 20));
        radiusPanel.add(lblradius);
        txtradius.setPreferredSize(new Dimension(90, 20));
        radiusPanel.add(txtradius);

        JPanel coordPanel = new JPanel();
        coordPanel.setOpaque(false);
        coordPanel.setPreferredSize(new Dimension(180, 24));
        JLabel lblcoord = new JLabel("coord");
        lblcoord.setPreferredSize(new Dimension(40, 20));
        txtx.setPreferredSize(new Dimension(60, 20));
        txty.setPreferredSize(new Dimension(60, 20));
        coordPanel.add(lblcoord);
        coordPanel.add(txtx);
        coordPanel.add(txty);

        final JPanel parmsPanel = new JPanel();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(radiusPanel);
        parmsPanel.add(coordPanel);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnCalc);
        parmsPanel.add(btnClear);

        final JPanel scatterPanel = new JPanel();
        scatterPanel.setOpaque(false);
        scatterPanel.add(lblImage);
        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().setBackground(new Color(200, 221, 242));
        getContentPane().add(parmsPanel);
        getContentPane().add(scatterPanel);

        lblImage.addMouseListener(new MouseAdapter()
        {
            @Override public void mousePressed(MouseEvent e)
            {
                //System.out.println(e.getX() + ", " + e.getY() + ", " + e.getButton() + ", " + e.getModifiers() + ", " + txtradius.getText());
                wx0 =  Double.parseDouble(txtradius.getText())*(2.0*e.getX() - image.getWidth())/image.getWidth();
                wy0 = -Double.parseDouble(txtradius.getText())*(2.0*e.getY() - image.getWidth())/image.getWidth();
                System.out.println("Mouse Clicked (" + wx0 + ", " + wy0 + ")");
                if (e.getX() > 0 && e.getX() < image.getWidth() && e.getY() > 0 && e.getY() < image.getHeight())
                {
                    image.setRGB(e.getX(), e.getY(), Color.BLUE.getRGB());
                    lblImage.repaint();
                }
                //w[0] = calc_xyz(wx0, wy0);                      // initiallize original coord at t = 0
                w[0] = back_transform_3D(5, wx0, wy0);
            }
        });

        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                //txtx.setText("" + Double.parseDouble(txtradius.getText())*(2*e.getX() - image.getWidth())/image.getWidth());
                txtx.setText(String.format("%.4f", Double.parseDouble(txtradius.getText())*(2*e.getX() - image.getWidth())/image.getWidth()));
                txty.setText(String.format("%.4f", -Double.parseDouble(txtradius.getText())*(2*e.getY() - image.getWidth())/image.getWidth()));
            }
        });

        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        btnCalc.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                //w[0] = calc_xyz(wx0, wy0);                      // re-initiallize original coord at t = 0
                //w[0] = w[N];
                btnCalc.setEnabled(false);
                DC.setColor(new Color(192, 128, 96));
                DC.drawLine(image.getWidth()/2 - 5, image.getHeight()/2, image.getWidth()/2 + 5, image.getHeight()/2);
                DC.drawLine(image.getWidth()/2, image.getHeight()/2 - 5, image.getWidth()/2, image.getHeight()/2 + 5);
                for (int i = 0; i < 3; i++)            // shoot from a fixed start w[0]
                    shoot(i);
                btnCalc.setEnabled(true);
                //w[0] = w[N];                          // copy by ref (DO NOT USE !!!)
                System.arraycopy(w[N], 0, w[0], 0, 3);  // arraycopy re-initiallize by value
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
        //double [] tempint = new double[N + 1];      // temporary integrand for Cotes
        //for (int i = 0; i < N + 1; i++)
        //    tempint[i] = Math.sin(i*2*Math.PI/N)*Math.sin(i*2*Math.PI/N);
        //System.out.println("test Cotes " + run_buckle.Cotes_4(tempint));
        //double [] tempint, vec;
        //for (int i = 0; i < N + 1; i++)
        //{
        //    tempint = new double[] {Math.sqrt(2.0/3.0)*Math.sin(i*2*Math.PI/N),
        //                            Math.sqrt(2.0/3.0)*Math.sin(i*2*Math.PI/N - Math.PI/3),
        //                            Math.sqrt(2.0/3.0)*Math.sin(i*2*Math.PI/N - 2*Math.PI/3)};
        //    vec = project_2D(tempint);
        //    //System.out.println("T_inv ," + tempint[0] + ", " + tempint[1] + ", " + tempint[2]);
        //    System.out.println("vec ," + vec[0] + ", " + vec[1] + ", " + vec[2]);
        //}
        //for (int i = 0; i < 3; i++)
        //for (int j = 0; j < 3; j++)
        //{
        //    double temp = 0;
        //    for (int k = 0; k < 3; k++)
        //        temp += T[i][k]*T_inv[k][j];
        //    System.out.println("T*T_inv," + i + ", " + j + ", " + temp);
        //}
    }

    private static void load_data()
    {
        //double [] u;                                            // unit vector in 'lateral' direction
        //double d_in;                                            // Lateral distance 'in'
        //double d_out;                                           // Lateral distance 'out'
        for (int i = 0; i < N + 1; i++)
            for (int j = 0; j < 3; j++)
                P[i][j] = 0;
        for (int i = 0; i < N + 1; i++)             // initiallize temporary w to calculate phi
            for (int j = 0; j < 3; j++)
                w[i][j] = 0;
        System.out.println("N_eps, " + N + ", " + eps + ", " + eta + ", " + tau);
        w[0][0] = Math.sqrt(2.0/3)*Math.sin(0.0);   // initiallize for calc of 'phi'
        w[0][1] = Math.sqrt(2.0/3)*Math.sin(0.0 - Math.PI/3);
        w[0][2] = Math.sqrt(2.0/3)*Math.sin(0.0 - Math.PI*2/3);
        for (int i = 0; i < N; i++)                 // calculate linear response 'phi'
            runge_kutta_Glass_iter(i, new double[] {0, 0, 0});
        for (int i = 0; i < N + 1; i++)             // copy w[] to phi[]
            System.arraycopy(w[i], 0, phi[i], 0, 3);
        //System.out.println("\ninitiallize phi, " + eps + ", " + eta + ", " + tau);
        //System.out.println("iter, phi[0], phi[1], phi[2]");
        //for (int i = 0; i < N + 1; i++)
        //    System.out.println(i + ", " + phi[i][0] + ", " + phi[i][1] + ", " + phi[i][2]);
        for (int i = 0; i < N + 1; i++)
            for (int j = 0; j < 3; j++)
                w[i][j] = 0;

        // update w first, then rhs (P[N][3])
        //w[0] = new double[] {0, -.5*10, -.5*10};              // initiallize at t = 0 TEST ONLY !!!
        w[0] = back_transform_3D(5, wx0, wy0);
        //w[0] = calc_xyz(wx0, wy0);                      // initiallize original coord at t = 0
        //w[0][0] += 5.0;
        //for (int loop = 0; loop < 1; loop++)        // increment start w[0]
        //{
        //    System.out.println("loop " + loop);
        //    for (int i = 0; i < 10; i++)            // shoot from a fixed start w[0]
        //        shoot(i);
            //u = calc_lateral();                     // lateral unit vector
            //d_in = w[0][0]*u[0] + w[0][1]*u[1] + w[0][2]*u[2];
            //d_out = w[N][0]*u[0] + w[N][1]*u[1] + w[N][2]*u[2];
        //    for (int i = 0; i < 3; i++)
        //        w[0][i] += d_in*u[i];
            //System.out.println("d_in_out, " + d_in + ", " + d_out + ", " + u[0] + ", " + u[1] + ", " + u[2]);
        //}
    }

    private static void shoot(int iT)
    {
        double [] rhs;
        double [] vec;                      // temporary transformed position vector
        double [] vecphi;                   // temporary transformed phi
        double r = Double.parseDouble(txtradius.getText());

        System.out.println("\nstart, " + iT + ", " + eps + ", " + eta + ", " + tau
                                       + ", " + (4.0 + eta*eps*eps) + ", " + 2.0*Math.PI*(1 + eps*eps*tau)/Math.sqrt(3));
        System.out.println("iter, w[0], w[1], w[2]");
        //System.out.println("0, " + w[0][0] + ", " + w[0][1] + ", " + w[0][2]);
        vec = project_2D(w[0]);
        System.out.println("0, " + vec[0] + ", " + vec[1] + ", " + vec[2]);
        if (vec[1] > -r && vec[1] < r && vec[2] > -r && vec[2] < r)
            image.setRGB((int) (image.getWidth()*(1 + vec[1]/r)/2), (int) (image.getHeight()*(1 - vec[2]/r)/2), Color.BLACK.getRGB());
        for (int i = 0; i < N; i++)
        {
            //rhs = new double[] {(-P[(i - 1 + N) % N][0] + 9*P[i][0] + 9*P[(i + 1) % N][0] - P[(i + 2) % N][0])/16,
            //                    (-P[(i - 1 + N) % N][1] + 9*P[i][1] + 9*P[(i + 1) % N][1] - P[(i + 2) % N][1])/16,
            //                    (-P[(i - 1 + N) % N][2] + 9*P[i][2] + 9*P[(i + 1) % N][2] - P[(i + 2) % N][2])/16};
            rhs = new double[] {(P[i][0] + P[i + 1][0])/2,  // this should include 4 data points
                                (P[i][1] + P[i + 1][1])/2,
                                (P[i][2] + P[i + 1][2])/2};
            runge_kutta_Glass_iter(i, rhs);
            //System.out.println((i + 1) + ", " + w[i + 1][0] + ", " + w[i + 1][1] + ", " + w[i + 1][2]);
            vec = project_2D(w[i + 1]);         // state at time 'i + 1'
            System.out.println((i + 1) + ", " + vec[0] + ", " + vec[1] + ", " + vec[2]);
            if (vec[1] > -r && vec[1] < r && vec[2] > -r && vec[2] < r)
                image.setRGB((int) (image.getWidth()*(1 + vec[1]/r)/2), (int) (image.getHeight()*(1 - vec[2]/r)/2), Color.BLACK.getRGB());
        }
        lblImage.repaint();

        // calculate <w, phi>

        double [] tempint = new double[N + 1];      // temporary integrand for Cotes
        for (int i = 0; i < N + 1; i++)
            tempint[i] = phi[i][0]*w[i][0] + phi[i][1]*w[i][1] + phi[i][2]*w[i][2]; // integrate phi[]*w[]
        System.out.println("\nend  , " + iT + ", " + eps + ", " + run_buckle.Cotes_4(tempint)
                                       + ", " + eta + ", " + tau
                                       + ", " + (4.0 + eta*eps*eps) + ", " + 2.0*Math.PI*(1 + eps*eps*tau)/Math.sqrt(3));
        for (int i = 0; i < N + 1; i++)
        {
            vec = project_2D(w[i]);                 // state at time 'i + 1'
            vecphi = project_2D(phi[i]);            // state at time 'i + 1'
            //System.out.println(i + ", " + vec[0] + ", " + vec[1] + ", " + vec[2] + ", " + vecphi[0] + ", " + vecphi[1] + ", " + vecphi[2]);
            tempint[i] = vec[1]*vecphi[1] + vec[2]*vecphi[2]; // integrate phi[]*w[]
        }
        System.out.println("<w, phi> 2-D, " + iT + ", " + run_buckle.Cotes_4(tempint));

        calc_2D_components();                       // check for internal symmetry

        //double z0out =   - w[N][0]/Math.sqrt(3) + w[N][1]/Math.sqrt(3) - w[N][2]/Math.sqrt(3);
        //double z1out = 2.0*w[N][0]/Math.sqrt(3) + w[N][1]/Math.sqrt(3) - w[N][2]/Math.sqrt(3);
        //double z2out =                          - w[N][1]              - w[N][2];
        //System.out.println("z0, z1, z2, z0out, z1out, z2out");
        //System.out.println(z0 + ", " + z1 + ", " + z2 + ", " + z0out + ", " + z1out + ", " + z2out);

        //System.out.println("w[0][0], w[0][1], w[0][2], w[N][0], w[N][1], w[N][2]");
        //System.out.println("[" + w[0][0] + ", " + w[0][1] + ", " + w[0][2] + ", " + w[N][0] + ", " + w[N][1] + ", " + w[N][2] + "],");
        //Point2D.Double pt2init = solve_boundary();      // initiallize w[][]
        //System.out.println("boundary solution, " + pt2init);

        // generate inhomogeneous rhs as a fxn of time: P(N + 1, 3)

        //System.out.println("\nshoot, " + iT + ", " + eps + ", " + eta + ", " + tau);
        //System.out.println("iter, P[0], P[1], P[2], (init P before Eq. 5.18");
        for (int i = 0; i < N + 1; i++)
        {
            // A0
            P[i][0] = eps*tau*( -w[i][0] - 2*w[i][2]);
            P[i][1] = eps*tau*(2*w[i][0] -   w[i][1]);
            P[i][2] = eps*tau*(2*w[i][1] -   w[i][2]);
            // A1
            P[i][0] += eps*eta*(-w[i][2] - tau*eps*(phi[i][2] + eps*w[i][2]))/2;    // re-defined eta
            P[i][1] += eps*eta*( w[i][0] + tau*eps*(phi[i][0] + eps*w[i][0]))/2;
            P[i][2] += eps*eta*( w[i][1] + tau*eps*(phi[i][1] + eps*w[i][1]))/2;
            // Q_R
            P[i][0] += -(1 + eps*eps*tau)*calc_Glass_G(eps*(phi[i][2] + eps*w[i][2]))/eps/eps/eps;
            P[i][1] +=  (1 + eps*eps*tau)*calc_Glass_G(eps*(phi[i][0] + eps*w[i][0]))/eps/eps/eps;
            P[i][2] +=  (1 + eps*eps*tau)*calc_Glass_G(eps*(phi[i][1] + eps*w[i][1]))/eps/eps/eps;
        }

        // calculate eta, tau (mu, T) Eq. 5.18

//        double [] tempint = new double[N + 1];      // temporary integrand for Cotes
//        for (int i = 0; i < N + 1; i++)
//            tempint[i] = P[i][0]*Math.sin(i*2*Math.PI/N) + P[i][1]*Math.sin(i*2*Math.PI/N - Math.PI/3) + P[i][2]*Math.sin(i*2*Math.PI/N - 2*Math.PI/3);
//        double v1 = -Math.sqrt(2.0/3)*run_buckle.Cotes_4(tempint);
//        for (int i = 0; i < N + 1; i++)
//            tempint[i] = P[i][0]*Math.cos(i*2*Math.PI/N) + P[i][1]*Math.cos(i*2*Math.PI/N - Math.PI/3) + P[i][2]*Math.cos(i*2*Math.PI/N - 2*Math.PI/3);
//        double v2 = -Math.sqrt(2.0/3)*run_buckle.Cotes_4(tempint);
//        eta = 4.0*v1;                           // page 9 of "Implementation of Glass Model"
//        tau = v2/Math.sqrt(3) - v1;

        double [] arg_P1_cos = new double[N + 1];      // temporary integrand for Cotes
        double [] arg_P1_sin = new double[N + 1];
        double [] arg_P2_cos = new double[N + 1];
        double [] arg_P2_sin = new double[N + 1];
        for (int i = 0; i < N + 1; i++)
        {
            vec = project_2D(P[i]);                     // projected forcing function at time 'i'
            arg_P1_cos[i] = vec[1]*Math.cos(i*2*Math.PI/N);
            arg_P1_sin[i] = vec[1]*Math.sin(i*2*Math.PI/N);
            arg_P2_cos[i] = vec[2]*Math.cos(i*2*Math.PI/N);
            arg_P2_sin[i] = vec[2]*Math.sin(i*2*Math.PI/N);
        }
        //double v1 = (run_buckle.Cotes_4(arg_P1_cos) - run_buckle.Cotes_4(arg_P2_sin))/beta;
        //double v2 = (run_buckle.Cotes_4(arg_P1_sin) + run_buckle.Cotes_4(arg_P2_cos))/beta;
        double v1 = (run_buckle.Cotes_4(arg_P1_cos) - run_buckle.Cotes_4(arg_P2_sin))/Math.sqrt(2);
        double v2 = (run_buckle.Cotes_4(arg_P1_sin) + run_buckle.Cotes_4(arg_P2_cos))/Math.sqrt(2);
        eta = -4.0*v2;                           // page 15 of "Implementation of Glass Model III"
        tau = -v1/beta + v2;

        // update P[i][j] to produce final rhs

        //System.out.println("\nshoot, " + iT + ", " + eps + ", " + eta + ", " + tau);
        //System.out.println("transformed P, " + iT + ", " + eps + ", " + eta + ", " + tau);
        //System.out.println("iter, P[0], P[1], P[2], (init P after Eq. 5.18");
        for (int i = 0; i < N + 1; i++)
        {
            P[i][0] = eps*(P[i][0] + tau*( -phi[i][0] - 2*phi[i][2]) - eta*phi[i][2]/2)/beta;    // re-defined eta
            P[i][1] = eps*(P[i][1] + tau*(2*phi[i][0] -   phi[i][1]) + eta*phi[i][0]/2)/beta;
            P[i][2] = eps*(P[i][2] + tau*(2*phi[i][1] -   phi[i][2]) + eta*phi[i][1]/2)/beta;
            //System.out.println(i + ", " + P[i][0] + ", " + P[i][1] + ", " + P[i][2]);
            // transformed P
            vec = project_2D(P[i]);                 // forcing function at time 'i'
            //System.out.println(i + ", " + vec[0] + ", " + vec[1] + ", " + vec[2]);
        }
        // solve boundary-value problem for w[i][j]
        //Point2D.Double pt2init = solve_boundary();      // initiallize w[][]
        //System.out.println("boundary solution, " + pt2init);
    }

    public static double [] project_2D(double [] coord)
    {
        // use the transform T_inv to project original 3-D coord
        // onto 2_D normalized plane for output
        // see "Implementation of Glass Model III)
        double [] ret = new double[3];
        for (int i = 0; i < 3; i++)
        {
            ret[i] = 0;
            for (int j = 0; j < 3; j++)
                ret[i] += T_inv[i][j]*coord[j];
        }
        return ret;
    }

    private static double [] back_transform_3D(double x, double y, double z)
    {
        // use the transform T to project back from a 2-D space (with z' fixed or 0)
        // to the original 3-D space
        // see "Implementation of Glass Model III)
        double [] ret = new double[3];
        for (int i = 0; i < 3; i++)
            ret[i] = T[i][0]*x + T[i][1]*y + T[i][2]*z;
        return ret;
    }

    private static void calc_2D_components()
    {
        // calculate 2-D projections of non-linear Glass_G
        // to check for internal symmetries

        double [] vec;
        double [] arg_w1_sin = new double[N + 1];
        double [] arg_w1_cos = new double[N + 1];      // temporary integrand for Cotes
        double [] arg_w2_sin = new double[N + 1];
        double [] arg_w2_cos = new double[N + 1];
        //System.out.println("i, phi_0, phi_1, phi_2");
        System.out.println("i, w_0, w_1, w_2");
        //System.out.println("i, Glass_0, Glass_1, Glass_2");
        for (int i = 0; i < N + 1; i++)
        {
            //vec = project_2D(new double[] {phi[i][0], phi[i][1], phi[i][2]});
            vec = project_2D(new double[] {w[i][0], w[i][1], w[i][2]});
            //vec = project_2D(new double[] {calc_Glass_G(eps*(phi[i][0] + eps*w[i][0])),
            //                               calc_Glass_G(eps*(phi[i][1] + eps*w[i][1])),
            //                               calc_Glass_G(eps*(phi[i][2] + eps*w[i][2]))});
            arg_w1_sin[i] = vec[1]*Math.sin(i*2*Math.PI/N);
            arg_w1_cos[i] = vec[1]*Math.cos(i*2*Math.PI/N);
            arg_w2_sin[i] = vec[2]*Math.sin(i*2*Math.PI/N);
            arg_w2_cos[i] = vec[2]*Math.cos(i*2*Math.PI/N);
//            System.out.println(i + ", " + vec[0] + ", " + vec[1] + ", " + vec[2]);
        }
        System.out.println("w1_2_sin_cos, " + run_buckle.Cotes_4(arg_w1_sin) + ", " +  + run_buckle.Cotes_4(arg_w1_cos)
                                     + ", " + run_buckle.Cotes_4(arg_w2_sin) + ", " +  + run_buckle.Cotes_4(arg_w2_cos)
                                     + ", " + (run_buckle.Cotes_4(arg_w1_sin) + run_buckle.Cotes_4(arg_w2_cos)));
    }

/*
    private static double [] calc_lateral()
    {
        // calculate a unit vector perpendicular to null vector and velocity (cross product)
        double [] unull = new double[] {1/Math.sqrt(3), -1/Math.sqrt(3), 1/Math.sqrt(3)};   // null vector
        double [] uvel  = new double[] {P[N][0] +  -w[N][0]           - 2*w[N][2]
                                       ,P[N][1] + 2*w[N][0] - w[N][1]
                                       ,P[N][2] +           2*w[N][1] - w[N][2]};           // velocity vector
        double norm = uvel[0]*uvel[0] + uvel[1]*uvel[1] + uvel[2]*uvel[2];
        for (int i = 0; i < 3; i++)
            uvel[i] = uvel[i]/Math.sqrt(norm);
        return new double[] {unull[1]*uvel[2] - unull[2]*uvel[1], unull[2]*uvel[0] - unull[0]*uvel[2], unull[0]*uvel[1] - unull[1]*uvel[0]};
    }

    private static Point2D.Double solve_boundary()
    {
        double [] tempint = new double[N + 1];      // temporary integrand for Cotes

        // produce a cyclic solution w[N + 1, 3]

        System.out.println("\nboundary solution, " + N + ", " + eps + ", " + eta + ", " + tau);
        for (int i = 0; i < N + 1; i++)
            tempint[i] = (-P[i][1] - P[i][2])*Math.sin(i*2*Math.PI/N);
        g2_sin = run_buckle.Cotes_4(tempint);
        for (int i = 0; i < N + 1; i++)
            tempint[i] = (-P[i][1] - P[i][2])*Math.cos(i*2*Math.PI/N);
        g2_cos = run_buckle.Cotes_4(tempint);
        for (int i = 0; i < N + 1; i++)
            tempint[i] = (2*P[i][0] + P[i][1] - P[i][2])*Math.sin(i*2*Math.PI/N);
        g1_sin = run_buckle.Cotes_4(tempint)/Math.sqrt(3);
        for (int i = 0; i < N + 1; i++)
            tempint[i] = (2*P[i][0] + P[i][1] - P[i][2])*Math.cos(i*2*Math.PI/N);
        g1_cos = run_buckle.Cotes_4(tempint)/Math.sqrt(3);
        System.out.println("g1_g2_sin_cos, " + g2_sin + ", " + g1_cos + ", " + g2_cos + ", " + g1_sin);
        System.out.println("cos_sin_coef , " + (g2_sin - g1_cos) + ", " + (g2_cos + g1_sin));
        return new Point2D.Double(0, 1);
    }
*/
    private static void runge_kutta_Glass_iter(int iter, double rhs[])
    {
        // Glass model as per Langford 1977
        // see p. 184, "Feedback-Inhibition Oscillator" iterative scheme
        // (x, y, z) = limit cycle
        // rhs = interpolated value of time-dependemt forcing function
        double delt = 2*Math.PI/N/beta;
        double x = w[iter][0];
        double y = w[iter][1];
        double z = w[iter][2];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;
        double m1, m2, m3, m4;

        k1 = delt*(-x - 2*z + P[iter][0]);
        l1 = delt*(2*x - y  + P[iter][1]);
        m1 = delt*(2*y - z  + P[iter][2]);

        k2 = delt*( -(x + k1/2) - 2*(z + m1/2) + rhs[0]);
        l2 = delt*(2*(x + k1/2) -   (y + l1/2) + rhs[1]);
        m2 = delt*(2*(y + l1/2) -   (z + m1/2) + rhs[2]);

        k3 = delt*( -(x + k2/2) - 2*(z + m2/2) + rhs[0]);
        l3 = delt*(2*(x + k2/2) -   (y + l2/2) + rhs[1]);
        m3 = delt*(2*(y + l2/2) -   (z + m2/2) + rhs[2]);

        k4 = delt*( -(x + k3)   - 2*(z + m3)   + P[iter + 1][0]);
        l4 = delt*(2*(x + k3)   -   (y + l3)   + P[iter + 1][1]);
        m4 = delt*(2*(y + l3)   -   (z + m3)   + P[iter + 1][2]);

        w[iter + 1][0] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        w[iter + 1][1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
        w[iter + 1][2] = z + (m1 + 2*m2 + 2*m3 + m4)/6;
    }

    private static double calc_Glass_G(double y)
    {
        // nonlinear response of Glass model: Langford_1977_Numerical_Solution
        double temp = Math.pow(1.0 + 2.0*y, 4 + eps*eps*eta);
        // subtract the term linear in mu (Eq. 6.2)
        //System.out.println("Glass " + y + ", " + (1.0 + 2*y) + ", " + (0.5*(temp - 1)/(temp + 1) - 2.0*y - 0.5*eps*eps*eta*y));
        return 0.5*(temp - 1)/(temp + 1) - 2.0*y - 0.5*eps*eps*eta*y;
    }

//    private static double [] calc_xyz(double tx, double ty)
//    {
//        // get normal coord from the transformed 2-D plane (transform T)
//        return new double[] {tx/Math.sqrt(3), tx*0.5/Math.sqrt(3) - ty*0.5, -tx*0.5/Math.sqrt(3) - ty*0.5};
//    }

    public static void main(String[] args)
    {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                load_data();                        // initiallize Glass data
                Glass_w dlg = new Glass_w();
            }
        });
    }
}
