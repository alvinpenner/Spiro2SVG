
package chua;

/*
 * this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Simul_5.java
 * see jdk file : \APP\Java\Demos\jdk_Demos\ButtonDemo.java
 *
 * simulate a Poincare map of degree 5 in the x'-y' plane, odd terms only
 * based on standard cubic model, Book IV, page 62
*/

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;
import java.util.Properties;

public class Chua_Simul_5 extends JDialog
{
    private static Properties simulProp = new Properties();
    private static final BufferedImage image = new BufferedImage(480, 480, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static JButton btnCalc = new JButton("Calc");
    private static JButton btnClear = new JButton("Clear");
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
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()}};
    private static JTextField txtstart = new JTextField();
    private static JTextField txtrange = new JTextField();
    private static double x0, y0;
    private static int iT = 0;

    public Chua_Simul_5()
    {
        setTitle("Chua System - Simulate x'-y' Scatter (5th degree)");
        setIconImage(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        setSize(780, 532);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        lblImage.setBorder(BorderFactory.createEtchedBorder());

        final JLabel[] lblarr = {new JLabel(C_or_g + " 50"),
                                 new JLabel(C_or_g + " 41"),
                                 new JLabel(C_or_g + " 32"),
                                 new JLabel(C_or_g + " 23"),
                                 new JLabel(C_or_g + " 14"),
                                 new JLabel(C_or_g + " 05")};
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
        parmsPanel.add(btnClear);
        parmsPanel.add(printPanel);
        parmsPanel.setMaximumSize(new Dimension(250, 3000));
        parmsPanel.setPreferredSize(new Dimension(250, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        final JPanel scatterPanel = new JPanel();
        scatterPanel.setOpaque(false);
        scatterPanel.add(lblImage);

        load_prefs();
        x0 = 0*Double.parseDouble(txtstart.getText());
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

    private static void refresh_graph()
    {
        PrintWriter fout = null;
        double alpha = Double.parseDouble(txtalpha.getText());
        double costheta = Math.cos(Double.parseDouble(txttheta.getText())*Math.PI/180);
        double sintheta = Math.sin(Double.parseDouble(txttheta.getText())*Math.PI/180);
        double a = Double.parseDouble(txta.getText());
        double b = Double.parseDouble(txtb.getText());
        double[][] Carr = new double[6][2];
        double range = Double.parseDouble(txtrange.getText());
        double x, y, xtemp, ytemp;
        int N = 100000;

        if (C_or_g.startsWith("C"))
            for (int i = 0; i < txtCarr.length; i++)            // default to Cxy input
            {
                Carr[i][0] = Double.parseDouble(txtCarr[i][0].getText());   // Cx
                Carr[i][1] = Double.parseDouble(txtCarr[i][1].getText());   // Cy
            }
        else                // override Cxy input with gij input (May 19/23 looseleaf)
        {
            final double[][] coeffR = new double[][] {{1,  0,-10,  0,  5,  0},   // for g50
                                                      {1,  0, -2,  0, -3,  0},
                                                      {1,  0,  2,  0,  1,  0},
                                                      {1,  0,  2,  0,  1,  0},
                                                      {1,  0, -2,  0, -3,  0},
                                                      {1,  0,-10,  0,  5,  0}};  // for g05
            final double[][] coeffI = new double[][] {{0,  5,  0,-10,  0,  1},
                                                      {0,  3,  0,  2,  0, -1},
                                                      {0,  1,  0,  2,  0,  1},
                                                      {0, -1,  0, -2,  0, -1},
                                                      {0, -3,  0, -2,  0,  1},
                                                      {0, -5,  0, 10,  0, -1}};
            System.out.println("Carr from gi");
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
                boolean fexists = new File("C:\\Windows\\Temp\\Chua_Simul_5_" + txtalpha.getText() + "_" + txttheta.getText() + "_" + txta.getText() + "_" + txtb.getText() + ".csv").exists();
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\Chua_Simul_5_" + txtalpha.getText() + "_" + txttheta.getText() + "_" + txta.getText() + "_" + txtb.getText() + ".csv", true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    fout.println("Cxi, " + Carr[0][0] + ", " + Carr[1][0] + ", " + Carr[2][0] + ", " + Carr[3][0] + ", " + Carr[4][0] + ", " + Carr[5][0]);
                    fout.println("Cyi, " + Carr[0][1] + ", " + Carr[1][1] + ", " + Carr[2][1] + ", " + Carr[3][1] + ", " + Carr[4][1] + ", " + Carr[5][1]);
                    fout.println("scatter hdr, NaN, " + (1 + alpha) + ", " + txttheta.getText());
                    fout.println("init x0 y0, 0, 0, NaN, NaN");
                    fout.println("iter       , x', y'");
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}
        image.setRGB(image.getWidth()/2, image.getHeight()/2, Color.BLACK.getRGB());
        for (int i = 0; i < N; i++)
        {
            xtemp = (1 + alpha)*x0 + (x0*x0 + y0*y0)*(a*x0 - b*y0);
            ytemp = (1 + alpha)*y0 + (x0*x0 + y0*y0)*(b*x0 + a*y0);
            xtemp += Carr[0][0]*x0*x0*x0*x0*x0 + Carr[1][0]*x0*x0*x0*x0*y0 + Carr[2][0]*x0*x0*x0*y0*y0 + Carr[3][0]*x0*x0*y0*y0*y0 + Carr[4][0]*x0*y0*y0*y0*y0 + Carr[5][0]*y0*y0*y0*y0*y0;
            ytemp += Carr[0][1]*x0*x0*x0*x0*x0 + Carr[1][1]*x0*x0*x0*x0*y0 + Carr[2][1]*x0*x0*x0*y0*y0 + Carr[3][1]*x0*x0*y0*y0*y0 + Carr[4][1]*x0*y0*y0*y0*y0 + Carr[5][1]*y0*y0*y0*y0*y0;
            x = costheta*xtemp - sintheta*ytemp;
            y = sintheta*xtemp + costheta*ytemp;
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
        int r = (int) (Math.sqrt(x0*x0 + y0*y0)/range*image.getWidth()/2);
        DC.setColor(Color.blue);
        DC.drawLine(0, image.getHeight()/2, image.getWidth(), image.getHeight()/2);
        DC.drawLine(image.getWidth()/2, 0, image.getWidth()/2, image.getHeight());
        DC.drawOval(image.getHeight()/2 - r, image.getWidth()/2 - r, 2*r, 2*r);
        lblImage.repaint();
        if (fout != null)
            fout.close();
    }

    protected static void load_prefs()
    {
        try                                         // recall simulation properties
        {
            if (new File(System.getProperty("user.home"), "ChuaSimul.ini").exists())
            {
                simulProp.load(new FileInputStream(new File(System.getProperty("user.home"), "ChuaSimul.ini")));
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
            {System.out.println("error reading ChuaSimul.ini : " + e);}
    }

    protected static void save_prefs()
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
            {simulProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "ChuaSimul.ini"), "Chua Simulate Degree 5 Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                Chua_Simul_5 dlg = new Chua_Simul_5();
            }
        });
    }
}
