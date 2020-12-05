
package rossler;

// plot y versus x
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import javax.swing.*;
import java.util.GregorianCalendar;

public class Rossler_y_vs_x extends JDialog
{
    private final static int N = 2000;            // total # of iterations 160000
    private static double[] pt3_old;
    protected Path2D.Double path1 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    //protected Path2D.Double path2 = new Path2D.Double(Path2D.WIND_NON_ZERO, N);
    protected Line2D.Double xaxis = new Line2D.Double(0, 0, 0, 0);
    protected Line2D.Double yaxis = new Line2D.Double(0, 0, 0, 0);
    private static JCheckBox printChk = new JCheckBox("  print  ");
    private static JLabel lblxrange;
    private static JLabel lblyrange;
    private static JPanel phasePanel = new Plot_Phase_Panel();

    public Rossler_y_vs_x(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        Main.type = "phase";
        setTitle(" Rossler System - phase y vs. x");
        setIconImage(img);
        setSize(600, 483);
        setLocationByPlatform(true);

        final JLabel[] lbl = {new JLabel("a"),
                              new JLabel("b"),
                              new JLabel("c"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.a)),
                                  new JTextField(Double.toString(Main.b)),
                                  new JTextField(Double.toString(Main.c)),
                                  new JTextField(Double.toString(Main.x0)),
                                  new JTextField(Double.toString(Main.y0)),
                                  new JTextField(Double.toString(Main.z0))};
        JPanel[] spacerPanel = new JPanel[3];
        JPanel[] dataPanel = new JPanel[lbl.length];
        JButton btnRun = new JButton("Run");

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

        printChk.setOpaque(false);
        JPanel printPanel = new JPanel();
        printPanel.setOpaque(false);
        printPanel.add(printChk);

        JPanel xrangePanel = new JPanel();
        xrangePanel.setOpaque(false);
        lblxrange = new JLabel("x = ");
        lblxrange.setPreferredSize(new Dimension(110, 18));
        xrangePanel.add(lblxrange);

        JPanel yrangePanel = new JPanel();
        yrangePanel.setOpaque(false);
        lblyrange = new JLabel("y = ");
        lblyrange.setPreferredSize(new Dimension(110, 18));
        yrangePanel.add(lblyrange);

        parmsPanel.removeAll();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(dataPanel[0]);
        parmsPanel.add(dataPanel[1]);
        parmsPanel.add(dataPanel[2]);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(printPanel);
        parmsPanel.add(spacerPanel[2]);
        parmsPanel.add(xrangePanel);
        parmsPanel.add(yrangePanel);

        parmsPanel.setMaximumSize(new Dimension(140, 3000));
        parmsPanel.setPreferredSize(new Dimension(140, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().add(parmsPanel);
        getContentPane().add(phasePanel);
        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        pt3_old = new double[] {Main.x0, Main.y0, Main.z0};

        btnRun.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                boolean changed = false;
                for (int i = 0; i < txt.length; i++)
                    if (txt[i].getText().isEmpty())
                        txt[i].setText("0");
                Main.a = Double.parseDouble(txt[0].getText());
                Main.b = Double.parseDouble(txt[1].getText());
                Main.c = Double.parseDouble(txt[2].getText());
                if (Main.x0 != Double.parseDouble(txt[3].getText())) changed = true;
                Main.x0 = Double.parseDouble(txt[3].getText());
                if (Main.y0 != Double.parseDouble(txt[4].getText())) changed = true;
                Main.y0 = Double.parseDouble(txt[4].getText());
                if (Main.z0 != Double.parseDouble(txt[5].getText())) changed = true;
                Main.z0 = Double.parseDouble(txt[5].getText());
                phase_space(changed);
                //System.out.println("btnRun " + phasePanel.getSize());
                //System.out.println("btnRun " + parmsPanel.getSize());
            }
        });

        addWindowListener(new WindowAdapter() {
            @Override public void windowClosing(WindowEvent ev)
                {Main.save_prefs();}
            });

        for (JTextField tx: txt)
            tx.addKeyListener(new KeyAdapter()
            {
                @Override public void keyPressed(KeyEvent e)
                {
                    if (e.getKeyCode() == KeyEvent.VK_ESCAPE)
                    {
                        Main.save_prefs();
                        Main.plot_y_vs_x.dispose();
                        Main.z_bifurcate = new Rossler_z_bifurcate(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
                    }
                }
            });
    }

    private void phase_space(boolean ch)
    {
        GregorianCalendar now = new GregorianCalendar();
        double delt = 0.1;
        double[] pt3 = new double[] {Main.x0, Main.y0, Main.z0};
        double xmin, xmax, ymin, ymax;

        System.out.println("Rossler PhaseSpace, " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + delt + ", " + N);
        if (!ch)                                            // re-use previous run
            System.arraycopy(pt3_old, 0, pt3, 0, pt3.length);
        if (printChk.isSelected())
        {
            System.out.println("i, x, y, z");
            System.out.println("0, " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
        }
        xmin = pt3[0];
        xmax = pt3[0];
        ymin = pt3[1];
        ymax = pt3[1];
        path1.reset();                                      // transient path
        path1.moveTo(pt3[0], pt3[1]);
        for (int i = 1; i <= N; i++)
        {
            Main.runge_kutta_rossler(pt3, delt);
            if (printChk.isSelected())
                System.out.println(i + ", " + pt3[0] + ", " + pt3[1] + ", " + pt3[2]);
            if (pt3[0] > xmax) xmax = pt3[0];
            if (pt3[0] < xmin) xmin = pt3[0];
            if (pt3[1] > ymax) ymax = pt3[1];
            if (pt3[1] < ymin) ymin = pt3[1];
            path1.lineTo(pt3[0], pt3[1]);
        }
        System.arraycopy(pt3, 0, pt3_old, 0, pt3.length);               // save
        AffineTransform at = new AffineTransform((phasePanel.getWidth() - 8)/(xmax - xmin),
                                          0, 0, -(phasePanel.getHeight() - 8)/(ymax - ymin),
                                                -xmin*phasePanel.getWidth()/(xmax - xmin),
                                                 ymax*phasePanel.getHeight()/(ymax - ymin));
        path1.transform(at);
        //System.out.println(xmin + ", " + xmax + ", " + ymin + ", " + ymax);
        lblxrange.setText("x = " + String.format("%.3f", xmin) + ", " + String.format("%.3f", xmax));
        lblyrange.setText("y = " + String.format("%.3f", ymin) + ", " + String.format("%.3f", ymax));
        xaxis = new Line2D.Double(0, ymax*phasePanel.getHeight()/(ymax - ymin), phasePanel.getWidth(), ymax*phasePanel.getHeight()/(ymax - ymin));
        yaxis = new Line2D.Double(-xmin*phasePanel.getWidth()/(xmax - xmin), 0, -xmin*phasePanel.getWidth()/(xmax - xmin), phasePanel.getHeight());
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
        g2.draw(Main.plot_y_vs_x.path1);
        //g2.setPaint(new Color(0, 0, 0));
        //g2.setStroke(new BasicStroke(2));
        //g2.draw(Main.plot_y_vs_x.path2);
        g2.setPaint(Color.BLUE);
        g2.draw(Main.plot_y_vs_x.xaxis);
        g2.draw(Main.plot_y_vs_x.yaxis);
    }
}
