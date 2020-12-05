
package rossler;

// bifurcation plot of max z versus c
// implement Runge-Kutta: see Froberg page 269
// see also: Numerical Recipes in C, page 713

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.GregorianCalendar;

public class Rossler_z_bifurcate extends JDialog
{
    private final static int N = 2000;            // total # of iterations 160000
    private static double[] pt3_old;

    public Rossler_z_bifurcate(Image img)
    {
        final JPanel parmsPanel = new JPanel();
        Main.type = "bifurcate";
        setTitle(" Rossler System - z bifurcate vs. c");
        setIconImage(img);
        setSize(600, 483);
        setLocationByPlatform(true);

        final JLabel[] lbl = {new JLabel("a"),
                              new JLabel("b"),
                              new JLabel("c start"),
                              new JLabel("c end"),
                              new JLabel("x0"),
                              new JLabel("y0"),
                              new JLabel("z0")};
        final JTextField[] txt = {new JTextField(Double.toString(Main.a)),
                                  new JTextField(Double.toString(Main.b)),
                                  new JTextField(Double.toString(Main.cstart)),
                                  new JTextField(Double.toString(Main.cend)),
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

        parmsPanel.removeAll();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(dataPanel[0]);
        parmsPanel.add(dataPanel[1]);
        parmsPanel.add(dataPanel[2]);
        parmsPanel.add(dataPanel[3]);
        parmsPanel.add(spacerPanel[0]);
        parmsPanel.add(dataPanel[4]);
        parmsPanel.add(dataPanel[5]);
        parmsPanel.add(dataPanel[6]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnRun);
        parmsPanel.add(spacerPanel[2]);

        parmsPanel.setMaximumSize(new Dimension(140, 3000));
        parmsPanel.setPreferredSize(new Dimension(140, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().add(parmsPanel);
        //getContentPane().add(phasePanel);
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
                Main.cstart = Double.parseDouble(txt[2].getText());
                Main.cend = Double.parseDouble(txt[3].getText());
                if (Main.x0 != Double.parseDouble(txt[4].getText())) changed = true;
                Main.x0 = Double.parseDouble(txt[4].getText());
                if (Main.y0 != Double.parseDouble(txt[5].getText())) changed = true;
                Main.y0 = Double.parseDouble(txt[5].getText());
                if (Main.z0 != Double.parseDouble(txt[6].getText())) changed = true;
                Main.z0 = Double.parseDouble(txt[6].getText());
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
                        Main.z_bifurcate.dispose();
                        Main.plot_y_vs_x = new Rossler_y_vs_x(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
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

        System.out.println("Rossler Bifurcate, " + now.getTime() + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + delt + ", " + N);
        if (!ch)                                            // re-use previous run
            System.arraycopy(pt3_old, 0, pt3, 0, pt3.length);
        xmin = pt3[0];
        xmax = pt3[0];
        ymin = pt3[1];
        ymax = pt3[1];
        for (int i = 1; i <= N; i++)
        {
            Main.runge_kutta_rossler(pt3, delt);
            if (pt3[0] > xmax) xmax = pt3[0];
            if (pt3[0] < xmin) xmin = pt3[0];
            if (pt3[1] > ymax) ymax = pt3[1];
            if (pt3[1] < ymin) ymin = pt3[1];
        }
        System.arraycopy(pt3, 0, pt3_old, 0, pt3.length);               // save
        //System.out.println(xmin + ", " + xmax + ", " + ymin + ", " + ymax);
    }
}
