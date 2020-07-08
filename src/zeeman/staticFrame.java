
// see: C:\APP\Java\CoreJava\v1ch07\DrawTest

package zeeman;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.io.*;
import java.util.Properties;

public class staticFrame extends JFrame
{
    protected static PrintWriter out = null;
    protected static boolean bBoundary = true;
    private static Properties pgmProp = new Properties();
//    private final JCheckBoxMenuItem ploty = new JCheckBoxMenuItem("Plot (y, theta)");
//    private final JCheckBoxMenuItem plotF = new JCheckBoxMenuItem("Plot F");
    private final JCheckBoxMenuItem writeData = new JCheckBoxMenuItem("Write to file");
    private staticComponent component = new staticComponent();
    private static JRadioButtonMenuItem[] distanceA = new JRadioButtonMenuItem[4];
    private static JRadioButtonMenuItem[] keyspeed = new JRadioButtonMenuItem[5];

    public staticFrame()
    {
        JMenu fileMenu = new JMenu("   File");
        fileMenu.setMnemonic('F');
        JMenuItem exitItem = fileMenu.add(new AbstractAction("Exit")
        {
            public void actionPerformed(ActionEvent event)
            {
                save_prefs();
                System.exit(0);
            }
        });
        exitItem.setAccelerator(KeyStroke.getKeyStroke("ctrl X"));

        JMenu configMenu = new JMenu("   Configure");
        configMenu.setMnemonic('C');
        JMenu distanceHdr = new JMenu("Fixed Distance ");
        ButtonGroup distanceGroup = new ButtonGroup();
        for (int i = 0; i < distanceA.length; i++)
        {
            distanceA[i] = new JRadioButtonMenuItem("" + (i + 3));
            distanceA[i].addActionListener(new AbstractAction()
            {
                public void actionPerformed(ActionEvent event)
                {
                    staticComponent.A = Double.parseDouble(event.getActionCommand());
                    component.setSize(component.getWidth() + 1, component.getHeight());
                    //System.out.println(component.A);
                }
            });
            distanceGroup.add(distanceA[i]);
            distanceHdr.add(distanceA[i]);
        }
        JMenu speedHdr = new JMenu("Key Speed ");
        ButtonGroup speedGroup = new ButtonGroup();
        for (int i = 0; i < keyspeed.length; i++)
        {
            keyspeed[i] = new JRadioButtonMenuItem("" + (i + 1));
            keyspeed[i].addActionListener(new AbstractAction()
            {
                public void actionPerformed(ActionEvent event)
                {
                    staticComponent.keyincr = Double.parseDouble(event.getActionCommand());
                }
            });
            speedGroup.add(keyspeed[i]);
            speedHdr.add(keyspeed[i]);
        }
        configMenu.add(distanceHdr);
        configMenu.add(speedHdr);

        JMenu prefsMenu = new JMenu("   Preferences");
        prefsMenu.setMnemonic('P');
        final JCheckBoxMenuItem ploty = new JCheckBoxMenuItem("Plot (y, theta)");
        ploty.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (ploty.isSelected())
                {
                    staticComponent.plt_y_dlg = new Plot_y_Dialog(getWidth(), getHeight(), Toolkit.getDefaultToolkit().getImage(main.class.getResource("images/icon.gif")), staticComponent.x, staticComponent.y);
                    component.addMouseListener(staticComponent.plt_y_pnl);
                    component.addMouseMotionListener(staticComponent.plt_y_pnl);
                    component.addKeyListener(staticComponent.plt_y_pnl);
                    staticComponent.plt_y_dlg.addWindowListener(new WindowAdapter() {
                        @Override public void windowClosing(WindowEvent ev)
                            {ploty.setSelected(!ploty.isSelected());}
                    });
                }
                else
                    staticComponent.plt_y_dlg.dispose();
            }
        });
        ploty.setAccelerator(KeyStroke.getKeyStroke("ctrl P"));

        final JCheckBoxMenuItem plotF = new JCheckBoxMenuItem("Plot F");
        plotF.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (plotF.isSelected())
                {
                    staticComponent.plt_F_dlg = new Plot_F_Dialog(Toolkit.getDefaultToolkit().getImage(main.class.getResource("images/icon.gif")), staticComponent.theta, staticComponent.A, (staticComponent.x - component.getWidth()/2)*staticComponent.A/0.4/component.getHeight(), (staticComponent.y - 0.4*component.getHeight())*staticComponent.A/0.4/component.getHeight());
                    component.addMouseListener(staticComponent.plt_F_pnl);
                    component.addMouseMotionListener(staticComponent.plt_F_pnl);
                    component.addKeyListener(staticComponent.plt_F_pnl);
                    staticComponent.plt_F_dlg.addWindowListener(new WindowAdapter() {
                        @Override public void windowClosing(WindowEvent ev)
                            {plotF.setSelected(!plotF.isSelected());}
                    });
                }
                else
                    staticComponent.plt_F_dlg.dispose();
            }
        });
        plotF.setAccelerator(KeyStroke.getKeyStroke("ctrl F"));

        writeData.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (writeData.isSelected())
                    write_static_data();
                else
                    out.close();
            }
        });
        writeData.setAccelerator(KeyStroke.getKeyStroke("ctrl W"));
        final JCheckBoxMenuItem showBoundary = new JCheckBoxMenuItem("Show Boundary");
        showBoundary.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                bBoundary = showBoundary.isSelected();
                repaint();
            }
        });
        final JCheckBoxMenuItem phase = new JCheckBoxMenuItem("Phase Space");
        phase.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (phase.isSelected())
                {
                    staticComponent.phase_dlg = new PhaseSpace(Toolkit.getDefaultToolkit().getImage(main.class.getResource("images/icon.gif")), (staticComponent.x - component.getWidth()/2)*staticComponent.A/0.4/component.getHeight(), (staticComponent.y - 0.4*component.getHeight())*staticComponent.A/0.4/component.getHeight());
                    staticComponent.phase_dlg.addWindowListener(new WindowAdapter() {
                        @Override public void windowClosing(WindowEvent ev)
                            {phase.setSelected(!phase.isSelected());}
                    });
                }
                else
                    staticComponent.phase_dlg.dispose();
            }
        });
        final JCheckBoxMenuItem bifurcate = new JCheckBoxMenuItem("Bifurcate Diagram");
        bifurcate.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (bifurcate.isSelected())
                {
                    staticComponent.bifurcate_dlg = new Bifurcate(Toolkit.getDefaultToolkit().getImage(main.class.getResource("images/icon.gif")), (staticComponent.x - component.getWidth()/2)*staticComponent.A/0.4/component.getHeight());
                    staticComponent.bifurcate_dlg.addWindowListener(new WindowAdapter() {
                        @Override public void windowClosing(WindowEvent ev)
                            {bifurcate.setSelected(!bifurcate.isSelected());}
                    });
                }
                else
                    staticComponent.bifurcate_dlg.dispose();
            }
        });
        prefsMenu.add(ploty);
        prefsMenu.add(plotF);
        prefsMenu.addSeparator();
        prefsMenu.add(writeData);
        prefsMenu.addSeparator();
        prefsMenu.add(showBoundary);
        showBoundary.setSelected(bBoundary);
        prefsMenu.addSeparator();
        prefsMenu.add(phase);
        prefsMenu.add(bifurcate);

        JMenu helpMenu = new JMenu("   Help");
        helpMenu.setMnemonic('H');
        helpMenu.add(new AbstractAction("System Info")
        {
            public void actionPerformed(ActionEvent event)
            {
                JOptionPane.showMessageDialog(staticFrame.this, main.getInfo(), " Zeeman Catastrophe Machine v" + main.VERSION_NO + " System Info ", JOptionPane.INFORMATION_MESSAGE);
            }
        });
        JMenuItem aboutItem = helpMenu.add(new AbstractAction("About")
        {
            public void actionPerformed(ActionEvent event)
            {
                ZeemanAbout.showDialog(staticFrame.this);
            }
        });
        aboutItem.setAccelerator(KeyStroke.getKeyStroke("ctrl A"));

        JMenuBar menuBar = new JMenuBar();
        menuBar.setLayout(new BoxLayout(menuBar, BoxLayout.X_AXIS));
        setJMenuBar(menuBar);
        menuBar.add(fileMenu);
        menuBar.add(configMenu);
        menuBar.add(prefsMenu);
        menuBar.add(helpMenu);

        addWindowListener(new WindowAdapter() {
            @Override public void windowClosing(WindowEvent ev) {
                save_prefs();
            }
        });

        final Toolkit kit = Toolkit.getDefaultToolkit();
        Dimension screenSize = kit.getScreenSize();
        setSize(screenSize.height/2, 9*screenSize.height/10);
        //setLocationByPlatform(true);
        setIconImage(kit.getImage(main.class.getResource("images/icon.gif")));
        setTitle("Zeeman Catastrophe Machine v" + main.VERSION_NO);
        add(component);
        try                                         // recall program properties
        {
            if (new File(System.getProperty("user.home"), "ZCMPrefs.ini").exists())
            {
                pgmProp.load(new FileInputStream(new File(System.getProperty("user.home"), "ZCMPrefs.ini")));
                staticComponent.A = Double.parseDouble(pgmProp.getProperty("initA", "4"));
                staticComponent.x = Double.parseDouble(pgmProp.getProperty("initx", "200"));
                staticComponent.y = Double.parseDouble(pgmProp.getProperty("inity", "400"));
                staticComponent.keyincr = Double.parseDouble(pgmProp.getProperty("keyincr", "1"));
                PhaseSpace.c = Double.parseDouble(pgmProp.getProperty("c", "1"));
                PhaseSpace.x0 = Double.parseDouble(pgmProp.getProperty("x0", "0"));
                PhaseSpace.w0 = Double.parseDouble(pgmProp.getProperty("w0", "0"));
                PhaseSpace.Tx = Double.parseDouble(pgmProp.getProperty("Tx", "1"));
                PhaseSpace.phi0 = Double.parseDouble(pgmProp.getProperty("phi0", "0"));
                PhaseSpace.NLimit = Integer.parseInt(pgmProp.getProperty("NLimit", "0"));
                Bifurcate.c = Double.parseDouble(pgmProp.getProperty("c", "1"));
                Bifurcate.x0 = Double.parseDouble(pgmProp.getProperty("x0", "0"));
                Bifurcate.ymin = Double.parseDouble(pgmProp.getProperty("ymin", "0"));
                Bifurcate.ymax = Double.parseDouble(pgmProp.getProperty("ymax", "2"));
                Bifurcate.w0 = Double.parseDouble(pgmProp.getProperty("w0", "0"));
                Bifurcate.Tx = Double.parseDouble(pgmProp.getProperty("Tx", "1"));
                Bifurcate.phi0 = Double.parseDouble(pgmProp.getProperty("phi0", "0"));
            }
            else
            {
                staticComponent.A = 4;
                staticComponent.x = getWidth()/2;
                staticComponent.y = 8.5*getHeight()/10;
                staticComponent.keyincr = 1;
                PhaseSpace.c = 1;
                PhaseSpace.x0 = 0;
                PhaseSpace.w0 = 0;
                PhaseSpace.Tx = 1;
                PhaseSpace.phi0 = 0;
                PhaseSpace.NLimit = 0;
                Bifurcate.c = 1;
                Bifurcate.x0 = 0;
                Bifurcate.ymin = 0;
                Bifurcate.ymax = 2;
                Bifurcate.w0 = 0;
                Bifurcate.Tx = 1;
                Bifurcate.phi0 = 0;
            }
        }
        catch (IOException e)
            {System.out.println("error reading ZCMPrefs.ini : " + e);}
        distanceA[(int) staticComponent.A - 3].setSelected(true);
        keyspeed[(int) staticComponent.keyincr - 1].setSelected(true);
    }

    private void write_static_data()
    {
        // generate a file consisting of (x, y, theta)
        if (new File(System.getProperty("user.home"), "ZCM_Output.csv").exists())
        {
            if (JOptionPane.showConfirmDialog(staticFrame.this, "The data file '" + System.getProperty("user.home") + System.getProperty("file.separator") + "ZCM_Output.csv' already exists. \n Do you wish to overwrite it ?", " Write to File ", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION)
            {
                writeData.setSelected(false);
                return;
            }
        }
        else
            if (JOptionPane.showConfirmDialog(staticFrame.this, "This will create the data file '" + System.getProperty("user.home") + System.getProperty("file.separator") + "ZCM_Output.csv' \n Continue ?", " Write to File ", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION)
            {
                writeData.setSelected(false);
                return;
            }
        try
        {
            FileWriter fw = new FileWriter(System.getProperty("user.home") + System.getProperty("file.separator") + "ZCM_Output.csv", false);
            out = new PrintWriter(fw);
            out.println("Zeeman Static Machine v" + main.VERSION_NO);
            out.println("A, x, y, theta");
        }
        catch (java.io.IOException e)
            {System.out.println("ZCM_Output.csv save error = " + e);}
    }

    private void save_prefs()
    {
        pgmProp.setProperty("initA", "" + staticComponent.A);
        if (staticComponent.phase_dlg != null && PhaseSpace.xa != 0 && PhaseSpace.y0 != 0)
        {
            pgmProp.setProperty("initx", "" + (0.4*component.getHeight()*PhaseSpace.xa/staticComponent.A + component.getWidth()/2));
            pgmProp.setProperty("inity", "" + 0.4*component.getHeight()*(1 + PhaseSpace.y0/staticComponent.A));
        }
        else if (staticComponent.bifurcate_dlg != null && Bifurcate.xa != 0 && Bifurcate.ymin != 0 && Bifurcate.ymax != 0)
        {
            pgmProp.setProperty("ymin", "" + Bifurcate.ymin);
            pgmProp.setProperty("ymax", "" + Bifurcate.ymax);
        }
        else
        {
            pgmProp.setProperty("initx", "" + staticComponent.x);
            pgmProp.setProperty("inity", "" + staticComponent.y);
        }
        pgmProp.setProperty("keyincr", "" + staticComponent.keyincr);
        pgmProp.setProperty("c", "" + PhaseSpace.c);
        pgmProp.setProperty("x0", "" + PhaseSpace.x0);
        pgmProp.setProperty("w0", "" + PhaseSpace.w0);
        pgmProp.setProperty("Tx", "" + PhaseSpace.Tx);
        pgmProp.setProperty("phi0", "" + PhaseSpace.phi0);
        pgmProp.setProperty("NLimit", "" + PhaseSpace.NLimit);
        try
            {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "ZCMPrefs.ini"), "Zeeman Catastrophe Machine v" + main.VERSION_NO + " Prefs");}
        catch (IOException e)
            {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
        if (out != null)
            out.close();
    }
}

class staticComponent extends JComponent
{
    protected static double A;                                     // fixed point A
    protected static Plot_y_Dialog plt_y_dlg;
    protected static Plot_y_Panel plt_y_pnl;
    protected static Plot_F_Dialog plt_F_dlg;
    protected static Plot_F_Panel plt_F_pnl;
    protected static PhaseSpace phase_dlg;
    protected static Bifurcate bifurcate_dlg;
    protected static double x, y, theta = Math.PI;
    protected static double keyincr;
    private double width, height = Double.NaN;
    private double radius;
    private int x1, x2, y1, y2;
    private Ellipse2D.Double drag;                          // variable point B

    public staticComponent()
    {
        addMouseListener(new MouseAdapter()
        {
            @Override public void mousePressed(MouseEvent e)
            {
                //System.out.println(e.getX() + ", " + e.getY() + ", " + e.getButton() + ", " + e.getModifiers());
                if (getCursor().getType() == Cursor.CROSSHAIR_CURSOR)
                {
                    x = e.getX();
                    y = e.getY();
                    theta = main.solve_for_critical(theta, A, (x - width/2)/radius, (y - 0.4*height)/radius);
                    if (plt_y_pnl != null)
                    {
                        plt_y_pnl.theta = 180*theta/Math.PI - 90;
                        if (plt_y_pnl.theta > 180)
                            plt_y_pnl.theta -= 360;
                        plt_y_pnl.theta += 180; //pltpnl.getHeight()/2;
                    }
                    //System.out.println(width + ", " + height + ", " + radius + ", " + A + ", " + x + ", " + y + ", " + 180*theta/Math.PI);
                    repaint();
                }
            }
            @Override public void mouseReleased(MouseEvent e)
            {
                setCursor(Cursor.getDefaultCursor());
                if (plt_y_pnl != null)
                    plt_y_pnl.cursor = Cursor.DEFAULT_CURSOR;
                if (plt_F_pnl != null)
                    plt_F_pnl.cursor = Cursor.DEFAULT_CURSOR;
            }
        });

        addMouseMotionListener(new MouseMotionListener()
        {
            public void mouseMoved(MouseEvent e)
            {
                if (drag.contains(new Point2D.Double(e.getX(), e.getY())))
                {
                    setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
                    if (plt_y_pnl != null)
                        plt_y_pnl.cursor = Cursor.CROSSHAIR_CURSOR;
                    if (plt_F_pnl != null)
                        plt_F_pnl.cursor = Cursor.CROSSHAIR_CURSOR;
                }
                else
                {
                    setCursor(Cursor.getDefaultCursor());
                    if (plt_y_pnl != null)
                        plt_y_pnl.cursor = Cursor.DEFAULT_CURSOR;
                    if (plt_F_pnl != null)
                        plt_F_pnl.cursor = Cursor.DEFAULT_CURSOR;
                }
            }

            public void mouseDragged(MouseEvent e)
            {
                if (getCursor().getType() == Cursor.CROSSHAIR_CURSOR)
                {
                    x = e.getX();
                    y = e.getY();
                    theta = main.solve_for_critical(theta, A, (x - width/2)/radius, (y - 0.4*height)/radius);
                    if (plt_y_pnl != null && plt_y_pnl.isShowing())
                    {
                        plt_y_pnl.theta = 180*theta/Math.PI - 90;
                        if (plt_y_pnl.theta > 180)
                            plt_y_pnl.theta -= 360;
                        plt_y_pnl.theta += 180; //pltpnl.getHeight()/2;
                        //System.out.println("dragged height = " + pltpnl.getSize().height);
                    }
                    if (plt_F_pnl != null && plt_F_pnl.isShowing())
                    {
                        plt_F_pnl.x = (x - width/2)/radius;
                        plt_F_pnl.y = (y - 0.4*height)/radius;
                        plt_F_pnl.theta = theta;
                    }
                    repaint();
                }
            }
        });

        addKeyListener(new KeyAdapter()
        {
            @Override public void keyPressed(KeyEvent e)
            {
                switch (e.getKeyCode())
                {
                case KeyEvent.VK_LEFT:
                    x -= keyincr;
                break;
                case KeyEvent.VK_RIGHT:
                    x += keyincr;
                break;
                case KeyEvent.VK_UP:
                    y -= keyincr;
                break;
                case KeyEvent.VK_DOWN:
                    y += keyincr;
                break;
                }
                theta = main.solve_for_critical(theta, A, (x - width/2)/radius, (y - 0.4*height)/radius);
                if (plt_y_pnl != null && plt_y_pnl.isShowing())
                {
                    plt_y_pnl.x = x;
                    plt_y_pnl.y = y;
                    plt_y_pnl.theta = 180*theta/Math.PI - 90;
                    if (plt_y_pnl.theta > 180)
                        plt_y_pnl.theta -= 360;
                    plt_y_pnl.theta += 180; //pltpnl.getHeight()/2;
                }
                if (plt_F_pnl != null && plt_F_pnl.isShowing())
                {
                    plt_F_pnl.x = (x - width/2)/radius;
                    plt_F_pnl.y = (y - 0.4*height)/radius;
                    plt_F_pnl.theta = theta;
                }
                repaint();
            }
        });

        addComponentListener(new ComponentAdapter() {
            @Override public void componentResized(ComponentEvent e) {
                //System.out.println("start resized " + x + ", " + y + ", " + width + ", " + height + ", " + getWidth() + ", " + getHeight());
                radius = 0.4*getHeight()/A;
                x += (getWidth() - width)/2;
                width = getWidth();
                if (Double.isNaN(height))
                {
                    x -= getWidth()/2;          // restore original x
                    theta = main.solve_for_critical(theta, A, (x - width/2)/radius, (y - 0.4*getHeight())/radius);
                }
                else
                {
                    x  = width/2 + (x - width/2)*getHeight()/height;
                    y *= getHeight()/height;
                }
                height = getHeight();
                setFont(new Font(Font.SERIF, Font.BOLD, (int) (2*radius)));
                //System.out.println("end   resized " + x + ", " + y + ", " + width + ", " + height + ", " + getWidth() + ", " + getHeight() + ", " + theta);
                grabFocus();
            }
        });
    }

    @Override protected void paintComponent(Graphics g)
    {
        //System.out.println("paint " + width + ", " + height + ", " + getWidth() + ", " + getHeight() + ", " + 180*theta/Math.PI);
        Graphics2D g2 = (Graphics2D) g;
        g2.fill(new Ellipse2D.Double(width/2 - 3.5, -3.5, 7, 7));                       // fixed point A
        Ellipse2D circle = new Ellipse2D.Double(width/2 - radius, 0.4*height - radius, 2*radius, 2*radius);
        Ellipse2D wheel = new Ellipse2D.Double(width/2 - 3.5 + radius*Math.sin(theta), 0.4*height - 3.5 - radius*Math.cos(theta), 7, 7);
        drag = new Ellipse2D.Double(x - 3.5, y - 3.5, 7, 7);   // variable point B
        Line2D lineA = new Line2D.Double(width/2, 0, width/2 + radius*Math.sin(theta), 0.4*height - radius*Math.cos(theta));
        Line2D lineB = new Line2D.Double(width/2 + radius*Math.sin(theta), 0.4*height - radius*Math.cos(theta), x, y);
        g2.setPaint(new Color(204, 204, 51));
        g2.fill(circle);
        g2.setStroke(new BasicStroke(3));
        g2.setPaint(new Color(160, 160, 40));
        g2.draw(circle);
        g2.setPaint(new Color(0, 0, 224));
        g2.draw(drag);
        g2.setStroke(new BasicStroke(1));
        g2.setPaint(new Color(0, 0, 0));
        g2.fill(wheel);
        g2.draw(lineA);
        g2.draw(lineB);
        if (staticFrame.bBoundary)
        {
            x1 = (int) (width/2 + radius*main.xbound[(int) A - 3][0]);
            y1 = (int) (0.4*height + radius*main.ybound[(int) A - 3][0]);
            for (int i = 1; i < main.xbound[(int) A - 3].length; i++)
            {
                x2 = (int) (width/2 + radius*main.xbound[(int) A - 3][i]);
                y2 = (int) (0.4*height + radius*main.ybound[(int) A - 3][i]);
                g2.drawLine(x1, y1, x2, y2);
                g2.drawLine((int) width - x1, y1, (int) width - x2, y2);
                x1 = x2;
                y1 = y2;
            }
        }
        g2.rotate(theta + Math.PI, width/2, 0.4*height);
        g2.drawString("Z", (float) (width/2 - 30*radius/48), (float) (0.4*height + 32*radius/48));
        if (staticFrame.out != null)
            staticFrame.out.println(A + ", " + x + ", " + y + ", " + theta*180/Math.PI);
    }
}
