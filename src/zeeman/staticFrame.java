
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
    private final JCheckBoxMenuItem plotData = new JCheckBoxMenuItem("Plot Data");
    private staticComponent component = new staticComponent();
    private static JRadioButtonMenuItem[] distanceA = new JRadioButtonMenuItem[4];

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
        ButtonGroup group = new ButtonGroup();
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
            group.add(distanceA[i]);
            distanceHdr.add(distanceA[i]);
        }
        configMenu.add(distanceHdr);

        JMenu prefsMenu = new JMenu("   Preferences");
        prefsMenu.setMnemonic('P');
//        final JCheckBoxMenuItem plotData = new JCheckBoxMenuItem("Plot Data");
        plotData.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                if (plotData.isSelected())
                {
                    staticComponent.pltdlg = new PlotDialog(getWidth(), getHeight(), Toolkit.getDefaultToolkit().getImage(main.class.getResource("images/icon.gif")), staticComponent.x, staticComponent.y);
                    component.addMouseListener(staticComponent.pltpnl);
                    component.addMouseMotionListener(staticComponent.pltpnl);
                    component.addKeyListener(staticComponent.pltpnl);
                    staticComponent.pltdlg.addWindowListener(new WindowAdapter() {
                        @Override public void windowClosing(WindowEvent ev)
                            {plotData.setSelected(!plotData.isSelected());}
                    });
                }
                else
                    staticComponent.pltdlg.dispose();
            }
        });
        plotData.setAccelerator(KeyStroke.getKeyStroke("ctrl P"));
        final JCheckBoxMenuItem writeData = new JCheckBoxMenuItem("Write to file");
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
        prefsMenu.add(plotData);
        prefsMenu.add(writeData);
        prefsMenu.add(showBoundary);
        showBoundary.setSelected(bBoundary);

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
            }
            else
            {
                staticComponent.A = 4;
                staticComponent.x = getWidth()/2;
                staticComponent.y = 8.5*getHeight()/10;
            }
        }
        catch (IOException e)
            {System.out.println("error reading ZCMPrefs.ini : " + e);}
        distanceA[(int) staticComponent.A - 3].setSelected(true);
    }

    private static void write_static_data()
    {
        // generate a file consisting of (x, y, theta)
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

    private static void save_prefs()
    {
        pgmProp.setProperty("initA", "" + staticComponent.A);
        pgmProp.setProperty("initx", "" + staticComponent.x);
        pgmProp.setProperty("inity", "" + staticComponent.y);
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
    protected static PlotDialog pltdlg;
    protected static PlotPanel pltpnl;
    protected static double x, y, theta = Math.PI;
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
                    if (pltpnl != null)
                    {
                        pltpnl.theta = 180*theta/Math.PI - 90;
                        if (pltpnl.theta > 180)
                            pltpnl.theta -= 360;
                        pltpnl.theta += 180; //pltpnl.getHeight()/2;
                    }
                    //System.out.println(width + ", " + height + ", " + radius + ", " + A + ", " + x + ", " + y + ", " + 180*theta/Math.PI);
                    repaint();
                }
            }
            @Override public void mouseReleased(MouseEvent e)
            {
                setCursor(Cursor.getDefaultCursor());
                if (pltpnl != null)
                    pltpnl.cursor = Cursor.DEFAULT_CURSOR;
            }
        });

        addMouseMotionListener(new MouseMotionListener()
        {
            public void mouseMoved(MouseEvent e)
            {
                if (drag.contains(new Point2D.Double(e.getX(), e.getY())))
                {
                    setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
                    if (pltpnl != null)
                        pltpnl.cursor = Cursor.CROSSHAIR_CURSOR;
                }
                else
                {
                    setCursor(Cursor.getDefaultCursor());
                    if (pltpnl != null)
                        pltpnl.cursor = Cursor.DEFAULT_CURSOR;
                }
            }

            public void mouseDragged(MouseEvent e)
            {
                if (getCursor().getType() == Cursor.CROSSHAIR_CURSOR)
                {
                    x = e.getX();
                    y = e.getY();
                    theta = main.solve_for_critical(theta, A, (x - width/2)/radius, (y - 0.4*height)/radius);
                    if (pltpnl != null)
                    {
                        pltpnl.theta = 180*theta/Math.PI - 90;
                        if (pltpnl.theta > 180)
                            pltpnl.theta -= 360;
                        pltpnl.theta += 180; //pltpnl.getHeight()/2;
                        //System.out.println("dragged height = " + pltpnl.getSize().height);
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
                    x--;
                break;
                case KeyEvent.VK_RIGHT:
                    x++;
                break;
                case KeyEvent.VK_UP:
                    y--;
                break;
                case KeyEvent.VK_DOWN:
                    y++;
                break;
                }
                theta = main.solve_for_critical(theta, A, (x - width/2)/radius, (y - 0.4*height)/radius);
                if (pltpnl != null)
                {
                    pltpnl.theta = 180*theta/Math.PI - 90;
                    if (pltpnl.theta > 180)
                        pltpnl.theta -= 360;
                    pltpnl.theta += 180; //pltpnl.getHeight()/2;
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
