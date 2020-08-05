
package zeeman;

// plot y versus x
// plot theta versus x
// see : C:\APP\Java\CoreJava\v1_Edition7\v1ch8\Sketch

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import javax.swing.*;

public class Plot_y_Dialog extends JDialog
{
    //private static final BufferedImage img = new BufferedImage(300, 300, BufferedImage.TYPE_3BYTE_BGR);
    //private static final Graphics2D DC = img.createGraphics();
    private JButton btnClear = new JButton("Clear");

    public Plot_y_Dialog(int width, int height, Image img, double xorg, double yorg)
    {
        setTitle(" Plot Static Zeeman Data - (y, theta) vs. x");
        setIconImage(img);
        setSize(width, height);
        setLocationByPlatform(true);
        setVisible(true);
        staticComponent.plt_y_pnl = new Plot_y_Panel(xorg, yorg);
        btnClear.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)          // btnOK
            {
                staticComponent.plt_y_pnl.linesth.clear();
                staticComponent.plt_y_pnl.linesy.clear();
                repaint();
            }
        });
        JPanel pnl = new JPanel();
        pnl.add(btnClear);
        pnl.setMaximumSize(new Dimension(5000, 10));
        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
        getContentPane().add(pnl);
        getContentPane().add(staticComponent.plt_y_pnl);
    }
}

class Plot_y_Panel extends JPanel implements MouseListener, MouseMotionListener, KeyListener
{
    protected int cursor;
    protected double theta, x, y;
    protected ArrayList<Line2D> linesy = new ArrayList<Line2D>();
    protected ArrayList<Line2D> linesth = new ArrayList<Line2D>();
    private Point2D lasty;
    private Point2D lastth;

    public Plot_y_Panel(double xorg, double yorg)
    {
        lasty = new Point2D.Double(xorg, yorg);
    }

    public void mousePressed(MouseEvent e)
    {
        if (cursor == Cursor.CROSSHAIR_CURSOR)
        {
            lasty = new Point2D.Double(e.getX(), e.getY());
            lastth = new Point2D.Double(e.getX(), theta);
        }
        //System.out.println("pressed  " + e.getX() + ", " + e.getY() + ", " + cursor);
    }

    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mouseClicked(MouseEvent e) {}
    public void mouseMoved(MouseEvent e) {}

    public void mouseDragged(MouseEvent e)
    {
        //System.out.println("dragged  " + e.getX() + ", " + e.getY() + ", " + cursor + ", " + 180*theta/Math.PI);
        if (cursor == Cursor.CROSSHAIR_CURSOR && isShowing())
            add(e.getX(), e.getY());
    }

    public void keyPressed(KeyEvent e)
    {
        if (!isShowing()) return;
        if (lastth == null)
            lastth = new Point2D.Double(lasty.getX(), theta);
        add ((int) x, (int) y);
    }

    public void keyReleased(KeyEvent e) {}
    public void keyTyped(KeyEvent e) {}

    public void add(int xnew, int ynew)
    {
        Point2D endy = new Point2D.Double(xnew, ynew);
        Point2D endth = new Point2D.Double(xnew, theta);
        Line2D line = new Line2D.Double(lasty, endy);
        linesy.add(line);
        line = new Line2D.Double(lastth, endth);
        linesth.add(line);
        repaint();
        lasty = endy;
        lastth = endth;
    }

    @Override public void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;

        for (Line2D l : linesy)
            g2.draw(l);
        g2.setPaint(new Color(0, 0, 255));
        for (Line2D l : linesth)
            g2.draw(l);
    }
}
