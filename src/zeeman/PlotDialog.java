
package zeeman;

// see : C:\APP\Java\CoreJava\v1_Edition7\v1ch8\Sketch

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
//import java.awt.image.BufferedImage;
import javax.swing.*;

public class PlotDialog extends JDialog
{
    //private static final BufferedImage img = new BufferedImage(300, 300, BufferedImage.TYPE_3BYTE_BGR);
    //private static final Graphics2D DC = img.createGraphics();

    public PlotDialog(int width, int height, Image img, double xorg, double yorg)
    {
        //System.out.println("PlotDialog " + xorg + ", " + yorg + ", " + thetaorg + ", " + height);
        //addWindowListener(new WindowAdapter() {
        //    @Override public void windowClosing(WindowEvent ev) {
        //        System.out.println("PlotDialog closing");
        //    }
        //});
        setTitle(" Plot Static Zeeman Data");
        setIconImage(img);
        setSize(width, height);
        setLocationByPlatform(true);
        setVisible(true);
        staticComponent.pltpnl = new PlotPanel(xorg, yorg);
        add(staticComponent.pltpnl);
    }
}

class PlotPanel extends JPanel implements MouseListener, MouseMotionListener, KeyListener
{
    protected int cursor;
    protected double theta;
    private ArrayList<Line2D> linesy = new ArrayList<Line2D>();
    private ArrayList<Line2D> linesth = new ArrayList<Line2D>();
    private Point2D lasty;
    private Point2D lastth;

    public PlotPanel(double xorg, double yorg)
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

    public void mouseReleased(MouseEvent e)
    {
        //System.out.println("released " + e.getX() + ", " + e.getY() + ", " + cursor);
    }

    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mouseClicked(MouseEvent e) {}

    public void mouseMoved(MouseEvent e) {}

    public void mouseDragged(MouseEvent e)
    {
        //System.out.println("dragged  " + e.getX() + ", " + e.getY() + ", " + cursor + ", " + 180*theta/Math.PI);
        if (cursor == Cursor.CROSSHAIR_CURSOR)
            add(e.getX(), e.getY());
    }

    public void keyPressed(KeyEvent e)
    {
        if (lastth == null)
            lastth = new Point2D.Double(lasty.getX(), theta);
        //System.out.println("PlotDialog keyPressed = " + e.getKeyCode() + ", " + lasty.getX());
        switch (e.getKeyCode())
        {
        case KeyEvent.VK_LEFT:
            add((int) lasty.getX() - 1, (int) lasty.getY());
        break;
        case KeyEvent.VK_RIGHT:
            add((int) lasty.getX() + 1, (int) lasty.getY());
        break;
        case KeyEvent.VK_UP:
            add((int) lasty.getX(), (int) lasty.getY() - 1);
        break;
        case KeyEvent.VK_DOWN:
            add((int) lasty.getX(), (int) lasty.getY() + 1);
        break;
        }
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
