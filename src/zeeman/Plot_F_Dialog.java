
package zeeman;

// plot F versus theta
// see : C:\APP\Java\CoreJava\v1_Edition7\v1ch8\Sketch

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import javax.swing.*;

public class Plot_F_Dialog extends JDialog
{
    public Plot_F_Dialog(Image img, double thetaorg, double A, double xorg, double yorg)
    {
        //System.out.println("Plot_F_Dialog " + xorg + ", " + yorg + ", " + thetaorg + ", " + height);
        setTitle(" Plot Static Zeeman Data - F vs. theta");
        setIconImage(img);
        setSize(360, 360);
        setLocationByPlatform(true);
        setVisible(true);
        staticComponent.plt_F_pnl = new Plot_F_Panel(thetaorg, A, xorg, yorg);
        add(staticComponent.plt_F_pnl);
    }
}

class Plot_F_Panel extends JPanel implements MouseListener, MouseMotionListener, KeyListener
{
    protected int cursor;
    protected double theta;
    protected double A, x, y;
    private double minF = 1000000000, maxF = 0;

    public Plot_F_Panel(double thetaorg, double Aorg, double xorg, double yorg)
    {
        double F;
        theta = thetaorg;
        A = Aorg;
        x = xorg;
        y = yorg;
        for (int i = 0; i < 360; i++)
        {
            F = main.calc_F(i*Math.PI/180, A, x, y);
            if (F > maxF) maxF = F;
            if (F < minF) minF = F;
        }
        //System.out.println(theta + ", " + A + ", " + x + ", " + y + ", " + main.calc_F(theta, A, x, y));
        //System.out.println(minF + ", " + maxF);
    }

    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mouseClicked(MouseEvent e) {}
    public void mouseMoved(MouseEvent e) {}

    public void mouseDragged(MouseEvent e)
    {
        if (cursor == Cursor.CROSSHAIR_CURSOR && isShowing())
            repaint();
    }

    public void keyPressed(KeyEvent e)
    {
        if (isShowing())
            repaint();
    }

    public void keyReleased(KeyEvent e) {}
    public void keyTyped(KeyEvent e) {}

    @Override public void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        double tempF, tempmin, tempmax;

        Path2D.Double path = new Path2D.Double(Path2D.WIND_NON_ZERO, getWidth());
        tempF = main.calc_F(0, A, x, y);
        path.moveTo(0, (getHeight() - 10)*(maxF - tempF)/(maxF - minF));
        tempmin = tempF;
        tempmax = tempF;
        for (int i = 1; i < getWidth(); i++)
        {
            tempF = main.calc_F(2*i*Math.PI/getWidth(), A, x, y);
            path.lineTo(i, (getHeight() - 10)*(maxF - tempF)/(maxF - minF));
            if (tempF > tempmax) tempmax = tempF;
            if (tempF < tempmin) tempmin = tempF;
        }
        g2.draw(path);
        g2.setPaint(new Color(0, 0, 255));
        g2.fill(new Ellipse2D.Double(theta*getWidth()/2/Math.PI - 5, (getHeight() - 10)*(maxF - main.calc_F(theta, A, x, y))/(maxF - minF) - 5, 9, 9));
        minF = tempmin;
        maxF = tempmax;
    }
}
