
package rossler;

/*
 * this is file : \Documents\NetBeansProjects\RosslerSystem\src\rossler\Scatter_Tau.java
 * see demo     : \Documents\NetBeansProjects\ButtonDemoProject\ButtonDemo.java
 * see jdk file : \APP\Java\Demos\jdk_Demos\ButtonDemo.java
 *
 * calculate the new angle in the Rossler scatter plot in the x'-y' plane, versus the old
 * scatter angle is measured once per period of the limit cycle.
 * this allows a measurement of the average increment in scatter angle
 * which can be used to calculate the modulation period of the N-S torus
*/ 

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Scatter_Tau extends JDialog
{
    private static final BufferedImage image = new BufferedImage(360, 360, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static JSlider slider_start;
    private static JSlider slider_end;
    private static JButton btnCalc = new JButton("Calc");
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static final String fDir = "\\APP\\Java\\RosslerSystem\\scatter_period\\";
    //private static final String fName = "scatter_angle_0.8493_0.6018_2.0_new";
    //private static final String fName = "scatter_angle_0.613613_0.6_1.25";
    private static final String fName = "scatter_angle_0.6154_0.6_1.25";
    //private static final String fName = "scatter_angle_0.8475_0.6_2.0";
    private static final JLabel lblfile = new JLabel("file = '" + fName + "'");
    private static double[] angles; //, times;
    private static int Nfinal = 0;
    private static double a, b, c, Period, delt, eig, angle, x0, y0;

    public Scatter_Tau()
    {
        setTitle("Rossler System - Modulation Period of x'-y' Scatter");
        setIconImage(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        setSize(530, 650);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        lblImage.setBorder(BorderFactory.createEtchedBorder());

        slider_start = new JSlider(JSlider.HORIZONTAL, 0, Nfinal, 0);
        slider_start.setMajorTickSpacing(Nfinal/10);
        slider_start.setMinorTickSpacing(1);
        slider_start.setPaintTicks(true);
        slider_start.setPaintLabels(true);
        slider_start.setBorder(BorderFactory.createEtchedBorder());

        slider_end = new JSlider(JSlider.HORIZONTAL, 0, Nfinal, Nfinal);
        slider_end.setMajorTickSpacing(Nfinal/10);
        slider_end.setMinorTickSpacing(1);
        slider_end.setPaintTicks(true);
        slider_end.setPaintLabels(true);
        slider_end.setBorder(BorderFactory.createEtchedBorder());

        final JPanel pnlfile = new JPanel();
        pnlfile.setOpaque(false);
        pnlfile.add(lblfile);
        final JPanel pnlstart = new JPanel();
        pnlstart.setOpaque(false);
        pnlstart.add(new JLabel("start"));
        final JPanel scatterPanel = new JPanel();
        scatterPanel.setOpaque(false);
        scatterPanel.add(lblImage);
        final JPanel spacerPanel = new JPanel();
        spacerPanel.setOpaque(false);
        final JPanel pnlend = new JPanel();
        pnlend.setOpaque(false);
        pnlend.add(new JLabel("end"));
        final JPanel pnlfooter = new JPanel();
        pnlfooter.setOpaque(false);
        pnlfooter.add(btnCalc);
        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
        getContentPane().setBackground(new Color(200, 221, 242));
        getContentPane().add(pnlfile);
        getContentPane().add(pnlstart);
        getContentPane().add(slider_start);
        getContentPane().add(spacerPanel);
        getContentPane().add(scatterPanel);
        getContentPane().add(pnlend);
        getContentPane().add(slider_end);
        getContentPane().add(pnlfooter);

        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        refresh_graph();

        slider_start.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                refresh_graph();
            }
        });

        slider_end.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                refresh_graph();
            }
        });
        btnCalc.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                calc_dist(slider_start.getValue(), slider_end.getValue());
            }
        });
    }

    private static void refresh_graph()
    {
        lblfile.setText("file = '" + fName + "' (" + slider_start.getValue() + ", " + slider_end.getValue() + ")");
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        DC.setColor(new Color(255, 128, 0));
        DC.drawLine(0, image.getHeight() - 1, image.getWidth() - 1, 0);
        DC.setColor(new Color(0, 64, 192));
        DC.drawLine(0, image.getHeight() - 1 - (int) angle, image.getHeight() - 1 - (int) angle, 0);
        DC.drawLine(image.getHeight() - 1 - (int) angle, image.getHeight() - 1, image.getHeight() - 1, image.getHeight() - 1 - (int) angle);
        for (int i = slider_start.getValue(); i < slider_end.getValue() - 1; i++)
        {
            //System.out.println(i + ", " + angles[i] + ", " + angles[i + 1]);
            image.setRGB((int) angles[i], image.getHeight() - 1 - (int) angles[i + 1], Color.BLACK.getRGB());
        }
        lblImage.repaint();
    }

    private static void calc_dist(int nstart, int nend)
    {
        // calculate the distribution frequency of the increment angles[i+1] - angles[i]

        int i, istart, iend;
        double[] diff = new double[nend - nstart - 1];
        int[] count = new int[360 + 90];             // number of points at each start angle (padded by 90)
        double[] average = new double[360 + 90];     // average increment at each start angle
        double av_inc = 0;
        double total_inc = 0;

        //System.out.println("org angles");
        //for (i = nstart; i <= nend; i++)
        //    System.out.println(i + ", " + angles[i]);
        for (i = 0; i < diff.length; i++)
        {
            diff[i] = (angles[nstart + i + 1] - angles[nstart + i] + 360) % 360;
            total_inc += diff[i];
        }
        //for (i = 0; i < diff.length; i++)
        //    System.out.println((nstart + i) + ", " + diff[i]);
        for (i = 0; i < 360; i++)
        {
            count[i] = 0;
            average[i] = 0;
        }
        for (i = 0; i < diff.length; i++)   // classify by integer of start angle
        {
            count[(int) angles[nstart + i]] += 1;
            average[(int) angles[nstart + i]] += diff[i];
        }
        for (i = 0; i < 360; i++)           // generate one data point at each start angle
            if (count[i] > 0)
                average[i] = average[i]/count[i];

        for (i = 0; i < 360; i++)
            if (count[i] > 0)
                System.out.println(i + ", " + count[i] + ", " + average[i]);

        //double tempsum = 0;
        //double countsum = 0;
        //for (i = 0; i < 360; i++)
        //{
        //    countsum += count[i];
        //    tempsum += count[i]*average[i];
        //}
        //System.out.println(countsum + ", " + tempsum/countsum);

        i = 0;
        while (count[i] == 0) i++;
        istart = i;                         // find nonzero start angle
        System.out.println("calc_dist, " + nstart + ", " + nend + ", " + istart + ", " + total_inc);
        if (istart > count.length - 360)
        {
            System.out.println("BAD data: istart too large : " + istart);
            return;
        }
        for (i = 0; i <= istart; i++)       // append blanks to end of data
        {
            count[360 + i] = count[i];
            average[360 + i] = average[i];
        }
        //for (i = 0; i <= 360 + istart; i++)
        //    System.out.println(i + ", " + count[i] + ", " + average[i]);

        i = istart;
        iend = istart;
        //System.out.println("\nstart loop");
        while (i < istart + 360)
        {
            while (count[iend + 1] == 0) iend++;        // find zeros
            //System.out.println();
            if (iend > i)                               // interpolate
                for (int j = i + 1; j < iend + 1; j++)
                    average[j] = average[i] + (j - i)*(average[iend + 1] - average[i])/(iend + 1 - i);
                    //System.out.println(j + ", " + (average[i] + (j - i)*(average[iend + 1] - average[i])/(iend + 1 - i)));
            iend = iend + 1;
            i = iend;
        }
        for (i = 0; i < istart; i++)                    // copy interpolated back to beginning
        {
            count[i] = count[360 + i];
            average[i] = average[360 + i];
        }
        //System.out.println("final");
        //for (i = 0; i < 360; i++)
        //    System.out.println(i + ", " + count[i] + ", " + average[i]);
        for (i = 0; i < 360; i++)
            av_inc += 1.0/average[i];
            //av_inc += average[i];
        av_inc = 360.0/av_inc;
        //av_inc = av_inc/360;
        //System.out.println(fName.substring(14) + ", " + a + ", " + b + ", " + c + ", " + Period + ", " + delt + ", " + eig + ", " + angle + ", " + av_inc);
        System.out.println(fName.substring(14) + ", " + a + ", " + b + ", " + c + ", " + Period + ", " + delt + ", " + eig + ", " + angle + ", " + total_inc/(nend - nstart - 1));
        //double diffsum = 0;
        //for (i = 0; i < diff.length; i++)
        //    diffsum += diff[i];
        //System.out.println("scatter_tau, " + a + ", " + b + ", " + c + ", " + Period + ", " + delt + ", " + eig + ", " + angle + ", " + av_inc + ", " + diffsum/diff.length);
    }

    private static void load_angles()
    {
        int i;
        int Nestimate = 0;              // preliminary over-estimate of number of lines
        String str = "";
        double x, y, r;
        System.out.println("load_angles : " + fName);
        try
        {
            BufferedReader istr = new BufferedReader(new FileReader(fDir + fName + ".csv"));
            try
            {
                while (istr.ready())                                // loop to end
                {
                    istr.readLine();
                    Nestimate++;
                }
                istr.close();
                angles = new double[Nestimate];
                //times = new double[Nestimate];
                //System.out.println("angles len = " + angles.length);

                istr = new BufferedReader(new FileReader(fDir + fName + ".csv"));
                for (i = 0; i < 5; i++)
                    str = istr.readLine();
                a = Double.parseDouble(str.split(",")[1]);
                b = Double.parseDouble(str.split(",")[2]);
                c = Double.parseDouble(str.split(",")[3]);
                str = istr.readLine();
                Period = Double.parseDouble(str.split(",")[1]);
                delt = Double.parseDouble(str.split(",")[2]);
                for (i = 0; i < 9; i++)
                    str = istr.readLine();
                eig = Double.parseDouble(str.split(",")[1]);
                angle = Double.parseDouble(str.split(",")[2]);
//                for (i = 0; i < 7; i++)
//                    str = istr.readLine();
                while (!istr.readLine().startsWith(" ,iT,x,y,z,x")) {}  // loop until after this line
                str = istr.readLine();                                  // should start with "init"
                x0 = Double.parseDouble(str.split(",")[5]);
                y0 = Double.parseDouble(str.split(",")[6]);
                System.out.println("init,        " + a + ", " + b + ", " + c + ", " + Period + ", " + delt + ", " + eig + ", " + angle + ", " + x0 + ", " + y0);
                System.out.println("scatter_tau, a, b, c, Period, delt, eig, angle, av_inc");
                while (istr.ready())                                    // should start with "z inter"
                {
                    str = istr.readLine();
                    x = Double.parseDouble(str.split(",")[5]);
                    y = Double.parseDouble(str.split(",")[6]);
                    r = Math.sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
                    //times[Nfinal] = Double.parseDouble(str.split(",")[7]);
                    angles[Nfinal] = (Math.atan2(y - y0, x - x0)*180.0/Math.PI + 360.0) % 360;
                    System.out.println(Nfinal + "," + str.split(",")[1] + "," + str.split(",")[2] + "," + str.split(",")[3] + "," + str.split(",")[4] + ", " + r + ", " + angles[Nfinal] + "," + str.split(",")[7]);
                    Nfinal++;
                }
            }
            catch (IOException e)
            {
                System.out.println("read error : " + e.getMessage());
                return;
            }
        }
        catch (FileNotFoundException e)
        {
            System.out.println("file not found : " + e.getMessage());
            return;
        }
        //for (i = 0; i < Nfinal; i++)
        //    System.out.println(i + ", " + angles[i] + ", " + times[i]);
        //calc_dist(255, 4600);
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                load_angles();
                Scatter_Tau dlg = new Scatter_Tau();
            }
        });
    }
}
