
package chua;

/*
 * this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Scatter_Tau.java
 * see jdk file : \APP\Java\Demos\jdk_Demos\ButtonDemo.java
 *
 * calculate the new angle in the Chua scatter plot in the x'-y' plane, versus the old
 * scatter angle is measured once per period of the limit cycle.
 * this allows a measurement of the average increment in scatter angle
 * which can be used to calculate the modulation period of the N-S torus
*/ 

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Chua_Scatter_Tau extends JDialog
{
    private static final BufferedImage image = new BufferedImage(360, 360, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static JSlider slider_start;
    private static JSlider slider_end;
    private static JButton btnCalc = new JButton("Calc");
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    //private static final JLabel lblPosn = new JLabel("posn");
    //private static final String fDir = "\\APP\\Java\\ChuaOscillator\\scatter_period\\";
    private static final String fDir = "\\APP\\Java\\ChuaOscillator\\Simul_4\\";
    //private static final String fName = "scatter_angle_99.98_530999600";
    //private static final String fDir = "\\Windows\\Temp\\";
    //private static final String fName = "Chua_scatter_530998000_0.021_0.0_99.98";
    //private static final String fName = "Chua_scatter_56640000_0.0010_0.0_99.994822_3";
    //private static final String fName = "Chua_scatter_384960000_0.0010_0.0_99.9948_800_1200";
    //private static final String fName = "Chua_scatter_180000000_0.0010_0.0_99.99481_100_1200";
    //private static final String fName = "Chua_scatter_91180000_0.0010_0.0_99.9948_0_500";
    //private static final String fName = "Chua_scatter_347999420_0.0010_0.0_99.98_0_600";
    //private static final String fName = "Chua_scatter_237600000_5.0E-4_0.0_99.9944_100_1200";
    //private static final String fName = "Chua_scatter_121440000_0.0010_0.0_99.994_100_1400";
    //private static final String fName = "Chua_scatter_320640000_0.0010_0.0_99.992_100_800";
    //private static final String fName = "Chua_scatter_530880000_0.0010_0.0_99.99_100_400";
    //private static final String fName = "Chua_scatter_757920000_0.0010_0.0_99.988_100_400";
    //private static final String fName = "Chua_scatter_997920000_0.0010_0.0_99.986_100_400";
    //private static final String fName = "Chua_scatter_1396320000_0.0010_0.0_99.983_100_400";
    //private static final String fName = "Chua_scatter_1958880000_0.0010_0.0_99.98_100_600";
    //private static final String fName = "Chua_scatter_478080000_0.0010_0.0_99.978_100_800";
    //private static final String fName = "Chua_scatter_478080000_5.0E-4_0.0_99.978_0_1200";
    //private static final String fName = "Chua_scatter_347998212_0.0050_0.0_99.98_0_800";
    //private static final String fName = "Chua_scatter_69598071_0.0020_0.0_99.9948_0_1100";
    //private static final String fName = "Chua_scatter_347997600_1.0E-4_0.0_99.98_0_1200";
    //private static final String fName = "Chua_scatter_347999425_0.0010_0.0_99.98_0_800";
    //private static final String fName = "Chua_scatter_1141920000_0.0010_0.0_99.977_100_1200";
    //private static final String fName = "Chua_scatter_1609920000_0.0010_0.0_99.977_100_1200";
    //private static final String fName = "Chua_Simul_5_0.00004_10.75_-5_5_0_3.5E10";
    //private static final String fName = "Chua_Simul_5_0.00004_10.75_-5_5_100000_-1500000000";
    //private static final String fName = "Chua_Simul_5_0.00004_10.75_-5_5_g0_10E8_0";
    //private static final String fName = "Chua_Simul_5_0.00004_10.75_-5_5_20.186E4_0_g";
    //private static final String fName = "Chua_Simul_5_0.00004_10.75_-5_5_0_0_C";
    //private static final String fName = "Chua_Simul_4_0.00004_10.75_-5_5_0_4.386E6_2";
    //private static final String fName = "Chua_Simul_4g_0.00004_10.75_-5_5_0_3.5E6";
    //private static final String fName = "Chua_Simul_4g_0.00004_10.75_-5_5_-7E11_0";
    //private static final String fName = "Chua_Simul_4g_0.00004_10.75_-5_5_0_4.465E6";
    //private static final String fName = "Chua_Simul_4_0.00008_10.762308_-5_-13.425_04y";
    //private static final String fName = "Chua_scatter_1252790000_0.0010_0.0_99.98";
    //private static final String fName = "Chua_Simul_4_0.0002_10.75_-5_-247";
    //private static final String fName = "Chua_Simul_4C_0.00004_46_-5_-5_0_1.67E7";
    private static final String fName = "Chua_Simul_4C_0.00004_30_-5_5_0_1.4160919E7";

    private static final JLabel lblfile = new JLabel("file = '" + fName + "'");
    private static final JLabel lblangle = new JLabel(" : angle = ");
    private static double[] angles; //, times;
    private static int Nfinal = 0;
    private static double alpha, beta, gamma, a, c, delt, Nhdr, eig, angle, x0, y0, xstat, ystat;
    private static double ymin = 0;
    private static double ymax = 60;

    public Chua_Scatter_Tau()
        {
        setTitle("Chua System - Modulation Period of x'-y' Scatter (yrange = " + ymin + "-" + ymax + ")");
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
        pnlend.add(lblangle);
        final JPanel pnlfooter = new JPanel();
        pnlfooter.setOpaque(false);
        pnlfooter.add(btnCalc);
        //pnlfooter.add(lblPosn);
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
                slider_start.setValue(slider_end.getValue() - 5000);
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
        lblImage.addMouseMotionListener(new MouseMotionAdapter()
        {
            @Override public void mouseMoved(MouseEvent e)
            {
                lblangle.setText(" :  angle = " + e.getX());
            }
        });
    }

    private static void refresh_graph()
    {
        double delta;
        lblfile.setText("file = '" + fName + "' (" + slider_start.getValue() + ", " + slider_end.getValue() + ")");
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        //DC.setColor(new Color(255, 128, 0));
        //DC.drawLine(0, image.getHeight() - 1, image.getWidth() - 1, 0);
        DC.setColor(Color.gray);
        DC.drawLine(0, image.getHeight() - 1 - (int) (-ymin*360/(ymax - ymin)), image.getWidth() - 1, image.getHeight() - 1 - (int) (-ymin*360/(ymax - ymin)));
        for (int i = 1; i < 12; i++)
            DC.drawLine(30*i, image.getHeight() - 1, 30*i, image.getHeight() - 8);
        DC.setColor(new Color(0, 64, 192));
        DC.drawLine(0, image.getHeight() - 1 - (int) ((angle - ymin)*360/(ymax - ymin)), image.getWidth() - 1, image.getHeight() - 1 - (int) ((angle - ymin)*360/(ymax - ymin)));

        for (int i = slider_start.getValue(); i < slider_end.getValue() - 1; i++)
        {
            delta = (angles[i + 1] - angles[i] + 540) % 360 - 180.0;
            if (delta < ymax)
                image.setRGB((int) angles[i], image.getHeight() - 1 - (int) ((delta - ymin)*360/(ymax - ymin)), Color.BLACK.getRGB());
            else
                System.out.println("BAD DATA, " + i + ", " + angles[i] + ", " + angles[i + 1] + ", " + delta);
        }
        lblImage.repaint();
    }

    private static void calc_dist(int nstart, int nend)
    {
        // calculate the distribution frequency of the increment angles[i+1] - angles[i]

        int i;          // , istart, iend;
        double[] diff = new double[nend - nstart - 1];
        int[] count = new int[360 + 90];             // number of points at each start angle (padded by 90)
        double[] average = new double[360 + 90];     // average increment at each start angle
        double total_inc = 0;
        double nwrap = 0, totalwrap = 0;

        //System.out.println("org angles, diff");
        //for (i = nstart; i <= nend; i++)
        //    System.out.println(i + ", " + angles[i]);
        for (i = 0; i < diff.length; i++)
        {
            //diff[i] = (angles[nstart + i + 1] - angles[nstart + i] + 360) % 360;
            diff[i] = (angles[nstart + i + 1] - angles[nstart + i] + 540) % 360 - 180.0;
            //System.out.println(i + ", " + angles[nstart + i] + ", " + diff[i]);
            if ((total_inc + diff[i]) % 360 < total_inc % 360)
            {
                nwrap = i;
                totalwrap = total_inc;
                //System.out.println(",,," + nwrap + ", " + totalwrap + ", " + totalwrap/nwrap);
            }
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
            if (count[i] > 0 && true)
                System.out.println(i + ", " + count[i] + ", " + average[i]);

        System.out.println("calc_dist, " + fName + ", " + (int) Nhdr + ", " + nstart + ", " + nend + ", " + alpha + ", " + beta + ", " + gamma + ", " + a + ", " + c + ", " + delt + ", " + eig + ", " + angle + ", " + total_inc/(nend - nstart - 1) + ", " + totalwrap/nwrap);
/*
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
*/
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
                System.out.println("angles len = " + angles.length);

                istr = new BufferedReader(new FileReader(fDir + fName + ".csv"));
                for (i = 0; i < 2; i++)
                    str = istr.readLine();
                alpha = Double.parseDouble(str.split(",")[1]);
                beta  = Double.parseDouble(str.split(",")[2]);
                gamma = Double.parseDouble(str.split(",")[3]);
                a     = Double.parseDouble(str.split(",")[4]);
                c     = Double.parseDouble(str.split(",")[5]);
                delt  = Double.parseDouble(str.split(",")[6]);
                str = istr.readLine();
                Nhdr  = Double.parseDouble(str.split(",")[1]);
                eig   = Double.parseDouble(str.split(",")[2]);
                angle = Double.parseDouble(str.split(",")[3]);
                //while (!istr.readLine().startsWith(" ,iT,x,y,z,x")) {}  // loop until after this line
                str = istr.readLine();                                    // should start with "init x0 y0"
                x0 = Double.parseDouble(str.split(",")[1]);
                y0 = Double.parseDouble(str.split(",")[2]);
                xstat = Double.parseDouble(str.split(",")[3]);
                ystat = Double.parseDouble(str.split(",")[4]);
                System.out.println("stat1 = " + (Math.atan2(-ystat - y0, -xstat - x0)*180.0/Math.PI + 360.0) % 360);
                System.out.println("stat2 = " + (Math.atan2(-y0, -x0)*180.0/Math.PI + 360.0) % 360);
                System.out.println("stat3 = " + (Math.atan2(ystat - y0, xstat - x0)*180.0/Math.PI + 360.0) % 360);
                str = istr.readLine();                                  // dummy
                System.out.println("scatter_tau, " + alpha + ", " + beta + ", " + gamma + ", " + a + ", " + c + ", " + delt + ", " + (int) Nhdr + ", " + eig + ", " + angle + ", " + x0 + ", " + y0 + ", " + xstat + ", " + ystat);
                System.out.println("Nfinal, x, y, r, angle");
                while (istr.ready())                                    // should start with "z inter"
                {
                    str = istr.readLine();
                    x = Double.parseDouble(str.split(",")[1]);
                    y = Double.parseDouble(str.split(",")[2]);
                    r = Math.sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
                    angles[Nfinal] = (Math.atan2(y - y0, x - x0)*180.0/Math.PI + 360.0) % 360;
                    //angles[Nfinal] = 360 - angles[Nfinal];                     // fix fix bug bug reverse direction temporarily
                    System.out.println(Nfinal + ", " + x + ", " + y + ", " + r + ", " + angles[Nfinal]);
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
                Chua_Scatter_Tau dlg = new Chua_Scatter_Tau();
            }
        });
    }
}
