
package rossler;

// slider for Euler angles phi, theta

import java.awt.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Euler_Slider extends JDialog
{
//    protected static JCheckBox reverseChk = new JCheckBox("reverse time");

    public Euler_Slider(Image img)
    {
        setTitle(" Rossler System - Euler Slider");
        setIconImage(img);
        setSize(650, 200);
        setLocationByPlatform(true);

        final JSlider slider_phi = new JSlider(JSlider.HORIZONTAL, 0, 180, (int) Main.project_phi);
        slider_phi.setBackground(new Color(200, 221, 242));
        slider_phi.setMajorTickSpacing(10);
        slider_phi.setMinorTickSpacing(1);
        slider_phi.setPaintTicks(true);
        slider_phi.setPaintLabels(true);
        slider_phi.setBorder(BorderFactory.createEtchedBorder());

        final JSlider slider_theta = new JSlider(JSlider.HORIZONTAL, -90, 90, (int) Main.project_theta);
        slider_theta.setBackground(new Color(200, 221, 242));
        slider_theta.setMajorTickSpacing(10);
        slider_theta.setMinorTickSpacing(1);
        slider_theta.setPaintTicks(true);
        slider_theta.setPaintLabels(true);
        slider_theta.setBorder(BorderFactory.createEtchedBorder());

        final JPanel spacerPanel1 = new JPanel();
        spacerPanel1.setBackground(new Color(200, 221, 242));
        spacerPanel1.add(new JLabel("phi"));
        final JPanel spacerPanel2 = new JPanel();
        spacerPanel2.setBackground(new Color(200, 221, 242));
        spacerPanel2.add(new JLabel("theta"));
        final JPanel spacerPanel3 = new JPanel();
        spacerPanel3.setBackground(new Color(200, 221, 242));

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
        getContentPane().add(spacerPanel1);
        getContentPane().add(slider_phi);
        getContentPane().add(spacerPanel2);
        getContentPane().add(slider_theta);
        getContentPane().add(spacerPanel3);
        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        slider_phi.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                Main.project_phi = slider_phi.getValue();
            }
        });
        slider_theta.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                Main.project_theta = slider_theta.getValue();
            }
        });
    }
}
