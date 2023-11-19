
package chua;

// slider for Euler angles phi, theta

import java.awt.*;
//import java.awt.geom.Point2D;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Chua_Euler_Slider extends JDialog
{
//    protected static JCheckBox reverseChk = new JCheckBox("reverse time");

    public Chua_Euler_Slider(Image img)
    {
        setTitle(" Chua System - Euler Slider");
        setIconImage(img);
        setSize(650, 300);
        setLocationByPlatform(true);

        final JSlider slider_phi = new JSlider(JSlider.HORIZONTAL, -180, 180, (int) Main.project_phi);
        slider_phi.setMajorTickSpacing(20);
        slider_phi.setMinorTickSpacing(2);
        slider_phi.setPaintTicks(true);
        slider_phi.setPaintLabels(true);
        slider_phi.setBorder(BorderFactory.createEtchedBorder());

        final JSlider slider_theta = new JSlider(JSlider.HORIZONTAL, 0, 180, (int) Main.project_theta);
        slider_theta.setMajorTickSpacing(10);
        slider_theta.setMinorTickSpacing(1);
        slider_theta.setPaintTicks(true);
        slider_theta.setPaintLabels(true);
        slider_theta.setBorder(BorderFactory.createEtchedBorder());

        final JSlider slider_psi = new JSlider(JSlider.HORIZONTAL, 0, 360, (int) Main.project_psi);
        slider_psi.setMajorTickSpacing(15);
        slider_psi.setMinorTickSpacing(5);
        slider_psi.setPaintTicks(true);
        slider_psi.setPaintLabels(true);
        slider_psi.setBorder(BorderFactory.createEtchedBorder());

        final JPanel spacerPanel1 = new JPanel();
        spacerPanel1.setOpaque(false);
        final JLabel lblphi = new JLabel("phi = " + slider_phi.getValue());
        final JPanel spacerPanel2 = new JPanel();
        spacerPanel2.setOpaque(false);
        final JLabel lbltheta = new JLabel("theta = " + slider_theta.getValue());
        final JPanel spacerPanel3 = new JPanel();
        spacerPanel3.setOpaque(false);
        final JLabel lblpsi = new JLabel("psi = " + slider_psi.getValue());
        final JPanel spacerPanel4 = new JPanel();
        spacerPanel4.setOpaque(false);

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
        getContentPane().setBackground(new Color(200, 221, 242));
        getContentPane().add(spacerPanel1);
        getContentPane().add(lblphi);
        getContentPane().add(slider_phi);
        getContentPane().add(spacerPanel2);
        getContentPane().add(lbltheta);
        getContentPane().add(slider_theta);
        getContentPane().add(spacerPanel3);
        getContentPane().add(lblpsi);
        getContentPane().add(slider_psi);
        getContentPane().add(spacerPanel4);
        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        slider_phi.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                Main.project_phi = slider_phi.getValue();
                lblphi.setText("phi = " + slider_phi.getValue());
            }
        });
        slider_theta.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                Main.project_theta = slider_theta.getValue();
                lbltheta.setText("theta = " + slider_theta.getValue());
            }
        });
        slider_psi.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                Main.project_psi = slider_psi.getValue();
                lblpsi.setText("psi = " + slider_psi.getValue());
            }
        });
        //slider_psi.setSnapToTicks(true);    // snap psi to nearest 5 degrees
    }
}
