
// Rossler System limit point at a = b = c/2.
// for the tangent phase space response equations,
// calculate the time at which we have a stationary point.
// see book Chaos II, and loose sheets April 20, 2021

package rossler;

// this is file : \Documents\NetBeansProjects\RosslerSystem\src\rossler\stationary_pt.java

public class stationary_pt
{
    private static double old_phase = 3.11;
    private static double old_det;

    public static void main (String[] args)
    {
        //double temp_phase;
        //for (int i = 14; i < 15; i++)
        {
            //regula_falsi(0.1*i + 0.012);
            regula_falsi(1.414);
            //calc_det(0.1*i, temp_phase);
        }
    }

    private static double regula_falsi (double a)
    {
        double temp_phase = old_phase + 2*Math.PI/1152;
        temp_phase = old_phase + 0.01; // Math.PI - 0.001;  // override the default
        double new_phase, new_det;

        old_det = calc_det(a, old_phase);
        do
        {
            new_phase = temp_phase;
            new_det = calc_det(a, new_phase);
            temp_phase = new_phase - new_det*(new_phase - old_phase)/(new_det - old_det);
            old_phase = new_phase;
            old_det = new_det;
        } while (Math.abs(new_phase - temp_phase) > 0.00000000001);
        return temp_phase;
    }

    private static double calc_det (double a, double phase)
    {
        double delx = Math.atan(Math.sqrt(2 - a*a)/a);
        //double dely = Math.atan(-a*Math.sqrt(2 - a*a)/(1 - a*a));
        double dely = 2*delx - Math.PI;
        double xdot = -Math.sqrt(2)*Math.sin(phase + delx);
        double ydot = -Math.sin(phase + dely);
        double zdot = -Math.sin(phase);
        double phi = Math.atan2(xdot, -ydot);
        double theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));
        double M00, M01, M10, M11;
        M00 = -a*Math.sin(phi)*Math.sin(phi);
        M01 = Math.cos(phi)/Math.sin(theta);
        M10 = -2.0*a*Math.sin(phi)*Math.cos(phi)*Math.cos(theta) - Math.cos(phi)/Math.sin(theta);
        M11 = - a*Math.cos(phi)*Math.cos(phi)*Math.cos(theta)*Math.cos(theta)
              + a*Math.sin(theta)*Math.sin(theta);
        System.out.println(a + ", " + phase + ", " + M00 + ", " + M01 + ", " + M10 + ", " + M11 + ", " + phi + ", " + theta + ", " + xdot + ", " + ydot + ", " + zdot
                             + ", " + (M00*M11 - M01*M10));
        return M00*M11 - M01*M10;
    }
}
