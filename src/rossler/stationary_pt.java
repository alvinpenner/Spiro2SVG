
// Rossler System limit point at a = b = c/2.
// for the tangent phase space response equations,
// calculate the time at which we have a stationary point.
// see book Chaos II, and loose sheets April 20, 2021

package rossler;

// this is file : \Documents\NetBeansProjects\RosslerSystem\src\rossler\stationary_pt.java
// see book Chaos III, p. 1

public class stationary_pt
{
    private static final int Period = 1152;
    private static final double a = 1;
    private static final double b = 1;
    private static final double c = 2;
    private static final double z0 = (c - Math.sqrt(c*c - 4*a*b))/2/a;
    private static final double w = 1; // Math.sqrt(2 - a*a);
    private static final double absx = Math.sqrt((z0 + 1)*(z0 + 1) + w*w*(c - a - a*z0)*(c - a - a*z0));
    private static final double absy = Math.sqrt((c - a*z0)*(c - a*z0) + w*w);
    private static final double absz = Math.sqrt(z0*z0*(a*a + w*w));
    private static final double delx = Math.atan2(w*(c - a - a*z0), -z0 - 1);
    private static final double dely = Math.atan2(w, c - a*z0);
    private static final double delz = Math.atan2(w, -a);
    private static double old_phase = 3.11;
    private static double old_det;

    public static void main (String[] args)
    {
        //absx = absx/absz;
        //absy = absy/absz;
        //absz = 1;
        //delx = delx - delz;
        //dely = dely - delz;
        //delz = 0;
        System.out.println("limit pt. a b c z0 w    , " + a + ", " + b + ", " + c + ", " + z0 + ", " + w);
        System.out.println("limit pt. absx absy absz, " + absx + ", " + absy + ", " + absz);
        System.out.println("limit pt. delx dely delz, " + delx + ", " + dely + ", " + delz);
        System.out.println();
        for (int i = 0; i < Period; i++)
        {
            //regula_falsi(0.1*i + 0.012);
            //regula_falsi ();
            //calc_det(i*2*Math.PI/Period);
        }
        //gen_limit_array();
        //one_shot_calc();
        regula_falsi_ba(1.020263121, -8.520878737);     // b/a as fxn of (eig, phi)
    }

    private static double regula_falsi ()
    {
        double temp_phase = old_phase + 2*Math.PI/Period;
        temp_phase = old_phase + 0.01; // Math.PI - 0.001;  // override the default
        double new_phase, new_det;

        old_det = calc_det(old_phase);
        do
        {
            new_phase = temp_phase;
            new_det = calc_det(new_phase);
            temp_phase = new_phase - new_det*(new_phase - old_phase)/(new_det - old_det);
            old_phase = new_phase;
            old_det = new_det;
        } while (Math.abs(new_phase - temp_phase) > 0.00000000001);
        return temp_phase;
    }

    private static double calc_det (double phase)
    {
        double xdot = -absx*Math.sin(phase + delx);
        double ydot = -absy*Math.sin(phase + dely);
        double zdot = -absz*Math.sin(phase + delz);
        double phi = Math.atan2(xdot, -ydot);
        double theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));
        double M00, M01, M10, M11;
        M00 = -a*Math.sin(phi)*Math.sin(phi);
        M01 = Math.cos(phi)/Math.sin(theta);
        M10 = -2.0*a*Math.sin(phi)*Math.cos(phi)*Math.cos(theta) - Math.cos(phi)/Math.sin(theta)
            - (z0 - 1)*Math.cos(phi)*Math.sin(theta);
        M11 = - a*Math.cos(phi)*Math.cos(phi)*Math.cos(theta)*Math.cos(theta)
              - (a*z0 - c)*Math.sin(theta)*Math.sin(theta)
              + (z0 - 1)*Math.sin(phi)*Math.sin(theta)*Math.cos(theta);
        System.out.println(a + ", " + phase + ", " + M00 + ", " + M01 + ", " + M10 + ", " + M11 + ", " + phi + ", " + theta + ", " + xdot + ", " + ydot + ", " + zdot
                             + ", " + (M00*M11 - M01*M10));
        return M00*M11 - M01*M10;
    }

    private static void gen_limit_array()
    {
        // generate data for Python program 'Rossler_limit_eigsh.py'
        // see book Chaos III, p. 1
        // these formulas assume that (x0, y0, z0) is a stationary point

        double[][][] Mout = new double[2][2][Period];
        double xdot, ydot, zdot, phi, theta;

        System.out.println("# Rossler limit point d.e. elements of M (Python)");
        System.out.println("a_in   = " + a);
        System.out.println("b_in   = " + b);
        System.out.println("c_in   = " + c);
        System.out.println("delt   = " + 2.0*Math.PI/w/Period);
//        System.out.println("i, xdot, ydot, zdot, phi, theta");

        for (int i = 0; i < Period; i++)                        // i = time index
        {
            xdot = -absx*Math.sin(2*Math.PI*i/Period + delx);
            ydot = -absy*Math.sin(2*Math.PI*i/Period + dely);
            zdot = -absz*Math.sin(2*Math.PI*i/Period + delz);
            phi = Math.atan2(xdot, -ydot);
            theta = Math.acos(zdot/Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot));
//            System.out.println(i + ", " + xdot + ", " + ydot + ", " + zdot + ", " + phi + ", " + theta);
            Mout[0][0][i] = -a*Math.sin(phi)*Math.sin(phi);
            Mout[0][1][i] = Math.cos(phi)/Math.sin(theta);
            Mout[1][0][i] = -2.0*a*Math.sin(phi)*Math.cos(phi)*Math.cos(theta) - Math.cos(phi)/Math.sin(theta)
                          - (z0 - 1)*Math.cos(phi)*Math.sin(theta);
            Mout[1][1][i] = - a*Math.cos(phi)*Math.cos(phi)*Math.cos(theta)*Math.cos(theta)
                            - (a*z0 - c)*Math.sin(theta)*Math.sin(theta)
                            + (z0 - 1)*Math.sin(phi)*Math.sin(theta)*Math.cos(theta);
        }
        for (int vari = 0; vari < 2; vari++)            // var = (dx'/dc, dy'/dc, dz'/dc)
            for (int varj = 0; varj < 2; varj++)
            {
                System.out.print("M" + vari + varj + " = np.array([" + Mout[vari][varj][0]);
                for (int i = 1; i < Period; i++)            // i = time index
                    System.out.print("," + Mout[vari][varj][i]);
                System.out.println("])");
            }
    }

    private static void one_shot_calc()
    {
        // these formulas are valid for any arbitrary (x, y, z) point
        // book Chaos II, p. 44
        double x = 6.359945076332767;
        double z = 0.4333762323920686;
        double phi = 2.5208774818236215;
        double theta = 1.3719733708820308;
        double Mout00 = -a*Math.sin(phi)*Math.sin(phi);
        double Mout01 = Math.cos(phi)/Math.sin(theta);
        double Mout10 = -2.0*a*Math.sin(phi)*Math.cos(phi)*Math.cos(theta) - Math.cos(phi)/Math.sin(theta)
                        - (z - 1)*Math.cos(phi)*Math.sin(theta);
        double Mout11 = - a*Math.cos(phi)*Math.cos(phi)*Math.cos(theta)*Math.cos(theta)
                        - (x - c)*Math.sin(theta)*Math.sin(theta)
                        + (z - 1)*Math.sin(phi)*Math.sin(theta)*Math.cos(theta);
        System.out.println("\nMout :\n" + x + ", " + z + ", " + Mout00 + ", " + Mout01 + ", " + Mout10 + ", " + Mout11);
    }

    private static void regula_falsi_ba(double eig, double phi_set)
    {
        // for a cubic model of a N_S torus, calculate b/a required to produce a phase shift phi

        double old_ba = 4.92115;
        double new_ba = 4.92124;
        double old_phi = calc_phi(old_ba, eig);
        double new_phi = calc_phi(new_ba, eig);
        double temp_ba;

        System.out.println("eig, phi, b/a");
        System.out.println(eig + ", " + old_phi + ", " + old_ba);
        System.out.println(eig + ", " + new_phi + ", " + new_ba);
        do
        {
            temp_ba = old_ba + (phi_set - old_phi)*(new_ba - old_ba)/(new_phi - old_phi);
            old_ba = new_ba;
            old_phi = new_phi;
            new_ba = temp_ba;
            new_phi = calc_phi(new_ba, eig);
            System.out.println(eig + ", " + new_phi + ", " + new_ba);
        } while (Math.abs(new_ba - old_ba) > 0.00000000001);
    }

    private static double calc_phi (double ba, double eig)
    {
        // see book Chaos IV - page 16
        // for a cubic model of a Neimark-Sacker torus
        // calculate phase shift phi as fxn of eig = (1 + alpha) and b/a
        double sqr = Math.sqrt(1 + ba*ba - ba*ba*eig*eig);
        return Math.atan(ba*(-eig + sqr)/(ba*ba*eig + sqr))*180/Math.PI;
    }
}
