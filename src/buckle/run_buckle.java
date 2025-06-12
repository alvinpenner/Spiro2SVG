
// solve the Euler rod-buckling problem
// as per Langford 1977: "Numerical Solution of Bifurcation Problems"
// see: \APP\Java\ChuaOscillator\Langford_1977_Numerical\Buckled_Rod_numerical.py - May 7, 2025
// use Runge-Kutta to solve a two-point boundary value problem
// with constraints on slope at endpoints
// try Bjorck p. 353, Eq. 8.3.20-8.3.21
// see also Stormer Method: Stoer+Bulirsch, p. 539

package buckle;

// this is file : \Documents\NetBeansProjects\ChuaOscillator\src\buckle\run_buckle.java

//import java.io.*;

public class run_buckle
{
    private static final int N = 64;                    // number of steps in (0, 1)
    private static double eps = 1.0;
    private static double lamb  = Math.PI*Math.PI;      // eigenvalue lambda
    private static double [] w = new double[N + 1];     // incremental response
    private static double [] wx = new double[N + 1];    // wx = dw/dx
    private static double [] theta = new double[N + 1]; // new eigenvector
    private static double [] phi = new double[N + 1];   // linear response
    private static double [] f = new double[N + 1];     // forcing function (rhs)

    public static void main (String[] args)
    {
        for (int i = 0; i < N + 1; i++)
        {
            w[i] = 0;
            wx[i] = 0;
            theta[i] = 0;
            phi[i] = Math.sqrt(2)*Math.cos(i*Math.PI/N);
            //System.out.println(i + ", " + w[i] + ", " + theta[i] + ", " + phi[i]);
        }
        System.out.println("N_eps_lambda, " + N + ", " + eps + ", " + lamb);
        //for (int i = 0; i < N + 1; i++)
        //    theta[i] = Math.sin(i*Math.PI/N)*Math.PI/2;
        //System.out.println("\n" + Cotes_4(theta));
        System.out.println("iter, w, wx, lambda");
        for (int i = 0; i < 10; i++)
            solve();
    }

    private static void solve()
    {
        // for 'rhs interpolate', see \ChuaOscillator\Langford_1977_Numerical\interpolate_cubic.xls
        double [] temp = new double[N + 1];
        for (int i = 0; i < N + 1; i++)
            theta[i] = eps*phi[i] + eps*eps*w[i];
        for (int i = 0; i < N + 1; i++)
            temp[i] = phi[i]*Math.sin(theta[i]);        // temporary vector
        lamb = eps*Math.PI*Math.PI/Cotes_4(temp);       // Simpson(temp)
        //System.out.println("\nw_in " + w_in + ", lambda = " + lamb);
        for (int i = 0; i < N + 1; i++)
            f[i] = (Math.PI*Math.PI*theta[i] - lamb*Math.sin(theta[i]))/eps/eps;    // forcing function
        //System.out.println("solve lambda, " + lamb);

        double w_old = w[0], w_new = w_old + 0.01, wx_old, wx_new, w_temp;
        //System.out.println("i, w, wx");
        w_old = -0.015929619;                                  // temporary override
        wx_old = shoot(w_old); //w_old);
        //System.out.println("0    " + ", " + w_old + ", " + wx_old);
        //wx_old = shoot(w_old + 0.01);
        //System.out.println("test0" + ", " + (w_old + 0.01) + ", " + wx_old);
        //for (int i = 0; i < 0; i++)
        //    if (Math.abs(wx_old) > 0.01)
        //{
        //    wx_new = shoot(w_new);
        //    System.out.println((i + 1) + ", " + w_new + ", " + wx_new + ", " + (wx_new - wx_old)/(w_new - w_old));
        //    w_temp = (w_old*wx_new - w_new*wx_old)/(wx_new - wx_old);
        //    w_old = w_new;
        //    wx_old = wx_new;
        //    w_new = w_temp;
        //}

        for (int i = 0; i < N + 1; i++)
            temp[i] = phi[i]*w[i];                     // integrate phi[]*w[]
        System.out.println("final, " + eps + ", " + w_old + ", " + wx_old + ", " + theta[0] + ", " + theta[16] + ", " + lamb + ", " + Cotes_4(temp));
    }

    private static double shoot(double w_in)
    {
        double [] rhs;

        w[0] = w_in;
        wx[0] = 0;
        //System.out.println("i, phi, theta, f, w, wx, phi*w");
        //System.out.println("0, " + phi[0] + ", " + theta[0] + ", " + f[0] + ", " + w[0] + ", " + wx[0] + ", " + phi[0]*w[0]);
        for (int i = 0; i < N; i++)         // calculate w[i+1], wx[i+1]
        {
            if (i == 0)
                rhs = new double[] {f[0], (3*f[0] + 6*f[1] - f[2])/8, f[1]};
            else if (i == N - 1)
                rhs = new double[] {f[N - 1], (-f[N - 2] + 6*f[N - 1] + 3*f[N])/8, f[N]};
            else
                rhs = new double[] {f[i], (-f[i - 1] + 9*f[i] + 9*f[i + 1] - f[i + 2])/16, f[i + 1]};
            runge_kutta_buckle2(i, rhs);
            //System.out.println((i + 1) + ", " + phi[i + 1] + ", " + theta[i + 1] + ", " + f[i + 1] + ", " + w[i + 1] + ", " + wx[i + 1] + ", " + phi[i + 1]*w[i + 1]);
        }
        return wx[N];
    }

    private static void runge_kutta_buckle2(int iter, double[] rhs)
    {
        // assume we have interpolated RHS = f(t)
        // f = f[t, t + delt/2, t + delt]
        double delt = 1.0/N;
        double x = w[iter];
        double y = wx[iter];
        double k1, k2, k3, k4;
        double l1, l2, l3, l4;

        k1 = delt*y;
        l1 = delt*(-Math.PI*Math.PI*x + rhs[0]);

        k2 = delt*(y + l1/2);
        l2 = delt*(-Math.PI*Math.PI*(x + k1/2) + rhs[1]);

        k3 = delt*(y + l2/2);
        l3 = delt*(-Math.PI*Math.PI*(x + k2/2) + rhs[1]);

        k4 = delt*(y + l3);
        l4 = delt*(-Math.PI*Math.PI*(x + k3) + rhs[2]);

        w[iter + 1] = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        wx[iter + 1] = y + (l1 + 2*l2 + 2*l3 + l4)/6;
    }

    private static double Simpson(double[] v)
    {
        // integrate vector v (assume N is even)
        double ret = 0;
        //System.out.println("v " + v.length);
        for (int i = 1; i < v.length; i += 2)
            ret += v[i - 1] + 4*v[i] + v[i + 1];
        return ret/(v.length - 1)/3;
    }

    private static double Cotes_4(double[] v)
    {
        // integrate vector v (assume N is multiple of 4)
        // see Froberg p.201, Table of Cote's numbers
        double ret = 0;
        //System.out.println("v " + v.length);
        for (int i = 2; i < v.length; i += 4)
            ret += 7*v[i - 2] + 32*v[i - 1] + 12*v[i] + 32*v[i + 1] + 7*v[i + 2];
        return ret*2/(v.length - 1)/45;
    }
}
