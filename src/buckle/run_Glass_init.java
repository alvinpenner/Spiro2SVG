
// solve the Glass - 'Hopf bifurcation' problem
// as per Langford 1977: "Numerical Solution of Bifurcation Problems"
// solve it as an initial-value problem
// assuming the parameters are already known as specified by Langford in Table 3

package buckle;

// this is file : \Documents\NetBeansProjects\ChuaOscillator\src\buckle\run_Glass_init.java

public class run_Glass_init
{
    private static final int N = 128;                           // number of steps in (0, 1)
    private static final double beta = Math.sqrt(3);            // eigenvalue at bifurc
    private static final double eps = 0.20;
    private static final double eta = 13.00720304;
    private static final double mu = 4 + eps*eps*eta;
    private static final double tau = -0.01311693;
    private static final double T = 2*Math.PI*(1 + eps*eps*tau)/beta;
    private static final double [] y_init = new double[] { 0, -0.14399230, -0.13934930};
    private static final double [][] A0 = new double[][] {{ -1,  0, -2},
                                                          {  2, -1,  0},
                                                          {  0,  2, -1}};
    private static final double [][] A1 = new double[][] {{  0,  0, -1},
                                                          {  1,  0,  0},
                                                          {  0,  1,  0}};
    private static double [][] phi = new double[N + 1][3];      // linear response (null vector)
    private static double [][] w = new double[N + 1][3];        // nonlinear incremental response

    public static void main (String[] args)
    {
        System.out.println("Glass_Hopf-Bifurcation System - initial-value solution");
        System.out.println("N_beta_eps, " + N + ", " + beta + ", " + eps);
        System.out.println("eta_mu,     " + eta + ", " + mu);
        System.out.println("tau_T,      " + tau + ", " + T);
        System.out.println("y_init,     " + y_init[0] + ", " + y_init[1] + ", " + y_init[2]);
        for (int i = 0; i < N + 1; i++)
        {
            phi[i][0] = Math.sqrt(2.0/3)*Math.sin(i*2*Math.PI/N);
            phi[i][1] = Math.sqrt(2.0/3)*Math.sin(i*2*Math.PI/N - Math.PI/3);
            phi[i][2] = Math.sqrt(2.0/3)*Math.sin(i*2*Math.PI/N - Math.PI*2/3);
            //System.out.println(i + ", " + phi[i][0] + ", " + phi[i][1] + ", " + phi[i][2]);
        }
        w[0] = new double[] {(y_init[0] - eps*phi[0][0])/eps/eps, (y_init[1] - eps*phi[0][1])/eps/eps, (y_init[2] - eps*phi[0][2])/eps/eps};
        //w[0] = new double[] {phi[0][0], phi[0][1], phi[0][2]};     // temporary over-ride to calc phi
        for (int i = 0; i < N; i++)
            runge_kutta_Glass_initial_value(i);
        for (int i = 0; i < N + 1; i++)
            System.out.println(i + ", " + w[i][0] + ", " + w[i][1] + ", " + w[i][2]);
    }

    private static void runge_kutta_Glass_initial_value(int iter)
    {
        // Glass model as per Langford 1977
        // solved as an initial value problem: F = F(t, x, y, z)
        // Note: the contribution of phi will be treated as a time-dependent forcing function
        double[] klm1 = F(iter      , w[iter][0]            , w[iter][1]            , w[iter][2]);
        double[] klm2 = F(iter + 0.5, w[iter][0] + klm1[0]/2, w[iter][1] + klm1[1]/2, w[iter][2] + klm1[2]/2);
        double[] klm3 = F(iter + 0.5, w[iter][0] + klm2[0]/2, w[iter][1] + klm2[1]/2, w[iter][2] + klm2[2]/2);
        double[] klm4 = F(iter + 1.0, w[iter][0] + klm3[0]  , w[iter][1] + klm3[1]  , w[iter][2] + klm3[2]);
        //double[] klm1 = transform(delt, A0, w[iter][0], w[iter][1], w[iter][2]);
        //double[] klm2 = transform(delt, A0, w[iter][0] + klm1[0]/2, w[iter][1] + klm1[1]/2, w[iter][2] + klm1[2]/2);
        //double[] klm3 = transform(delt, A0, w[iter][0] + klm2[0]/2, w[iter][1] + klm2[1]/2, w[iter][2] + klm2[2]/2);
        //double[] klm4 = transform(delt, A0, w[iter][0] + klm3[0], w[iter][1] + klm3[1], w[iter][2] + klm3[2]);
        for (int i = 0; i < 3; i++)
            w[iter + 1][i] = w[iter][i] + (klm1[i] + 2*klm2[i] + 2*klm3[i] + klm4[i])/6;
    }

    private static double [] F(double t, double x, double y, double z)      // Froberg, p. 269
    {
        double delt = 2*Math.PI/N/beta;
        double [] ret = new double[] {0, 0, 0};                             // implement Glass Eq. 5.19
        double phix = Math.sqrt(2.0/3)*Math.sin(t*2*Math.PI/N);             // phi_x(t)
        double phiy = Math.sqrt(2.0/3)*Math.sin(t*2*Math.PI/N - Math.PI/3);
        double phiz = Math.sqrt(2.0/3)*Math.sin(t*2*Math.PI/N - Math.PI*2/3);
        //System.out.println(t + ", " + x + ", " + y + ", " + z);
        transform(ret, delt, A0, x, y, z);                                  // term A0*w
        transform(ret, delt*eps*eta, A1, phix, phiy, phiz);                 // term eta*A1*phi
        transform(ret, delt*eps*tau, A0, phix, phiy, phiz);                 // term tau*A0*phi
        transform(ret, delt*eps*eps*eta, A1, x, y, z);                      // term eps*eps*eta*A1*w
        transform(ret, delt*eps*eps*tau, A0, x, y, z);                      // term eps*eps*tau*A0*w
        transform(ret, delt*eps*eps*eps*eta*tau, A1, phix, phiy, phiz);     // term eps*eps*eps*eta*tau*A1*phi
        transform(ret, delt*eps*eps*eps*eps*eta*tau, A1, x, y, z);          // term eps*eps*eps*eps*eta*tau*A1*w
        // add nonlinear term (1 + eps*eps*tau)/eps/eps*Q
        ret[0] += -delt*(1 + eps*eps*tau)*calc_Glass_G(eps*(phiz + eps*z))/eps/eps;
        ret[1] +=  delt*(1 + eps*eps*tau)*calc_Glass_G(eps*(phix + eps*x))/eps/eps;
        ret[2] +=  delt*(1 + eps*eps*tau)*calc_Glass_G(eps*(phiy + eps*y))/eps/eps;
        return ret;
    }

    private static void transform(double [] in, double scale, double [][] trans, double x, double y, double z)
    {
        for (int i = 0; i < 3; i++)
            in[i] = in[i] + scale*(trans[i][0]*x + trans[i][1]*y + trans[i][2]*z);
    }

    private static double calc_Glass_G(double y)
    {
        // nonlinear response of Glass model: Langford_1977_Numerical_Solution
        double temp = Math.pow(1.0 + 2.0*y, 4 + eps*eps*eta);
        // subtract the term linear in mu (Eq. 6.2)
        return 0.5*(temp - 1)/(temp + 1) - 2.0*y - 0.5*eps*eps*eta*y;
    }
}
