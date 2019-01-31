
package components;

// fit a quartic polynomial in 2D to a set of 3 critical points.
// Input data consists of : c, d1, d2, rms*rms*180.0*180.0/2.0, eig0, eigangle.
// which is measured at 3 critical points, of which one is a saddle point.
// Input data is from a BezierCubic ODF fit.
// At each point we impose 5 constraints: F, dF/dx, dF/dy, eig0, eigangle.
// the constraints F, dF/dx = 0, dF/dy = 0 are self-explanatory;
// the constraints on eigenvectors are imposed using the linear equations:
// M00*cos(theta) - M01*sin(theta) =  eig*cos(theta)
// M01*cos(theta) - M11*sin(theta) = -eig*sin(theta)
// where M is the second-order response matrix.
// we are ignoring the larger eigenvalue since we only wish to map the path of the minimum F
// joining the saddle point to the two bracketing local minima.
// The curve has 15 d.f.: 1, x, y, x^2, x y, y^2, x^3, x^2 y, x y^2, y^3, x^4, x^3 y, x^2 y^2, x y^3, y^4
// The matrix equation is: M*soln = v

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\fitquartic.java

public class fitquartic
{
    private static final double[] bot = new double[] {3.596, 57.37560224314887, 30.967774777857894, 9.551748694876165E-6, 1.8540526058874462E-10, 1.0688249862172572};
    private static final double[] sad = new double[] {3.596, 57.385613334213964, 30.949532747053127, 9.55174870814764E-6, -1.8218471002640069E-10, 1.068936197685124};
    private static final double[] top = new double[] {3.596, 57.934460259245746, 29.942019686141055, 9.550589977320092E-6, 1.048041579736203E-8, 1.0750307025263028};
    private static double[][] M;
    private static double[] v = new double[] {bot[3], 0, 0, bot[4]*Math.cos(bot[5]), -bot[4]*Math.sin(bot[5]),
                                              sad[3], 0, 0, sad[4]*Math.cos(sad[5]), -sad[4]*Math.sin(sad[5]),
                                              top[3], 0, 0, top[4]*Math.cos(top[5]), -top[4]*Math.sin(top[5])};
    private static double[] soln;

    public static void main (String[] args)
    {
        System.out.println("fitquartic: c = , " + bot[0]);
        System.out.println("         , d1, d2, F, eig, eigangle");
        System.out.println("bottom = , " + bot[1] + ", " + bot[2] + ", " + bot[3] + ", " + bot[4] + ", " + bot[5]);
        System.out.println("saddle = , " + sad[1] + ", " + sad[2] + ", " + sad[3] + ", " + sad[4] + ", " + sad[5]);
        System.out.println("top    = , " + top[1] + ", " + top[2] + ", " + top[3] + ", " + top[4] + ", " + top[5]);
        solve_quartic();
    }

    private static void solve_quartic()
    {
        System.out.println("\nsolving at c = " + bot[0]);
        M = new double[][] {
            {1, bot[1], bot[2], bot[1]*bot[1], bot[1]*bot[2], bot[2]*bot[2], bot[1]*bot[1]*bot[1], bot[1]*bot[1]*bot[2], bot[1]*bot[2]*bot[2], bot[2]*bot[2]*bot[2], bot[1]*bot[1]*bot[1]*bot[1], bot[1]*bot[1]*bot[1]*bot[2], bot[1]*bot[1]*bot[2]*bot[2], bot[1]*bot[2]*bot[2]*bot[2], bot[2]*bot[2]*bot[2]*bot[2]},
            {0,      1,      0,      2*bot[1],        bot[2],             0,      3*bot[1]*bot[1],      2*bot[1]*bot[2],        bot[2]*bot[2],                    0,      4*bot[1]*bot[1]*bot[1],      3*bot[1]*bot[1]*bot[2],      2*bot[1]*bot[2]*bot[2],        bot[2]*bot[2]*bot[2],                           0},
            {0,      0,      1,             0,        bot[1],      2*bot[2],                    0,        bot[1]*bot[1],      2*bot[1]*bot[2],      3*bot[2]*bot[2],                           0,        bot[1]*bot[1]*bot[1],      2*bot[1]*bot[1]*bot[2],      3*bot[1]*bot[2]*bot[2],      4*bot[1]*bot[1]*bot[1]},
            {0,      0,      0, 2*Math.cos(bot[5]), -Math.sin(bot[5]),                   0, 6*bot[1]*Math.cos(bot[5]), 2*bot[2]*Math.cos(bot[5]) - 2*bot[1]*Math.sin(bot[5]),                             -2*bot[2]*Math.sin(bot[5]),                          0, 12*bot[1]*bot[1]*Math.cos(bot[5]), 6*bot[1]*bot[2]*Math.cos(bot[5]) - 3*bot[1]*bot[1]*Math.sin(bot[5]), 2*bot[2]*bot[2]*Math.cos(bot[5]) - 4*bot[1]*bot[2]*Math.sin(bot[5]),                                    -3*bot[2]*bot[2]*Math.sin(bot[5]),                                  0},
            {0,      0,      0,                  0,  Math.cos(bot[5]), -2*Math.sin(bot[5]),                         0, 2*bot[1]*Math.cos(bot[5])                            ,  2*bot[2]*Math.cos(bot[5]) - 2*bot[1]*Math.sin(bot[5]), -6*bot[2]*Math.sin(bot[5]),                                 0, 3*bot[1]*bot[1]*Math.cos(bot[5])                                   , 4*bot[1]*bot[2]*Math.cos(bot[5]) - 2*bot[1]*bot[1]*Math.sin(bot[5]),  3*bot[2]*bot[2]*Math.cos(bot[5]) - 6*bot[1]*bot[2]*Math.sin(bot[5]), -12*bot[2]*bot[2]*Math.sin(bot[5])},
            {1, sad[1], sad[2], sad[1]*sad[1], sad[1]*sad[2], sad[2]*sad[2], sad[1]*sad[1]*sad[1], sad[1]*sad[1]*sad[2], sad[1]*sad[2]*sad[2], sad[2]*sad[2]*sad[2], sad[1]*sad[1]*sad[1]*sad[1], sad[1]*sad[1]*sad[1]*sad[2], sad[1]*sad[1]*sad[2]*sad[2], sad[1]*sad[2]*sad[2]*sad[2], sad[2]*sad[2]*sad[2]*sad[2]},
            {0,      1,      0,      2*sad[1],        sad[2],             0,      3*sad[1]*sad[1],      2*sad[1]*sad[2],        sad[2]*sad[2],                    0,      4*sad[1]*sad[1]*sad[1],      3*sad[1]*sad[1]*sad[2],      2*sad[1]*sad[2]*sad[2],        sad[2]*sad[2]*sad[2],                           0},
            {0,      0,      1,             0,        sad[1],      2*sad[2],                    0,        sad[1]*sad[1],      2*sad[1]*sad[2],      3*sad[2]*sad[2],                           0,        sad[1]*sad[1]*sad[1],      2*sad[1]*sad[1]*sad[2],      3*sad[1]*sad[2]*sad[2],      4*sad[1]*sad[1]*sad[1]},
            {0,      0,      0, 2*Math.cos(sad[5]), -Math.sin(sad[5]),                   0, 6*sad[1]*Math.cos(sad[5]), 2*sad[2]*Math.cos(sad[5]) - 2*sad[1]*Math.sin(sad[5]),                             -2*sad[2]*Math.sin(sad[5]),                          0, 12*sad[1]*sad[1]*Math.cos(sad[5]), 6*sad[1]*sad[2]*Math.cos(sad[5]) - 3*sad[1]*sad[1]*Math.sin(sad[5]), 2*sad[2]*sad[2]*Math.cos(sad[5]) - 4*sad[1]*sad[2]*Math.sin(sad[5]),                                    -3*sad[2]*sad[2]*Math.sin(sad[5]),                                  0},
            {0,      0,      0,                  0,  Math.cos(sad[5]), -2*Math.sin(sad[5]),                         0, 2*sad[1]*Math.cos(sad[5])                            ,  2*sad[2]*Math.cos(sad[5]) - 2*sad[1]*Math.sin(sad[5]), -6*sad[2]*Math.sin(sad[5]),                                 0, 3*sad[1]*sad[1]*Math.cos(sad[5])                                   , 4*sad[1]*sad[2]*Math.cos(sad[5]) - 2*sad[1]*sad[1]*Math.sin(sad[5]),  3*sad[2]*sad[2]*Math.cos(sad[5]) - 6*sad[1]*sad[2]*Math.sin(sad[5]), -12*sad[2]*sad[2]*Math.sin(sad[5])},
            {1, top[1], top[2], top[1]*top[1], top[1]*top[2], top[2]*top[2], top[1]*top[1]*top[1], top[1]*top[1]*top[2], top[1]*top[2]*top[2], top[2]*top[2]*top[2], top[1]*top[1]*top[1]*top[1], top[1]*top[1]*top[1]*top[2], top[1]*top[1]*top[2]*top[2], top[1]*top[2]*top[2]*top[2], top[2]*top[2]*top[2]*top[2]},
            {0,      1,      0,      2*top[1],        top[2],             0,      3*top[1]*top[1],      2*top[1]*top[2],        top[2]*top[2],                    0,      4*top[1]*top[1]*top[1],      3*top[1]*top[1]*top[2],      2*top[1]*top[2]*top[2],        top[2]*top[2]*top[2],                           0},
            {0,      0,      1,             0,        top[1],      2*top[2],                    0,        top[1]*top[1],      2*top[1]*top[2],      3*top[2]*top[2],                           0,        top[1]*top[1]*top[1],      2*top[1]*top[1]*top[2],      3*top[1]*top[2]*top[2],      4*top[1]*top[1]*top[1]},
            {0,      0,      0, 2*Math.cos(top[5]), -Math.sin(top[5]),                   0, 6*top[1]*Math.cos(top[5]), 2*top[2]*Math.cos(top[5]) - 2*top[1]*Math.sin(top[5]),                             -2*top[2]*Math.sin(top[5]),                          0, 12*top[1]*top[1]*Math.cos(top[5]), 6*top[1]*top[2]*Math.cos(top[5]) - 3*top[1]*top[1]*Math.sin(top[5]), 2*top[2]*top[2]*Math.cos(top[5]) - 4*top[1]*top[2]*Math.sin(top[5]),                                    -3*top[2]*top[2]*Math.sin(top[5]),                                  0},
            {0,      0,      0,                  0,  Math.cos(top[5]), -2*Math.sin(top[5]),                         0, 2*top[1]*Math.cos(top[5])                            ,  2*top[2]*Math.cos(top[5]) - 2*top[1]*Math.sin(top[5]), -6*top[2]*Math.sin(top[5]),                                 0, 3*top[1]*top[1]*Math.cos(top[5])                                   , 4*top[1]*top[2]*Math.cos(top[5]) - 2*top[1]*top[1]*Math.sin(top[5]),  3*top[2]*top[2]*Math.cos(top[5]) - 6*top[1]*top[2]*Math.sin(top[5]), -12*top[2]*top[2]*Math.sin(top[5])}};
        for (int i = 0; i < M.length; i++)
        {
            for (int j = 0; j < M[0].length - 1; j++)
                System.out.print(M[i][j] + ", ");
            System.out.println(M[i][M[0].length - 1]);
        }
        soln = BSpline5.gaussj(M, v);
        System.out.print("\nsoln = , ");
        for (int i = 0; i < soln.length - 1; i++)
            System.out.print(soln[i] + ", ");
        System.out.println(soln[soln.length - 1]);
        System.out.print("\nv = , ");
        for (int i = 0; i < v.length - 1; i++)
            System.out.print(v[i] + ", ");
        System.out.println(v[v.length - 1]);
    }
}
