
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
    private static final double[] bot = new double[] {8.5, 14.831255905699996, 92.22256068144186, 7.555722926464988E-4, 2.0004236297398864E-5, 0.808149959774216, 0.015992107606605212};
    private static final double[] sad = new double[] {8.5, 34.447600065539774, 66.04531557578703, 0.0021065998157774246, -1.14717469323726E-6, 1.0383595814460744, 0.017449175786289933};
    private static final double[] top = new double[] {8.5, 38.01733799460117, 59.67559639591379, 0.002096645845061304, 1.0890957850267363E-6, 1.0810032247751613, 0.01812772511340366};
    private static double[][] M;
    private static double[] v;
    private static double[] soln;

    public static void main (String[] args)
    {
        System.out.println("         , c, d1, d2, F, eig0, eigangle0, eig1");
        System.out.println("input bottom = , " + bot[0] + ", " + bot[1] + ", " + bot[2] + ", " + bot[3] + ", " + bot[4] + ", " + bot[5] + ", " + bot[6]);
        System.out.println("input saddle = , " + sad[0] + ", " + sad[1] + ", " + sad[2] + ", " + sad[3] + ", " + sad[4] + ", " + sad[5] + ", " + sad[6]);
        System.out.println("input top    = , " + top[0] + ", " + top[1] + ", " + top[2] + ", " + top[3] + ", " + top[4] + ", " + top[5] + ", " + top[6]);
        bot[1] -= sad[1];                           // relocate the origin for precision
        bot[2] -= sad[2];
        top[1] -= sad[1];
        top[2] -= sad[2];
        sad[1] = 0;
        sad[2] = 0;

        //sad[4] = sad[6];                            //use eig1 instead of eig0
        //sad[5] -= Math.PI/2;
        v = new double[] {bot[3], 0, 0, bot[4]*Math.cos(bot[5]), -bot[4]*Math.sin(bot[5]),
//                          sad[3], 0, 0, sad[4]*Math.cos(sad[5]), -sad[4]*Math.sin(sad[5]),
                          sad[3], 0, 0, bot[6]*Math.cos(bot[5] - Math.PI/2), -top[6]*Math.sin(top[5] - Math.PI/2),
                          top[3], 0, 0, top[4]*Math.cos(top[5]), -top[4]*Math.sin(top[5])};
        System.out.print("\nv = , ");
        for (int i = 0; i < v.length - 1; i++)
            System.out.print(v[i] + ", ");
        System.out.println(v[v.length - 1]);
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
            //{0,      0,      0, 2*Math.cos(sad[5]), -Math.sin(sad[5]),                   0, 6*sad[1]*Math.cos(sad[5]), 2*sad[2]*Math.cos(sad[5]) - 2*sad[1]*Math.sin(sad[5]),                             -2*sad[2]*Math.sin(sad[5]),                          0, 12*sad[1]*sad[1]*Math.cos(sad[5]), 6*sad[1]*sad[2]*Math.cos(sad[5]) - 3*sad[1]*sad[1]*Math.sin(sad[5]), 2*sad[2]*sad[2]*Math.cos(sad[5]) - 4*sad[1]*sad[2]*Math.sin(sad[5]),                                    -3*sad[2]*sad[2]*Math.sin(sad[5]),                                  0},
            {0,      0,      0, 2*Math.cos(bot[5] - Math.PI/2), -Math.sin(bot[5] - Math.PI/2),                   0, 6*bot[1]*Math.cos(bot[5] - Math.PI/2), 2*bot[2]*Math.cos(bot[5] - Math.PI/2) - 2*bot[1]*Math.sin(bot[5] - Math.PI/2),                             -2*bot[2]*Math.sin(bot[5] - Math.PI/2),                          0, 12*bot[1]*bot[1]*Math.cos(bot[5] - Math.PI/2), 6*bot[1]*bot[2]*Math.cos(bot[5] - Math.PI/2) - 3*bot[1]*bot[1]*Math.sin(bot[5] - Math.PI/2), 2*bot[2]*bot[2]*Math.cos(bot[5] - Math.PI/2) - 4*bot[1]*bot[2]*Math.sin(bot[5] - Math.PI/2),                                    -3*bot[2]*bot[2]*Math.sin(bot[5] - Math.PI/2),                                  0},
            //{0,      0,      0,                  0,  Math.cos(sad[5]), -2*Math.sin(sad[5]),                         0, 2*sad[1]*Math.cos(sad[5])                            ,  2*sad[2]*Math.cos(sad[5]) - 2*sad[1]*Math.sin(sad[5]), -6*sad[2]*Math.sin(sad[5]),                                 0, 3*sad[1]*sad[1]*Math.cos(sad[5])                                   , 4*sad[1]*sad[2]*Math.cos(sad[5]) - 2*sad[1]*sad[1]*Math.sin(sad[5]),  3*sad[2]*sad[2]*Math.cos(sad[5]) - 6*sad[1]*sad[2]*Math.sin(sad[5]), -12*sad[2]*sad[2]*Math.sin(sad[5])},
            {0,      0,      0,                  0,  Math.cos(top[5] - Math.PI/2), -2*Math.sin(top[5] - Math.PI/2),                         0, 2*top[1]*Math.cos(top[5] - Math.PI/2)                            ,  2*top[2]*Math.cos(top[5] - Math.PI/2) - 2*top[1]*Math.sin(top[5] - Math.PI/2), -6*top[2]*Math.sin(top[5] - Math.PI/2),                                 0, 3*top[1]*top[1]*Math.cos(top[5] - Math.PI/2)                                   , 4*top[1]*top[2]*Math.cos(top[5] - Math.PI/2) - 2*top[1]*top[1]*Math.sin(top[5] - Math.PI/2),  3*top[2]*top[2]*Math.cos(top[5] - Math.PI/2) - 6*top[1]*top[2]*Math.sin(top[5] - Math.PI/2), -12*top[2]*top[2]*Math.sin(top[5] - Math.PI/2)},
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
        verify_F();
    }

    private static void verify_F()
    {
        System.out.println("\n             , c, dFdx, dFdy, F, eig0, eigangle0, eig1");
        // re-calculate the properties of F at the 3 critical points
        double[][] curv = new double[][] {{calc_d2Fdx2(bot[1], bot[2]), calc_d2Fdxdy(bot[1], bot[2])}, {calc_d2Fdxdy(bot[1], bot[2]), calc_d2Fdy2(bot[1], bot[2])}};
        double eig0 = (curv[0][0] + curv[1][1] - Math.sqrt((curv[0][0] - curv[1][1])*(curv[0][0] - curv[1][1]) + 4*curv[0][1]*curv[0][1]))/2;
        double eig1 = (curv[0][0] + curv[1][1] + Math.sqrt((curv[0][0] - curv[1][1])*(curv[0][0] - curv[1][1]) + 4*curv[0][1]*curv[0][1]))/2;
        System.out.println("output bot  =, " + bot[0] + ", " + calc_dFdx(bot[1], bot[2]) + ", " + calc_dFdy(bot[1], bot[2]) + ", " + calc_F(bot[1], bot[2])
                               + ", " + eig0 + ", " + Math.atan((curv[0][0] - eig0)/curv[0][1]) + ", " + eig1);
        curv = new double[][] {{calc_d2Fdx2(sad[1], sad[2]), calc_d2Fdxdy(sad[1], sad[2])}, {calc_d2Fdxdy(sad[1], sad[2]), calc_d2Fdy2(sad[1], sad[2])}};
        eig0 = (curv[0][0] + curv[1][1] - Math.sqrt((curv[0][0] - curv[1][1])*(curv[0][0] - curv[1][1]) + 4*curv[0][1]*curv[0][1]))/2;
        eig1 = (curv[0][0] + curv[1][1] + Math.sqrt((curv[0][0] - curv[1][1])*(curv[0][0] - curv[1][1]) + 4*curv[0][1]*curv[0][1]))/2;
        System.out.println("output sad  =, " + sad[0] + ", " + calc_dFdx(sad[1], sad[2]) + ", " + calc_dFdy(sad[1], sad[2]) + ", " + calc_F(sad[1], sad[2])
                               + ", " + eig0 + ", " + Math.atan((curv[0][0] - eig0)/curv[0][1]) + ", " + eig1);
        curv = new double[][] {{calc_d2Fdx2(top[1], top[2]), calc_d2Fdxdy(top[1], top[2])}, {calc_d2Fdxdy(top[1], top[2]), calc_d2Fdy2(top[1], top[2])}};
        eig0 = (curv[0][0] + curv[1][1] - Math.sqrt((curv[0][0] - curv[1][1])*(curv[0][0] - curv[1][1]) + 4*curv[0][1]*curv[0][1]))/2;
        eig1 = (curv[0][0] + curv[1][1] + Math.sqrt((curv[0][0] - curv[1][1])*(curv[0][0] - curv[1][1]) + 4*curv[0][1]*curv[0][1]))/2;
        System.out.println("output top  =, " + top[0] + ", " + calc_dFdx(top[1], top[2]) + ", " + calc_dFdy(top[1], top[2]) + ", " + calc_F(top[1], top[2])
                               + ", " + eig0 + ", " + Math.atan((curv[0][0] - eig0)/curv[0][1]) + ", " + eig1);
    }

    private static double calc_F(double x, double y)
    {
        // F = (1, x, y, x^2, x y, y^2, x^3, x^2 y, x y^2, y^3, x^4, x^3 y, x^2 y^2, x y^3, y^4)
        return soln[0] + soln[1]*x + soln[2]*y + soln[3]*x*x + soln[4]*x*y + soln[5]*y*y
             + soln[6]*x*x*x + soln[7]*x*x*y + soln[8]*x*y*y + soln[9]*y*y*y
             + soln[10]*x*x*x*x + soln[11]*x*x*x*y + soln[12]*x*x*y*y + soln[13]*x*y*y*y + soln[14]*y*y*y*y;
    }

    private static double calc_dFdx(double x, double y)
    {
        return soln[1] + 2*soln[3]*x + soln[4]*y
             + 3*soln[6]*x*x + 2*soln[7]*x*y + soln[8]*y*y
             + 4*soln[10]*x*x*x + 3*soln[11]*x*x*y + 2*soln[12]*x*y*y + soln[13]*y*y*y;
    }

    private static double calc_dFdy(double x, double y)
    {
        return soln[2] + soln[4]*x + 2*soln[5]*y
             + soln[7]*x*x + 2*soln[8]*x*y + 3*soln[9]*y*y
             + soln[11]*x*x*x + 2*soln[12]*x*x*y + 3*soln[13]*x*y*y + 4*soln[14]*y*y*y;
    }

    private static double calc_d2Fdx2(double x, double y)
    {
        return 2*soln[3]
             + 6*soln[6]*x + 2*soln[7]*y
             + 12*soln[10]*x*x + 6*soln[11]*x*y + 2*soln[12]*y*y;
    }

    private static double calc_d2Fdxdy(double x, double y)
    {
        return soln[4]
             + 2*soln[7]*x + 2*soln[8]*y
             + 3*soln[11]*x*x + 4*soln[12]*x*y + 3*soln[13]*y*y;
    }

    private static double calc_d2Fdy2(double x, double y)
    {
        return 2*soln[5]
             + 2*soln[8]*x + 6*soln[9]*y
             + 2*soln[12]*x*x + 6*soln[13]*x*y + 12*soln[14]*y*y;
    }
}
