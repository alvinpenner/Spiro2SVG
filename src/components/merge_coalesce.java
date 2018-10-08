
package components;

// model two local minima plus a saddle point using a quartic equation
// error functional F = x^4/4 + a*x^2/2 + b*x
// see Spiro2SVG Book6, Oct 2018, page 48
// simulate straight-line motion in the (a, b) plane and calculate
// critical points x, F(x), and curvatures
// a_touch is the point at which we tangentially touch the 'merge' boundary

// this is file : \Documents\NetBeansProjects\MyDemo\src\components\merge_coalesce.java

public class merge_coalesce
{
    private static final double a_touch = -3;
    private static final double b_touch = 2*Math.sqrt(-a_touch*a_touch*a_touch/27);
    private static final double dadb = -Math.sqrt(-3/a_touch);
    private static final double a_start = -5; // 2*a_touch;
    private static final double a_end = 0;
    private static final double b_start = b_touch + (a_start - a_touch)/dadb;
    private static final double b_end = b_touch + (a_end - a_touch)/dadb;
    private static final double TOL = 0.0000000001;

    public static void main (String[] args)
    {
        final int N = 100;
        double a, b;
        double[] x;                                         // critical points x

        System.out.println("from ," + b_start + ", " + a_start + ", to," + b_end + ", " + a_end + ", @," + b_touch + ", " + a_touch);
        System.out.println("i,  b,   a,   x0,  x1,  x2,  F0,  F1,  F2,  F''0, F''1, F''2");
        for (int i = 0; i <= N; i++)
        {
            a = a_start + i*(a_end - a_start)/N;
            b = b_start + i*(b_end - b_start)/N;
            x = solve_cubic(a, b);
            System.out.println(i + ", " + (float) b + ", " + (float) a + ", " + x[0] + ", " + x[1] + ", " + x[2]
                                 + ", " + getF(a, b, x[0]) + ", " + getF(a, b, x[1]) + ", " + getF(a, b, x[2])
                                 + ", " + getF2prime(a, x[0])+ ", " + getF2prime(a, x[1])+ ", " + getF2prime(a, x[2]));
            //System.out.println("     " + getFprime(a, b, x[0])+ ", " + getFprime(a, b, x[1])+ ", " + getFprime(a, b, x[2]));
        }
        //scan_quartic();
    }

    private static void scan_quartic()
    {
        final int N = 10;          // number of values of (a, b)
        double a, b, x;
        double xlo = -8;
        double xhi = 10;
        int xsteps = 100;           // number of values of x

        System.out.println("\nscan F from ," + b_start + ", " + a_start + ", to," + b_end + ", " + a_end + ", @," + b_touch + ", " + a_touch);
        for (int i = 0; i <= xsteps; i++)
        {
            x = xlo + i*(xhi - xlo)/xsteps;
            System.out.print((float) x + ", ");
            for (int j = 0; j <= N; j++)
            {
                a = a_start + j*(a_end - a_start)/N;
                b = b_start + j*(b_end - b_start)/N;
                System.out.print((float) getF(a, b, x) + ", ");
            }
            System.out.println();
        }
    }

    private static double[] solve_cubic(double cua, double cub)
    {
        // see Math CRC book, page 392
        double[] xret = new double[] {Double.NaN, Double.NaN, Double.NaN};
        double cud = cub*cub/4 + cua*cua*cua/27;

        if (cud < TOL)
        {
            double myphi = Math.acos(-cub/2/Math.sqrt(-cua*cua*cua/27));
            xret[0] = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3);
            xret[1] = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 2*Math.PI/3);
            xret[2] = 2*Math.sqrt(-cua/3)*Math.cos(myphi/3 + 4*Math.PI/3);
        }
        else
            xret[0] = Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud));
        return xret;
    }

    private static double getF(double a, double b, double x)
    {
        return x*x*x*x/4 + a*x*x/2 + b*x;
    }

    private static double getFprime(double a, double b, double x)
    {
        return x*x*x + a*x + b;
    }

    private static double getF2prime(double a, double x)
    {
        return 3*x*x + a;
    }
}
