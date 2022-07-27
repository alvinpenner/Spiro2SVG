
// calculate quadratic and cubic coefficients of a Rossler torus
// assuming that the linear response has been obtained from 'fit_linear_response()'
// use perturbation theory as per Iooss, Bouc, and my book Chaos IV, p. 27

package rossler;

// this is file : \Documents\NetBeansProjects\RosslerSystem\src\rossler\Perturb_Torus.java

import java.awt.geom.Point2D;

public final class Perturb_Torus
{
    private static double[][][] S = new double[Main.final_Period + 1][3][3];
    private static double[][] xyz = new double[Main.final_Period + 1][3];       // limit cycle
    private static double[][] S_inv = new double[3][3];
    private static double[][] coeff = new double[3][7];     // i = (fx, fy, fz), j = (Cx20, Cx11, Cx02, Cx30, Cx21, Cx12, Cx03)
    private static double fx12 = 9999, fy12 = 9999, fz12 = 9999;

    protected static void calc_coeff()
    {
        // calculate first-order response, S matrix, during one cycle
        // assume that the output has already been made uniform, by running 'fit_linear_response()'

        if (Main.final_Re_V21 == 0 || Main.final_Im_V21 == 0 || Main.project_phi == 0 || Main.project_theta == 0
        ||  Main.final_Period == 0 || Main.skew_transform)          // || eig == 0 || angle == 0)
        {
            System.out.println("Bad data in 'Perturb_Torus.calc_coeff()'");
            return;
        }

        // perform one uniformized pass to save S

        double[] pt6 = new double[6];
        Main.skew_transform = true;                 // make the linear response "uniform"
        double xdot = -Main.final_y - Main.final_z;
        double ydot =  Main.final_x + Main.a*Main.final_y;
        double zdot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
        //System.out.println("vel , " + xdot + ", " + ydot + ", " + zdot + ", " + v);

        System.out.println("\nPython output - Rossler S matrix");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Neimark-Sacker - fit S-matrix perturbation theory response\\n\\");
        System.out.println("incr        , 1.0, " + "\\n\\");
        System.out.println("a_b_c       , " + Main.a + ", " + Main.b + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt , " + Main.final_Period + ", " + Main.final_delt + ",\\n\\");
        System.out.println("x_y_z       , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi, " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("Re_V21_Im_V21, " + Main.final_Re_V21 + ", " + Main.final_Im_V21 + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',z'\"");
        System.out.print("data = np.array([");

        for (int i = 0; i < 3; i++)         // initiallize dx/du, dy/du, or dz/du
        {
            if (i == 0)                                                     // initial (dxdu', dydu', dzdu')
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, Main.invert_from_xp_yp(1, 0, 0, "x"), Main.invert_from_xp_yp(1, 0, 0, "y"), Main.invert_from_xp_yp(1, 0, 0, "z")};
            else if (i == 1)
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, Main.invert_from_xp_yp(0, 1, 0, "x"), Main.invert_from_xp_yp(0, 1, 0, "y"), Main.invert_from_xp_yp(0, 1, 0, "z")};
            else
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, xdot/v, ydot/v, zdot/v};
            //System.out.println("axis_" + i);
            //pt2 = Main.project_2D(pt6[3], pt6[4], pt6[5]);   // final (dxdu', dydu', dzdu')
            //System.out.println("0" + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt6[3], pt6[4], pt6[5]));
            if (i == 0)
            {
                xyz[0][0] = pt6[0];     // save limit cyclr
                xyz[0][1] = pt6[1];
                xyz[0][2] = pt6[2];
            }
            S[0][0][i] = pt6[3];        // new state (x, y, z) is given as a column vector
            S[0][1][i] = pt6[4];
            S[0][2][i] = pt6[5];
            for (int k = 1; k <= Main.final_Period; k++)                     // loop through one cycle
            {
                Main.runge_kutta_rossler6_ddu3(pt6, Main.final_delt, Main.c);
                //pt2 = Main.project_2D(pt6[3], pt6[4], pt6[5]);
                //S[k][i][0] = pt2.x;
                //S[k][i][1] = pt2.y;
                //S[k][i][2] = Main.project_zp(pt6[3], pt6[4], pt6[5]);
                S[k][0][i] = pt6[3];
                S[k][1][i] = pt6[4];
                S[k][2][i] = pt6[5];
                if (i == 0)
                {
                    xyz[k][0] = pt6[0];     // save limit cyclr
                    xyz[k][1] = pt6[1];
                    xyz[k][2] = pt6[2];
                }
                //System.out.println((k + 1) + ", " + pt6[0] + ", " + pt6[1] + ", " + pt6[2] + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt6[3], pt6[4], pt6[5]));
                //System.out.println(i + ", " + k + ", " + pt6[3] + ", " + pt6[4] + ", " + pt6[5]);     // original S matrix
            }
        }
        //System.out.println("\nraw xyz");
        //for (int k = 0; k < S.length; k++)
        //    System.out.println(k + ", " + xyz[k][0] + ", " + xyz[k][1] + ", " + xyz[k][2]);

        //System.out.println("\nS matrix");
        //for (int i = 0; i < 3; i++)
        //{
        //    System.out.println("axis_" + i);
        //    for (int k = 0; k < S.length; k++)
        //    {
        //        Point2D.Double pt2 = Main.project_2D(S[k][0][i], S[k][1][i], S[k][2][i]);
        //        System.out.println(i + ", " + k + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(S[k][0][i], S[k][1][i], S[k][2][i]));
        //        //System.out.println(i + ", " + k + ", " + S[k][0][i] + ", " + S[k][1][i] + ", " + S[k][2][i]);
        //    }
        //}

        //System.out.println("\nS matrix (org)");
        //for (int i = 0; i < 3; i++)
        //    System.out.println(S[Main.final_Period][i][0] + ", " + S[Main.final_Period][i][1] + ", " + S[Main.final_Period][i][2]);
/*
        System.out.println("\nS matrix (primed)");
        int col = 1;                // temporary code
        Point2D.Double pt2 = Main.project_2D(S[Main.final_Period][0][col], S[Main.final_Period][1][col], S[Main.final_Period][2][col]);
        System.out.println(pt2.x);
        System.out.println(pt2.y);
        System.out.println(Main.project_zp(S[Main.final_Period][0][col], S[Main.final_Period][1][col], S[Main.final_Period][2][col]));
        System.out.println();       // temporary code
        pt2 = Main.project_2D(xdot, ydot, zdot);
        System.out.println("xyz'dot = " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(xdot, ydot, zdot));
        System.out.println("v     = " + v);
*/
        //System.out.println("gen_quadratic_Cxy");
        for (int i = -2; i < 3; i++)                            // increment x' by i
            for (int j = -2; j < 3; j++)                        // increment y' by j
                gen_time_sync_Cxy(i, j);
        System.out.println("])");

        gen_quadratic_extrapolate_Cxy();        // enhance the t_sync result (numerical curve-fit)
        gen_quadratic_response();               // calculate analytical partial derivatives
        //gen_quadratic_Cxy(0, 0);              // temporary test code only!
        //integrate_state_i(1);
        //differentiate_state_i(2, 1000);
        //partial_integrate(2, 3, 2001);

        //int k = Main.final_Period;                              // print S (projected) at time = k
        //System.out.println("\nS matrix at time t = " + k);
        //for (int i = 0; i < 3; i++)
        //{
        //    Point2D.Double pt2 = Main.project_2D(S[k][i][0], S[k][i][1], S[k][i][2]);
        //    System.out.println(i + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(S[k][i][0], S[k][i][1], S[k][i][2]));
        //}

        // Define S0, the transform from the primed system to the unprimed, at t = 0
/*
        Point2D.Double pt2;
        System.out.println("\nS0 matrix (at t = 0)");
        for (int i = 0; i < 3; i++)
        {
            System.out.print(i + ", " + S[0][i][0] + ", " + S[0][i][1] + ", " + S[0][i][2]);
            if (i == 0)
            {
                pt2 = Main.project_2D(1, 0, 0);
                System.out.println(", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(1, 0, 0));
            }
            else if (i == 1)
            {
                pt2 = Main.project_2D(0, 1, 0);
                System.out.println(", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(0, 1, 0));
            }
            else if (i == 2)
            {
                pt2 = Main.project_2D(0, 0, 1);
                System.out.println(", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(0, 0, 1));
            }
        }
*/
    }

    private static void differentiate_state_i(int istate, int tk)   // test code only
    {
        System.out.println(istate + ", " + tk + ", " + S[tk][0][istate] + ", " + S[tk][1][istate] + ", " + S[tk][2][istate]);
        System.out.println(istate + ", " + (tk + 1) + ", " + S[tk + 1][0][istate] + ", " + S[tk + 1][1][istate] + ", " + S[tk + 1][2][istate]);
        System.out.println("diff,, " + (S[tk + 1][0][istate] - S[tk - 1][0][istate])/2/Main.final_delt + ", " + (S[tk + 1][1][istate] - S[tk - 1][1][istate])/2/Main.final_delt + ", " + (S[tk + 1][2][istate] - S[tk - 1][2][istate])/2/Main.final_delt);
        System.out.println("anal,, " + (-S[tk][1][istate] - S[tk][2][istate])
                              + ", " + (S[tk][0][istate] + Main.a*S[tk][1][istate])
                              + ", " + (xyz[tk][2]*S[tk][0][istate] + (xyz[tk][0] - Main.c)*S[tk][2][istate]));
    }

    private static void integrate_state_i(int istate)               // test code only
    {
        // test the integration of the linear system in a single pure state i
        double[] accum_vec = new double[3];
        double[] state_vec = new double[3];
        Point2D.Double pt2;

        pt2 = Main.project_2D(S[0][0][istate], S[0][1][istate], S[0][2][istate]);   // final (dxdu', dydu', dzdu')
        System.out.println("start, " + istate + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(S[0][0][istate], S[0][1][istate], S[0][2][istate]));
        for (int i = 0; i < 3; i++)
            accum_vec[i] = S[0][i][istate];
        for (int k = 0; k < S.length; k++)                  // scan one limit cycle
        {
            state_vec[0] = -S[k][1][istate] - S[k][2][istate];
            state_vec[1] =  S[k][0][istate] + Main.a*S[k][1][istate];
            state_vec[2] =  xyz[k][2]*S[k][0][istate] + (xyz[k][0] - Main.c)*S[k][2][istate];
            if (k == 0 || k == Main.final_Period)
                for (int i = 0; i < 3; i++)
                    state_vec[i] /= 2.0;                    // correct the endpoints
            for (int i = 0; i < 3; i++)
                accum_vec[i] += state_vec[i]*Main.final_delt;
        }
        pt2 = Main.project_2D(accum_vec[0], accum_vec[1], accum_vec[2]);   // final (dxdu', dydu', dzdu')
        System.out.println("end  , " + istate + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(accum_vec[0], accum_vec[1], accum_vec[2]));
    }

    private static void partial_integrate(int ix, int iy, int tk)
    {
        // partial integration up to time tk

        double[] accum_vec = new double[3];
        double[] state_vec = new double[3];
        double[] M_vec = new double[3];

        accum_vec[0] = ix;                              // initial state, x'
        accum_vec[1] = iy;                              // y' coord
        accum_vec[2] = 0;                               // z' coord
        for (int k = 0; k < tk; k++)                  // scan up to time tk
        {
            // linear combination of dx and dy response (col 0 and col 1)

            for (int i = 0; i < 3; i++)
                state_vec[i] = ix*S[k][i][0] + iy*S[k][i][1];   // new state

            // form non-linear response M = (0, 0, dz*dx)

            M_vec[0] = 0;
            M_vec[1] = 0;
            M_vec[2] = state_vec[0]*state_vec[2];

            invert_S(k);    // calculate S_inverse
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    accum_vec[i] += S_inv[i][j]*M_vec[j]*Main.final_delt;
        }
        System.out.println("delt        ,,,," + Main.final_delt);
        System.out.println("Mz(tk)      ,,,,0,0," + M_vec[2]);

        for (int i = 0; i < 3; i++)                     // project to end of cycle
        {
            state_vec[i] = 0;
            for (int j = 0; j < 3; j++)
                state_vec[i] += S[tk][i][j]*accum_vec[j];
        }
        System.out.println("final state ," + tk + ", " + ix + ", " + iy + ", " + state_vec[0] + ", " + state_vec[1] + ", " + state_vec[2]);
        // calculate Jacobian*final_state
        System.out.println("Jac*final_state ,,,," + (-state_vec[1] - state_vec[2])
                                           + ", " + (state_vec[0] + Main.a*state_vec[1])
                                           + ", " + (xyz[tk][2]*state_vec[0] + (xyz[tk][0] - Main.c)*state_vec[2]));
    }

    private static void gen_time_sync_Cxy(int ix, int iy)
    {
        // calculate second-order response, based on S matrix, during one cycle
        // note that this is an integration to a specific time endpoint, not a space-endpoint
        // ix and iy are coordinates in the primed coord system (assuming z' = 0)

        double[] state_y1 = new double[3];
        double[][] state_y2 = new double[Main.final_Period + 1][3];
        double[] state_y3 = new double[3];
        double[] state_final = new double[3];
        double[] y_incr = new double[3];
        double[] y_accum = new double[3];
        double[] y_accum_half = new double[3];     // accumulate only half at an endpoint
        double[] M_vec = new double[3];
        Point2D.Double pt2;
        double zp;
        int i;

        // quadratic component

        for (int k = 0; k < S.length; k++)                  // scan one limit cycle
        {
            // linear combination of dx' and dy' response (col 0 and col 1)

            for (i = 0; i < 3; i++)
                state_y1[i] = ix*S[k][i][0] + iy*S[k][i][1];    // perturbed state (first-order)

            // form non-linear response M = (0, 0, dz*dx)

            M_vec[0] = 0;
            M_vec[1] = 0;
            M_vec[2] = state_y1[0]*state_y1[2];
            invert_S(k);                                            // calculate S_inverse
            for (i = 0; i < 3; i++)
            {
                y_incr[i] = 0;
                for (int j = 0; j < 3; j++)
                    y_incr[i] += S_inv[i][j]*M_vec[j];
                if (k == 0)
                {
                    y_accum_half[i] = 0;                            // state is unchanged
                    y_accum[i] = y_incr[i]/2;                       // use only half of first point
                }
                else
                {
                    y_accum_half[i] = y_accum[i] + y_incr[i]/2;     // use only for endpoint
                    y_accum[i] += y_incr[i];                        // use for next time increment
                }
            }
            for (i = 0; i < 3; i++)                                 // project to end of cycle
            {
                state_y2[k][i] = 0;
                for (int j = 0; j < 3; j++)
                    state_y2[k][i] += S[k][i][j]*y_accum_half[j];
                state_y2[k][i] *= Main.final_delt;
            }
            if (ix == 1 && iy == 1)
                System.out.println("state_y2 , " + ix + ", " + iy + ", " + k + ", " + state_y2[k][0] + ", " + state_y2[k][1] + ", " + state_y2[k][2]);
        }

        // cubic component

        for (int k = 0; k < S.length; k++)                          // scan one limit cycle
        {
            for (i = 0; i < 3; i++)
                state_y1[i] = ix*S[k][i][0] + iy*S[k][i][1];        // perturbed state (first-order)

            // form non-linear response M = (0, 0, dz*dx)

            M_vec[0] = 0;
            M_vec[1] = 0;
            M_vec[2] = state_y1[0]*state_y2[k][2] + state_y1[2]*state_y2[k][0];
            invert_S(k);                                            // calculate S_inverse
            for (i = 0; i < 3; i++)
            {
                y_incr[i] = 0;
                for (int j = 0; j < 3; j++)
                    y_incr[i] += S_inv[i][j]*M_vec[j];
                if (k == 0)
                {
                    y_accum_half[i] = 0;                            // state is unchanged
                    y_accum[i] = y_incr[i]/2;                       // use only half of first point
                }
                else
                {
                    y_accum_half[i] = y_accum[i] + y_incr[i]/2;     // use only for endpoint
                    y_accum[i] += y_incr[i];                        // use for next time increment
                }
            }
        }
        for (i = 0; i < 3; i++)                                     // project to end of cycle
        {
            state_y3[i] = 0;
            for (int j = 0; j < 3; j++)
                state_y3[i] += S[Main.final_Period][i][j]*y_accum_half[j];
            state_y3[i] *= Main.final_delt;
        }

        // summary

        for (i = 0; i < 3; i++)
            state_final[i] = ix*S[Main.final_Period][i][0] + iy*S[Main.final_Period][i][1] + state_y2[Main.final_Period][i] + state_y3[i];
        pt2 = Main.project_2D(state_final[0], state_final[1], state_final[2]);   // final (dxdu', dydu', dzdu')
        zp  = Main.project_zp(state_final[0], state_final[1], state_final[2]);
        System.out.print("[" + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp + "]");
        if (ix < 2 || iy < 2)
            System.out.println(",");

        // calculate t_sync version of quadratic coeff (Cx20, Cx11, Cx02)

        pt2 = Main.project_2D(state_y2[Main.final_Period][0], state_y2[Main.final_Period][1], state_y2[Main.final_Period][2]);
        zp  = Main.project_zp(state_y2[Main.final_Period][0], state_y2[Main.final_Period][1], state_y2[Main.final_Period][2]);
        if (ix == 1 && iy == 0)         // Cx20
        {
            //System.out.println("unprojected coeff (col 0) =, " + state_y2[Main.final_Period][0] + ", " + state_y2[Main.final_Period][1] + ", " + state_y2[Main.final_Period][2]);
            coeff[0][0] = pt2.x;
            coeff[1][0] = pt2.y;
            coeff[2][0] = zp;
        }
        if (ix == 0 && iy == 1)         // Cx02
        {
            //System.out.println("unprojected coeff (col 2) =, " + state_y2[Main.final_Period][0] + ", " + state_y2[Main.final_Period][1] + ", " + state_y2[Main.final_Period][2]);
            coeff[0][2] = pt2.x;
            coeff[1][2] = pt2.y;
            coeff[2][2] = zp;
        }
        if (ix == 1 && iy == 1)         // Cx11 (this must be chronologically the last to execute)
        {
            //System.out.println("unprojected coeff (col 1 RAW) =, " + state_y2[Main.final_Period][0] + ", " + state_y2[Main.final_Period][1] + ", " + state_y2[Main.final_Period][2]);
            coeff[0][1] = pt2.x - coeff[0][0] - coeff[0][2];
            coeff[1][1] = pt2.y - coeff[1][0] - coeff[1][2];
            coeff[2][1] = zp - coeff[2][0] - coeff[2][2];
        }

        // calculate t_sync version of cubic coeff (Cx30, Cx21, Cx12, Cx03) (see Book IV, p.41)

        pt2 = Main.project_2D(state_y3[0], state_y3[1], state_y3[2]);
        zp  = Main.project_zp(state_y3[0], state_y3[1], state_y3[2]);
        if (ix == 1 && iy == 0)         // Cx30
        {
            coeff[0][3] = pt2.x;
            coeff[1][3] = pt2.y;
            coeff[2][3] = zp;
        }
        if (ix == 0 && iy == 1)         // Cx03
        {
            coeff[0][6] = pt2.x;
            coeff[1][6] = pt2.y;
            coeff[2][6] = zp;
        }
        if (ix == 1 && iy == 2)
        {
            fx12 = pt2.x;               // fx12 (temporary storage)
            fy12 = pt2.y;
            fz12 = zp;
        }
        if (ix == 2 && iy == 1)         // Cx21 (this must be chronologically the last to execute)
        {
            coeff[0][4] = (-15*coeff[0][3] + 2*pt2.x - 1*fx12 + 6*coeff[0][6])/6;   // Cx21
            coeff[1][4] = (-15*coeff[1][3] + 2*pt2.y - 1*fy12 + 6*coeff[1][6])/6;
            coeff[2][4] = (-15*coeff[2][3] + 2*zp - 1*fz12 + 6*coeff[2][6])/6;
            coeff[0][5] = (  6*coeff[0][3] - 1*pt2.x + 2*fx12 - 15*coeff[0][6])/6;  // Cx12
            coeff[1][5] = (  6*coeff[1][3] - 1*pt2.y + 2*fy12 - 15*coeff[1][6])/6;
            coeff[2][5] = (  6*coeff[2][3] - 1*zp + 2*fz12 - 15*coeff[2][6])/6;
        }
    }

    private static void gen_quadratic_extrapolate_Cxy()
    {
        // assume we have calculated second-order response, based on S matrix, during one cycle
        // then enhance the theory with a quadratic extrapolation
        // from t_sync to z_sync

        double incr = 0.001; // 0.001;
        double x0dot = -Main.final_y - Main.final_z;
        double y0dot =  Main.final_x + Main.a*Main.final_y;
        double z0dot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        double zpnew;
        //pt2new = Main.project_2D(x0dot, y0dot, z0dot);
        //zpnew = Main.project_zp(x0dot, y0dot, z0dot);
        //System.out.print("test velocity projection ," + pt2new.x + ", " + pt2new.y + ", " + zpnew + ", " + Math.sqrt(x0dot*x0dot + y0dot*y0dot + z0dot*z0dot));

        System.out.println("\nPython output - Rossler S matrix (enhanced)");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Neimark-Sacker - fit S-matrix perturbation theory (enhanced)\\n\\");
        System.out.println("incr        , " + incr + ",\\n\\");
        System.out.println("a_b_c       , " + Main.a + ", " + Main.b + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt , " + Main.final_Period + ", " + Main.final_delt + ",\\n\\");
        System.out.println("x_y_z       , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi, " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("Re_V21_Im_V21, " + Main.final_Re_V21 + ", " + Main.final_Im_V21 + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',z'\"");
        System.out.print("data = np.array([");

        for (int i = -2; i < 3; i++)                            // increment x' by i
            for (int j = -2; j < 3; j++)                        // increment y' by j
            {
                double xp = incr*incr*(coeff[0][0]*i*i + coeff[0][1]*i*j + coeff[0][2]*j*j)    // primed coordinates
                          + incr*incr*incr*(coeff[0][3]*i*i*i + coeff[0][4]*i*i*j + coeff[0][5]*i*j*j + coeff[0][6]*j*j*j);
                double yp = incr*incr*(coeff[1][0]*i*i + coeff[1][1]*i*j + coeff[1][2]*j*j)
                          + incr*incr*incr*(coeff[1][3]*i*i*i + coeff[1][4]*i*i*j + coeff[1][5]*i*j*j + coeff[1][6]*j*j*j);
                double zp = incr*incr*(coeff[2][0]*i*i + coeff[2][1]*i*j + coeff[2][2]*j*j)
                          + incr*incr*incr*(coeff[2][3]*i*i*i + coeff[2][4]*i*i*j + coeff[2][5]*i*j*j + coeff[2][6]*j*j*j);
                double x = Main.final_x + S[Main.final_Period][0][0]*incr*i + S[Main.final_Period][0][1]*incr*j + Main.invert_from_xp_yp(xp, yp, zp, "x");     // original coordinates
                double y = Main.final_y + S[Main.final_Period][1][0]*incr*i + S[Main.final_Period][1][1]*incr*j + Main.invert_from_xp_yp(xp, yp, zp, "y");
                double z = Main.final_z + S[Main.final_Period][2][0]*incr*i + S[Main.final_Period][2][1]*incr*j + Main.invert_from_xp_yp(xp, yp, zp, "z");
                double xdot = -y - z;
                double ydot =  x + Main.a*y;
                double zdot =  Main.b + z*(x - Main.c);
                double x2dot = -ydot - zdot;
                double y2dot =  xdot + Main.a*ydot;
                double z2dot =  zdot*(x - Main.c) + z*xdot;
                double x3dot = -y2dot - z2dot;
                double y3dot =  x2dot + Main.a*y2dot;
                double z3dot =  z2dot*(x - Main.c) + 2*zdot*xdot + z*x2dot;
                double A0 = (x0dot*x3dot + y0dot*y3dot + z0dot*z3dot)/6;
                double A = (x0dot*x2dot + y0dot*y2dot + z0dot*z2dot)/2;
                double B = x0dot*xdot + y0dot*ydot + z0dot*zdot;
                double C = x0dot*(x - Main.final_x) + y0dot*(y - Main.final_y) + z0dot*(z - Main.final_z);
                double alpha = (-B + Math.sqrt(B*B - 4*A*C))/2/A;
                alpha = Main.solve_cubic(A/A0, B/A0, C/A0);             // over-ride quadratic solution
                //alpha = 0;
                //alpha = -C/B; // *(1 + A*C/B/B);                      // fix fix bug bug
                //System.out.println("vel = " + i + ", " + j + ", " + x + ", " + y + ", " + z + ", " + xdot + ", " + ydot + ", " + zdot);
                //System.out.println("alpha = " + i + ", " + j + ", " + (-C/B) + ", " + (-C/B*(1 + A*C/B/B)) + ", " + alpha);
                //System.out.println("dot = ," + i + ", " + j + ", " + (x0dot*(x - Main.final_x + alpha*xdot + alpha*alpha*x2dot/2 + alpha*alpha*alpha*x3dot/6)
                //                                                   +  y0dot*(y - Main.final_y + alpha*ydot + alpha*alpha*y2dot/2 + alpha*alpha*alpha*y3dot/6)
                //                                                   +  z0dot*(z - Main.final_z + alpha*zdot + alpha*alpha*z2dot/2 + alpha*alpha*alpha*z3dot/6)));
                pt2new = Main.project_2D(x + alpha*xdot + alpha*alpha*x2dot/2 + alpha*alpha*alpha*x3dot/6, y + alpha*ydot + alpha*alpha*y2dot/2 + alpha*alpha*alpha*y3dot/6, z + alpha*zdot + alpha*alpha*z2dot/2 + alpha*alpha*alpha*z3dot/6);
                zpnew = Main.project_zp (x + alpha*xdot + alpha*alpha*x2dot/2 + alpha*alpha*alpha*x3dot/6, y + alpha*ydot + alpha*alpha*y2dot/2 + alpha*alpha*alpha*y3dot/6, z + alpha*zdot + alpha*alpha*z2dot/2 + alpha*alpha*alpha*z3dot/6);
                //pt2new = Main.project_2D(x + alpha*xdot + alpha*alpha*x2dot/2, y + alpha*ydot + alpha*alpha*y2dot/2, z + alpha*zdot + alpha*alpha*z2dot/2);
                //zpnew = Main.project_zp (x + alpha*xdot + alpha*alpha*x2dot/2, y + alpha*ydot + alpha*alpha*y2dot/2, z + alpha*zdot + alpha*alpha*z2dot/2);
                //pt2new = Main.project_2D(x + alpha*xdot + 0*alpha*alpha*x2dot/2, y + alpha*ydot + 0*alpha*alpha*y2dot/2, z + alpha*zdot + 0*alpha*alpha*z2dot/2);
                //zpnew = Main.project_zp (x + alpha*xdot + 0*alpha*alpha*x2dot/2, y + alpha*ydot + 0*alpha*alpha*y2dot/2, z + alpha*zdot + 0*alpha*alpha*z2dot/2);
                //System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + (Main.final_Period + alpha/Main.final_delt) + "]");
                //System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + alpha + "]");
                System.out.print("[" + i + ", " + j + ", " + pt2new.x + ", " + pt2new.y + ", " + zpnew + "]");  // USE THIS!!!
                if (i < 2 || j < 2)
                    System.out.println(",");
            }
        System.out.println("])");
    }

    private static void gen_quadratic_response()
    {
        // assume we have calculated second-order (t-sync) response, based on S matrix, during one cycle
        // then define the partial derivatives of the output (x, y, z) response wrt input (x', y')

        double x0dot = -Main.final_y - Main.final_z;
        double y0dot =  Main.final_x + Main.a*Main.final_y;
        double z0dot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        double x02dot = -y0dot - z0dot;
        double y02dot =  x0dot + Main.a*y0dot;
        double z02dot =  z0dot*(Main.final_x - Main.c) + Main.final_z*x0dot;
        double x03dot = -y02dot - z02dot;
        double y03dot =  x02dot + Main.a*y02dot;
        double z03dot =  z02dot*(Main.final_x - Main.c) + 2*z0dot*x0dot + Main.final_z*x02dot;
        double[][] dPdi = new double[][] {{S[Main.final_Period][0][0], S[Main.final_Period][0][1]},                     // dP/di (redundant)
                                           {S[Main.final_Period][1][0], S[Main.final_Period][1][1]},                    // used to calculate dPdot/di
                                           {S[Main.final_Period][2][0], S[Main.final_Period][2][1]}};                   // row = (x,y,z), col = (i,j)
        double[][] dPdotdi = new double[][] {{-dPdi[1][0] - dPdi[2][0]       , -dPdi[1][1] - dPdi[2][1]},               // dPdot/di
                                             { dPdi[0][0] + Main.a*dPdi[1][0],  dPdi[0][1] + Main.a*dPdi[1][1]},
                                             { dPdi[2][0]*(Main.final_x - Main.c) + Main.final_z*dPdi[0][0], dPdi[2][1]*(Main.final_x - Main.c) + Main.final_z*dPdi[0][1]}};
        double[][] dP2dotdi = new double[][] {{-dPdotdi[1][0] - dPdotdi[2][0]       , -dPdotdi[1][1] - dPdotdi[2][1]},  // dP2dot/di
                                              { dPdotdi[0][0] + Main.a*dPdotdi[1][0],  dPdotdi[0][1] + Main.a*dPdotdi[1][1]},
                                              { dPdotdi[2][0]*(Main.final_x - Main.c) + z0dot*dPdi[0][0] + Main.final_z*dPdotdi[0][0] + dPdi[2][0]*x0dot,
                                                dPdotdi[2][1]*(Main.final_x - Main.c) + z0dot*dPdi[0][1] + Main.final_z*dPdotdi[0][1] + dPdi[2][1]*x0dot}};
        double[][] d2Pdidj = new double[3][3];                                  // row = (x,y,z), col = (ii,ij,jj)
        double[][] d2Pdotdidj = new double[3][3];                               // row = (x,y,z), col = (ii,ij,jj)
        double PdotPdot = x0dot*x0dot + y0dot*y0dot + z0dot*z0dot;
        double PdotP2dot = x0dot*x02dot + y0dot*y02dot + z0dot*z02dot;
        double[] dalphadi = new double[] {-(x0dot*dPdi[0][0] + y0dot*dPdi[1][0] + z0dot*dPdi[2][0])/PdotPdot,
                                          -(x0dot*dPdi[0][1] + y0dot*dPdi[1][1] + z0dot*dPdi[2][1])/PdotPdot};
        double[] d2alphadidj = new double[3];
        int i;

        // summarize quadratic coeff

        System.out.println("\ndPdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dPdi[i][0] + ", " + dPdi[i][1]);

        System.out.println("\nd2Pdidj :");          // convert quadratic coeff from primed to org units
        for (i = 0; i < 3; i++)
        {
            d2Pdidj[0][i] = Main.invert_from_xp_yp(coeff[0][i], coeff[1][i], coeff[2][i], "x");
            d2Pdidj[1][i] = Main.invert_from_xp_yp(coeff[0][i], coeff[1][i], coeff[2][i], "y");
            d2Pdidj[2][i] = Main.invert_from_xp_yp(coeff[0][i], coeff[1][i], coeff[2][i], "z");
            d2Pdotdidj[0][i] = -d2Pdidj[1][i] - d2Pdidj[2][i];                  // time derivative
            d2Pdotdidj[1][i] =  d2Pdidj[0][i] + Main.a*d2Pdidj[1][i];           // see Chaos IV, p. 40
            d2Pdotdidj[2][i] =  Main.final_z*d2Pdidj[0][i] + (Main.final_x - Main.c)*d2Pdidj[2][i];
        }
        d2Pdotdidj[2][0] += S[Main.final_Period][0][0]*S[Main.final_Period][2][0];  // add quadratic M(S, S) term
        d2Pdotdidj[2][1] += S[Main.final_Period][0][0]*S[Main.final_Period][2][1] + S[Main.final_Period][0][1]*S[Main.final_Period][2][0];  // (1,1) - (0,1) - (1,0)
        d2Pdotdidj[2][2] += S[Main.final_Period][0][1]*S[Main.final_Period][2][1];
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + d2Pdidj[i][0] + ", " + d2Pdidj[i][1] + ", " + d2Pdidj[i][2]);
        System.out.println("\nd2Pdotdidj :");          // convert quadratic coeff from primed to org units
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + d2Pdotdidj[i][0] + ", " + d2Pdotdidj[i][1] + ", " + d2Pdotdidj[i][2]);

        System.out.println("\ndalphadi :");
        System.out.println(-dalphadi[0] + ", " + -dalphadi[1]);
        System.out.println("\nd2alphadidj :");
        d2alphadidj[0] = 0.5*(2.0*(x0dot*d2Pdidj[0][0] + y0dot*d2Pdidj[1][0] + z0dot*d2Pdidj[2][0])
                         + 2*dalphadi[0]*(x0dot*dPdotdi[0][0] + y0dot*dPdotdi[1][0] + z0dot*dPdotdi[2][0])
                         + dalphadi[0]*dalphadi[0]*PdotP2dot)/PdotPdot;
        d2alphadidj[1] = ((x0dot*d2Pdidj[0][1] + y0dot*d2Pdidj[1][1] + z0dot*d2Pdidj[2][1])
                         + dalphadi[0]*(x0dot*dPdotdi[0][1] + y0dot*dPdotdi[1][1] + z0dot*dPdotdi[2][1])
                         + dalphadi[1]*(x0dot*dPdotdi[0][0] + y0dot*dPdotdi[1][0] + z0dot*dPdotdi[2][0])
                         + dalphadi[0]*dalphadi[1]*PdotP2dot)/PdotPdot;
        d2alphadidj[2] = 0.5*(2.0*(x0dot*d2Pdidj[0][2] + y0dot*d2Pdidj[1][2] + z0dot*d2Pdidj[2][2])
                         + 2*dalphadi[1]*(x0dot*dPdotdi[0][1] + y0dot*dPdotdi[1][1] + z0dot*dPdotdi[2][1])
                         + dalphadi[1]*dalphadi[1]*PdotP2dot)/PdotPdot;
        System.out.println(d2alphadidj[0] + ", " + d2alphadidj[1] + ", " + d2alphadidj[2]);

        System.out.println("\ndPdotdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dPdotdi[i][0] + ", " + dPdotdi[i][1]);

        System.out.println("\ndP2dotdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dP2dotdi[i][0] + ", " + dP2dotdi[i][1]);

        System.out.println("\nCxy quadratic/cubic coeff (Type 1):");
        for (i = 0; i < 3; i++)                     // Type 1: proportional to P0
            System.out.println(i + ", " + coeff[i][0] + ", " + coeff[i][1] + ", " + coeff[i][2] + ", " + coeff[i][3] + ", " + coeff[i][4] + ", " + coeff[i][5] + ", " + coeff[i][6]);

        double[][] coeff2 = new double[3][3];       // i = (fx, fy, fz), j = (Cxx, Cxy< Cyy)
        System.out.println("\nCxy quadratic coeff (Type 2):");
        for (i = 0; i < 3; i++)                     // Type 2: proportional to dalpha*dP0dot
        {
            coeff2[i][0] = dalphadi[0]*dPdotdi[i][0];
            coeff2[i][1] = dalphadi[0]*dPdotdi[i][1] + dalphadi[1]*dPdotdi[i][0];
            coeff2[i][2] = dalphadi[1]*dPdotdi[i][1];
            //System.out.println(i + ", " + coeff2[i][0] + ", " + coeff2[i][1] + ", " + coeff2[i][2]);
        }
        double[][] coeff2proj = new double[3][3];   // project coeff2 to primed coord
        Point2D.Double pt2temp;
        double zptemp;
        for (i = 0; i < 3; i++)
        {
            pt2temp = Main.project_2D(coeff2[0][i], coeff2[1][i], coeff2[2][i]);
            zptemp = Main.project_zp (coeff2[0][i], coeff2[1][i], coeff2[2][i]);
            coeff2proj[0][i] = pt2temp.x;
            coeff2proj[1][i] = pt2temp.y;
            coeff2proj[2][i] = zptemp;
        }
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + coeff2proj[i][0] + ", " + coeff2proj[i][1] + ", " + coeff2proj[i][2]);

        System.out.println("\nCxy quadratic coeff (Type 3):");
        pt2temp = Main.project_2D(x02dot, y02dot, z02dot);
        zptemp = Main.project_zp (x02dot, y02dot, z02dot);
        System.out.println("x, " + pt2temp.x*dalphadi[0]*dalphadi[0]/2 + ", " + pt2temp.x*dalphadi[0]*dalphadi[1] + ", " + pt2temp.x*dalphadi[1]*dalphadi[1]/2);
        System.out.println("y, " + pt2temp.y*dalphadi[0]*dalphadi[0]/2 + ", " + pt2temp.y*dalphadi[0]*dalphadi[1] + ", " + pt2temp.y*dalphadi[1]*dalphadi[1]/2);
        System.out.println("z, " + zptemp*dalphadi[0]*dalphadi[0]/2 + ", " + zptemp*dalphadi[0]*dalphadi[1] + ", " + zptemp*dalphadi[1]*dalphadi[1]/2);

        System.out.println("\nCxy coeff (final):");
        Point2D.Double pt2eig = Main.project_2D(S[Main.final_Period][0][0], S[Main.final_Period][1][0], S[Main.final_Period][2][0]);
        System.out.print(Main.project_psi + ", " + Main.a + ", " + Main.b + ", " + Main.c + ", " + Main.final_Period + ", " + Main.final_delt);
        //System.out.print(",0," + pt2eig.x + ",0, " + (coeff[0][0] + tempout[0][0] + pt2temp.x*vec[0]) + ", " + (coeff[0][1] + tempout[0][1] + pt2temp.x*vec[1]) + ", " + (coeff[0][2] + tempout[0][2] + pt2temp.x*vec[2]) + ",0,0,0,0");
        //System.out.println(",0," + pt2eig.y + ",0, " + (coeff[1][0] + tempout[1][0] + pt2temp.y*vec[0]) + ", " + (coeff[1][1] + tempout[1][1] + pt2temp.y*vec[1]) + ", " + (coeff[1][2] + tempout[1][2] + pt2temp.y*vec[2]) + ",0,0,0,0");
        //System.out.print(",0," + pt2eig.x + ",0, " + coeff[0][0] + ", " + coeff[0][1] + ", " + coeff[0][2] + ", " + coeff[0][3] + ", " + coeff[0][4] + ", " + coeff[0][5] + ", " + coeff[0][6]);    // temporary fudge
        //System.out.println(",0," + pt2eig.y + ",0, " + coeff[1][0] + ", " + coeff[1][1] + ", " + coeff[1][2] + ", " + coeff[1][3] + ", " + coeff[1][4] + ", " + coeff[1][5] + ", " + coeff[1][6]);  // temporary fudge
        System.out.print(",0," + pt2eig.x + ",0, " + (coeff[0][0] + coeff2proj[0][0] + pt2temp.x*dalphadi[0]*dalphadi[0]/2) + ", " + (coeff[0][1] + coeff2proj[0][1] + pt2temp.x*dalphadi[0]*dalphadi[1]) + ", " + (coeff[0][2] + coeff2proj[0][2] + pt2temp.x*dalphadi[1]*dalphadi[1]/2) + ", " + coeff[0][3] + ", " + coeff[0][4] + ", " + coeff[0][5] + ", " + coeff[0][6]);
        System.out.println(",0," + pt2eig.y + ",0, " + (coeff[1][0] + coeff2proj[1][0] + pt2temp.y*dalphadi[0]*dalphadi[0]/2) + ", " + (coeff[1][1] + coeff2proj[1][1] + pt2temp.y*dalphadi[0]*dalphadi[1]) + ", " + (coeff[1][2] + coeff2proj[1][2] + pt2temp.y*dalphadi[1]*dalphadi[1]/2) + ", " + coeff[1][3] + ", " + coeff[1][4] + ", " + coeff[1][5] + ", " + coeff[1][6]);
    }

    private static void invert_S(int row)
    {
        double det = S[row][0][0]*S[row][1][1]*S[row][2][2] + S[row][0][1]*S[row][1][2]*S[row][2][0] + S[row][0][2]*S[row][1][0]*S[row][2][1]
                   - S[row][0][2]*S[row][1][1]*S[row][2][0] - S[row][0][0]*S[row][1][2]*S[row][2][1] - S[row][0][1]*S[row][1][0]*S[row][2][2];
        S_inv[0][0] =  (S[row][1][1]*S[row][2][2] - S[row][1][2]*S[row][2][1])/det;
        S_inv[1][0] = -(S[row][1][0]*S[row][2][2] - S[row][1][2]*S[row][2][0])/det;
        S_inv[2][0] =  (S[row][1][0]*S[row][2][1] - S[row][1][1]*S[row][2][0])/det;

        S_inv[0][1] = -(S[row][0][1]*S[row][2][2] - S[row][0][2]*S[row][2][1])/det;
        S_inv[1][1] =  (S[row][0][0]*S[row][2][2] - S[row][0][2]*S[row][2][0])/det;
        S_inv[2][1] = -(S[row][0][0]*S[row][2][1] - S[row][0][1]*S[row][2][0])/det;

        S_inv[0][2] =  (S[row][0][1]*S[row][1][2] - S[row][0][2]*S[row][1][1])/det;
        S_inv[1][2] = -(S[row][0][0]*S[row][1][2] - S[row][0][2]*S[row][1][0])/det;
        S_inv[2][2] =  (S[row][0][0]*S[row][1][1] - S[row][0][1]*S[row][1][0])/det;

        //System.out.println("\nS S_inv ," + row + ", " + det);
        //for (int i = 0; i < 3; i++)
        //    System.out.println(S[row][i][0] + ", " + S[row][i][1] + ", " + S[row][i][2] + ", " + S_inv[i][0] + ", " + S_inv[i][1] + ", " + S_inv[i][2]);
    }
}
