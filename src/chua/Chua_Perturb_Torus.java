
// calculate quadratic and cubic coefficients of a Chua Oscillator torus
// assuming that the linear response has been obtained from 'fit_linear_response()'
// use perturbation theory as per Iooss, Bouc, and my book Chaos IV, p. 27

package chua;

// this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Perturb_Torus.java

import java.awt.geom.Point2D;

public final class Chua_Perturb_Torus
{
    private static double[][][] S = new double[Main.final_Period + 1][3][3];
    private static double[][] xyz = new double[Main.final_Period + 1][3];       // limit cycle
    private static double[][] S_inv = new double[3][3];
    private static double[][] coeff = new double[3][7];     // i = (fx, fy, fz), j = (Cx20, Cx11, Cx02, Cx30, Cx21, Cx12, Cx03)
    private static double fx12 = 9999, fy12 = 9999, fz12 = 9999;
    private static int Nrect = 2;                           // # of coordinates > 0 (for rectangular grid)

    protected static void calc_coeff()
    {
        // calculate first-order response, S matrix, during one cycle
        // assume that the output has already been made uniform, by running 'fit_linear_response()'

        if (Chua_y_vs_x.first_order_hdr == null)
        {
            System.out.println("Bad data in 'Perturb_Torus.calc_coeff()' : 'first_order_hdr' is not initialized");
            return;
        }

        // perform one uniformized pass to save S

        double[] pt6 = new double[6];
        Main.skew_transform = true;                 // make the linear response "uniform"
        double xdot = Main.calc_xdot(Main.final_x, Main.final_y, Main.final_z);
        double ydot = Main.calc_ydot(Main.final_x, Main.final_y, Main.final_z);
        double zdot = Main.calc_zdot(Main.final_x, Main.final_y, Main.final_z);
        double v = Math.sqrt(xdot*xdot + ydot*ydot + zdot*zdot);
        //System.out.println("vel , " + xdot + ", " + ydot + ", " + zdot + ", " + v);

        System.out.println("\nPython output - Chua S matrix");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Chua - Neimark-Sacker - fit S-matrix perturbation theory response\\n\\");
        System.out.println("incr_iT_eig_angle, 1.0, " + Chua_y_vs_x.first_order_hdr + ", \\n\\");
        System.out.println("alpha_beta_gamma , " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt      , " + Main.final_Period + ", " + Main.final_delt + ",\\n\\");
        System.out.println("x_y_z            , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi    , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("Re_V21_Im_V21    , " + Main.final_Re_V21 + ", " + Main.final_Im_V21 + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',tc\"");
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
                Main.runge_kutta_chua6_ddu3(pt6, Main.final_delt);
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
        //System.out.println("\nlimit cycle xyz det S");
        //for (int k = 0; k < S.length; k++)
        //    System.out.println(k + ", " + xyz[k][0] + ", " + xyz[k][1] + ", " + xyz[k][2] + ", " + det_S(k));

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

        //int ttemp = Main.final_Period;                                  // temporary time index
        //System.out.println("\nS matrix (org) @ time = " + ttemp);
        //for (int i = 0; i < 3; i++)
        //    System.out.println(S[ttemp][i][0] + ", " + S[ttemp][i][1] + ", " + S[ttemp][i][2]);

        //System.out.println("\nS matrix (projected) @ time = " + ttemp);
        //for (int i = 0; i < 3; i++)
        //{
        //    Point2D.Double pt2 = Main.project_2D(S[ttemp][0][i], S[ttemp][1][i], S[ttemp][2][i]);
        //    System.out.println(pt2.x + ", " + pt2.y + ", " + Main.project_zp(S[ttemp][0][i], S[ttemp][1][i], S[ttemp][2][i]));
        //}
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
        for (int i = -Nrect; i < Nrect + 1; i++)                            // increment x' by i
            for (int j = -Nrect; j < Nrect + 1; j++)                        // increment y' by j
                gen_time_sync_Cxy(i, j);
        System.out.println("])");

        //gen_TOF_extrapolate_Cxy();              // enhance the t_sync result (numerical curve-fit) (NOT necessary)
        gen_cubic_response();                   // calculate analytical partial derivatives
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
                              + ", " + (S[tk][0][istate] + Main.alpha*S[tk][1][istate])
                              + ", " + (xyz[tk][2]*S[tk][0][istate] + (xyz[tk][0] - Main.gamma)*S[tk][2][istate]));
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
            state_vec[1] =  S[k][0][istate] + Main.alpha*S[k][1][istate];
            state_vec[2] =  xyz[k][2]*S[k][0][istate] + (xyz[k][0] - Main.gamma)*S[k][2][istate];
            if (k == 0 || k == Main.final_Period)
                for (int i = 0; i < 3; i++)
                    state_vec[i] /= 2.0;                    // correct the endpoints
            for (int i = 0; i < 3; i++)
                accum_vec[i] += state_vec[i]*Main.final_delt;
        }
        pt2 = Main.project_2D(accum_vec[0], accum_vec[1], accum_vec[2]);   // final (dxdu', dydu', dzdu')
        System.out.println("end  , " + istate + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(accum_vec[0], accum_vec[1], accum_vec[2]));
    }

    private static void partial_integrate(int ix, int iy, int tk)   // test code only
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
                                           + ", " + (state_vec[0] + Main.alpha*state_vec[1])
                                           + ", " + (xyz[tk][2]*state_vec[0] + (xyz[tk][0] - Main.gamma)*state_vec[2]));
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

        // second-order perturbation (cubic coeff)

        for (int k = 0; k < S.length; k++)                      // scan one limit cycle
        {
            // linear combination of dx' and dy' response (col 0 and col 1)

            for (i = 0; i < 3; i++)
                state_y1[i] = ix*S[k][i][0] + iy*S[k][i][1];    // perturbed state (first-order)

            // form non-linear response M = (-alpha*a*x*x*x, 0, 0)

            M_vec[0] = -Main.alpha*Main.a*state_y1[0]*state_y1[0]*state_y1[0];
            M_vec[1] = 0;
            M_vec[2] = 0;
            invert_S(k);                                        // calculate S_inverse
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
            if (ix == 0 && iy == 1 && !true)
            {
                //pt2 = Main.project_2D(ix*S[k][0][0] + iy*S[k][0][1], ix*S[k][1][0] + iy*S[k][1][1], ix*S[k][2][0] + iy*S[k][2][1]);
                pt2 = Main.project_2D(state_y2[k][0], state_y2[k][1], state_y2[k][2]);
                System.out.println("state_y2 , " + ix + ", " + iy + ", " + k + ", " + pt2.x + ", " + pt2.y);
            }
        }

        // cubic component

        for (int k = 0; k < S.length; k++)                          // scan one limit cycle
        {
            for (i = 0; i < 3; i++)
                state_y1[i] = ix*S[k][i][0] + iy*S[k][i][1];    // perturbed state (first-order)

            // form non-linear response M = (0, 0, dz*dx)

            M_vec[0] = -Main.alpha*Main.a*(state_y1[0] + state_y2[k][0])*(state_y1[0] + state_y2[k][0])*(state_y1[0] + state_y2[k][0]);
            M_vec[1] = 0;
            M_vec[2] = 0;
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
        //for (i = 0; i < 3; i++)
        //    System.out.println("temporary out = " + i + ", " + ix + ", " + iy + ", " + state_y1[i] + ", " + state_y2[Main.final_Period][i] + ", " + state_y3[i]);

        // summary

        for (i = 0; i < 3; i++)
            //state_final[i] = ix*S[Main.final_Period][i][0] + iy*S[Main.final_Period][i][1] + state_y3[i];   // test CODE
            state_final[i] = ix*S[Main.final_Period][i][0] + iy*S[Main.final_Period][i][1] + state_y2[Main.final_Period][i];
            //state_final[i] = ix*S[Main.final_Period][i][0] + iy*S[Main.final_Period][i][1] + state_y2[Main.final_Period][i] + state_y3[i];
        pt2 = Main.project_2D(state_final[0], state_final[1], state_final[2]);   // final (dxdu', dydu', dzdu')
        zp  = Main.project_zp(state_final[0], state_final[1], state_final[2]);
        System.out.print("[" + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp + "]");
        if (ix < Nrect || iy < Nrect)
            System.out.println(",");

        // calculate t_sync version of quadratic coeff (Cx20, Cx11, Cx02) DISABLED !!!!!!!

        pt2 = Main.project_2D(state_y2[Main.final_Period][0], state_y2[Main.final_Period][1], state_y2[Main.final_Period][2]);
        zp  = Main.project_zp(state_y2[Main.final_Period][0], state_y2[Main.final_Period][1], state_y2[Main.final_Period][2]);
        if (ix == 1 && iy == 0)         // Cx20
        {
            //System.out.println("unprojected coeff (col 0) =, " + state_y2[Main.final_Period][0] + ", " + state_y2[Main.final_Period][1] + ", " + state_y2[Main.final_Period][2]);
            coeff[0][0] = 0*pt2.x;
            coeff[1][0] = 0*pt2.y;
            coeff[2][0] = 0*zp;
        }
        if (ix == 0 && iy == 1)         // Cx02
        {
            //System.out.println("unprojected coeff (col 2) =, " + state_y2[Main.final_Period][0] + ", " + state_y2[Main.final_Period][1] + ", " + state_y2[Main.final_Period][2]);
            coeff[0][2] = 0*pt2.x;
            coeff[1][2] = 0*pt2.y;
            coeff[2][2] = 0*zp;
        }
        if (ix == 1 && iy == 1)         // Cx11 (this must be chronologically the last to execute)
        {
            //System.out.println("unprojected coeff (col 1 RAW) =, " + state_y2[Main.final_Period][0] + ", " + state_y2[Main.final_Period][1] + ", " + state_y2[Main.final_Period][2]);
            coeff[0][1] = 0; // pt2.x - coeff[0][0] - coeff[0][2];
            coeff[1][1] = 0; // pt2.y - coeff[1][0] - coeff[1][2];
            coeff[2][1] = 0; // zp - coeff[2][0] - coeff[2][2];
        }

        // calculate t_sync version of cubic coeff (Cx30, Cx21, Cx12, Cx03) (see Book IV, p.41)

        // use second-order perturbation theory, which generates a cubic function
        pt2 = Main.project_2D(state_y2[Main.final_Period][0], state_y2[Main.final_Period][1], state_y2[Main.final_Period][2]);
        zp  = Main.project_zp(state_y2[Main.final_Period][0], state_y2[Main.final_Period][1], state_y2[Main.final_Period][2]);
        //pt2 = Main.project_2D(state_y3[0], state_y3[1], state_y3[2]);
        //zp  = Main.project_zp(state_y3[0], state_y3[1], state_y3[2]);
        if (ix == 1 && iy == 0)         // Cx30
        {
            coeff[0][3] = pt2.x;
            coeff[1][3] = pt2.y;
            coeff[2][3] = zp;
            //System.out.println("raw, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
        }
        if (ix == 0 && iy == 1)         // Cx03
        {
            coeff[0][6] = pt2.x;
            coeff[1][6] = pt2.y;
            coeff[2][6] = zp;
            //System.out.println("raw, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
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
        //if (ix == -1 && iy == 1)
        //    System.out.println("raw, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
        //if (ix ==  1 && iy == 1)
        //    System.out.println("raw, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
        //if (ix ==  1 && iy == -1)
        //    System.out.println("raw, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
    }
/*
    private static void gen_TOF_extrapolate_Cxy()
    {
        // assume we have calculated second-order response, based on S matrix, during one cycle
        // then enhance the theory with a quadratic extrapolation
        // from t_sync to z_sync
        // This is NOT necessary, theoretically, just an intermediate double-check

        double incr = 0.001; // 0.001;
        double x0dot = -Main.final_y - Main.final_z;
        double y0dot =  Main.final_x + Main.a*Main.final_y;
        double z0dot =  Main.b + Main.final_z*(Main.final_x - Main.c);
        Point2D.Double pt2new;          // projected, new, (x', y') using Euler angles
        double zpnew;
        //pt2new = Main.project_2D(x0dot, y0dot, z0dot);
        //zpnew = Main.project_zp(x0dot, y0dot, z0dot);
        //System.out.print("test velocity projection ," + pt2new.x + ", " + pt2new.y + ", " + zpnew + ", " + Math.sqrt(x0dot*x0dot + y0dot*y0dot + z0dot*z0dot));

        System.out.println("\nPython output - Chua S matrix (enhanced)");
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
*/
    private static void gen_cubic_response()
    {
        // assume we have calculated second-order (t-sync) response, based on S matrix, during one cycle
        // then define the partial derivatives of the output (x, y, z) response wrt input (x', y')

        double x0dot = Main.calc_xdot(Main.final_x, Main.final_y, Main.final_z);
        double y0dot = Main.calc_ydot(Main.final_x, Main.final_y, Main.final_z);
        double z0dot = Main.calc_zdot(Main.final_x, Main.final_y, Main.final_z);
        double[] P2dot = new double[] {-y0dot - z0dot, x0dot + Main.alpha*y0dot, z0dot*(Main.final_x - Main.gamma) + Main.final_z*x0dot};
        double[] P3dot = new double[] {-P2dot[1] - P2dot[2], P2dot[0] + Main.alpha*P2dot[1], P2dot[2]*(Main.final_x - Main.gamma) + 2*z0dot*x0dot + Main.final_z*P2dot[0]};
        double[][] dPdi = new double[][] {{S[Main.final_Period][0][0], S[Main.final_Period][0][1]},                     // dP/di (redundant)
                                          {S[Main.final_Period][1][0], S[Main.final_Period][1][1]},                     // used to calculate dPdot/di
                                          {S[Main.final_Period][2][0], S[Main.final_Period][2][1]}};                    // row = (x,y,z), col = (i,j)
        double[][] dPdotdi = new double[][] {{-dPdi[1][0] - dPdi[2][0]       , -dPdi[1][1] - dPdi[2][1]},               // dPdot/di
                                             { dPdi[0][0] + Main.alpha*dPdi[1][0],  dPdi[0][1] + Main.alpha*dPdi[1][1]},
                                             { dPdi[2][0]*(Main.final_x - Main.gamma) + Main.final_z*dPdi[0][0], dPdi[2][1]*(Main.final_x - Main.gamma) + Main.final_z*dPdi[0][1]}};
        double[][] dP2dotdi = new double[][] {{-dPdotdi[1][0] - dPdotdi[2][0]       , -dPdotdi[1][1] - dPdotdi[2][1]},  // dP2dot/di
                                              { dPdotdi[0][0] + Main.alpha*dPdotdi[1][0],  dPdotdi[0][1] + Main.alpha*dPdotdi[1][1]},
                                              { dPdotdi[2][0]*(Main.final_x - Main.gamma) + z0dot*dPdi[0][0] + Main.final_z*dPdotdi[0][0] + dPdi[2][0]*x0dot,
                                                dPdotdi[2][1]*(Main.final_x - Main.gamma) + z0dot*dPdi[0][1] + Main.final_z*dPdotdi[0][1] + dPdi[2][1]*x0dot}};
        double[][] d2Pdidj = new double[3][3];                                  // row = (x,y,z), col = (ii,ij,jj)
        double[][] d2Pdotdidj = new double[3][3];                               // row = (x,y,z), col = (ii,ij,jj)
        double PdotPdot = x0dot*x0dot + y0dot*y0dot + z0dot*z0dot;
        double PdotP2dot = x0dot*P2dot[0] + y0dot*P2dot[1] + z0dot*P2dot[2];
        double[] dalphadi = new double[] {-(x0dot*dPdi[0][0] + y0dot*dPdi[1][0] + z0dot*dPdi[2][0])/PdotPdot,
                                          -(x0dot*dPdi[0][1] + y0dot*dPdi[1][1] + z0dot*dPdi[2][1])/PdotPdot};
        double[] d2alphadidj = new double[3];
        double[][] collect_Pdot = new double[3][7];     // i = (fx, fy, fz), j = (Cx20, Cx11, Cx02, Cx30, Cx21, Cx12, Cx03)
        double[][] collect_P2dot = new double[3][7];    // terms proportional to P2dot
        double[][] collect_P3dot = new double[3][7];    // terms proportional to P3dot
        double[][] collect_Palldot = new double[3][7];
        int i;

        // summarize quadratic coeff

        System.out.println("\ndPdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dPdi[i][0] + ", " + dPdi[i][1]);

        System.out.println("\nd2Pdidj :");          // convert quadratic coeff from primed to org units
        for (i = 0; i < 3; i++)
        {
            d2Pdidj[0][i] = Main.invert_from_xp_yp(coeff[0][i], coeff[1][i], coeff[2][i], "x");     // terms proportional to P
            d2Pdidj[1][i] = Main.invert_from_xp_yp(coeff[0][i], coeff[1][i], coeff[2][i], "y");
            d2Pdidj[2][i] = Main.invert_from_xp_yp(coeff[0][i], coeff[1][i], coeff[2][i], "z");
            d2Pdotdidj[0][i] = -d2Pdidj[1][i] - d2Pdidj[2][i];                                      // Pdot: time derivative
            d2Pdotdidj[1][i] =  d2Pdidj[0][i] + Main.alpha*d2Pdidj[1][i];                           // see Chaos IV, p. 40
            d2Pdotdidj[2][i] =  Main.final_z*d2Pdidj[0][i] + (Main.final_x - Main.gamma)*d2Pdidj[2][i];
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
        System.out.println(dalphadi[0] + ", " + dalphadi[1]);
        System.out.println("\nd2alphadidj :");
        d2alphadidj[0] = - 0.5*(2.0*(x0dot*d2Pdidj[0][0] + y0dot*d2Pdidj[1][0] + z0dot*d2Pdidj[2][0])
                         + 2*dalphadi[0]*(x0dot*dPdotdi[0][0] + y0dot*dPdotdi[1][0] + z0dot*dPdotdi[2][0])
                         + dalphadi[0]*dalphadi[0]*PdotP2dot)/PdotPdot;
        d2alphadidj[1] = - ((x0dot*d2Pdidj[0][1] + y0dot*d2Pdidj[1][1] + z0dot*d2Pdidj[2][1])
                         + dalphadi[0]*(x0dot*dPdotdi[0][1] + y0dot*dPdotdi[1][1] + z0dot*dPdotdi[2][1])
                         + dalphadi[1]*(x0dot*dPdotdi[0][0] + y0dot*dPdotdi[1][0] + z0dot*dPdotdi[2][0])
                         + dalphadi[0]*dalphadi[1]*PdotP2dot)/PdotPdot;
        d2alphadidj[2] = - 0.5*(2.0*(x0dot*d2Pdidj[0][2] + y0dot*d2Pdidj[1][2] + z0dot*d2Pdidj[2][2])
                         + 2*dalphadi[1]*(x0dot*dPdotdi[0][1] + y0dot*dPdotdi[1][1] + z0dot*dPdotdi[2][1])
                         + dalphadi[1]*dalphadi[1]*PdotP2dot)/PdotPdot;
        System.out.println(d2alphadidj[0] + ", " + d2alphadidj[1] + ", " + d2alphadidj[2]);

        System.out.println("\ndPdotdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dPdotdi[i][0] + ", " + dPdotdi[i][1]);

        System.out.println("\ndP2dotdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dP2dotdi[i][0] + ", " + dP2dotdi[i][1]);

        System.out.println("\nCxy cubic coeff (Type t-sync):");
        for (i = 0; i < 3; i++)                     // proportional to P
            System.out.println(i + ", " + coeff[i][0] + ", " + coeff[i][1] + ", " + coeff[i][2] + ", " + coeff[i][3] + ", " + coeff[i][4] + ", " + coeff[i][5] + ", " + coeff[i][6]);
        System.out.println("Cxy cubic coeff (t-sync) (summary):");
        System.out.println(Chua_y_vs_x.first_order_hdr.split(",")[0] + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Main.final_Period + ", " + Main.final_delt
                       + ", 0, 0, 0, " + coeff[0][0] + ", " + coeff[0][1] + ", " + coeff[0][2] + ", " + coeff[0][3] + ", " + coeff[0][4] + ", " + coeff[0][5] + ", " + coeff[0][6]
                       + ", 0, 0, 0, " + coeff[1][0] + ", " + coeff[1][1] + ", " + coeff[1][2] + ", " + coeff[1][3] + ", " + coeff[1][4] + ", " + coeff[1][5] + ", " + coeff[1][6]);
/*      z-sync calc (temporarily disabled)
        for (i = 0; i < 3; i++)                                                 // proportional to Pdot
        {
            collect_Pdot[i][0] =                             dalphadi[0]*dPdotdi[i][0];
            collect_Pdot[i][1] = dalphadi[0]*dPdotdi[i][1] + dalphadi[1]*dPdotdi[i][0];
            collect_Pdot[i][2] = dalphadi[1]*dPdotdi[i][1];
            collect_Pdot[i][3] = dalphadi[0]*d2Pdotdidj[i][0]                                                               + d2alphadidj[0]*dPdotdi[i][0];
            collect_Pdot[i][4] = dalphadi[0]*d2Pdotdidj[i][1] + dalphadi[1]*d2Pdotdidj[i][0] + d2alphadidj[0]*dPdotdi[i][1] + d2alphadidj[1]*dPdotdi[i][0];
            collect_Pdot[i][5] = dalphadi[0]*d2Pdotdidj[i][2] + dalphadi[1]*d2Pdotdidj[i][1] + d2alphadidj[1]*dPdotdi[i][1] + d2alphadidj[2]*dPdotdi[i][0];
            collect_Pdot[i][6] =                                dalphadi[1]*d2Pdotdidj[i][2] + d2alphadidj[2]*dPdotdi[i][1];
        }
        for (i = 0; i < 3; i++)                                                 // proportional to P2dot
        {
            collect_P2dot[i][0] = dalphadi[0]*dalphadi[0]*P2dot[i]/2;
            collect_P2dot[i][1] = dalphadi[0]*dalphadi[1]*P2dot[i];
            collect_P2dot[i][2] = dalphadi[1]*dalphadi[1]*P2dot[i]/2;
            collect_P2dot[i][3] =  dalphadi[0]*d2alphadidj[0]*P2dot[i]                               + dalphadi[0]*dalphadi[0]*dP2dotdi[i][0]/2;
            collect_P2dot[i][4] = (dalphadi[0]*d2alphadidj[1] + dalphadi[1]*d2alphadidj[0])*P2dot[i] + dalphadi[0]*dalphadi[0]*dP2dotdi[i][1]/2 + dalphadi[0]*dalphadi[1]*dP2dotdi[i][0];
            collect_P2dot[i][5] = (dalphadi[0]*d2alphadidj[2] + dalphadi[1]*d2alphadidj[1])*P2dot[i] + dalphadi[0]*dalphadi[1]*dP2dotdi[i][1]   + dalphadi[1]*dalphadi[1]*dP2dotdi[i][0]/2;
            collect_P2dot[i][6] =                               dalphadi[1]*d2alphadidj[2]*P2dot[i]                                             + dalphadi[1]*dalphadi[1]*dP2dotdi[i][1]/2;
        }
        for (i = 0; i < 3; i++)                                                 // proportional to P3dot
        {
            collect_P3dot[i][0] = 0;
            collect_P3dot[i][1] = 0;
            collect_P3dot[i][2] = 0;
            collect_P3dot[i][3] = dalphadi[0]*dalphadi[0]*dalphadi[0]*P3dot[i]/6;
            collect_P3dot[i][4] = dalphadi[0]*dalphadi[0]*dalphadi[1]*P3dot[i]/2;
            collect_P3dot[i][5] = dalphadi[0]*dalphadi[1]*dalphadi[1]*P3dot[i]/2;
            collect_P3dot[i][6] = dalphadi[1]*dalphadi[1]*dalphadi[1]*P3dot[i]/6;
        }
        for (i = 0; i < 3; i++)
            for (int j = 0; j < 7; j++)
                collect_Palldot[i][j] = collect_Pdot[i][j] + collect_P2dot[i][j] + collect_P3dot[i][j];
        System.out.println("Cxy coeff (final):");
        Point2D.Double pt2eig = Main.project_2D(S[Main.final_Period][0][0], S[Main.final_Period][1][0], S[Main.final_Period][2][0]);

        System.out.print(Chua_y_vs_x.first_order_hdr.split(",")[0] + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Main.final_Period + ", " + Main.final_delt);
        System.out.print(",0," + pt2eig.x + ",0");
        Point2D.Double pt2temp;
        for (i = 0; i < 7; i++)
        {
            pt2temp = Main.project_2D(collect_Palldot[0][i], collect_Palldot[1][i], collect_Palldot[2][i]);
            System.out.print(", " + (coeff[0][i] + pt2temp.x));
        }
        System.out.print(",0," + pt2eig.y + ",0");
        for (i = 0; i < 7; i++)
        {
            pt2temp = Main.project_2D(collect_Palldot[0][i], collect_Palldot[1][i], collect_Palldot[2][i]);
            System.out.print(", " + (coeff[1][i] + pt2temp.y));
        }
        System.out.println();
 */
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

    private static double det_S(int row)
    {
        return S[row][0][0]*S[row][1][1]*S[row][2][2] + S[row][0][1]*S[row][1][2]*S[row][2][0] + S[row][0][2]*S[row][1][0]*S[row][2][1]
             - S[row][0][2]*S[row][1][1]*S[row][2][0] - S[row][0][0]*S[row][1][2]*S[row][2][1] - S[row][0][1]*S[row][1][0]*S[row][2][2];
    }
}
