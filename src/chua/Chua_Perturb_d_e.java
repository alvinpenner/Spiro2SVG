
// calculate linear and cubic coefficients of a Chua Oscillator torus
// assuming that the linear response has been obtained from 'fit_linear_response()'
// use perturbation theory cast into the form of a d.e., as per my book Chaos IV, p. 65

package chua;

// this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Perturb_d_e.java

import java.awt.geom.Point2D;

public final class Chua_Perturb_d_e
{
    private static final int Nrect = 2;           // # of coordinates > 0 (for rectangular grid)
    private static double fx12 = 9999, fy12 = 9999, fz12 = 9999;
    private static double[][][] S = new double[Main.final_Period + 1][3][3];
    private static double[][] coeff = new double[3][7];     // i = (fx, fy, fz), j = (Cx20, Cx11, Cx02, Cx30, Cx21, Cx12, Cx03)

    protected static void calc_coeff_d_e()
    {
        // calculate linear and cubic response, using a d.e. form of perturbation theory in one cycle
        // assume that the output has already been made uniform, by running 'fit_linear_response()'

        if (Chua_y_vs_x.first_order_hdr == null)
        {
            System.out.println("Bad data in 'Perturb_d_e.calc_coeff()' : 'first_order_hdr' is not initialized");
            return;
        }

        Main.skew_transform = true;                 // make the linear response "uniform"

        System.out.println("\nPython output - Chua_Perturb_d_e");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Chua - Neimark-Sacker - S-matrix perturbation theory (d.e.)\\n\\");
        System.out.println("incr_iT_eig_angle, 1.0, " + Chua_y_vs_x.first_order_hdr + ", \\n\\");
        System.out.println("alpha_beta_gamma , " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ",\\n\\");
        System.out.println("Period_delt      , " + Main.final_Period + ", " + Main.final_delt + ",\\n\\");
        System.out.println("x_y_z            , " + Main.final_x + ", " + Main.final_y + ", " + Main.final_z + ",\\n\\");
        System.out.println("phi_theta_psi    , " + Main.project_phi + ", " + Main.project_theta + ", " + Main.project_psi + ",\\n\\");
        System.out.println("Re_V21_Im_V21    , " + Main.final_Re_V21 + ", " + Main.final_Im_V21 + ",\\n\\");
        System.out.println("\\n\\");
        System.out.println("i,j,x',y',tc\"");

        System.out.print("data = np.array([");
        for (int i = -Nrect; i < Nrect + 1; i++)                            // increment x' by i
            for (int j = -Nrect; j < Nrect + 1; j++)                        // increment y' by j
                gen_time_sync_d_e(i, j);
        System.out.println("])");

        gen_cubic_response();                   // calculate analytical partial derivatives
    }

    private static void gen_time_sync_d_e(int ix, int iy)
    {
        Point2D.Double pt2;
        double zp;
        double[] pt9 = new double[] {Main.final_x, Main.final_y, Main.final_z,
                                     Main.invert_from_xp_yp(ix, iy, 0, "x"), Main.invert_from_xp_yp(ix, iy, 0, "y"), Main.invert_from_xp_yp(ix, iy, 0, "z"),
                                     Main.invert_from_xp_yp(ix, iy, 0, "x"), Main.invert_from_xp_yp(ix, iy, 0, "y"), Main.invert_from_xp_yp(ix, iy, 0, "z")};
//                                     0, 0, 0};          // fix fix bug bug

        for (int k = 1; k <= Main.final_Period; k++)                     // loop through one cycle
        {
            Main.runge_kutta_chua9_ddu6(pt9, Main.final_delt);
        }
        pt2 = Main.project_2D(pt9[3], pt9[4], pt9[5]);   // final (dxdu', dydu', dzdu')
        System.out.print("[" + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + Main.project_zp(pt9[3], pt9[4], pt9[5]) + "]");
        if (ix < Nrect || iy < Nrect)
            System.out.println(",");

        // calculate t_sync version of quadratic coeff (Cx20, Cx11, Cx02) DISABLED !!!!!!!

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                coeff[i][j] = 0;

        // calculate t_sync version of cubic coeff (Cx30, Cx21, Cx12, Cx03) (see Book IV, p.41)
        // use second-order perturbation theory, which generates a cubic function

        //if (ix == 1 && iy == 0 || ix == 0 && iy == 1)
        //    System.out.println("org    , " + ix + ", " + iy + ", " + pt9[0] + ", " + pt9[1] + ", " + pt9[2]
        //                                                    + ", " + pt9[3] + ", " + pt9[4] + ", " + pt9[5]
        //                                                    + ", " + pt9[6] + ", " + pt9[7] + ", " + pt9[8]);
        pt2 = Main.project_2D(pt9[6], pt9[7], pt9[8]);
        zp  = Main.project_zp(pt9[6], pt9[7], pt9[8]);
        if (ix == 1 && iy == 0)         // Cx30
        {
            coeff[0][3] = pt2.x;
            coeff[1][3] = pt2.y;
            coeff[2][3] = zp;
            System.out.println("project, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
        }
        if (ix == 0 && iy == 1)         // Cx03
        {
            coeff[0][6] = pt2.x;
            coeff[1][6] = pt2.y;
            coeff[2][6] = zp;
            System.out.println("project, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
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
            //invert_S(k);                                        // calculate S_inverse
            for (i = 0; i < 3; i++)
            {
                y_incr[i] = 0;
                //for (int j = 0; j < 3; j++)
                //    y_incr[i] += S_inv[i][j]*M_vec[j];
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
            //invert_S(k);                                            // calculate S_inverse
            for (i = 0; i < 3; i++)
            {
                y_incr[i] = 0;
                //for (int j = 0; j < 3; j++)
                //    y_incr[i] += S_inv[i][j]*M_vec[j];
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
            //fx12 = pt2.x;               // fx12 (temporary storage)
            //fy12 = pt2.y;
            //fz12 = zp;
        }
        if (ix == 2 && iy == 1)         // Cx21 (this must be chronologically the last to execute)
        {
            //coeff[0][4] = (-15*coeff[0][3] + 2*pt2.x - 1*fx12 + 6*coeff[0][6])/6;   // Cx21
            //coeff[1][4] = (-15*coeff[1][3] + 2*pt2.y - 1*fy12 + 6*coeff[1][6])/6;
            //coeff[2][4] = (-15*coeff[2][3] + 2*zp - 1*fz12 + 6*coeff[2][6])/6;
            //coeff[0][5] = (  6*coeff[0][3] - 1*pt2.x + 2*fx12 - 15*coeff[0][6])/6;  // Cx12
            //coeff[1][5] = (  6*coeff[1][3] - 1*pt2.y + 2*fy12 - 15*coeff[1][6])/6;
            //coeff[2][5] = (  6*coeff[2][3] - 1*zp + 2*fz12 - 15*coeff[2][6])/6;
        }
    }

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
}
