
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
    private static double[][] coeff = new double[3][7];     // i = (fx, fy, fz), j = (Cx20, Cx11, Cx02, Cx30, Cx21, Cx12, Cx03)

    protected static void calc_perturb_d_e()
    {
        // calculate linear and cubic response, using a d.e. form of perturbation theory in one cycle
        // assume that the output has already been made uniform, by running 'fit_linear_response()'

        if (Chua_y_vs_x.first_order_hdr == null)
        {
            System.out.println("Bad data in 'Chua_Perturb_d_e.calc_perturn_d_e()' : 'first_order_hdr' is not initialized");
            return;
        }

        Main.skew_transform = true;                 // make the linear response "uniform"

        System.out.println("\nPython output - Chua_Perturb_d_e");
        System.out.println("cubic_hdr = \"\\n\\");
        System.out.println("Chua - N-S - perturbation theory (d.e./quadratic/cubic)\\n\\");
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
                perturb_time_sync_d_e(i, j);
        System.out.println("])");

        gen_perturb_response();                   // calculate analytical partial derivatives
    }

    private static void perturb_time_sync_d_e(int ix, int iy)
    {
        Point2D.Double pt2;
        double zp;
        double[] pt12 = new double[] {Main.final_x, Main.final_y, Main.final_z,
                                      Main.invert_from_xp_yp(ix, iy, 0, "x"), Main.invert_from_xp_yp(ix, iy, 0, "y"), Main.invert_from_xp_yp(ix, iy, 0, "z"),
                                      0, 0, 0, 0, 0, 0};

        //if (ix ==1 || iy == 1)                                  // TEST code
        //{
        //    System.out.println("\nk , " + ix + ", " + iy);
        //    System.out.println("0, " + 0 + ", " + 0);           // start pt in (x', y') plane
        //}
        for (int k = 1; k <= Main.final_Period; k++)            // loop through one cycle
        {
            Main.runge_kutta_chua12(pt12, Main.final_delt);  // KEEP code
            //if (k % 10 == 0 && (ix ==1 || iy == 1))             // TEST code
            //{
            //    pt2 = Main.project_2D(pt9[6], pt9[7], pt9[8]);  // cubic response only
            //    System.out.println(k + ", " + pt2.x + ", " + pt2.y);
            //}
        }
        pt2 = Main.project_2D(pt12[3] + pt12[6] + pt12[9], pt12[4] + pt12[7] + pt12[10], pt12[5] + pt12[8] + pt12[11]);   // linear/quadratic/cubic response wrt du'
        zp = Main.project_zp(pt12[3] + pt12[6] + pt12[9], pt12[4] + pt12[7] + pt12[10], pt12[5] + pt12[8] + pt12[11]);
        System.out.print("[" + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp + "]");
        if (ix < Nrect || iy < Nrect)
            System.out.println(",");

        // calculate t_sync version of quadratic coeff (Cx20, Cx11, Cx02)

        pt2 = Main.project_2D(pt12[6], pt12[7], pt12[8]);   // quadratic response only
        zp  = Main.project_zp(pt12[6], pt12[7], pt12[8]);
        if (ix == 1 && iy == 0)         // Cx20
        {
            coeff[0][0] = pt2.x;
            coeff[1][0] = pt2.y;
            coeff[2][0] = zp;
        }
        if (ix == 0 && iy == 1)         // Cx02
        {
            coeff[0][2] = pt2.x;
            coeff[1][2] = pt2.y;
            coeff[2][2] = zp;
        }
        if (ix == 1 && iy == 1)         // Cx11 (this must be chronologically the last to execute)
        {
            coeff[0][1] = pt2.x - coeff[0][0] - coeff[0][2];
            coeff[1][1] = pt2.y - coeff[1][0] - coeff[1][2];
            coeff[2][1] = zp - coeff[2][0] - coeff[2][2];
        }

        // calculate t_sync version of cubic coeff (Cx30, Cx21, Cx12, Cx03) (see Book IV, p.41)

        pt2 = Main.project_2D(pt12[9], pt12[10], pt12[11]);  // cubic-response only
        zp  = Main.project_zp(pt12[9], pt12[10], pt12[11]);
        if (ix == 1 && iy == 0)         // Cx30
        {
            coeff[0][3] = pt2.x;
            coeff[1][3] = pt2.y;
            coeff[2][3] = zp;
            //System.out.println("project, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
        }
        if (ix == 0 && iy == 1)         // Cx03
        {
            coeff[0][6] = pt2.x;
            coeff[1][6] = pt2.y;
            coeff[2][6] = zp;
            //System.out.println("project, " + ix + ", " + iy + ", " + pt2.x + ", " + pt2.y + ", " + zp);
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

    private static void gen_perturb_response()
    {
        // assume we have calculated second-order (t-sync) response, based on linearized d.e., during one cycle
        // then define the partial derivatives of the output (x, y, z) response wrt input (x', y')
        // collect terms to produce z-sync response

        double[] pt6 = new double[6];
        double[] Pdot = new double[] {Main.calc_xdot(Main.final_x, Main.final_y, Main.final_z), Main.calc_ydot(Main.final_x, Main.final_y, Main.final_z), Main.calc_zdot(Main.final_x, Main.final_y, Main.final_z)};
        double[] P2dot = new double[] {Main.calc_x2dot(Main.final_x, Main.final_y, Main.final_z), Main.calc_y2dot(Main.final_x, Main.final_y, Main.final_z), Main.calc_z2dot(Main.final_x, Main.final_y, Main.final_z)};
        double[] P3dot = new double[] {Main.calc_x3dot(Main.final_x, Main.final_y, Main.final_z), Main.calc_y3dot(Main.final_x, Main.final_y, Main.final_z), Main.calc_z3dot(Main.final_x, Main.final_y, Main.final_z)};

        double[][] dPdi = new double[3][2];
        for (int i = 0; i < 2; i++)         // initiallize dx/du, dy/du, or dz/du
        {
            if (i == 0)                                                 // initial (dxdu', dydu', dzdu')
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, Main.invert_from_xp_yp(1, 0, 0, "x"), Main.invert_from_xp_yp(1, 0, 0, "y"), Main.invert_from_xp_yp(1, 0, 0, "z")};
            else if (i == 1)
                pt6 = new double[] {Main.final_x, Main.final_y, Main.final_z, Main.invert_from_xp_yp(0, 1, 0, "x"), Main.invert_from_xp_yp(0, 1, 0, "y"), Main.invert_from_xp_yp(0, 1, 0, "z")};
//            for (int k = 1; k <= 300; k++)                // test code ONLY !!!!!!!
            for (int k = 1; k <= Main.final_Period; k++)                // loop through one cycle KEEP !!!
                Main.runge_kutta_chua6_ddu3(pt6, Main.final_delt);
            dPdi[0][i] = pt6[3];                                        // dP/di
            dPdi[1][i] = pt6[4];                                        // used to calculate dPdot/di
            dPdi[2][i] = pt6[5];                                        // row = (x,y,z), col = (i,j)
        }
        double[][] dPdotdi = new double[][] {{ Main.alpha*(dPdi[1][0] - 3*Main.a*Main.final_x*Main.final_x*dPdi[0][0] - Main.c*dPdi[0][0]),
                                               Main.alpha*(dPdi[1][1] - 3*Main.a*Main.final_x*Main.final_x*dPdi[0][1] - Main.c*dPdi[0][1])},   // dPdot/di
                                             { dPdi[0][0] - dPdi[1][0] + dPdi[2][0], dPdi[0][1] - dPdi[1][1] + dPdi[2][1]},
                                             {-Main.beta*dPdi[1][0] - Main.gamma*dPdi[2][0], -Main.beta*dPdi[1][1] - Main.gamma*dPdi[2][1]}};
        double[][] dP2dotdi = new double[][] {{ Main.alpha*(dPdotdi[1][0] - 3*Main.a*Main.final_x*Main.final_x*dPdotdi[0][0] - Main.c*dPdotdi[0][0] - 6*Main.a*Main.final_x*Pdot[0]*dPdi[0][0]),
                                                Main.alpha*(dPdotdi[1][1] - 3*Main.a*Main.final_x*Main.final_x*dPdotdi[0][1] - Main.c*dPdotdi[0][1] - 6*Main.a*Main.final_x*Pdot[0]*dPdi[0][1])},   // dP2dot/di
                                              { dPdotdi[0][0] - dPdotdi[1][0] + dPdotdi[2][0], dPdotdi[0][1] - dPdotdi[1][1] + dPdotdi[2][1]},
                                              {-Main.beta*dPdotdi[1][0] - Main.gamma*dPdotdi[2][0], -Main.beta*dPdotdi[1][1] - Main.gamma*dPdotdi[2][1]}};
        double[][] d2Pdidj = new double[3][3];                                  // row = (x,y,z), col = (ii,ij,jj)
        double[][] d2Pdotdidj = new double[3][3];                               // row = (x,y,z), col = (ii,ij,jj)
        double PdotPdot = Pdot[0]*Pdot[0] + Pdot[1]*Pdot[1] + Pdot[2]*Pdot[2];
        double PdotP2dot = Pdot[0]*P2dot[0] + Pdot[1]*P2dot[1] + Pdot[2]*P2dot[2];
        double[] dalphadi = new double[] {-(Pdot[0]*dPdi[0][0] + Pdot[1]*dPdi[1][0] + Pdot[2]*dPdi[2][0])/PdotPdot,
                                          -(Pdot[0]*dPdi[0][1] + Pdot[1]*dPdi[1][1] + Pdot[2]*dPdi[2][1])/PdotPdot};
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
            d2Pdotdidj[0][i] =  Main.alpha*(d2Pdidj[1][i] - 3*Main.a*Main.final_x*Main.final_x*d2Pdidj[0][i] - Main.c*d2Pdidj[0][i]);
            d2Pdotdidj[1][i] =  d2Pdidj[0][i] - d2Pdidj[1][i] + d2Pdidj[2][i];                      // Pdot: time derivative
            d2Pdotdidj[2][i] = -Main.beta*d2Pdidj[1][i] - Main.gamma*d2Pdidj[2][i];                 // see Chaos IV, p. 50
        }
        d2Pdotdidj[0][0] += - 3*Main.alpha*Main.a*Main.final_x*dPdi[0][0]*dPdi[0][0];               // add nonlinear M(S, S, S) term
        d2Pdotdidj[0][1] += - 6*Main.alpha*Main.a*Main.final_x*dPdi[0][0]*dPdi[0][1];               // Note that this is a
        d2Pdotdidj[0][2] += - 3*Main.alpha*Main.a*Main.final_x*dPdi[0][1]*dPdi[0][1];               // non-Taylor expansion
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + d2Pdidj[i][0] + ", " + d2Pdidj[i][1] + ", " + d2Pdidj[i][2]);
        System.out.println("\nd2Pdotdidj :");          // convert quadratic coeff from primed to org units
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + d2Pdotdidj[i][0] + ", " + d2Pdotdidj[i][1] + ", " + d2Pdotdidj[i][2]);
        System.out.println("\ndalphadi :");
        System.out.println(dalphadi[0] + ", " + dalphadi[1]);
        //System.out.println(dalphadi[0]/Main.final_delt + ", " + dalphadi[1]/Main.final_delt);
        System.out.println("\nd2alphadidj :");
        d2alphadidj[0] = - 0.5*(2.0*(Pdot[0]*d2Pdidj[0][0] + Pdot[1]*d2Pdidj[1][0] + Pdot[2]*d2Pdidj[2][0])
                         + 2*dalphadi[0]*(Pdot[0]*dPdotdi[0][0] + Pdot[1]*dPdotdi[1][0] + Pdot[2]*dPdotdi[2][0])
                         + dalphadi[0]*dalphadi[0]*PdotP2dot)/PdotPdot;
        d2alphadidj[1] = - ((Pdot[0]*d2Pdidj[0][1] + Pdot[1]*d2Pdidj[1][1] + Pdot[2]*d2Pdidj[2][1])
                         + dalphadi[0]*(Pdot[0]*dPdotdi[0][1] + Pdot[1]*dPdotdi[1][1] + Pdot[2]*dPdotdi[2][1])
                         + dalphadi[1]*(Pdot[0]*dPdotdi[0][0] + Pdot[1]*dPdotdi[1][0] + Pdot[2]*dPdotdi[2][0])
                         + dalphadi[0]*dalphadi[1]*PdotP2dot)/PdotPdot;
        d2alphadidj[2] = - 0.5*(2.0*(Pdot[0]*d2Pdidj[0][2] + Pdot[1]*d2Pdidj[1][2] + Pdot[2]*d2Pdidj[2][2])
                         + 2*dalphadi[1]*(Pdot[0]*dPdotdi[0][1] + Pdot[1]*dPdotdi[1][1] + Pdot[2]*dPdotdi[2][1])
                         + dalphadi[1]*dalphadi[1]*PdotP2dot)/PdotPdot;
        System.out.println(d2alphadidj[0] + ", " + d2alphadidj[1] + ", " + d2alphadidj[2]);
        //System.out.println(d2alphadidj[0]/Main.final_delt + ", " + d2alphadidj[1]/Main.final_delt + ", " + d2alphadidj[2]/Main.final_delt);

        System.out.println("\ndPdotdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dPdotdi[i][0] + ", " + dPdotdi[i][1]);

        System.out.println("\ndP2dotdi :");
        for (i = 0; i < 3; i++)
            System.out.println(i + ", " + dP2dotdi[i][0] + ", " + dP2dotdi[i][1]);

        Point2D.Double pt2eig = Main.project_2D(dPdi[0][0], dPdi[1][0], dPdi[2][0]);
        System.out.println("\nCxy cubic coeff (Type t-sync):");
        for (i = 0; i < 3; i++)                     // proportional to P
            System.out.println(i + ", " + coeff[i][0] + ", " + coeff[i][1] + ", " + coeff[i][2] + ", " + coeff[i][3] + ", " + coeff[i][4] + ", " + coeff[i][5] + ", " + coeff[i][6]);
        System.out.println("\nCxy coeff (t-sync) (intermediate):");
        System.out.println(Chua_y_vs_x.first_order_hdr.split(",")[0] + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Main.final_Period + ", " + Main.final_delt
                       + ", 0, " + pt2eig.x + ", 0, " + coeff[0][0] + ", " + coeff[0][1] + ", " + coeff[0][2] + ", " + coeff[0][3] + ", " + coeff[0][4] + ", " + coeff[0][5] + ", " + coeff[0][6]
                       + ", 0, " + pt2eig.y + ", 0, " + coeff[1][0] + ", " + coeff[1][1] + ", " + coeff[1][2] + ", " + coeff[1][3] + ", " + coeff[1][4] + ", " + coeff[1][5] + ", " + coeff[1][6]);

//      z-sync calc

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

        System.out.println("\nCxy coeff (z-sync) (final):");
        System.out.print(Chua_y_vs_x.first_order_hdr.split(",")[0] + ", " + Main.alpha + ", " + Main.beta + ", " + Main.gamma + ", " + Main.a + ", " + Main.c + ", " + Main.final_Period + ", " + Main.final_delt);
        System.out.print(", 0, " + pt2eig.x + ", 0");
        Point2D.Double pt2temp;
        for (i = 0; i < 7; i++)
        {
            pt2temp = Main.project_2D(collect_Palldot[0][i], collect_Palldot[1][i], collect_Palldot[2][i]);
            System.out.print(", " + (coeff[0][i] + pt2temp.x));
        }
        System.out.print(", 0, " + pt2eig.y + ", 0");
        for (i = 0; i < 7; i++)
        {
            pt2temp = Main.project_2D(collect_Palldot[0][i], collect_Palldot[1][i], collect_Palldot[2][i]);
            System.out.print(", " + (coeff[1][i] + pt2temp.y));
        }
        System.out.println();
    }
}
