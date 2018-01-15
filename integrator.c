/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Mościbrodzka
 *
 */

#include "functions.h"
#include "parameters.h"
#include "constants.h"

// Updates the vector y (containing position/velocity) by one RK4 step.
void rk4_step(real *y, real dt){
        // Array containing all "update elements" (4 times Nelements because RK4)
        real dx[DIM * 2 * 4];
        // Create a copy of the "y vector" that can be shifted for the
        // separate function calls made by RK4
        real yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

        // fvector contains f(yshift), as applied to yshift (the 'current' y
        // during RK steps). It is used to compute the 'k-coefficients' (dx)
        real fvector[DIM * 2];

        // Compute the RK4 update coefficients ('K_n' in lit., 'dx' here)
        int i, q;
        real weights[4] = {0.5, 0.5, 1., 0.}; // Weights used for updating y
        for (q = 0; q < 4; q++) {
                f_geodesic(yshift, fvector); // Apply function f to current y to obtain fvector
                for (i = 0; i < DIM * 2; i++) {
                        dx[q * DIM * 2 + i] = dt * fvector[i]; // Use fvector to fill dx
                        yshift[i] = y[i] + dx[q * DIM * 2 + i] * weights[q]; // Update y
                }
        }

        // Update the y-vector (light ray)
        for (i = 0; i < DIM * 2; i++)
                y[i] = y[i] + 1. / 6. * (dx[0 * DIM * 2 + i] +
                                         dx[1 * DIM * 2 + i] * 2. +
                                         dx[2 * DIM * 2 + i] * 2. +
                                         dx[3 * DIM * 2 + i]);
}

// Updates the vector y (containing position/velocity) by one RK2 step.
void rk2_step(real *y, real dt){
        // Array containing all "update elements" (2 times Nelements because RK2)
        real dx[DIM * 2 * 2];
        // Create a copy of the "y vector" that can be shifted for the
        // separate function calls made by RK2
        real yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

        // fvector contains f(yshift), as applied to yshift (the 'current' y
        // during RK steps). It is used to compute the 'k-coefficients' (dx)
        real fvector[DIM * 2];

        // Compute the RK2 update coefficients ('K_n' in lit., 'dx' here)
        int i, q;
        real weights[2] = {0.5, 0.}; // Weights used for updating y
        for (q = 0; q < 2; q++) {
                //  f(yshift, fvector); // Apply function f to current y to obtain fvector
                f_geodesic(yshift, fvector);
                for (i = 0; i < DIM * 2; i++) {
                        dx[q * DIM * 2 + i] = dt * fvector[i]; // Use fvector to update dx
                        yshift[i] = y[i] + dx[q * DIM * 2 + i] * weights[q]; // Update y
                }
        }

        // Update the y-vector (light ray)
        for (i = 0; i < DIM * 2; i++)
                y[i] = y[i] + dx[1 * DIM * 2 + i];  // y_n+1 = y_n + k2 + O(h^3)
}

// Updates the vector y (containing position/velocity) using 'velocity Verlet'
// Ref: Dolence et al 2009 eqn 14a-14d
void verlet_step(real *y, void (*f)(real*, real*),real dl){
        // Create a copy of the "y vector" that can be shifted for the
        // separate function calls made by RK2
        real yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};
        // fvector contains f(yshift), as applied to yshift (the 'current' y
        // during RK steps). It is used to compute the 'k-coefficients' (dx)
        real fvector[DIM * 2];

        // Temporary acceleration vector
        real A_u_temp[DIM];

        // Step 1: Compute A_u(lambda) (Preparation for Eq 14a)
        //    f(yshift, fvector); // fvector now contains A_u(lambda)
        f_geodesic(yshift, fvector);
        // Step 2: Compute X_u(lambda + dlambda) and the temporary four-velocity
        // (Eq 14a, 14b)
        int i;
        LOOP_i {
                yshift[i] += dl * yshift[i + DIM] + 0.5 * dl * dl * fvector[i + DIM];
                yshift[i + DIM] = yshift[i + DIM] + fvector[i+ DIM] * dl;
                A_u_temp[i] = fvector[i + DIM]; // STORE A_u(lambda)
        }

        // Step 3: Compute A_u(lambda + dlambda) (Eq 14c)
        f(yshift, fvector); // fvector now contains A_u(lambda + dl)

        // Step 4: Compute new velocity (Eq 14d)
        LOOP_i {
                y[i] = yshift[i]; // X_u(l + dl)
                y[i + DIM] += 0.5 * (A_u_temp[i] + fvector[i + DIM]) * dl; // A_u(l+dl)
        }
}

// Returns an appropriate stepsize dlambda, which depends on position & velocity
// Ref. DOLENCE & MOSCIBRODZKA 2009
real stepsize(real X_u[4], real U_u[4]){
        real SMALL =  1.e-20;

        real dlx1  = STEPSIZE / (fabs(U_u[1]) + SMALL*SMALL);
        real dlx2  = STEPSIZE * fmin( X_u[2], 1. - X_u[2])/ (fabs(U_u[2]) + SMALL*SMALL); //TODO watch out with MKS raptor convention. 1. - X2 -> M_PI - X2
        real dlx3  = STEPSIZE / (fabs(U_u[3]) + SMALL*SMALL);

        real idlx1 = 1. / (fabs(dlx1) + SMALL*SMALL);
        real idlx2 = 1. / (fabs(dlx2) + SMALL*SMALL);
        real idlx3 = 1. / (fabs(dlx3) + SMALL*SMALL);

        return -(1. / (idlx1 + idlx2 + idlx3) + SMALL);
}

// The function to be used by the integrator for GR geodesic calculations.
// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_geodesic(real *y, real *fvector){
        // Create variable (on the stack) for the connection
        real gamma_udd[4][4][4];
        // Initialize position, four-velocity, and four-acceleration vectors based
        // on values of y
        real X_u[4] = {y[0], y[1], y[2], y[3]}; // X
        real U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
        real A_u[4] = {0.,   0.,   0.,   0.  }; // d^2X/dLambda^2

        // Obtain the Christoffel symbols at the current location
        connection_udd(X_u,gamma_udd);
        //connection_num_udd(X_u, gamma_udd);

        // Compute 4-acceleration using the geodesic equation
        int i, j, k; // Einstein summation over indices v and w
        LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];
        LOOP_i {
                fvector[i]       = U_u[i];
                fvector[i + DIM] = A_u[i];
        }
}

// Integrate the null geodesic defined by "photon_u"
void integrate_geodesic(int icur,int x, int y, real intensityfield2[maxsize][num_indices],real *frequencies, real **** p,real t,real Xcam[4],real Ucam[4]){
        int f,i, q;
        real t_init = 0.;
        real dlambda_adaptive;
        int theta_turns = 0;
        real thetadot_prev;
        real X_u[4], k_u[4];
        real tau[num_indices];
        real alpha,beta;
        for( f = 0; f < num_indices; f++)
                tau[f]=0.0;

#if (GEOD)
        struct stat st = {0};
        char geod_folder[256]="";
        sprintf(geod_folder, "geod");

        if (stat(geod_folder, &st) == -1) {
                mkdir(geod_folder, 0700);
        }

        char geod_filename[100];
        sprintf(geod_filename,"%s/geodesic_%d_%d.dat",geod_folder,x,y);
        FILE *geod = fopen(geod_filename,"w");
#endif

        real lightpath[15];
        for( i=0; i<15; i++)
                lightpath[i]=0;
        int steps = 0;
        real photon_u[8];
        for( i=0; i<8; i++)
                photon_u[i]=0;

        // Create initial ray conditions
#if (LINEAR_IMPACT_CAM)
        real stepx = CAM_SIZE_X / (real) IMG_WIDTH;
        real stepy = CAM_SIZE_Y / (real) IMG_HEIGHT;
        alpha = -CAM_SIZE_X * 0.5 + (x + 0.5) * stepx;
        beta  = -CAM_SIZE_Y * 0.5 + (y + 0.5) * stepy;
#elif (LOG_IMPACT_CAM)
        real stepx = (real)CAM_SIZE_X / (real) IMG_WIDTH;
        real stepy = (real)CAM_SIZE_Y / (real) IMG_HEIGHT;
        real r_i = exp(log(20.)*(real)(x+0.5)/((real)IMG_WIDTH)) - 1.;
        real theta_i = 2.0*M_PI  * (real)(y+0.5)/((real)IMG_HEIGHT);

        alpha = r_i * cos(theta_i);
        beta  = r_i * sin(theta_i);
#endif
        initialize_photon(alpha, beta, photon_u, t_init);

        LOOP_i {
                X_u[i] = photon_u[i];
        }


#if (GEOD)
        fprintf(geod,"%e %e %e\n",X_u[1],X_u[2],X_u[3]);
#endif

        // Current r-coordinate
        real r_current = logscale ? exp(photon_u[1]) : photon_u[1];
        // Reset lambda and steps
        real lambda = 0.;
        steps = 0;
        int TERMINATE = 0; // Termination condition for ray

        // Trace light ray until it reaches the event horizon or the outer
        // cutoff, or steps > max_steps


#if ( metric == BL || metric == MBL)
        // Stop condition for BL coords
        while (r_current > cutoff_inner && r_current < cutoff_outer &&
               steps < max_steps && !TERMINATE) { // && photon_u[0] < t_final){
#elif (metric == KS || metric == MKS)
        // Stop condition for KS coords
        while ( r_current < cutoff_outer && r_current > cutoff_inner &&
                steps < max_steps && !TERMINATE) { // 2.26 for Neuton star 3 solar masses
#endif

                // Current photon position/wave vector
                LOOP_i {
                        X_u[i] = photon_u[i];
                        k_u[i] = photon_u[i + 4];
                }
#if (GEOD)
//            real theta = M_PI * X_u[2] + 0.5 * (1. - hslope) * sin(2. * M_PI * X_u[2]);
                fprintf(geod,"%e %e %e\n",X_u[1],X_u[2],X_u[3]); //r_current*sin(theta)*cos(X_u[3]),r_current*sin(theta)*sin(X_u[3]),r_current*cos(theta));
#endif
                // Possibly terminate ray to eliminate higher order images
                if (thetadot_prev * photon_u[6] < 0. && steps > 2)
                        theta_turns += 1;
                thetadot_prev = photon_u[6];
                if((beta < 0. && theta_turns > max_order) || (beta > 0. && theta_turns > (max_order + 1)))
                        TERMINATE = 1;

                // Compute adaptive step size
                //        dlambda_adaptive = -STEPSIZE;

                dlambda_adaptive = stepsize(X_u, k_u);
                // Enter current position/velocity/dlambda into lightpath
                for (q = 0; q < 8; q++)
                        lightpath[ q] = photon_u[q];
                lightpath[ 8] = fabs(dlambda_adaptive);

                // Advance ray/particle
#if (int_method == RK4)
                rk4_step(photon_u, dlambda_adaptive);
#elif (int_method == VER)
                verlet_step(photon_u, &f_geodesic, dlambda_adaptive);
#elif (int_method == RK2)
                rk2_step(photon_u, dlambda_adaptive);
#endif

                LOOP_i {X_u[i] = photon_u[i];
                        k_u[i] = photon_u[i + 4]; }

                if(X_u[1] > stopx[1] &&  k_u[1] < 0) {
                        break;
                }

                r_current = logscale ? exp(photon_u[1]) : photon_u[1];

#if (RAD_TRANS)
                if(X_u[1] < stopx[1] && X_u[1] > startx[1] && X_u[2] > startx[2] && X_u[2] < stopx[2])
                        radiative_transfer(X_u,k_u,lightpath[8], frequencies,icur,intensityfield2,tau,p);
#endif

                lambda += fabs(dlambda_adaptive);
                r_current = logscale ? exp(photon_u[1]) : photon_u[1];
                steps = steps + 1;

        }
#if (GEOD)
        fclose(geod);
#endif
}
