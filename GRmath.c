/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Mościbrodzka
 *
 * Indices 0,1,2,3 correspond to t,r,theta,phi
 *
 * Sign convention: (-,+,+,+)
 *
 */

#include "functions.h"
#include "parameters.h"

// Lowers the index of the contravariant vector V_u, storing the results in a
// covariant one (V_d), based on the metric at position X_u
void lower_index(const real X_u[4], real g_dd[4][4], real V_u[4], real V_d[4]){
        // Obtain the covariant metric g_dd at X_u
        //real g_dd[4][4];
        //metric_dd(X_u, g_dd);
        // Initialize V_d
        V_d[0] = 0.; V_d[1] = 0.; V_d[2] = 0.; V_d[3] = 0.;

        // Lower the index of X_u
        int i, j; // Einstein summation over index j
        LOOP_ij V_d[i] += g_dd[i][j] * V_u[j];
}

// Lowers two indices on a rank (2, 0) tensor: T_uu -> T_dd at location X_u.
void lower_two_indices(real N_uu[4][4], real N_dd[4][4], real X_u[4]){
        real g_dd[4][4];
        int i, j, k, l;
        LOOP_ij N_dd[i][j] = 0.;
        metric_dd(X_u, g_dd);

        LOOP_ijkl N_dd[i][j] += g_dd[i][k] * g_dd[j][l] * N_uu[k][l];
}

// Lowers the index of a contravariant vector V_u in BL coordinates.
void BL_lower_index(const real X_u[4], real V_u[4], real V_d[4]){
        real r       = logscale ? exp(X_u[1]) : X_u[1];
        real rfactor = logscale ? r : 1.;
        real theta = X_u[2];
        real sint  = sin(theta);
        real cost  = cos(theta);
        real sigma = r * r + a * a * cost * cost;
        real delta = r * r + a * a - 2. * r;
        real A_    = (r * r + a * a) * (r * r + a * a) - delta * a * a *
                     sint * sint;

        // Covariant metric elements
        real g_dd_00 = -(1. - 2. * r / sigma);
        real g_dd_11 = sigma / delta * rfactor * rfactor;
        real g_dd_22 = sigma;
        real g_dd_33 = A_ / sigma * sint * sint;
        real g_dd_03 = -2. * a * r * sint * sint / sigma;

        V_d[0] = g_dd_00 * V_u[0] + g_dd_03 * V_u[3];
        V_d[1] = g_dd_11 * V_u[1];
        V_d[2] = g_dd_22 * V_u[2];
        V_d[3] = g_dd_33 * V_u[3] + g_dd_03 * V_u[0];
}

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
void raise_index(const real X_u[4], real V_d[4], real V_u[4]){
        // Obtain the contravariant metric g_uu at X_u
        real g_uu[4][4];
        metric_uu(X_u, g_uu);

        // Initialize V_u
        V_u[0] = 0.; V_u[1] = 0.; V_u[2] = 0.; V_u[3] = 0.;

        // Raise the index of X_d
        int i, j; // Einstein summation over index j
        LOOP_ij V_u[i] += g_uu[i][j] * V_d[j];
}

// Adjusts y[4] = U_u[0] = U^t so that y describes a lightray/null geodesic.
// This function works for all metrics.
// **NOTE** THIS FUNCTION IS CORRECT BUT NOT USED IN RAPTOR AT THE MOMENT
void normalize_null(real X_u[4], real U_u[4]){
        // Obtain the covariant metric at X_u
        real g_dd[4][4];
        metric_dd(X_u, g_dd);


        // Now we get a quadratic equation for U_u_t:
        real aa = g_dd[0][0];
        real bb = 2. * (g_dd[0][1] * U_u[1] + g_dd[0][2] * U_u[2] +
                        g_dd[0][3] * U_u[3]);
        real cc = g_dd[1][1] * U_u[1] * U_u[1] + g_dd[2][2] * U_u[2] * U_u[2] +
                  g_dd[3][3] * U_u[3] * U_u[3] + 2. *
                  (g_dd[1][2] * U_u[1] * U_u[2] + g_dd[1][3] * U_u[1] * U_u[3] +
                   g_dd[2][3] * U_u[2] * U_u[3]);

        // Two solutions, two directions for the ray
        real U_u_t_1 = (-bb + sqrt(bb * bb - 4. * aa * cc)) / (2. * aa);
        //real U_u_t_2 = -b - sqrt(b * b - 4. * a * c) / (2. * a);

        U_u[0] = U_u_t_1;
}

// Returns the norm of U_u, which is the scalar g_dd[a][b] * U_u[a] * U_u[b]
//MO is this just a dot product?why such a weird name?
real four_velocity_norm(real X_u[4], real U_u[4]){
        // Obtain the covariant metric at X_u
        real g_dd[4][4];
        metric_dd(X_u, g_dd);

        // Compute the norm
        real norm = 0.;
        int i, j; // Einstein summation over indices i and j
        LOOP_ij norm += g_dd[i][j] * U_u[i] * U_u[j];

        return norm;
}

real inner_product(real *X_u, real *A_u, real *B_u){
        // Obtain the covariant metric at X_u
        real g_dd[4][4];
        metric_dd(X_u, g_dd);

        // Compute the dot produt
        real dotproduct = 0.;
        int i, j; // Einstein summation over indices i and j
        LOOP_ij dotproduct += g_dd[i][j] * A_u[i] * B_u[j];

        return dotproduct;
}


// Transform a PHOTON (contravariant position and velocity vectors)
// from BL to KS coordinates
void BL_to_KS_u(real *BLphoton_u, real *KSphoton_u){
        real trans[4][4];
        real X_u[4], U_u[4];

        int i, j;
        LOOP_i {
                X_u[i] = BLphoton_u[i];
                U_u[i] = BLphoton_u[i+4];
        }

        // Construct BL -> MKS matrix
        LOOP_ij trans[i][j] = 0.;
        LOOP_i trans[i][i] = 1.;

        // Note that r and theta are identical in BL and KS.
        // See McKinney & Gammie (2004)
        real r_current2    = logscale ? exp(BLphoton_u[1]) : BLphoton_u[1];
        real delta_current = r_current2 * r_current2 - 2. * r_current2 +
                             a * a;
        real rfactor = logscale ? r_current2 : 1.;
        trans[0][1] = 2. * r_current2 / delta_current * rfactor;
        trans[3][1] = a / delta_current * rfactor;

        // Do the transformation
        real U_u_dummy[4], X_u_dummy[4];
        LOOP_i {
                U_u_dummy[i] = U_u[i];
                X_u_dummy[i] = X_u[i];
                U_u[i] = 0.;
                X_u[i] = 0.;
        }

        // Transform the wave vector
        LOOP_ij U_u[i] += trans[i][j] * U_u_dummy[j];

        real rplus = 1. + sqrt(1. - a * a);
        real rmin = 1. - sqrt(1. - a * a);

        // Transform t and phi for the position vector
        X_u[1] = X_u_dummy[1];
        X_u[2] = X_u_dummy[2];
        X_u[0] = X_u_dummy[0] + (log(delta_current) + 1. / sqrt(1. - a*a) * log((r_current2 - rplus) / (r_current2 - rmin)));
        X_u[3] = X_u_dummy[3] + (a / (2. * sqrt(1. - a * a)) * log((r_current2 - rplus) / (r_current2 - rmin)));

        // Put result in photon variable
        LOOP_i {
                KSphoton_u[i] = X_u[i];
                KSphoton_u[i+4] = U_u[i];
        }
}

// Transform a contravariant vector from KS to BL coordinates
void KS_to_BL_u(real *KSphoton_u, real *BLphoton_u){
        real trans[4][4];
        real X_u[4], U_u[4];

        int i, j;
        LOOP_i {
                X_u[i] = KSphoton_u[i];
                U_u[i] = KSphoton_u[i+4];
        }

        // Construct BL -> MKS matrix
        LOOP_ij trans[i][j] = 0.;
        LOOP_i trans[i][i] = 1.;

        // Note that r and theta are identical in BL and KS.
        // See McKinney & Gammie (2004)
        real r_current2    = logscale ? exp(KSphoton_u[1]) : KSphoton_u[1];
        real delta_current = r_current2 * r_current2 - 2. * r_current2 +
                             a * a;
        real rfactor = logscale ? r_current2 : 1.;
        trans[0][1] = -(2. * r_current2 / delta_current) * rfactor;
        trans[3][1] = -(a / delta_current) * rfactor;

        // Do the transformation
        real U_u_dummy[4], X_u_dummy[4];
        LOOP_i {
                U_u_dummy[i] = U_u[i];
                X_u_dummy[i] = X_u[i];
                U_u[i] = 0.;
                X_u[i] = 0.;
        }

        // Transform the wave vector using the BL->KS matrix given in literature
        LOOP_ij U_u[i] += trans[i][j] * U_u_dummy[j];

        real rplus = 1. + sqrt(1. - a * a);
        real rmin = 1. - sqrt(1. - a * a);

        // Transform t and phi for the position vector (transforms differently!)
        X_u[1] = X_u_dummy[1];
        X_u[2] = X_u_dummy[2];
        X_u[0] = X_u_dummy[0] - (log(delta_current) + 1. / sqrt(1. - a*a) * log((r_current2 - rplus) / (r_current2 - rmin)));
        X_u[3] = X_u_dummy[3] - (a / (2. * sqrt(1. - a * a)) * log((r_current2 - rplus) / (r_current2 - rmin)));

        // Put transformed photon in BLphoton_u variable
        LOOP_i {
                BLphoton_u[i] = X_u[i];
                BLphoton_u[i+4] = U_u[i];
        }
}
