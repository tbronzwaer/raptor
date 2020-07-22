/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 * Indices 0,1,2,3 correspond to t,r,theta,phi
 *
 * Sign convention: (-,+,+,+)
 *
 */
#include <math.h>
#include "functions.h"
#include "parameters.h"

// Lowers the index of the contravariant vector V_u, storing the results in a
// covariant one (V_d), based on the metric at position X_u
void lower_index(const double X_u[4], double V_u[4], double V_d[4]){
    // Obtain the covariant metric g_dd at X_u
    double g_dd[4][4];
    metric_dd(X_u, g_dd);

    // Initialize V_d
    V_d[0] = 0.; V_d[1] = 0.; V_d[2] = 0.; V_d[3] = 0.;

    // Lower the index of X_u
    int i, j; // Einstein summation over index j
    LOOP_ij V_d[i] += g_dd[i][j] * V_u[j];
}

// Lowers two indices on a rank (2, 0) tensor: T_uu -> T_dd at location X_u.
void lower_two_indices(double N_uu[4][4], double N_dd[4][4], double X_u[4]){
    double g_dd[4][4];
    int i, j, k, l;
    LOOP_ij N_dd[i][j] = 0.;
    metric_dd(X_u, g_dd);

    LOOP_ijkl N_dd[i][j] += g_dd[i][k] * g_dd[j][l] * N_uu[k][l];
}

// Lowers the index of a contravariant vector V_u in BL coordinates.
void BL_lower_index(const double X_u[4], double V_u[4], double V_d[4]){
    double r       = logscale ? exp(X_u[1]) : X_u[1];
    double rfactor = logscale ? r : 1.;
    double theta = X_u[2];
    double sint  = sin(theta);
    double cost  = cos(theta);
    double sigma = r * r + a * a * cost * cost;
    double delta = r * r + a * a - 2. * r;
    double A_    = (r * r + a * a) * (r * r + a * a) - delta * a * a *
                   sint * sint;

    // Covariant metric elements
    double g_dd_00 = -(1. - 2. * r / sigma);
    double g_dd_11 = sigma / delta * rfactor * rfactor;
    double g_dd_22 = sigma ;
    double g_dd_33 = A_ / sigma * sint * sint;
    double g_dd_03 = -2. * a * r * sint * sint / sigma;

    V_d[0] = g_dd_00 * V_u[0] + g_dd_03 * V_u[3];
    V_d[1] = g_dd_11 * V_u[1];
    V_d[2] = g_dd_22 * V_u[2];
    V_d[3] = g_dd_33 * V_u[3] + g_dd_03 * V_u[0];
}

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
void raise_index(const double X_u[4], double V_d[4], double V_u[4]){
    // Obtain the contravariant metric g_uu at X_u
    double g_uu[4][4];
    metric_uu(X_u, g_uu);

    // Initialize V_u
    V_u[0] = 0.; V_u[1] = 0.; V_u[2] = 0.; V_u[3] = 0.;

    // Raise the index of X_d
    int i, j; // Einstein summation over index j
    LOOP_ij V_u[i] += g_uu[i][j] * V_d[j];
}

// Adjusts k_u[0] = k^t so that k_u describes a lightray/null geodesic.
// This function works for all metrics.
void normalize_null(double X_u[4], double k_u[4]){
    // Obtain the covariant metric at X_u
    double g_dd[4][4];
    int i, j;
    LOOP_ij g_dd[i][j] = 0.;
    metric_dd(X_u, g_dd);

    // Now we get a quadratic equation for k_u_t:
    double aa = g_dd[0][0];
    double bb = 2. * (g_dd[0][1] * k_u[1] + g_dd[0][2] * k_u[2] +
                      g_dd[0][3] * k_u[3]);
    double cc = g_dd[1][1] * k_u[1] * k_u[1] + g_dd[2][2] * k_u[2] * k_u[2] +
                g_dd[3][3] * k_u[3] * k_u[3] + 2. *
                (g_dd[1][2] * k_u[1] * k_u[2] + g_dd[1][3] * k_u[1] * k_u[3] +
                 g_dd[2][3] * k_u[2] * k_u[3]);

    // Two solutions, two directions for the ray
    double k_u_t_1 = -bb + sqrt(bb * bb - 4. * aa * cc) / (2. * aa);

    k_u[0] = k_u_t_1;

    double k_u1 = k_u[1];
    double k_u2 = k_u[2];
    double k_u3 = k_u[3];

    double Betacap = -(g_dd[0][1] * k_u1)/(g_dd[0][0]) - (g_dd[0][2] * k_u2)/(g_dd[0][0]) - (g_dd[0][3] * k_u3)/(g_dd[0][0]);

    double gammaZiri = -(g_dd[1][1] * k_u1 * k_u1) / (g_dd[0][0]) - (g_dd[1][2] * k_u1 * k_u2) / (g_dd[0][0]) - (g_dd[1][3] * k_u1 * k_u3) / (g_dd[0][0])
                       -(g_dd[2][1] * k_u2 * k_u1) / (g_dd[0][0]) - (g_dd[2][2] * k_u2 * k_u2) / (g_dd[0][0]) - (g_dd[2][3] * k_u2 * k_u3) / (g_dd[0][0])
                       -(g_dd[3][1] * k_u3 * k_u1) / (g_dd[0][0]) - (g_dd[3][2] * k_u3 * k_u2) / (g_dd[0][0]) - (g_dd[3][3] * k_u3 * k_u3) / (g_dd[0][0]);

    // Temporal component of wave vector
    double k_u0 = Betacap + sqrt(Betacap * Betacap + gammaZiri);

    k_u[0] = k_u0;
}
// Returns the norm of U_u, which is the scalar g_dd[a][b] * U_u[a] * U_u[b]
//MO is this just a dot product?why such a weird name?
double four_velocity_norm(double X_u[4], double U_u[4]){
    // Obtain the covariant metric at X_u
    double g_dd[4][4];
    metric_dd(X_u, g_dd);

    // Compute the norm
    double norm = 0.;
    int i, j; // Einstein summation over indices i and j
    LOOP_ij norm += g_dd[i][j] * U_u[i] * U_u[j];

    return norm;
}

double inner_product(double *X_u, double *A_u, double *B_u){
    // Obtain the covariant metric at X_u
    double g_dd[4][4];
    metric_dd(X_u, g_dd);

    // Compute the dot produt
    double dotproduct = 0.;
    int i, j; // Einstein summation over indices i and j
    LOOP_ij dotproduct += g_dd[i][j] * A_u[i] * B_u[j];

    return dotproduct;
}

// This is a temporary function for debugging purpose:
// It takes the HARM "MKS" convention and transforms to a vector
// using the RAPTOR "MKS" convention.
void HARMMKS_to_TBMKS(double *HARM_MKS_vector_u, double *TB_MKS_vector_u){

}

// Transform a PHOTON (contravariant position and velocity vectors)
// from BL to KS coordinates
void BL_to_KS_u(double *BLphoton_u, double *KSphoton_u){
    double trans[4][4];
    double X_u[4], U_u[4];

    int i, j;
    LOOP_i{
        X_u[i] = BLphoton_u[i];
        U_u[i] = BLphoton_u[i+4];
    }

    // Construct BL -> MKS matrix
    LOOP_ij trans[i][j] = 0.;
    LOOP_i  trans[i][i] = 1.;

    // Note that r and theta are identical in BL and KS.
    // See McKinney & Gammie (2004)
    double r_current2    = logscale ? exp(BLphoton_u[1]) : BLphoton_u[1];
    double delta_current = r_current2 * r_current2 - 2. * r_current2 +
                           a * a;
    double rfactor = logscale ? r_current2 : 1.;
    trans[0][1] = 2. * r_current2 / delta_current * rfactor;
    trans[3][1] = a / delta_current * rfactor;

    // Do the transformation
    double U_u_dummy[4], X_u_dummy[4];
    LOOP_i{
        U_u_dummy[i] = U_u[i];
        X_u_dummy[i] = X_u[i];
        U_u[i] = 0.;
        X_u[i] = 0.;
    }

    // Transform the wave vector
    LOOP_ij U_u[i] += trans[i][j] * U_u_dummy[j];

    double rplus = 1. + sqrt(1. - a * a);
    double rmin = 1. - sqrt(1. - a * a);

    // Transform t and phi for the position vector
    X_u[1] = X_u_dummy[1];
    X_u[2] = X_u_dummy[2];
    X_u[0] = X_u_dummy[0] + (log(delta_current) + 1. / sqrt(1. - a*a) * log((r_current2 - rplus) / (r_current2 - rmin)));
    X_u[3] = X_u_dummy[3] + (a / (2. * sqrt(1. - a * a)) * log((r_current2 - rplus) / (r_current2 - rmin)));

    // Put result in photon variable
    LOOP_i{
        KSphoton_u[i] = X_u[i];
        KSphoton_u[i+4] = U_u[i];
    }
}

// Transform a contravariant vector from KS to BL coordinates
void KS_to_BL_u(double *KSphoton_u, double *BLphoton_u){
    double trans[4][4];
    double X_u[4], U_u[4];

    int i, j;
    LOOP_i{
        X_u[i] = KSphoton_u[i];
        U_u[i] = KSphoton_u[i+4];
    }

    // Construct BL -> MKS matrix
    LOOP_ij trans[i][j] = 0.;
    LOOP_i  trans[i][i] = 1.;

    // Note that r and theta are identical in BL and KS.
    // See McKinney & Gammie (2004)
    double r_current2    = logscale ? exp(KSphoton_u[1]) : KSphoton_u[1];
    double delta_current = r_current2 * r_current2 - 2. * r_current2 +
                           a * a;
    double rfactor = logscale ? r_current2 : 1.;
    trans[0][1] = -(2. * r_current2 / delta_current) * rfactor;
    trans[3][1] = -(a / delta_current) * rfactor;

    // Do the transformation
    double U_u_dummy[4], X_u_dummy[4];
    LOOP_i{
        U_u_dummy[i] = U_u[i];
        X_u_dummy[i] = X_u[i];
        U_u[i] = 0.;
        X_u[i] = 0.;
    }

    // Transform the wave vector using the BL->KS matrix given in literature
    LOOP_ij U_u[i] += trans[i][j] * U_u_dummy[j];

    double rplus = 1. + sqrt(1. - a * a);
    double rmin = 1. - sqrt(1. - a * a);

    // Transform t and phi for the position vector (transforms differently!)
    X_u[1] = X_u_dummy[1];
    X_u[2] = X_u_dummy[2];
    X_u[0] = X_u_dummy[0] - (log(delta_current) + 1. / sqrt(1. - a*a) * log((r_current2 - rplus) / (r_current2 - rmin)));
    X_u[3] = X_u_dummy[3] - (a / (2. * sqrt(1. - a * a)) * log((r_current2 - rplus) / (r_current2 - rmin)));

    // Put transformed photon in BLphoton_u variable
    LOOP_i{
        BLphoton_u[i] = X_u[i];
        BLphoton_u[i+4] = U_u[i];
    }
}
