/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Mościbrodzka
 *
 * A list of all RAPTOR functions.
 */


#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include "parameters.h"

void init_model();
#pragma acc routine (get_fluid_params)
int get_fluid_params(real X[NDIM], real *Ne,
                     real *Thetae, real *B, real *beta, real * Bern, real *B_u, real Ucon[NDIM], int *IN_VOLUME, real **** p,real gcov[NDIM][NDIM],real gcon[NDIM][NDIM]);
// GRMATH.C
///////////

// Lowers the index of the contravariant vector V_u, storing the results in a
// covariant one (V_d), based on the metric at position X_u
#pragma acc routine(lower_index)
void lower_index(const real X_u[4],real g_dd[4][4], real V_u[4], real V_d[4]);

// Lowers two indices of a rank (2, 0) tensor
void lower_two_indices(real N_uu[4][4], real N_dd[4][4], real X_u[4]);

// Lowers the index of a contravariant vector V_u in BL coordinates.
void BL_lower_index(const real X_u[4], real V_u[4], real V_d[4]);

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
void raise_index(const real X_u[4], real V_d[4], real V_u[4]);

// Adjusts y[4] = U_u[0] so that y describes a lightray/null geodesic
void normalize_null(real X_u[4], real U_u[4]);

// Returns the norm of U_u, which is the scalar g_dd[a][b] * U_u[a] * U_u[b]
real four_velocity_norm(real X_u[4], real U_u[4]);

// Returns the inner product of vectors A and B, i.e. A_u B_d
#pragma acc routine (inner_product)
real inner_product(real *X_u, real *A_u, real *B_u);

// Transform a contravariant vector from BL to KS coordinates
#pragma acc routine (BL_to_KS_u)
void BL_to_KS_u(real *BLphoton_u, real *KSphoton_u);

// Transform a contravariant vector from KS to BL coordinates
#pragma acc routine (KS_to_BL_u)
void KS_to_BL_u(real *KSphoton_u, real *BLphoton_u);

// CORE.C
///////////////
void output_files(real ** intesityfield,real spectrum[num_indices],real frequencies[num_indices]);

void calculate_image( real ** intensityfield, real spectrum[num_indices],real frequencies[num_indices]);

void read_model( char *argv[]);

// INTEGRATOR.C
///////////////

// Returns an appropriate stepsize dlambda, which depends on position & velocity
#pragma acc routine (stepsize)
real stepsize(real X_u[4], real U_u[4]);

// Updates the vector y (containing position/velocity) by one RK4 step.
#pragma acc routine (rk4_step)
void rk4_step(real *y,real dt);

// Updates the vector y (containing position/velocity) by one RK2 step.
#pragma acc routine (rk2_step)
void rk2_step(real *y, real dt);

// Updates the vector y (containing position/velocity) by one Verlet step.
#pragma acc routine (verlet_step)
void verlet_step(real *y, void (*f)(real*, real*), real dl);
// The function to be used by the integrator - it solves the geodesic eqn.
// y contains the 4-position and the 4-velocity for one lightray/particle
#pragma acc routine (f_geodesic)
void f_geodesic(real *y, real *fvector);

// Integrate the null geodesic specified by alpha and beta, store results
// in lightpath
#pragma acc routine (integrate_geodesic)
void integrate_geodesic(int icur,int x, int y, real intensityfield2[maxsize][num_indices],real *frequencies, real **** p,real t,real Xcam[4],real Ucam[4]);

// METRIC.C
///////////

// Computes the metric at location X
#pragma acc routine (metric_dd)
void metric_dd(const real X_u[4], real g_dd[4][4]);

// Computes the inverse metric at location X
#pragma acc routine (metric_uu)
void metric_uu(const real X_u[4], real g_uu[4][4]);

// Computes the Christoffel symbols at location X numerically (general metric)
void connection_num_udd(const real X_u[4], real gamma_udd[4][4][4]);

// Computes the Christoffel symbols at location X based on an exact metric
#pragma acc routine (connection_udd)
void connection_udd(const real X_u[4], real gamma_udd[4][4][4]);

// This function initializes a single 'superphoton' or light ray.
#pragma acc routine (initialize_photon)
void initialize_photon(real alpha, real beta, real photon_u[8], real t_init);


// RADIATIVE_TRANSFER.C
///////////////////////

void read_in_table(char* filename);//, real ****j_nu_data, real ****alpha_nu_data);
void init_memory_kappa();//real ****j_nu_data,real ****alpha_nu_data);
real interpolate_scalar_4d(real ****A, real nuratio, real kappa, real Thetae, real theta);
real interpolate_scalar_2d(real **var, real nuratio, real Thetae);

real absorption_coeff_kappa(real **** alpha_nu_data,real nu,real Ne, real Thetae, real B,real beta, real theta);
real emission_coeff_kappa(real **** j_nu_data,real nu,real Ne, real Thetae, real B,real beta, real theta);

// Return emission coefficient j_nu for kappa distribution function
//#pragma acc routine(emission_coeff_kappa_FIT)
real emission_coeff_kappa_FIT(real nu,real Ne, real Thetae, real B,real beta, real theta);

// Return absorption coefficient for kappa distribution function
//#pragma acc routine(absorption_coeff_kappa_FIT)
real absorption_coeff_kappa_FIT(real nu, real Ne, real Thetae, real B,real beta, real theta);

// Return emission coefficient j_nu for thermal synchrotron radiation
#pragma acc routine (emission_coeff_THSYNCH)
real emission_coeff_THSYNCH(real B_, real theta, real THETA_e_, real nu_plasma, real n_e);
#pragma acc routine (emission_coeff_THSYNCHAV)
// Return emission coefficient for angle-averaged thermal synchrotron radiation
real emission_coeff_THSYNCHAV(real B_, real THETA_e_, real nu_plasma, real n_e);

// Return emission coefficient for thermal free-free radiation
real emission_coeff_FFTHERMAL(real nu, real n_e, real T);

// Simple emissivity model: orbiting Gaussian hotspot (see Dexter 2009)
real emissivity_hotspot(real *X_u);

// Simple emissivity model: thin disk line emission (see Dexter 2009)
real emissivity_thindisk(real *X_u);

// Return absorption coefficient a_nu
#pragma acc routine (absorption_coeff_TH)
real absorption_coeff_TH(real j_nu, real nu, real THETA_e);

// Return the photon frequency in the co-moving frame of the plasma
#pragma acc routine (freq_in_plasma_frame)
real freq_in_plasma_frame(real Uplasma_u[4], real k_d[4]);

// Planck function
#pragma acc routine (planck_function)
real planck_function(real nu, real THETA_e);

// Angle between k_u and B_u in the plasma frame
#pragma acc routine (pitch_angle)
real pitch_angle(real *X_u, real *k_u, real *B_u, real *Uplasma_u,real B);

// Perform radiative transfer along the ray stored in "lightpath"
#pragma acc routine (radiative_transfer)
real radiative_transfer(real *X_u, real *k_u,real dl_current, real *frequency,int icur,real intensity[maxsize][num_indices], real *tau,real **** p);

// Backward transfer
real backward_transfer(real alpha, real beta, real *photon_u, int *steps);

// UTILITIES.C
//////////////

void set_constants(real t);

// Write the array "intensityfield" (scaled by "scalefactor") to the file "imgfile"
void write_image(FILE *imgfile, real **intensityfield,int f, real scalefactor);

// Write the arrays "intensityfield" (scaled by "scalefactor") and "lambdafield" to a VTK file
void write_VTK_image(FILE *imgfile, real *intensityfield, real *lambdafield, real scalefactor);

 void print_time(int start);


// VR checker bord sphere functions

#pragma acc routine (color_sphere_outer)
real color_sphere_outer(real X_u[NDIM],int icur,real intensity[maxsize][num_indices],real tau);

#pragma acc routine (color_sphere)
real color_sphere(real X_u[NDIM],int icur,real intensity[maxsize][num_indices]);

#pragma acc routine (intersect)
int intersect(real X_u[NDIM]);

#pragma acc routine (color)
real color(real X_u[NDIM]);

#pragma acc routine (color_outer)
real color_outer(real X_u[NDIM]);

// RCARRY FUNCTIONS
real genrandrcarry();
void initrcarry(int seed_);

// NEWTON-RAPHSON FUNCTIONS
real f_Xg2(real Xg2, real Xr2,real lr,int init);
real f_primed_Xg2(real Xg2,real lr,int init);
real NR_stepX(real Xg2_0, real Xr2,real lr,int init);
real f_Ug2(real Ug2, real Ug1, real Ur2, real Xg2, real lr,int init);
real f_primed_Ug2(real Ug2, real Xg2, real lr,int init);
real NR_stepU(real Ug2_0,real Ug1, real Ur2, real Xg2, real lr,int init);
real Xg2_approx_rand(real Xr2, real lr,int init2);
real Ug2_approx_rand(real Ur2, real Ug1,real Xg2, real lr,int init2);


#endif // FUNCTIONS_H
