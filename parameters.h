/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka
 *
 * RAPTOR uses cgs units for light transport calculations.
 * Entries marked [BL] are only applicable to Boyer-Lindquist coordinates.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef double real;


#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <png.h>
#include <omp.h>
#include <time.h>
#include "constants.h"

// OpenACC or OMP
////////////////////

#define ACC (0)
#define OMP (1)
#define COMP OMP


#include <math.h>

//#include <accel.h>
//#include <accelmath.h>
// OUTPUT SWITCHES
////////////////////

#define VTKFILE (0)
#define IMGFILE (1)
#define SPECFILE (1)

#define RAD_TRANS (1)

//GEOD FILE, print geodesics
#define GEOD (0)

// CAMERA SWITCHES
////////////////////
#define maxsize 1000000

#define LOG_IMPACT_CAM (0)
#define LINEAR_IMPACT_CAM (1)

//VR parameter
extern real Ucam[4];
extern real Xcam[4];
extern real tcam;

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

// GLOBAL VARIABLES
////////////////////
extern real L_unit;
extern real T_unit;
extern real RHO_unit;
extern real U_unit;
extern real B_unit;
extern real Ne_unit;
extern real Thetae_unit;

#define NDIM  (4)
#define NPRIM    (8)
#define N1 (256)
#define N2 (256)
#define N3 (1)

extern real R0, Rin, Rh, Rout, Rms;
extern real a;
extern real hslope;
extern real startx[NDIM], stopx[NDIM], dx[NDIM];
extern real dlE, lE0;
extern real gam;
extern real dMsim;
extern real ****p;


// METRIC PARAMETERS
////////////////////

// These are used for geodesic integration; BH mass normalized

//coordinate and metric choices
#define CAR      (0)        //  Minkowski
#define BL       (1)        // Boyer-Lindquist,               x1=r, x2=th, and x3=phi
#define MBL      (2)        // modified Boyer-Lindquist, x1=log(r), x2=th, and x3=phi
#define KS       (3)        // Kerr-Schild,                   x1=r, x2=th, and x3=phi
#define MKS      (4)        // Proper MKS coords


#define TT      (0)
#define RR      (1)
#define TH      (2)
#define PH      (3)

extern real a;
extern real R0;

// Metric
#define metric   (MKS)
#if(metric == CAR || metric == BL || metric == KS)
#define logscale (0)    // Standard BL/KS coordinates; no logarithmic radius
#elif(metric == MBL || metric == MKS)
#define logscale (1)    // Modified BL/KS coordinates; logarithmic radius
#endif

// MODEL PARAMETERS
///////////////////

// These are used for light transport computations; BH now has a specific mass

// GRMHD data file
extern char GRMHD_FILE[256];
extern char OUTPUT_FILE[256];
extern int SPHERICAL_ACC;
extern char TEMP_MODEL[100];
extern int ABSORPTION;

#define LIGHT_TRANSPORT     (1) // Toggle light transport calculation on/off for integration debugging

//	#define RT_OUTER_CUTOFF     (2000.)//1.01*rcam) // Outer boundary of radiative transfer computation
#define RT_OUTER_CUTOFF     (1.01*rcam) // Outer boundary of radiative transfer computation

// Black hole mass
extern double MBH;

//sets all units.
extern real M_UNIT;

extern real R_HIGH;
extern real R_LOW;
#define source_dist    (2.6228263e22) // Distance to Sgr A* (cm)
//#define source_dist    (5.061e25) // Distance to M87 (cm)

// These are for simple analytical model - move to different .c file!!
#define n_e0           (4.5e6)   // Electron density normalization
#define B0             (100.)    // Magnetic field scalar
#define THETA_e_0      (80.)     // Dimensionless temperature scalar
#define nblobs         (1)       // number of blobs?

// Dark matter spike parameters
#define qfactor      (0.1) // M_DM/M_BH, see Lacroix & Silk 2013

// OBSERVER PARAMETERS
//////////////////////

extern real CAM_FREQ;
extern real TIME_INIT;
extern real INCLINATION;

// SED parameters
extern int  FREQS_PER_DEC;
extern real FREQ_MIN;
extern real FREQ_MAX;


//	#define rcam         (15.)//(500.)    // Camera distance from the sing.(units of Rg)
#define rcam         (1.e4)
extern int IMG_WIDTH;
extern int IMG_HEIGHT;
extern real CAM_SIZE_X;
extern real CAM_SIZE_Y;
#define max_order    (100)       // Maximimum order of lensed image computed (0 = direct only)

// INTEGRATOR PARAMETERS
////////////////////////
extern real STEPSIZE;

#define delta_num    (1.e-7)    // Used for numerical derivatives
#define max_steps    (1e5)   // Maximum number of integration steps

#define cutoff_outer (rcam*1.01)    // Outer cutoff, near flat spacetime, in M
#define horizon_marg (1e-1)     // Stop tracing at this distance from E.H. [BL]
#define VER       (1)        //
#define RK4       (2)        //
#define RK2      (3)
#define int_method   (RK2)     // method of integration 2=Verlet, 4-RK4

// MACROS
/////////

#define DIM 4
#define LOOP_i    for(i = 0; i < DIM; i++)
#define LOOP_ij   for(i = 0; i < DIM; i++) \
for(j = 0; j < DIM; j++)

#define LOOP_kl    for(k = 0; k < DIM; k++) \
for(l = 0; l < DIM; l++)


#define LOOP_ijk2   for( i=NDIM; i--; )\
for( j=NDIM; j--; )\
for( k=NDIM; k--; )
#define LOOP_ijk    for(i = 0; i < DIM; i++) \
for(j = 0; j < DIM; j++) \
for(k = 0; k < DIM; k++)
#define LOOP_ijkl for(i = 0; i < DIM; i++) \
for(j = 0; j < DIM; j++) \
for(k = 0; k < DIM; k++) \
for(l = 0; l < DIM; l++)


#define num_indices 1

#define N_kappa 11
#define N_theta 11
#define N_theta_e 21
#define N_nuratio 56

#pragma acc copyin(IMG_WIDTH,IMG_HEIGHT,CAM_SIZE_X,CAM_SIZE_Y,cutoff_inner,ABSORPTION,MBH,INCLINATION,STEPSIZE,SPHERICAL_ACC,a,R_LOW,R_HIGH,hslope)

#endif // PARAMETERS_H
