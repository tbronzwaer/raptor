/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Moscibrodzka
 *
 * RAPTOR uses cgs units for light transport calculations.
 * Entries marked [BL] are only applicable to Boyer-Lindquist coordinates.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

//test
//double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;


// METRIC PARAMETERS
////////////////////

// These are used for geodesic integration; BH mass normalized

//coordinate and metric choices
#define CAR      (0)        //  Minkowski
#define BL       (1)        // Boyer-Lindquist,               x1=r, x2=th, and x3=phi
#define MBL      (2)        // modified Boyer-Lindquist, x1=log(r), x2=th, and x3=phi
#define KS       (3)        // Kerr-Schild,                   x1=r, x2=th, and x3=phi
#define MKS      (4)        // modified Kerr-Schild,     x1=log(r), x2=th, and x3=phi
#define MKS2     (9)        // Proper MKS coords
#define DM       (5)
//#define a        (0.0)    // Black hole spin (range: [0,1]), only for BL,MBL,KS,MKS
//#define R0       (0.0)      // Parameter for MKS coordinates.
extern double a;
extern double R0;

// Metric
#define metric   (MKS2)
#if(metric == BL || metric == KS || metric == DM)
    #define logscale (0)    // Standard BL/KS coordinates; no logarithmic radius
#elif(metric == MBL || metric == MKS || metric == MKS2)
    #define logscale (1)    // Modified BL/KS coordinates; logarithmic radius
#endif

// MODEL PARAMETERS
///////////////////

// These are used for light transport computations; BH now has a specific mass

// GRMHD data file
char GRMHD_FILE[256];
char OUTPUT_FILE[256];
//#define GRMHD_FILE "grmhd/spherical_accretion/dump0000"
//#define GRMHD_FILE "grmhd/michael/dump234"
//#define GRMHD_FILE "grmhd/bhac/dump500"
//#define GRMHD_FILE "grmhd/dump019"

int SPHERICAL_ACC;
//#define SPHERICAL_ACCRETION (1) // Indicate whether or not to use the spherical accretion model for the B field
char TEMP_MODEL[100];
//#define DISK_MODEL          (0) // DOES NOTHING RIGHT NOW, TODO DEPRECATED

//#define MICHAEL_MODEL       (0) // Michael Janssen dumps model
int ABSORPTION;
//#define ABSORPTION          (1)
#define LIGHT_TRANSPORT     (1) // Toggle light transport calculation on/off for integration debugging

#define RT_OUTER_CUTOFF     (40.) // Outer boundary of radiative transfer computation

// Black hole mass
double MBH;
//#define MBH           (8.95e39)  // 4.5e6 solar masses (Sgr A*, dump019)  (g)
//#define MBH           (6.962e42) // M87, Michael's run
//#define MBH (1.989e33) // Solar mass BH; for spherical accretion (g)

//#define M_unit  (1.e19) // 1.e19 = Sgr A*
double M_UNIT;
//#define M_unit    (1.e15) // For Jordy's big BHAC run, dump019
//#define M_unit    (2.3e28) // M87, for Michael's run
//#define M_unit  (4.543e-7) // For spherical accretion

//#define source_dist    (2.47e22) // Distance to Sgr A* (cm)
#define source_dist    (5.061e25) // Distance to M87 (cm)

// These are for simple analytical model - move to different .c file!!
#define n_e0           (4.5e6)   // Electron density normalization
#define B0             (100.)    // Magnetic field scalar
#define THETA_e_0      (80.)     // Dimensionless temperature scalar
#define nblobs         (1)       // number of blobs?

// Dark matter spike parameters
#define qfactor      (0.1) // M_DM/M_BH, see Lacroix & Silk 2013

// OBSERVER PARAMETERS
//////////////////////

double CAM_FREQ;
double TIME_INIT;
double INCLINATION;

// SED parameters
int    FREQS_PER_DEC;
double FREQ_MIN;
double FREQ_MAX;

#define rcam         (1e4)//(500.)    // Camera distance from the sing.(units of Rg)

int IMG_WIDTH;
int IMG_HEIGHT;
double CAM_SIZE_X;
double CAM_SIZE_Y;
#define max_order    (100)       // Maximimum order of lensed image computed (0 = direct only)

// INTEGRATOR PARAMETERS
////////////////////////

#define delta_num    (1.e-6)    // Used for numerical derivatives
#define max_steps    (1e7)   // Maximum number of integration steps

double STEPSIZE;
#define cutoff_outer (1.1 * rcam)    // Outer cutoff, near flat spacetime, in M
#define horizon_marg (1.e-5)     // Stop tracing at this distance from E.H. [BL]
#define VER       (1)        //
#define RK4       (2)        //
#define int_method   (2)     // method of integration 2=Verlet, 4-RK4

// MACROS
/////////

#define DIM 4
#define LOOP_i    for(i = 0; i < DIM; i++)
#define LOOP_ij   for(i = 0; i < DIM; i++) \
                  for(j = 0; j < DIM; j++)
#define LOOP_kl   for(k = 0; k < DIM; k++) \
                  for(l = 0; l < DIM; l++)
#define LOOP_ijk  for(i = 0; i < DIM; i++) \
                  for(j = 0; j < DIM; j++) \
                  for(k = 0; k < DIM; k++)
#define LOOP_ijkl for(i = 0; i < DIM; i++) \
                  for(j = 0; j < DIM; j++) \
                  for(k = 0; k < DIM; k++) \
                  for(l = 0; l < DIM; l++)

#endif // PARAMETERS_H
