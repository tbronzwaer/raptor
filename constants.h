/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 * A list of physical constants used by RAPTOR.
 * CGS units are used.
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

// PHYSICAL CONSTANTS
#define ELECTRON_CHARGE    (4.80320425e-10)
#define ELECTRON_MASS      (9.1093829e-28)
#define PROTON_MASS        (1.6726219e-24)
#define BOLTZMANN_CONSTANT (1.3806488e-16)
#define SPEED_OF_LIGHT     (2.99792458e10)
#define PLANCK_CONSTANT    (6.62606885e-27)
#define MPCL2              (0.0015033)
#define GGRAV              (6.674e-8)
#define MSUN               (1.989e33)
#define MPoME              (PROTON_MASS / ELECTRON_MASS)
#define M_PI           3.14159265358979323846


// Constants that must be evaluated at startup
// (They depend on spin and other user-supplied parameters)
double R_GRAV; // Gravitational radius
double R_ISCO; // Innermost stable circular orbit
double CUTOFF_INNER; // Inner integration boundary
double JANSKY_FACTOR; // Factor to scale image output

#endif // CONSTANTS_H
