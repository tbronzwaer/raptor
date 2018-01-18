/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Mościbrodzka
 *
 * A list of physical constants used by RAPTOR.
 * CGS units are used.
 */

#include "parameters.h"

#ifndef CONSTANTS_H
#define CONSTANTS_H

// PHYSICAL CONSTANTS
#define ELECTRON_CHARGE    (4.8032068e-10)
#define ELECTRON_MASS      (9.1093826e-28)
#define PROTON_MASS        (1.67262171e-24)
#define BOLTZMANN_CONSTANT (1.3806505e-16)
#define SPEED_OF_LIGHT     (2.99792458e10)
#define PLANCK_CONSTANT    (6.6260693e-27)
#define MPCL2              (0.0015033)
#define GGRAV              (6.6742e-8)
#define MSUN               (1.989e33)
#define MPoME              (1836.15) //(PROTON_MASS / ELECTRON_MASS)
#define M_PI               (3.14159265358979323846)
#define YEAR               (31536000.)    /* No. of seconds in year */
#define e2_c               (7.695589048e-30)

// Constants that must be evaluated at startup
// (They depend on spin and other user-supplied parameters)
real R_GRAV; // Gravitational radius
real R_ISCO; // Innermost stable circular orbit
real cutoff_inner; // Inner integration boundary
real JANSKY_FACTOR; // Factor to scale image output

#endif // CONSTANTS_H
