/*
 * raptor_harm_model.h
 *
 * Please note that most of the code in this file was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 */

#ifndef RAPTOR_HARM_MODEL_H
#define RAPTOR_HARM_MODEL_H

#include <stdio.h>
#include "parameters.h"
#include "functions.h"


int n1, n2, n3;

double Ladv, dMact;

/* mnemonics for primitive vars; conserved vars */
#define KRHO     0
#define UU      1
#define U1      2
#define U2      3
#define U3      4
#define B1      5
#define B2      6
#define B3      7

#define TP_OVER_TE    (0.)

#pragma acc copyin(startx,stopx,dx,Thetae_unit,Ne_unit,B_unit,U_unit,RHO_unit,T_unit,L_unit,R0,Rin,Rh,Rout,Rms,a,hslope)

// RAPTOR_HARM_MODEL.C
//////////////////////

// See grmonty paper by Dolence et al.
// HARM model internal utilities
void set_units(real);
void init_harm_data(char *fname);
void init_storage();
#pragma acc routine(interp_scalar)
real interp_scalar_2D(real ***var, int i, int j,int k, real coeff[4]);
#pragma acc routine(Xtoij)
void Xtoij(real *X, int *i, int *j, real *del);

#pragma acc routine (lower)
void lower(real *ucon, real Gcov[NDIM][NDIM], real *ucov);

#endif // RAPTOR_HARM_MODEL_H
