/*
 * raptor_harm_model.c
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

#define NDIM	4
#define NPRIM	8

double ****econ;
double ****ecov;
double ***bcon;
double ***bcov;
double ***ucon;
double ***ucov;
double ****p;
double **ne;
double **thetae;
double **b;

int N1, N2, N3;
int n_within_horizon;

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

#define TP_OVER_TE	(0.)

/* some coordinate parameters */
double R0, Rin, Rh, Rout, Rms;
double a;
double hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];
double dlE, lE0;
double gam;
double dMsim;

//double M_unit;
double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;
double Ne_unit;
double Thetae_unit;

#endif // RAPTOR_HARM_MODEL_H
