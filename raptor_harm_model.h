/*
 * raptor_harm_model.h
 *
 * Please note that most of the code in this file was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 *
 * A copy of the GRMONTY license is included below.
 */

/***********************************************************************************
    Copyright 2013 Joshua C. Dolence, Charles F. Gammie, Monika Mo\'scibrodzka,
                   and Po Kin Leung

                        GRMONTY  version 1.0   (released February 1, 2013)

    This file is part of GRMONTY.  GRMONTY v1.0 is a program that calculates the
    emergent spectrum from a model using a Monte Carlo technique.

    This version of GRMONTY is configured to use input files from the HARM code
    available on the same site.   It assumes that the source is a plasma near a
    black hole described by Kerr-Schild coordinates that radiates via thermal
    synchrotron and inverse compton scattering.

    You are morally obligated to cite the following paper in any
    scientific literature that results from use of any part of GRMONTY:

    Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009,
        Astrophysical Journal Supplement, 184, 387

    Further, we strongly encourage you to obtain the latest version of
    GRMONTY directly from our distribution website:
    http://rainman.astro.illinois.edu/codelib/

    GRMONTY is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    GRMONTY is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GRMONTY; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

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
#pragma acc routine(interp_scalar_2D)
real interp_scalar_2D(real ***var, int i, int j,int k, real coeff[4]);
#pragma acc routine(Xtoij)
void Xtoij(real *X, int *i, int *j, real *del);

#pragma acc routine (lower)
void lower(real *ucon, real Gcov[NDIM][NDIM], real *ucov);

#endif // RAPTOR_HARM_MODEL_H
