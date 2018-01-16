/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Mościbrodzka
 *
 */

//#include <math.h>
#include "functions.h"
#include "constants.h"
#include "parameters.h"
#include <stdio.h>
#include "gsl/gsl_sf_hyperg.h"
#include "raptor_harm_model.h"

#define ME  (9.10956e-28)
#define mc2 (8.18726e-07)
#define kb  (1.38e-16)
#define hpl (6.6262e-27)
#define CL  (2.99792458e10)
#define keV (1.602e-9)
#define alphaf  (7.29735e-3)
#define h__mc2  (8.09e-21)
#define SIGMATH (0.665245873e-24)

// Compute the photon frequency in the plasma frame:
real freq_in_plasma_frame(real Uplasma_u[4], real k_d[4]){
        real nu_plasmaframe = 0.;
        int i;
        LOOP_i nu_plasmaframe += Uplasma_u[i] * k_d[i];
        nu_plasmaframe *= -(ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
                          PLANCK_CONSTANT;

        if(nu_plasmaframe<=0.0) {
                nu_plasmaframe = 1.0;
        }
        //    if(nu_plasmaframe<0){

        //    nu_plasmaframe = 1e-5;
        //  }
        return nu_plasmaframe;
}

#pragma acc routine(acos_imp)
float acos_imp(float x) {
        float negate = 1; //x / fabs(x);
        x = abs(x);
        float ret = -0.0187293;
        ret = ret * x;
        ret = ret + 0.0742610;
        ret = ret * x;
        ret = ret - 0.2121144;
        ret = ret * x;
        ret = ret + 1.5707288;
        ret = ret * sqrt(1.0-x);
        ret = ret - 2 * negate * ret;

        return negate*3.14159265358979 + ret;
}

// See eqn 73 in Dexter 2016
real pitch_angle(real *X_u, real *k_d, real *B_u, real *Uplasma_u, real B){
        // Compute and clamp result (result can slightly exceed domain of acos due to numerics)
        real result;

        if(B==0.)
                return 0;
        real k = fabs(k_d[0] * Uplasma_u[0]+k_d[1] * Uplasma_u[1]+ k_d[2] * Uplasma_u[2] + k_d[3] * Uplasma_u[3]);
        result = (k_d[0] * B_u[0] + k_d[1] * B_u[1] + k_d[2] * B_u[2] +
                  k_d[3] * B_u[3]) / (k * B /B_unit);
        if (fabs(result) > 1.)
                result /= fabs(result);


        return result;

}

real radiative_transfer(real *X_u, real *k_u,real dl_current, real *frequency,int icur,real intensity[maxsize][num_indices], real *tau,real **** p){
        int IN_VOLUME=0, path_counter;
        real j_nu      = 0.;
        real B, THETA_e, pitch_ang, beta, nu_p, n_e, nu_p2, Bern;
        int i,f;
        real k_d[4],k_u_f[4], dl_current_f;
        real B_u[4], Uplasma_u[4];
        real Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm
        real g_dd[4][4],g_uu[4][4];

        real a_nu = 0.;

        metric_dd(X_u,g_dd);
        metric_uu(X_u,g_uu);

        if(get_fluid_params(X_u, &n_e, &THETA_e, &B, &beta, &Bern, B_u, Uplasma_u, &IN_VOLUME,p,g_dd,g_uu)) {
                // Obtain pitch angle: still no units (geometric)
                lower_index(X_u,g_dd,k_u,k_d);
                pitch_ang = pitch_angle(X_u, k_d, B_u, Uplasma_u,B);

                // CGS UNITS USED FROM HERE ON OUT
                //////////////////////////////////
                for(f = 0; f < num_indices; f++) {
                        if(tau[f] < log(1000.) ) {

                                real Icurrent = intensity[icur][f];
                                // Scale the wave vector to correct energy
                                LOOP_i k_u_f[i] = k_u[i] * PLANCK_CONSTANT * frequency[f] /
                                                  (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

                                //lower the index of the wavevector
                                lower_index(X_u,g_dd,k_u_f,k_d);

                                // Compute the photon frequency in the plasma frame:
                                nu_p = freq_in_plasma_frame(Uplasma_u, k_d);
                                nu_p2 = nu_p * nu_p;

                                // Obtain emission coefficient in current plasma conditions
                                emission_coeff(B,pitch_ang, THETA_e, nu_p, n_e,&j_nu,&a_nu);
                                dl_current_f = dl_current *  Rg / ( frequency[f]);
                                real dtau  = (nu_p * a_nu * dl_current_f);
                                tau[f] += dtau;
                                real dI =  (j_nu/nu_p2) * exp(-tau[f]) * dl_current_f; //j_nu*exp(-tau[f]);

                                if( tau[f] < log(1000.) && !(tau[f]<0.0)) {
                                        Icurrent+= dI;
                                }
                                intensity[icur][f]=Icurrent;
                        }
                }
                /*
                   if(0 && (j_nu != j_nu || isnan(j_nu))){
                   printf("NaN emissivity! X_u[2], theta = %+.15e %+.15e\n", X_u[2],M_PI*X_u[2] + 2*0.35/M_PI * sin(2*M_PI*X_u[2])*atan2(2.*(log(50.)-X_u[1]),1.));
                   printf("NaN emissivity! expX_u[1], innercutoff = %+.15e %e\n", exp(X_u[1]),cutoff_inner);
                   printf("NaN emissivity! pitch_angle = %+.15e\n", pitch_ang);
                   printf("NaN emissivity! B = %+.15e\n", B);
                   printf("NaN emissivity! THETAe = %+.15e\n", THETA_e);
                   printf("NaN emissivity! nu_plasmaframe = %+.15e\n", nu_p);
                   printf("NaN emissivity! n_e = %+.15e\n", n_e);
                   printf("NaN emissivity! j_nu = %+.15e\n", j_nu);
                   printf("NaN emissivity! tau = %+.15e\n", tau);
                   printf("NaN emissivity! U dot U = %+.15e\n\n", inner_product(X_u, Uplasma_u, Uplasma_u));
                   printf("NaN emissivity! k dot k = %+.15e\n\n", inner_product(X_u, k_u, k_u));

                   }
                 */

        }
        return 1;

}
