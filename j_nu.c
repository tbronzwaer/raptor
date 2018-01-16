/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar,  Monika Mościbrodzka
 *
 */

#include "functions.h"
#include "constants.h"
#include "parameters.h"
#include <stdio.h>
#include "gsl/gsl_sf_hyperg.h"
#include <gsl/gsl_sf_bessel.h>

#define e2_ch (1.161410e-03)

//non thermal emission


void emission_coeff(real B, real theta, real THETA_e, real nu_plasma, real n_e, real *j_nu, real *a_nu){

        *j_nu = emission_coeff_THSYNCH(B, theta,  THETA_e, nu_plasma, n_e);

        *a_nu = absorption_coeff_TH(*j_nu, nu_plasma, THETA_e);

}

real emission_coeff_kappa_FIT(real nu,real Ne, real Thetae, real B,real beta, real theta){
        //emissivity for the kappa distribution function, see Pandya et al. 2016
        real nuc, sth, nus, x,  w, X_kappa, factor;
        real J_low, J_high, J_s,kappa=0.;

        kappa=3.5;
        w = Thetae;
        nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
        sth = sqrt(1-theta*theta); //sin(theta);

        factor = (Ne * pow(ELECTRON_CHARGE, 2.) * nuc * sth)/SPEED_OF_LIGHT;

        nus = nuc * sth * pow( w * kappa,2);
        if(nu==1.0 || fabs(theta)==1 )
                return 0;

        X_kappa = nu / nus;
        J_low = pow(X_kappa,1./3.) * 4 * M_PI * tgamma(kappa - 4./3.) /(pow(3,7./3.) * tgamma(kappa - 2.));
        J_high = pow(X_kappa,-(kappa-2)/2.) * pow(3,(kappa-1)/2.) * (kappa - 2.)*(kappa - 1.)/4. * tgamma(kappa/4. - 1./3.) * tgamma(kappa/4. + 4./3.);
        x = 3 * pow(kappa,-3./2.);

        J_s = pow( ( pow(J_low,-x) + pow(J_high,-x) ), -1./x );

        if(J_s!=J_s) {
                printf("B %e\n",B);
                printf("J_s %e %e %e %e %e %e %e %e %e\n", X_kappa,B,nuc, theta, nus,J_low,J_high,x,kappa);
        }
        return ( J_s*factor);

}

real absorption_coeff_kappa_FIT(real nu, real Ne, real Thetae, real B, real beta, real theta){
        // absortion for the kappa distribution function, see Pandya et al. 2016
        real nuc, sth, nus, x, w, X_kappa, factor,kappa=0.;
        real A_low, A_high, A_s;

        w = Thetae;
        nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
        sth =  sqrt(1-theta*theta);
        kappa=3.5;

        factor = Ne * ELECTRON_CHARGE /(  B *sth);
        if(nu==1.0 || fabs(theta)==1.)
                return 0;
        nus =  nuc * sth *pow( w * kappa,2);

        X_kappa = nu / nus;

        //identity to be able to calculate a hypergeometric function, from the code symphony by Pandya et al. 2016

        real a = kappa - 1./3.;
        real b = kappa + 1.;
        real c = kappa + 2./3.;
        real z = -kappa*w;
        real hyp2F1;
        if(fabs(z)==1.)
                return 0;

        if(fabs(z)>1.) {
                hyp2F1 = pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
                         / (tgamma(b)*tgamma(c-a))
                         * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
                         + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
                         / (tgamma(a) * tgamma(c-b))
                         * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
        }

        else
                hyp2F1 = gsl_sf_hyperg_2F1(a,b,c,z);

        A_low = pow(X_kappa, -5./3.) * pow(3,1./6.) * (10./41.) * pow(2* M_PI,2) / pow(w*kappa, 16./3. - kappa) * (kappa - 2.)*(kappa - 1.) * kappa / (3.*kappa -1.) * tgamma(5./3.) * hyp2F1;
        A_high = pow(X_kappa, -(3.+kappa)/2.) * (pow(M_PI,5./2.)/3.) * ((kappa - 2. )* (kappa - 1.) * kappa /pow(w*kappa,5.)) * (2* tgamma(2.+ kappa/2.)/(2.+kappa) -1.) * (pow(3./kappa,19./4.) + 3./5.);

        x = pow(-7./4. +8. * kappa / 5., -43./50.);

        A_s = pow( ( pow(A_low,-x) + pow(A_high,-x) ), -1./x );

        return (factor*A_s);
}




// Return emissivity j_nu which depends on local plasma parameters
// Ref. Dolence & Moscibrodzka 2009
real emission_coeff_THSYNCH(real B, real theta, real THETA_e, real nu_plasma, real n_e){
        real sth = sqrt(1-theta*theta);
        real nu_c = ELECTRON_CHARGE * B /
                    (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
        real nu_s = (2. / 9.) * nu_c * THETA_e * THETA_e * sth;
        real X    =nu_plasma /( nu_s);
        real f    = pow(pow(X, 0.5) + pow(2., 11. / 12.) * pow(X, 1. / 6.), 2.);
        real K2;
        if(nu_plasma==1) {
                fprintf(stderr,"problems nu plasma equals 1");
                return 0;
        }

        if(THETA_e>1e-4)
                K2 = gsl_sf_bessel_Kn(2, 1. / THETA_e);
        else
                return 0;
        real j_nu = n_e * sqrt(2.) * M_PI *
                    nu_s / (3. * K2) * f *
                    exp(-pow(X,1./3.));
        return j_nu * e2_c;

}

// Return emission constant j_nu as described in Dexter (2009) (the geokerr paper)
real emission_coeff_THSYNCHAV(real B, real THETA_e, real nu_p, real n_e){
        real nu_c = 3. * ELECTRON_CHARGE * B * THETA_e * THETA_e /
                    (4. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
        real x_M = nu_p / nu_c;
        real I = 4.0505 / pow(x_M, 1. / 6.) * (1. + 0.4 / pow(x_M, 1. / 4.) +
                                               0.5316 / sqrt(x_M)) * exp(-1.8899 * pow(x_M, 1. / 3.));
        real j_nu = nu_p * n_e /
                    (2. * sqrt(3.)) * 1. / (THETA_e * THETA_e) * I;

        return j_nu;
}

// Return emission coefficient for thermal free-free radiation
real emission_coeff_FFTHERMAL(real nu, real n_e, real T){
        real n_i = n_e; // Assume neutral hydrogen plasma
        real Z = 1.; // For H, Z = 1
        real g_ff = 1.; // Gaunt factor

        real j_nu = 5.44e-39 * (Z * Z / sqrt(T)) * n_i * n_e * g_ff *
                    exp(-PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T));

        return j_nu;
}

// Return emissivity for the thin disk model described in Dexter 2009
real emissivity_thindisk(real *X_u){
        real r = logscale ? exp(X_u[1]) : X_u[1];
        return 1 / r / r;
}

// Planck function
real planck_function(real nu, real THETA_e){
        real T = THETA_e * ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT / BOLTZMANN_CONSTANT;
        return 2. * nu * nu * nu /(SPEED_OF_LIGHT * SPEED_OF_LIGHT)*1./ (exp(PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T)) - 1.);
}

// Return absorption coefficient - assume local thermodynamical equilibrium so that Kirchoff's Law applies
real absorption_coeff_TH(real j_nu, real nu, real THETA_e){
        //    if(nu==1)
        //        return 0;

        real B_nu = planck_function(nu, THETA_e); // Planck function
        return j_nu *e2_ch/ B_nu/e2_c;
}
