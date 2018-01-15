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

//kappa distribution function

#define ME  (9.10956e-28)
#define mc2 (8.18726e-07)
#define kb  (1.38e-16)
#define hpl (6.6262e-27)
#define CL  (2.99792458e10)
#define keV (1.602e-9)
#define alphaf  (7.29735e-3)
#define h__mc2  (8.09e-21)
#define SIGMATH (0.665245873e-24)

#define e2_ch (1.161410e-03)



//non thermal emission
real emission_coeff_kappa_FIT(real nu,real Ne, real Thetae, real B,real beta, real theta){
        //emissivity for the kappa distribution function, see Pandya et al. 2016
        real nuc, sth, nus, x,  w, X_kappa, factor;
        real J_low, J_high, J_s,kappa=0.;
        //  real Rhigh = 2.01;
        //  real Rlow = 2.01;
        real b2 =  beta*beta/0.1;
        // printf("%g\n", b2);
        //  kappa = Rhigh * b2/(1+b2) + Rlow / (1 + b2);
        kappa=3.5;
        w = Thetae; //sqrt(  2./9./kappa *Thetae * Thetae);
        nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
        sth = sqrt(1-theta*theta); //sin(theta);

        factor = (Ne * pow(ELECTRON_CHARGE, 2.) * nuc * sth)/SPEED_OF_LIGHT;

        //      fprintf(stderr,"sinth %g\n", sth);
        nus = nuc * sth * pow( w * kappa,2);
        if(nu==1.0 || fabs(theta)==1 )
                return 0;
        // if (nu > 1.e14 * nus || Thetae > 400. || Thetae < .1)
        //     return (0.000);
        X_kappa = nu / nus;
        //      fprintf(stderr, "X_kappa %g\n", X_kappa);
        J_low = pow(X_kappa,1./3.) * 4 * M_PI * tgamma(kappa - 4./3.) /(pow(3,7./3.) * tgamma(kappa - 2.));
        //      fprintf(stderr, "J_low %g\n", J_low);
        J_high = pow(X_kappa,-(kappa-2)/2.) * pow(3,(kappa-1)/2.) * (kappa - 2.)*(kappa - 1.)/4. * tgamma(kappa/4. - 1./3.) * tgamma(kappa/4. + 4./3.);
        //      fprintf(stderr, "J_high %g\n", J_high );
        x = 3 * pow(kappa,-3./2.);

        J_s = pow( ( pow(J_low,-x) + pow(J_high,-x) ), -1./x );
        //      fprintf(stderr,"J_s %g\n", J_s * factor);
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
        // real Rhigh = 5.;
        // real Rlow = 5.;
        real b2 =  beta*beta;
        //  kappa = Rhigh * b2/(1+b2) + Rlow / (1 + b2);
        w = Thetae; //sqrt(2./9./kappa *Thetae * Thetae);
        nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
        sth =  sqrt(1-theta*theta); //sin(theta);
        kappa=3.5;

        factor = Ne * ELECTRON_CHARGE /(  B *sth);
        if(nu==1.0 || fabs(theta)==1.)
                return 0;
        nus =  nuc * sth *pow( w * kappa,2);
        // if (nu > 1.e14 * nus || Thetae > 400. || Thetae < .1)
        //      return (0.0000);
        X_kappa = nu / nus;
//     gsl_set_error_handler_off();
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
        //  if(hyp2F1!=hyp2F1)
//	printf("%e %e %e\n",z,kappa,w);
        A_low = pow(X_kappa, -5./3.) * pow(3,1./6.) * (10./41.) * pow(2* M_PI,2) / pow(w*kappa, 16./3. - kappa) * (kappa - 2.)*(kappa - 1.) * kappa / (3.*kappa -1.) * tgamma(5./3.) * hyp2F1;
        //      fprintf(stderr, "A_low %g\n", A_low);
        A_high = pow(X_kappa, -(3.+kappa)/2.) * (pow(M_PI,5./2.)/3.) * ((kappa - 2. )* (kappa - 1.) * kappa /pow(w*kappa,5.)) * (2* tgamma(2.+ kappa/2.)/(2.+kappa) -1.) * (pow(3./kappa,19./4.) + 3./5.);

        //      fprintf(stderr, "A_high %g\n", A_high);

        x = pow(-7./4. +8. * kappa / 5., -43./50.);

        //      fprintf(stderr,"factor %g\n",factor);

        A_s = pow( ( pow(A_low,-x) + pow(A_high,-x) ), -1./x );
        //      fprintf(stderr, "A_s %g\n", A_s*factor);
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
        /*
           if(THETA_e >0)
            K2=2*THETA_e*THETA_e;
           else
            K2 =0.;//
         */
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
        //    printf("%lf %lf %lf j nu %.40lf\n",x_M, I,nu_c,j_nu);
        //    j_nu = THETA_e * (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / BOLTZMANN_CONSTANT;
        //    j_nu = n_e;

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

// Return emissivity for the simple Gaussian hot spot model discussed in Dexter 2009.
real emissivity_hotspot(real *X_u){
        real xspot[4];

        real r  = logscale ? exp(X_u[1]) : X_u[1];

        real Rspot = 0.5; // Spot size (radius)

        real r_spot  = 6.0; // Spot orbit radius
        real th_spot = 0.5 * M_PI;
        real r32     = pow(r_spot, 1.5);
        real omega   = 1. / (r32+a);
        real P       = 2. * M_PI/omega; // Period of the spot on a Keplerian orbit[M]

        //spot currrent position
        xspot[0] = X_u[0]; //current coordinate time
        xspot[1] = logscale ? log(r_spot) : r_spot;
        xspot[2] = th_spot;                      //equator 0.5*pi
        xspot[3] = fmod(X_u[0] / P, 1.) * 2. * M_PI + M_PI; //spot current phi at t=X[0]

        // Pseudo-Cartesian coordinates
        real xc = sqrt(r * r + a * a) * cos(X_u[3]);
        real yc = sqrt(r * r + a * a) * sin(X_u[3]);
        real zc = exp(X_u[1]) * cos(X_u[2]);

        real xs = sqrt(r_spot * r_spot + a * a) * cos(xspot[3]);
        real ys = sqrt(r_spot * r_spot + a * a) * sin(xspot[3]);
        real zs = r_spot * cos(xspot[2]);

        //distance^2 between photon position and spot center
        real xx = fabs(pow(xc - xs, 2) + pow(yc - ys, 2) + pow(zc - zs, 2));

        if (xx <= 4.)
                return exp(-(xx) / 2. / Rspot / Rspot);
        return 0.;
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
