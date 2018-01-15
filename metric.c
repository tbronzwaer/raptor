/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka
 *
 */
#include "parameters.h"
#include "raptor_harm_model.h"
#include "functions.h"
#include "constants.h"
#include "gsl/gsl_linalg.h"


// Returns the covariant metric g_dd at location X_u
void metric_dd(const real X_u[4], real g_dd[4][4]){
        // Initialize: set all elements to 0
        int i, j;
        LOOP_ij g_dd[i][j] = 0.;
        //memset(g_dd, 0, 16 * sizeof(real));

#if (metric == CAR) //Minkowski metric
        g_dd[0][0] = -1.;
        for (i = 1; i < DIM; i++) g_dd[i][i] = 1.;

#elif (metric == BL || metric == MBL) // (modified) Boyer-Lindquist coordinates

        real r       = logscale ? exp(X_u[1]) : X_u[1];
        real rfactor = logscale ? r : 1.;
        real theta = X_u[2];
        real sint  = sin(theta);
        real cost  = cos(theta);
        real sigma = r * r + a * a * cost * cost;
        real delta = r * r + a * a - 2. * r;
        real A_    = (r * r + a * a) * (r * r + a * a) - delta * a * a *
                     sint * sint;

        // Covariant metric elements
        g_dd[0][0] = -(1. - 2. * r / sigma);
        g_dd[1][1] = sigma / delta * rfactor * rfactor;
        g_dd[2][2] = sigma;
        g_dd[3][3] = A_ / sigma * sint * sint;
        g_dd[0][3] = -2. * a * r * sint * sint / sigma;
        g_dd[3][0] = g_dd[0][3];

#elif (metric == KS) // (Modified) Kerr-Schild metric

        real r   = logscale ? exp(X_u[1]) + R0 : X_u[1];
        real rfactor = logscale ? r : 1.;
        //  real theta = X_u[2];
        real theta = M_PI * X_u[2] + 0.5 * (1. - hslope) * sin(2. * M_PI * X_u[2]);
        real dthdx2  = M_PI * (1. + (1. - hslope) * cos(2. * M_PI * X_u[2]));
        real dthdx22 = dthdx2*dthdx2;
        real cth = cos(theta);
        real sth = sin(theta);
        real cth2 = cth * cth;
        real sth2 = sth * sth;
        real r2 = r*r;
        real a2 =a*a;
        real rho2 = r2 + a2 * cth2;
        real z = 2. * r / rho2;
        real delta = r2 - 2. * r + a2;
        real sigma2 = (r2 + a2) * (r2 + a2) - a2 * delta * sth2;
        real omega = sqrt(sigma2) * sth / sqrt(rho2);

        // Covariant metric elements
        g_dd[0][0] = (z - 1.);
        g_dd[0][1] = z * rfactor;
        g_dd[0][3] = -z * a * sth2;
        g_dd[1][0] = g_dd[0][1];
        g_dd[1][1] = (z + 1.) * rfactor * rfactor;
        g_dd[1][3] = -a * (z + 1.) * sth2 * rfactor;
        g_dd[2][2] = rho2 *dthdx22;
        g_dd[3][0] = g_dd[0][3];
        g_dd[3][1] = g_dd[1][3];
        g_dd[3][3] = omega * omega;


#elif (metric == DM) // Modified Schwarzschild metric with dark matter present)

        real r   = logscale ? exp(X_u[1]) + R0 : X_u[1];
        real theta = X_u[2];
        real Rfactor = 1.e3 * 2.; // 1.e3 * Rs
        real gammafactor = 7. / 3.; // Slope of the NFW profile
        real Afactor = 1. / (1. - 2. / r * (1. + qfactor * pow(r / Rfactor, 3. - gammafactor)));
        real Bfactor = 1. / Afactor;
        real sth = sin(theta);
        real sth2 = sth * sth;

        g_dd[0][0] = -Bfactor;
        g_dd[1][1] = Afactor;
        g_dd[2][2] = r * r;
        g_dd[3][3] = r * r * sth2;

#elif (metric == MKS) // Proper MKS coords

        // TAKEN FROM GRMONTY
        real r     = exp(X_u[1]) + R0;
        real theta = M_PI * X_u[2] + 0.5 * (1. - hslope) * sin(2. * M_PI * X_u[2]);
        /*
           if (fabs(theta) < 1e-10) {
           if ((theta) >= 0)
           theta = 1e-10;
           if ((theta) < 0)
           theta = -1e-10;
           }
           if (fabs(M_PI - (theta)) < 1e-10) {
           if ((theta) >= M_PI)
           theta = M_PI + 1e-10;
           if ((theta) < M_PI)
           theta = M_PI - 1e-10;
           }
         */
        real sinth = sin(theta);
        real sin2th = sinth * sinth;
        real costh = cos(theta);
        real tfac, rfac, hfac, pfac;
        real rho2 = r * r + a * a * costh * costh;

        tfac = 1.;
        rfac = r - R0;
        hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X_u[2]);
        pfac = 1.;

        g_dd[0][0] = (-1. + 2. * r / rho2) * tfac * tfac;
        g_dd[0][1] = (2. * r / rho2) * tfac * rfac;
        g_dd[0][3] = (-2. * a * r * sin2th / rho2) * tfac * pfac;

        g_dd[1][0] = g_dd[0][1];
        g_dd[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
        g_dd[1][3] = (-a * sin2th * (1. + 2. * r / rho2)) * rfac * pfac;

        g_dd[2][2] = rho2 * hfac * hfac;

        g_dd[3][0] = g_dd[0][3];
        g_dd[3][1] = g_dd[1][3];
        g_dd[3][3] =
                sin2th * (rho2 + a * a * sin2th * (1. + 2. * r / rho2)) * pfac * pfac;



#endif
}


// Returns the contravariant metric g_uu at location X
void metric_uu(const real X_u[4], real g_uu[4][4]){
        // Initialize: set all elements to 0
        int i, j;
        LOOP_ij g_uu[i][j] = 0.;
#if (metric == CAR) //Minkowski metric

        g_uu[0][0] = -1.;
        for (i = 1; i < 4; i++) g_uu[i][i] = 1.;

#elif (metric == BL || metric == MBL) // (modified) Boyer-Lindqust coordinates

        real r       = logscale ? exp(X_u[1]) : X_u[1];
        real rfactor = logscale ? r : 1;
        real theta   = X_u[2];
        real sint    = sin(theta);
        real cost    = cos(theta);
        real sigma   = r * r + a * a * cost * cost;
        real delta   = r * r + a * a - 2. * r;
        real A_      = (r * r + a * a) * (r * r + a * a) - delta * a * a *
                       sint * sint;

        // Contravariant metric elements
        g_uu[0][0] = -A_ / (sigma * delta);
        g_uu[1][1] = delta / sigma / rfactor / rfactor;
        g_uu[2][2] = 1. / sigma;
        g_uu[3][3] = (delta - a * a * sint * sint) /
                     (sigma * delta * sint * sint);
        g_uu[0][3] = -2. * a * r / (sigma * delta);
        g_uu[3][0] = g_uu[0][3];

#elif (metric == KS )     // (modified) Kerr-Schild metric

        real r   = logscale ? exp(X_u[1]) + R0 : X_u[1];
        real rfactor = logscale ? r : 1.;
        //  real theta = X_u[2];
        real theta = M_PI * X_u[2] + 0.5 * (1. - hslope) * sin(2. * M_PI * X_u[2]);
        real dthdx2  = M_PI * (1. + (1. - hslope) * cos(2. * M_PI * X_u[2]));
        real dthdx22 = dthdx2*dthdx2;
        real cth = cos(theta);
        real sth = sin(theta);
        real sth2 = sth * sth;
        real cth2 = cth * cth;
        real r2 = r*r;
        real a2 =a*a;
        real rho2 = r2 + a2 * cth2;
        real z = 2. * r / rho2;
        real delta = r2 - 2. * r + a2;
        real sigma2 = (r2 + a2) * (r2 + a2) - a2 * delta * sth2;
        real omega = sqrt(sigma2) * sth / sqrt(rho2);
        real ztilde = omega * omega - (z + 1.) * a2 * sth2 * sth2;

        // Contravariant metric elements
        g_uu[0][0] = -(z + 1.);
        g_uu[0][1] = z / rfactor;
        g_uu[1][0] = g_uu[0][1];
        g_uu[1][1] = (z * z * a2 * sth2 * sth2 - (z - 1.) * omega * omega) / (ztilde * rfactor * rfactor);
        g_uu[1][3] = a * sth2 / (ztilde * rfactor);
        g_uu[3][1] = g_uu[1][3];
        g_uu[2][2] = 1. / (rho2*dthdx22);
        g_uu[3][3] = 1. / (ztilde);

#elif (metric == MKS) // Proper MKS units

        // TAKEN FROM GRMONTY
        real r     = exp(X_u[1]) + R0;
        real theta = M_PI * X_u[2] + 0.5 * (1. - hslope) * sin(2. * M_PI * X_u[2]);

        real sinth = sin(theta);
        real sin2th = sinth * sinth;
        real costh = cos(theta);

        real irho2 = 1. / (r * r + a * a * costh * costh);

        real hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X_u[2]);

        g_uu[0][0] = -1. - 2. * r * irho2;
        g_uu[0][1] = 2. * irho2;

        g_uu[1][0] = g_uu[0][1];
        g_uu[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
        g_uu[1][3] = a * irho2 / r;

        g_uu[2][2] = irho2 / (hfac * hfac);

        g_uu[3][1] = g_uu[1][3];
        g_uu[3][3] = irho2 / (sin2th);

#elif (metric == DM) // Modified Schwarzschild metric with dark matter present)

        real r   = logscale ? exp(X_u[1]) + R0 : X_u[1];
        real theta = X_u[2];
        real Rfactor = 1.e3 * 2.; // 1.e3 * Rs
        real gammafactor = 7. / 3.; // Slope of the NFW profile
        real Afactor = 1. / (1. - 2. / r * (1. + qfactor * pow(r / Rfactor, 3. - gammafactor)));
        real Bfactor = 1. / Afactor;
        real sth = sin(theta);
        real sth2 = sth * sth;

        g_uu[0][0] = 1. / (-Bfactor);
        g_uu[1][1] = 1. / (Afactor);
        g_uu[2][2] = 1. / (r * r);
        g_uu[3][3] = 1. / (r * r * sth2);


#endif
}

// Returns the partial derivative of metric element g_dd[alpha][beta] w.r.t.
// direction 'dir'
real part_deriv_metric_dd(const real X_u[4], int dir, int alpha, int beta){
        // Temporary variable X_u_temp will be modified, leaving original X_u intact
        real X_u_temp[4] = {X_u[0], X_u[1], X_u[2], X_u[3]};
        real g_dd[4][4];

        // Get metric_dd at X + delta(dir), g_ab1[][]...
        X_u_temp[dir] += delta_num;
        metric_dd(X_u_temp, g_dd);
        real plusdelta = g_dd[alpha][beta];

        // ...and X - delta(dir), g_ab2[][]
        X_u_temp[dir] -= 2. * delta_num;
        metric_dd(X_u_temp, g_dd);
        real minusdelta = g_dd[alpha][beta];

        // Return the numerical derivative (central difference rule)
        return (plusdelta - minusdelta) / (2. * delta_num);
}

// Computes the Christoffel symbols at location X numerically
// (Requires the metric to be specified everywhere!)
void connection_num_udd(const real X_u[4], real gamma_udd[4][4][4]){
        // CASE 0: Minkowski metric
        //////////////////////
        int i, j, k;
        LOOP_ijk gamma_udd[i][j][k] = 0.;

        // Obtain metric at current position (contravariant form)
        real g_uu[4][4];
        metric_uu(X_u, g_uu);

        // Solve the Christoffel connection equation
        int alpha;
        for (alpha = 0; alpha < 4; alpha++)
                LOOP_ijk // Index summation over k
                        gamma_udd[alpha][i][j] += 0.5 * g_uu[alpha][k] *
                                                  (part_deriv_metric_dd(X_u, j, k, i) +
                                                   part_deriv_metric_dd(X_u, i, k, j) -
                                                   part_deriv_metric_dd(X_u, k, i, j));
}

// Computes the Christoffel symbols at location X based on an exact metric
// In the form gamma_udd (contr, cov, cov indices)
void connection_udd(const real X_u[4], real gamma[4][4][4]){
        // CASE 0: Minkowski metric
        //////////////////////
        int i, j, k;
//    for(k=0;k<NDIM*NDIM*NDIM;k++){ gamma[0][0][k] = 0.;}
        LOOP_ijk gamma[i][j][k]=0.0;

#if (metric == CAR)

        // Do nothing, connection vanishes

#elif (metric == BL || metric == MBL)  // (modified) Boyer-Lindquist coordinates

        // Compute relevant/often used quantities
        // Depending on whether logscale is active, we interpret X[1] as r or
        // log(r) and use the metric/connection correction factor 'rfactor'
        real r       = logscale ? exp(X_u[1]) : X_u[1];
        real rfactor = logscale ? r : 1;
        real theta   = X_u[2];
        real sint    = sin(theta);
        real cost    = cos(theta);
        real sigma   = r * r + a * a * cost * cost;
        real delta   = r * r + a * a - 2. * r;
        real A       = (r * r + a * a) * (r * r + a * a) - delta * a * a *
                       sint * sint;
        real sigma3  = sigma * sigma * sigma;

        // Unique, non-zero connection elements
        gamma[1][0][0] =  delta / sigma3 * (2. * r * r - sigma) / rfactor;
        gamma[2][0][0] = -2. * a * a * r * sint * cost / sigma3;
        gamma[1][1][1] =  logscale ? ((1. - r) / delta + r / sigma + 1. / r) *
                         r : (1. - r) / delta + r / sigma;
        gamma[2][1][1] =  a * a * sint * cost / (sigma * delta) * rfactor *
                         rfactor;
        gamma[1][2][2] = -r * delta / sigma / rfactor;
        gamma[2][2][2] = -a * a * sint * cost / sigma;
        gamma[1][3][3] = -delta * sint * sint / sigma * (r - a * a * sint *
                                                         sint / (sigma * sigma) * (2. * r * r - sigma)) /
                         rfactor;
        gamma[2][3][3] = -sint * cost / sigma3 * ((r * r + a * a) *
                                                  A - sigma * delta * a * a * sint * sint);
        gamma[0][0][1] =  (r * r + a * a) / (sigma * sigma * delta) *
                         (2. * r * r - sigma) * rfactor;
        gamma[3][0][1] =  a / (sigma * sigma * delta) * (2. * r * r - sigma)
                         * rfactor;
        gamma[0][0][2] = -2. * a * a * r * sint * cost / (sigma * sigma);
        gamma[3][0][2] = -2. * a * r * cost / (sigma * sigma * sint);
        gamma[1][0][3] = -a * delta * sint * sint / sigma3 *
                         (2. * r * r - sigma) / rfactor;
        gamma[2][0][3] =  2. * a * r * (r * r + a * a) * sint * cost /
                         sigma3;
        gamma[1][1][2] = -a * a * sint * cost / sigma;
        gamma[2][1][2] =  r / sigma * rfactor;
        gamma[0][1][3] = -a * sint * sint / (sigma * delta) * (2. * r *
                                                               r / sigma * (r * r + a * a) + r * r - a * a) *
                         rfactor;
        gamma[3][1][3] =  (r / sigma - a * a * sint * sint / (sigma * delta) *
                           (r - 1. + 2. * r * r / sigma)) * rfactor;
        gamma[0][2][3] =  2. * a * a * a * r * sint * sint * sint * cost /
                         (sigma * sigma);
        gamma[3][2][3] =  cost / sint * (1. + 2. * a * a * r * sint * sint /
                                         (sigma * sigma));

        // Take symmetries into account
        gamma[0][1][0] = gamma[0][0][1];
        gamma[3][1][0] = gamma[3][0][1];
        gamma[0][2][0] = gamma[0][0][2];
        gamma[3][2][0] = gamma[3][0][2];
        gamma[1][3][0] = gamma[1][0][3];
        gamma[2][3][0] = gamma[2][0][3];
        gamma[1][2][1] = gamma[1][1][2];
        gamma[2][2][1] = gamma[2][1][2];
        gamma[0][3][1] = gamma[0][1][3];
        gamma[3][3][1] = gamma[3][1][3];
        gamma[0][3][2] = gamma[0][2][3];
        gamma[3][3][2] = gamma[3][2][3];

        // Kerr-Schild coordinates
#elif (metric == KS)

        real r1 = logscale ? exp(X_u[1]) : X_u[1];
        real r2 = r1*r1;
        real r3 = r2*r1;
        real r4 = r3*r1;

        real th = X_u[2];
        real dthdx2 = 1.0;
        real d2thdx22 = 0.0;
        real dthdx22 = dthdx2*dthdx2;
        real sth=sin(th);
        real cth=cos(th);
        real sth2 = sth*sth;
        real r1sth2 = r1*sth2;
        real sth4 = sth2*sth2;
        real cth2 = cth*cth;
        real cth4 = cth2*cth2;
        real s2th = 2.*sth*cth;
        real c2th = 2.*cth2 - 1.;

        real a2 = a*a;
        real a2sth2 = a2*sth2;
        real a2cth2 = a2*cth2;
        real a3 = a2*a;
        real a4 = a3*a;

        real rho2 = r2 + a2cth2;
        real rho22 = rho2*rho2;
        real rho23 = rho22*rho2;
        real irho2 = 1./rho2;
        real irho22 = irho2*irho2;
        real irho23 = irho22*irho2;
        real irho23_dthdx2 = irho23/dthdx2;

        real fac1 = r2 - a2cth2;
        real fac1_rho23 = fac1*irho23;
        real fac2 = a2 + 2.*r2 + a2*c2th;
        real fac3 = a2 + r1*(-2. + r1);

        real fac4 = r2 + a2 * cth2;
        real fac5 = r1 * (r1 + 2.);

        gamma[0][0][0] = 2.*r1*fac1_rho23;
        gamma[0][0][1] = (r2 - a2 * cth2) * (r1 * (2. + r1) + a2 * cth2) / pow(r2 + a2 * cth2, 3.); //
        gamma[0][0][2] = -a2*r1*s2th*dthdx2*irho22;
        gamma[0][0][3] = -2.*a*r1sth2*fac1_rho23;

        gamma[0][1][0] = gamma[0][0][1];
        gamma[0][1][1] = 2. * (r2 - a2 * cth2) * (r1 + r2 + a2 * cth2) / pow(r2 + a2 * cth2, 3.); //
        gamma[0][1][2] = -(2. * a2 * r1 * cth * sth) / pow(r2 + a2 * cth2, 2.); //
        gamma[0][1][3] = a * (-r2 + a2 * cth2) * (r1 * (2. + r1) + a2 * cth2) *sth2 / pow(r2 + a2 * cth2, 3.); //

        gamma[0][2][0] = gamma[0][0][2];
        gamma[0][2][1] = gamma[0][1][2];
        gamma[0][2][2] = -2.*r2*dthdx22*irho2;
        gamma[0][2][3] = a3*r1sth2*s2th*dthdx2*irho22;

        gamma[0][3][0] = gamma[0][0][3];
        gamma[0][3][1] = gamma[0][1][3];
        gamma[0][3][2] = gamma[0][2][3];
        gamma[0][3][3] = 2.*r1sth2*(-r1*rho22 + a2sth2*fac1)*irho23;

        gamma[1][0][0] = (r2 - a2 * cth2) * (0.5 * (a2 + r1 * (r1 - 2.)) * (a4 + 2. * r4 + a2 * r1 * (3. * r1 - 2.) + a2 * (a2 + r1 * (2. + r1)) * cos(2. * th)) - a2 * ((r1 - 2.) * r1 + a2 * cth2) * fac3 * sth2) / (pow(r2 + a2 * cth2, 3.) * (pow(a2 + r2, 2.) - a2 * (r1 * (2. + r1) + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //
        gamma[1][0][1] = 2. * (-a2 + 4. * r1 + a2 * cos(2. * th)) * (a2 - 2. * r2 + a2 * cos(2. * th)) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //
        gamma[1][0][2] = 0.; //
        gamma[1][0][3] = a * (a2 - 2. * r2 + a2 * cos(2. * th)) * sth2 * ((a2 + (r1 - 2.) * r1) * (a4 + 2. * r4 + a2 * r1 * (3. * r1 - 2.) + a2 * (a2 + fac5) * cos(2. * th)) - a2 * (a2 + 2. * (r1 - 2.) * r1 + a2 * cos(2. * th)) * fac3 * sth2) /
                         (4. * pow(r2 + a2 * cth2, 3.) * (pow(a2 + r2, 2.) - a2 * (fac5 + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //

        gamma[1][1][0] = gamma[1][0][1];
        gamma[1][1][1] = 4. * (a2 - 2. * r2 + a2 * cos(2. * th)) * (fac5 + a2 * cos(2. * th)) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

        gamma[1][1][2] = -a2 * sin(2. * th) / (a2 + 2. * r2 + a2 * cos(2. * th)); //
        gamma[1][1][3] = a * (a4 - 8. * a2 * r1 + 3. * a4 * r1 - 4. * a2 * r2 + 16. * r3 + 8. * a2 * r3 + 8. * r4 * r1 + 4. * a2 * r1 * (-2. + a2 + r1 + 2. * r2) * cos(2. * th) + a4 * (r1 - 1.) * cos(4. * th)) * sth2 / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

        gamma[1][2][0] = gamma[1][0][2];
        gamma[1][2][1] = gamma[1][1][2];
        gamma[1][2][2] = -(r1 * ((a2 + r1 * (r1 - 2.)) * (a4 + 2. * r4 + a2 * r1 * (3. * r1 - 2.) + a2 * (a2 + fac5) * cos(2. * th)) - 2. * a2 * (r1 * (r1 - 2.) + a2 * cth2) * fac3 * sth2)) / (2. * fac4 * (pow(a2 + r2, 2.) - a2 * (fac5 + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //
        gamma[1][2][3] = 0.; //

        gamma[1][3][0] = gamma[1][0][3];
        gamma[1][3][1] = gamma[1][1][3];
        gamma[1][3][2] = gamma[1][2][3];
        gamma[1][3][3] = -(a2 + r1 * (r1 - 2.)) * (8. * r4 * r1 + 4. * a2 * r2 * (2. * r1 - 1.) + a4 * (1. + 3. * r1) + 4. * a2 * r1 * (a2 + r1 + 2. * r2) * cos(2. * th) + a4 * (r1 - 1.) * cos(4. * th)) * sth2 / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

        gamma[2][0][0] = -a2*r1*s2th*irho23_dthdx2;
        gamma[2][0][1] = -2. * a2 * r1 * cth * sth / pow(r2 + a2 * cth2, 3.); //
        gamma[2][0][2] = 0.0;
        gamma[2][0][3] = a*r1*(a2+r2)*s2th*irho23_dthdx2;

        gamma[2][1][0] = gamma[2][0][1];
        gamma[2][1][1] = -2. * a2 * r1 * cth * sth / pow(r2 + a2 * cth2, 3.); //
        gamma[2][1][2] = r1 / (r2 + a2 * cth2); //
        gamma[2][1][3] = a * cth * sth * (r3 * (2. + r1) + 2. * a2 * r1 * (1. + r1) * cth2 + a4 * cth4 + 2. * a2 * r1 * sth2) / pow(r2 + a2 * cth2, 3.); //

        gamma[2][2][0] = gamma[2][0][2];
        gamma[2][2][1] = gamma[2][1][2];
        gamma[2][2][2] = -a2*cth*sth*dthdx2*irho2 + d2thdx22/dthdx2;
        gamma[2][2][3] = 0.0;

        gamma[2][3][0] = gamma[2][0][3];
        gamma[2][3][1] = gamma[2][1][3];
        gamma[2][3][2] = gamma[2][2][3];
        gamma[2][3][3] = -cth * sth * (rho23 + a2sth2 * rho2 * (r1 * (4. + r1) +
                                                                a2cth2) + 2. * r1 * a4 * sth4) * irho23_dthdx2;

        gamma[3][0][0] = a * fac1_rho23;
        gamma[3][0][1] = (a * r2 - a3 * cth2) / ((r2 + a2 * cth2) * (pow(a2 + r2, 2.) - a2 * (fac5 + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //
        gamma[3][0][2] = -2. * a * r1 * cth * dthdx2 / (sth * rho22);
        gamma[3][0][3] = -a2sth2 * fac1_rho23;

        gamma[3][1][0] = gamma[3][0][1];
        gamma[3][1][1] = 8. * (a * r2 - a3 * cth2) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //
        gamma[3][1][2] = (a * (fac5 + a2 * cth2) * (1. / tan(th))) / (-pow(a2 + r2, 2.) + a2 * (fac5 + a2 * cth2) * sth2 + a2 * fac3 * sth2); //
        gamma[3][1][3] = (8. * r4 * r1 + 4. * a2 * r2 * (2. * r1 - 1.) + a4 * (1. + 3. * r1) + 4. * a2 * r1 * (a2 + r1 + 2. * r2) * cos(2. * th) + a4 * (r1 - 1.) * cos(4. * th)) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

        gamma[3][2][0] = gamma[3][0][2];
        gamma[3][2][1] = gamma[3][1][2];
        gamma[3][2][2] = -a * r1 * dthdx22 * irho2;
        gamma[3][2][3] = dthdx2 * (0.25 * fac2 * fac2 * cth / sth + a2 * r1 *
                                   s2th) * irho22;

        gamma[3][3][0] = gamma[3][0][3];
        gamma[3][3][1] = gamma[3][1][3];
        gamma[3][3][2] = gamma[3][2][3];
        gamma[3][3][3] = (-a*r1sth2*rho22 + a3*sth4*fac1)*irho23;


#elif (metric == MKS)     //  (Modified) Kerr-Schild metric

        real r           = R0 + exp(X_u[1]);
        real r2          = r * r;
        real rprime      = exp(X_u[1]);
        real rprime2     = rprime * rprime;
        real rprimeprime = rprime;
        real theta       = M_PI * X_u[2] + 0.5 * (1. - hslope) * sin(2. * M_PI * X_u[2]);
        real thetaprime  = M_PI * (1. + (1. - hslope) * cos(2. * M_PI * X_u[2]));
        real thetaprime2 = thetaprime * thetaprime;
        real thetaprimeprime = -2. * M_PI * M_PI * (1. - hslope) * sin(2. * M_PI * X_u[2]);
        real costh       = cos(theta);
        real cos2th      = costh * costh;
        real sinth       = sin(theta);

        real sin2th      = sinth * sinth;
        real sintwoth    = sin(2. * theta);
        real cotth       = 1. / tan(theta);
        real a2          = a * a;

        real Sigma       = r2 + a2 * cos2th;
        real Sigma2      = Sigma * Sigma;
        real Delta       = r2 - 2. * r + a2;
        real A           = Sigma * Delta + 2. * r * (r2 + a2);
        real Sigmabar    = 1. / Sigma;
        real Sigmabar2   = Sigmabar * Sigmabar;
        real Sigmabar3   = Sigmabar * Sigmabar * Sigmabar;
        real B           = (2. * r2 - Sigma) * Sigmabar3;
        real C           = r * Sigmabar - a2 * B * sin2th;

        // Gamma[t][mu][nu]
        gamma[0][0][0] = 2. * r * B;
        gamma[0][0][1] = (Sigma + 2. * r) * B * rprime;
        gamma[0][1][0] = gamma[0][0][1];
        gamma[0][0][2] = -a2 * r * sintwoth * Sigmabar2 * thetaprime;
        gamma[0][2][0] = gamma[0][0][2];
        gamma[0][0][3] = -2. * a * r * B * sin2th;
        gamma[0][3][0] = gamma[0][0][3];
        gamma[0][1][1] = 2. * (Sigma + r) * B * rprime2;
        gamma[0][1][2] = gamma[0][0][2] * rprime;
        gamma[0][2][1] = gamma[0][1][2];
        gamma[0][1][3] = -a * sin2th * gamma[0][0][1];
        gamma[0][3][1] = gamma[0][1][3];
        gamma[0][2][2] = -2. * r2 * Sigmabar * thetaprime2;
        gamma[0][2][3] = -a * sin2th * gamma[0][0][2];
        gamma[0][3][2] = gamma[0][2][3];
        gamma[0][3][3] = -2. * r * C * sin2th;

        // Gamma[r][mu][nu]
        gamma[1][0][0] = Delta * B / rprime;
        gamma[1][0][1] = (Delta - Sigma) * B;
        gamma[1][1][0] = gamma[1][0][1];
        gamma[1][2][0] = 0;
        gamma[1][0][2] = gamma[1][2][0];
        gamma[1][0][3] = -a * Delta * B * sin2th / rprime;
        gamma[1][3][0] = gamma[1][0][3];
        gamma[1][1][1] = rprimeprime / rprime - (2. * Sigma - Delta) * B * rprime;
        gamma[1][1][2] = -a2 * sinth * costh * Sigmabar * thetaprime;
        gamma[1][3][2] = 0;
        gamma[1][2][1] = gamma[1][1][2];
        gamma[1][2][3] = gamma[1][3][2];
        gamma[1][1][3] = a * (r * Sigmabar + (Sigma - Delta) * B) * sin2th;
        gamma[1][3][1] = gamma[1][1][3];
        gamma[1][2][2] = -r * Delta * Sigmabar /rprime * thetaprime2;
        gamma[1][3][3] = -Delta * C * sin2th / rprime;

        // Gamma[theta][mu][nu]
        gamma[2][0][0] = Sigmabar * gamma[0][0][2] / thetaprime2;
        gamma[2][0][1] = gamma[2][0][0] * rprime;
        gamma[2][1][0] = gamma[2][0][1];
        gamma[2][2][0] = 0;
        gamma[2][0][2] = gamma[2][2][0];
        gamma[2][0][3] = a * r * (r2 + a2) * sintwoth * Sigmabar3 / thetaprime;
        gamma[2][3][0] = gamma[2][0][3];
        gamma[2][1][1] = gamma[2][0][0] * rprime2;
        gamma[2][1][2] = r * Sigmabar * rprime;
        gamma[2][2][1] = gamma[2][1][2];
        gamma[2][1][3] = a * sinth * costh * (A + Sigma * (Sigma - Delta)) * Sigmabar3 * rprime / thetaprime;
        gamma[2][3][1] = gamma[2][1][3];
        gamma[2][3][2] = 0;
        gamma[2][2][3] = gamma[2][3][2];
        gamma[2][2][2] = thetaprimeprime / thetaprime + gamma[1][1][2];
        gamma[2][3][3] = -sinth * costh * (Delta * Sigma2 + 2. * r * (r2 + a2) * (r2 + a2)) * Sigmabar3 / thetaprime;

        // Gamma[phi][mu][nu]
        gamma[3][0][0] = a * B;
        gamma[3][0][1] = gamma[3][0][0] * rprime;
        gamma[3][1][0] = gamma[3][0][1];
        gamma[3][0][2] = -2. * a * r * cotth * Sigmabar2 * thetaprime;
        gamma[3][2][0] = gamma[3][0][2];
        gamma[3][0][3] = -a2 * B * sin2th;
        gamma[3][3][0] = gamma[3][0][3];
        gamma[3][1][1] = gamma[3][0][0] * rprime2;
        gamma[3][1][2] = -a * (Sigma + 2. * r) * cotth * Sigmabar2 * rprime * thetaprime;
        gamma[3][2][1] = gamma[3][1][2];
        gamma[3][1][3] = C * rprime;
        gamma[3][3][1] = gamma[3][1][3];
        gamma[3][2][2] = -a * r * Sigmabar * thetaprime2;
        gamma[3][2][3] = (cotth + a2 * r * sintwoth * Sigmabar2) * thetaprime;
        gamma[3][3][2] = gamma[3][2][3];
        gamma[3][3][3] = -a * C * sin2th;
#endif


}

// Initialize a contravariant photon wave vector based on impact parameters
// alpha and beta
// Ref. Cunningham & Bardeen 1973
// We want to pick E, L, Q based on impact params alpha, beta
// Then construct k_u using E, L, Q
// The photons all start at the camera location
void initialize_photon(real alpha, real beta, real photon_u[8], real t_init){

        //    beta *= -1.;
        real mu0 = cos(INCLINATION / 180. * M_PI);
        real Xcam_u[4] = {t_init, logscale ? log(rcam) : rcam, acos(mu0), 90./180. * M_PI};
        real En = 1.;
        real E2 = En * En;
        real ll = -alpha * sqrt(1. - mu0 * mu0);
        real qq = beta * beta + mu0 * mu0 * (alpha * alpha - 1.);
        real L1 = ll * En;
        real Q1 = qq * E2;
        real k_d[4], k_u[4];
        real sinc2 = sin(Xcam_u[2]) * sin(Xcam_u[2]);
        real cosc2 = cos(Xcam_u[2]) * cos(Xcam_u[2]);

        // Covariant wave vector entries are known:
        k_d[0] = -En;
        k_d[3] = L1;
        k_d[2] = sign(beta)*sqrt( fabs(Q1 - L1 * L1 * (cosc2/sinc2) + E2*cosc2));

        // Construct contravariant wave vector k_u using the BL metric
        real r       = logscale ? (rcam) : rcam;
        real rfactor = logscale ? r : 1.;

        real theta = acos(mu0);
        real sint  = sin(theta);
        real cost  = cos(theta);
        real sigma = r * r + a * a * cost * cost;
        real delta = r * r + a * a - 2. * r;
        real A_    = (r * r + a * a) * (r * r + a * a) - delta * a * a *
                     sint * sint;

        real g_dd_11 = sigma / delta * rfactor * rfactor;
        real g_uu_00 = -A_ / (sigma * delta);
        real g_uu_03 = -2. * a * r / (sigma * delta);
        real g_uu_33 = (delta - a * a * sint * sint) /
                       (sigma * delta * sint * sint);
        real g_uu_22 = 1. / sigma;

        k_u[0] = g_uu_00 * k_d[0] + g_uu_03 * k_d[3];
        k_u[3] = g_uu_33 * k_d[3] + g_uu_03 * k_d[0];
        k_u[2] = g_uu_22 * k_d[2];
        k_u[1] = sqrt((-k_u[0]*k_d[0]-k_u[2]*k_d[2]-k_u[3]*k_d[3]) / g_dd_11);

        // Normalize the photon wavevector with cam_freq in Hz
        int i;
        //    LOOP_i k_u[i] *= PLANCK_CONSTANT * cam_freq /
        //                   (ELECTRON_MASS * SPEED_OF_LIGHT*SPEED_OF_LIGHT);
        // Place wave vector into "photon_u"
        photon_u[0] = Xcam_u[0];
        photon_u[1] = Xcam_u[1];
        photon_u[2] = Xcam_u[2];
        photon_u[3] = Xcam_u[3];
        photon_u[4] = k_u[0];
        photon_u[5] = k_u[1];
        photon_u[6] = k_u[2];
        photon_u[7] = k_u[3];

        //    printf("\nFirst we build k in BL coords");
        //    LOOP_i printf("\n%+.15e", photon_u[i]);
        //    LOOP_i printf("\n%+.15e", photon_u[i+4]);

        // Convert k_u to the coordinate system that is currently used
#if (metric == KS  || metric == MKS)
        real KSphoton_u[8];
        BL_to_KS_u(photon_u, KSphoton_u);
        LOOP_i {
                photon_u[i] = KSphoton_u[i];
                photon_u[i+4] = KSphoton_u[i+4];
        }
#endif

#if (metric == MKS )
        photon_u[2] = Xg2_approx_rand(photon_u[2],photon_u[1],1); // We only transform theta - r is already exponential and R0 = 0
        photon_u[6] = Ug2_approx_rand(photon_u[6], photon_u[5],photon_u[2],photon_u[1],1); // We only transform theta - r is already exponential and R0 = 0
#endif

}

// Initialize photon using a simple Euclidean virtual camera consisting of eye point + img plane
void initialize_photon_perspective(real alpha, real beta, real photon_u[8], real t_init){

        real camdist = 140.;

        real x = (logscale ? log(camdist) : camdist); //rcam;
        real y = 0.;
        real z = 0.;

        real plane_dist = 30.;

        real xprime = (logscale ? log(camdist - plane_dist) : camdist - plane_dist); //rcam - 10.;
        real yprime = alpha;
        real zprime = beta;

        real ux = xprime - x;
        real uy = yprime - y;
        real uz = zprime - z;

        ux /= sqrt (ux * ux + uy * uy + uz * uz);
        uy /= sqrt (ux * ux + uy * uy + uz * uz);
        uz /= sqrt (ux * ux + uy * uy + uz * uz);

        real t, r, theta, phi, ut, ur, utheta, uphi;

        t     = 0;
        r     = sqrt(x*x + y*y + z*z);
        theta = acos(z / sqrt(x*x + y*y + z*z)); //atan(sqrt(x*x + y*y) / z); // abs() hack fix! why is theta negative sometimes?
        phi   = asin(y / sqrt(x*x + y*y)); //atan(y / x);

        real rprime     = sqrt(xprime*xprime + yprime*yprime + zprime*zprime);
        real thetaprime = acos(zprime / sqrt(xprime*xprime + yprime*yprime + zprime*zprime)); //atan(sqrt(x*x + y*y) / z); // abs() hack fix! why is theta negative sometimes?
        real phiprime   = asin(yprime / sqrt(xprime*xprime + yprime*yprime)); //atan(y / x);

        ut     = 0.;
        ur     = rprime - r; //ux * sin(theta) * cos(phi)  + uy * sin(theta) * sin(phi) + uz * cos(theta);
        utheta = thetaprime - theta; //ux * cos(theta) * cos(phi) + uy * cos(theta) * sin(phi) - uz * sin(theta);;
        uphi   = phiprime - phi; //-ux * sin(phi) + uy * cos(phi);

        real X_u[4];
        real k_u[4];
        X_u[0] = t;
        X_u[1] = r;
        X_u[2] = theta;
        X_u[3] = phi;
        k_u[0] = ut;
        k_u[1] = ur;
        k_u[2] = utheta;
        k_u[3] = uphi;

        normalize_null(X_u, k_u);

        photon_u[0] = X_u[0];
        photon_u[1] = X_u[1];
        photon_u[2] = X_u[2];
        photon_u[3] = X_u[3];
        photon_u[4] = k_u[0];
        photon_u[5] = -k_u[1];
        photon_u[6] = k_u[2];
        photon_u[7] = k_u[3];

}
