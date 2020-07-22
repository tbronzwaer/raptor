/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 */

#include <math.h>
#include "functions.h"
#include "constants.h"
#include "parameters.h"
#include <stdio.h>

// Obtain the electron density, dimensionless temperature, and B field strength
void get_plasma_parameters(double X_u[4], double *n_e, double *THETA_e_,
                           double *B, double Uplasma_u[4]){
    *n_e      = electron_density(X_u, Uplasma_u);
    *THETA_e_ = THETA_e(X_u);
    *B        = B_fun(X_u, *n_e);
}

void hotspot_velocity(double *X_u, double *Uplasma_u){
    double omega, r2, r32, rho2, cth2, sth2;//, u_t, l;
    double g_tt, g_tp, g_pp;
    double Ucon_bl[DIM];
    double m1oA2, AA;

    double r  = 6.;//logscale ? exp(X_u[1]) : X_u[1];
    double th = X_u[2];

    //backgound density
    rho2 = r * r;

    //Keplerian velocity of plasma
    sth2  = sin(th) * sin(th);
    cth2  = cos(th) * cos(th);
    r2    = r * r;
    r32   = pow(r, 1.5);
    rho2  = r2 + a * a * cth2;
    g_tt  = -(1. - 2. * r / rho2);
    g_tp  = -2. * a * r * sth2 / rho2;
    g_pp  = (r2 + a * a + 2 * a * a * r * sth2 / rho2) * sth2;
    omega = 1. / (r32 + a);
    m1oA2 = g_tt + 2. * omega * g_tp + omega * omega * g_pp;
    AA    = sqrt(-1. / m1oA2);

    // Formula for Keplerian velocity in BL metric, notice that has to be transformedinto MBL,Ks or MKS
    Ucon_bl[0] = AA;
    Ucon_bl[1] = 0.;
    Ucon_bl[2] = 0.;
    Ucon_bl[3] = AA * omega;

    int i;

    // In case of KS/MKS coords, convert this four-vector into corresponding coordinate system
    #if(metric == KS || metric == MKS)
    #endif

    // Put plasma velocity into u_u[4]
    LOOP_i Uplasma_u[i] = Ucon_bl[i];
}

// Return electron density at position X_u
double electron_density(double X_u[4], double Uplasma_u[4]){
//    double z2;
    double omega, r2, r32, rho2, cth2, sth2;//, u_t, l;
    double g_tt, g_tp, g_pp;
    double Ucon_bl[DIM];
    double m1oA2, AA;
//    double ne_background;

    double spot_ne, spot_ne0, Rspot;
    double xspot[DIM];//curent position of a spot center in KS'
    double th_spot, r_spot;
    double xx;//, ux;
    double P;// X_cov[DIM], xspot_cov[DIM];
    double xc, yc, zc;
    double xs, ys, zs;

    double r  = logscale ? exp(X_u[1]) : X_u[1];
    double th = X_u[2];

    //backgound density
//    z2   = r * r * sin(th) * sin(th);
    rho2 = r * r;
//    ne_background = n_e0 * pow(r, -1.1) * exp(z2 / 2. /rho2);

    //Keplerian velocity of plasma
    sth2  = sin(th) * sin(th);
    cth2  = cos(th) * cos(th);
    r2    = r * r;
    r32   = pow(r, 1.5);
    rho2  = r2 + a * a * cth2;
    g_tt  = -(1. - 2. * r / rho2);
    g_tp  = -2. * a * r * sth2 / rho2;
    g_pp  = (r2 + a * a + 2 * a * a * r * sth2 / rho2) * sth2;
    omega = 1. / (r32 + a);
    m1oA2 = g_tt + 2. * omega * g_tp + omega * omega * g_pp;
    AA    = sqrt(-1. / m1oA2);

    // Formula for Keplerian velocity in BL metric, notice that has to be transformedinto MBL,Ks or MKS
    Ucon_bl[0] = AA;
    Ucon_bl[1] = 0.;
    Ucon_bl[2] = 0.;
    Ucon_bl[3] = AA * omega;

    int i;

    // In case of KS/MKS coords, convert this four-vector into corresponding coordinate system
    #if(metric == KS || metric == MKS)
        double BLcoords[8], KScoords[8];
        LOOP_i{
            BLcoords[i] = X_u[i];
            BLcoords[i+4] = Ucon_bl[i];
        }
        BL_to_KS_u(BLcoords, KScoords);
        LOOP_i Ucon_bl[i] = KScoords[i+4];
    #endif

    // Put plasma velocity into u_u[4]
    LOOP_i Uplasma_u[i] = Ucon_bl[i];

    //double norm = four_velocity_norm(X_u, Uplasma_u);

//if (norm > -0.999)
//printf("\nNorm = %.15e", norm);

    // HOT SPOT PARAMETERS
    //////////////////////
    // period of a spot = 2pi(rspot^1.5+a) [M]
    r_spot  = 6.0;//5.2650; // for a=0.9
    th_spot = 0.5 * M_PI;
    r32     = pow(r_spot, 1.5);
    omega   = 1. / (r32+a);
    P       = 2. * M_PI/omega; // Period of the spot on a Keplerian orbit[M]
    spot_ne = 0.;

    int qq;
    for (qq = 0; qq < nblobs; qq++){
        //spot currrent position
        xspot[0] = X_u[0]; //current coordinate time
        xspot[1] = logscale ? log(r_spot) : r_spot;
        xspot[2] = th_spot;                          //equator 0.5*pi
        xspot[3] = fmod(X_u[0] / P, 1.) * 2. * M_PI + (double) qq /
                   (double) nblobs * 2. * M_PI; //spot current phi at t=X[0]

        // Pseudo-Cartesian coordinates
        xc = sqrt(r * r + a * a) * cos(X_u[3]);
        yc = sqrt(r * r + a * a) * sin(X_u[3]);
        zc = exp(X_u[1]) * cos(X_u[2]);

        xs = sqrt(r_spot * r_spot + a * a) * cos(xspot[3]);
        ys = sqrt(r_spot * r_spot + a * a) * sin(xspot[3]);
        zs = r_spot * cos(xspot[2]);

        //distance^2 between photon position and spot center
        xx = fabs(pow(xc - xs, 2) + pow(yc - ys, 2) + pow(zc - zs, 2));

        spot_ne0 = 3.e6;
        Rspot    = 0.5; //spot size
        // this should be equivalent to the hot spot by Broderick
        spot_ne += spot_ne0 * exp(-(xx) / 2. / Rspot / Rspot);
    }

    return 0.000000000001 + spot_ne;
}

// Return the magnetic field intensity at position X_u
// B strength defined through beta plasma parameter
double B_fun(double X_u[4], double ne){
    double r    = logscale ? exp(X_u[1]) : X_u[1];
    double beta = 10.; //plasma parameter
    return sqrt(8. * M_PI / beta * ne * MPCL2 * 2. / 12. / r);
}

// Return the value for Te at position X_u
double THETA_e(double X_u[4]){
    double r = logscale ? exp(X_u[1]) : X_u[1];
    return THETA_e_0 * pow(r, -0.84);
}
