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

#include "raptor_harm_model.h"
#include "constants.h"
#include "functions.h"
#include "parameters.h"

void init_model()
{
        /* find dimensional quantities from black hole
           mass and its accretion rate */
        set_units(M_UNIT);

        fprintf(stderr, "\nHARM2D SIMULATION DATA\n");
        fprintf(stderr, "Reading simulation data...\n");

        init_harm_data(GRMHD_FILE);
}


void init_harm_data(char *fname)
{
        FILE *fp;
        real x[4];

        int i, j, k;

        /* header variables not used except locally */
        real t, tf, cour, DTd, DTl, DTi, dt;
        int nstep, DTr, dump_cnt, image_cnt, rdump_cnt, lim, failed;
        real r, h, divb, vmin, vmax, m,gdet;
        real Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
        fprintf(stderr, "Use slice %s\n", fname);
        fp = fopen(fname, "r");

        if (fp == NULL) {
                fprintf(stderr, "can't open sim data file\n");
                exit(1);
        } else {
                fprintf(stderr, "successfully opened %s\n", fname);
        }

        /* get standard HARM header */
        fscanf(fp, "%lf ", &t);
        fscanf(fp, "%d ", &n1);
        fscanf(fp, "%d ", &n2);
        fscanf(fp, "%lf ", &startx[1]);
        fscanf(fp, "%lf ", &startx[2]);
        fscanf(fp, "%lf ", &dx[1]);
        fscanf(fp, "%lf ", &dx[2]);
        fscanf(fp, "%lf ", &tf);
        fscanf(fp, "%d ", &nstep);
        fscanf(fp, "%lf ", &a);
        fscanf(fp, "%lf ", &gam);
        fscanf(fp, "%lf ", &cour);
        fscanf(fp, "%lf ", &DTd);
        fscanf(fp, "%lf ", &DTl);
        fscanf(fp, "%lf ", &DTi);
        fscanf(fp, "%d ", &DTr);
        fscanf(fp, "%d ", &dump_cnt);
        fscanf(fp, "%d ", &image_cnt);
        fscanf(fp, "%d ", &rdump_cnt);
        fscanf(fp, "%lf ", &dt);
        fscanf(fp, "%d ", &lim);
        fscanf(fp, "%d ", &failed);
        fscanf(fp, "%lf ", &Rin);
        fscanf(fp, "%lf ", &Rout);
        fscanf(fp, "%lf ", &hslope);
        fscanf(fp, "%lf ", &R0);
        if(n1!=N1 || n2!=N2) {
                printf("FATAL ERROR, n1!=N1 or n2!=N2\n");
                exit;
        }


        // nominal non-zero values for axisymmetric simulations
        startx[0] = 0.;
        startx[2]=startx[2]; // *M_PI;
        startx[3] = 0.;
        dx[2]=dx[2]; // *M_PI;
        stopx[0] = 1.;
        stopx[1] = startx[1] + N1 * dx[1];
        stopx[2] = startx[2] + N2 * dx[2];
        stopx[3] = 2. * M_PI;

        fprintf(stderr, "Sim range x1, x2:  %g %g, %g %g\n", startx[1],
                stopx[1], startx[2], stopx[2]);

        dx[0] = 1.;
        dx[3] = 2. * M_PI;

        // Allocate storage for all model size dependent variables
        fprintf(stderr,"memory alloc...");
        init_storage();
        fprintf(stderr,"done\n");


        Thetae_unit = (gam - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + TP_OVER_TE);

        dMact = 0.;
        Ladv = 0.;

        for (k = 0; k < N1 * N2; k++) {
                j = k % N2;
                i = (k - j) / N2;
                fscanf(fp, "%lf %lf %lf %lf", &x[1], &x[2], &r, &h);

                fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
                       &p[KRHO][i][j][0],
                       &p[UU][i][j][0],
                       &p[U1][i][j][0],
                       &p[U2][i][j][0],
                       &p[U3][i][j][0],
                       &p[B1][i][j][0], &p[B2][i][j][0], &p[B3][i][j][0]);
                //#if (metric==MKS )d
                fscanf(fp, "%lf", &divb);

                fscanf(fp, "%lf %lf %lf %lf",
                       &Ucon[0], &Ucon[1], &Ucon[2], &Ucon[3]);
                fscanf(fp, "%lf %lf %lf %lf", &Ucov[0],
                       &Ucov[1], &Ucov[2], &Ucov[3]);
                fscanf(fp, "%lf %lf %lf %lf", &Bcon[0],
                       &Bcon[1], &Bcon[2], &Bcon[3]);
                fscanf(fp, "%lf %lf %lf %lf", &Bcov[0],
                       &Bcov[1], &Bcov[2], &Bcov[3]);
                fscanf(fp, "%lf ", &vmin);
                fscanf(fp, "%lf ", &vmax);
                fscanf(fp, "%lf ", &vmin);
                fscanf(fp, "%lf ", &vmax);
                //#endif
                fscanf(fp, "%lf\n", &gdet);

                // Check accretion rate
                if (i <= 20)
                        dMact += gdet * p[KRHO][i][j][0] * Ucon[1];
                if (i >= 20 && i < 40)
                        Ladv += gdet * p[UU][i][j][0] * Ucon[1] * Ucov[0];
        }

        dMact *= dx[3] * dx[2];
        dMact /= 21.;
        Ladv *= dx[3] * dx[2];
        Ladv /= 21.;
        fprintf(stderr, "dMact: %g Macc: %e [Msun/year], Ladv: %g\n", dMact, dMact * M_UNIT / T_unit / (MSUN / YEAR),Ladv); //0.05279
        fprintf(stderr, "Done reading data\n\n");
        //p[0:NPRIM][0:N1][0:N2]
#pragma acc copyin(Ne_unit,B_unit,U_unit,Thetae_unit,RHO_unit,R0,Rin,Rh,Rout,Rms,hslope,dMact,Ladv)

}

// Get the flud parameters - IN THE PLASMA FRAME?
int get_fluid_params(real X[NDIM], real *Ne,
                     real *Thetae, real *B, real *beta, real * Bern, real Bcon[NDIM], real Ucon[NDIM], int *IN_VOLUME,real **** p,real gcov[NDIM][NDIM],real gcon[NDIM][NDIM])
{

        int i, j,k;
        real del[NDIM];
        real rho, uu;
        real Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
        real coeff[4];
        real Ucov[NDIM], Bcov[NDIM];

        metric_uu(X, gcon);
        metric_dd(X, gcov);

        *IN_VOLUME = 1;

        real smalll = 1.e-20;

        //convert position to index
        Xtoij(X, &i, &j, del);
        k=0; //2D DATA

        //interpolation coefficients
        coeff[0] = (1. - del[1]) * (1. - del[2]);
        coeff[1] = (1. - del[1]) * del[2];
        coeff[2] = del[1] * (1. - del[2]);
        coeff[3] = del[1] * del[2];

        //inteprolate density and internal energy
        rho = interp_scalar_2D(p[KRHO], i, j,k, coeff);
        uu = interp_scalar_2D(p[UU], i, j,k, coeff);

        //Number density
        *Ne = rho * Ne_unit;

        //interpolate primitive B and V
        Bp[1] = interp_scalar_2D(p[B1], i, j,k, coeff);
        Bp[2] = interp_scalar_2D(p[B2], i, j,k, coeff);
        Bp[3] = interp_scalar_2D(p[B3], i, j,k, coeff);

        Vcon[1] = interp_scalar_2D(p[U1], i, j,k, coeff);
        Vcon[2] = interp_scalar_2D(p[U2], i, j,k, coeff);
        Vcon[3] = interp_scalar_2D(p[U3], i, j,k, coeff);

        //Reconstruction of the four-velocity
        VdotV = 0.;
        for (i = 1; i < NDIM; i++)
                for (j = 1; j < NDIM; j++)
                        VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
        Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
        Ucon[0] = -Vfac * gcon[0][0];
        for (i = 1; i < NDIM; i++)
                Ucon[i] = Vcon[i] - Vfac * gcon[0][i];


        lower_index(X, gcov, Ucon, Ucov);

        // Get B and Bcov
        UdotBp = 0.;
        for (i = 1; i < NDIM; i++)
                UdotBp += Ucov[i] * Bp[i];
        Bcon[0] = UdotBp;
        for (i = 1; i < NDIM; i++)
                Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];
        lower_index(X, gcov, Bcon, Bcov);

        //caclulate magnetic field strength
        *B =sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
                 Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]) * B_unit + smalll;
        real Bsq = (Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
                    Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]);

        // plasma beta model for temperature
        real beta_trans=1.;
        *beta = uu*(4./3.-1.) / (0.5 *(Bsq + smalll) *beta_trans);
        real b2=pow(uu*(4./3.-1.) / (0.5 *(Bsq + smalll) *beta_trans),2);

        real trat = R_HIGH * b2/(1. + b2) + R_LOW /(1. + b2);

        real two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat +2.)) + 4./3.);

        Thetae_unit = (two_temp_gam - 1.) * (PROTON_MASS/ELECTRON_MASS) / (1. + trat);

        *Thetae = ((uu)/(rho))*Thetae_unit;

        if(*Thetae < 0 ) {
                *Thetae=smalll;
        }

//some plasma quantaties that could be usefull, lorentz factor, beta (velocity), bernoulli factor
        real lor = 1./sqrt(-gcon[0][0])*Ucon[0];
        real betaf= sqrt(lor*lor-1.)/lor;
        real gam = 4./3.;
        *Bern = -(1.+ uu/rho*gam)*Ucov[0];

// Cant trust high magnetized regions.
        if(Bsq/rho>1.0) {
                *Ne=smalll;
                *Thetae = smalll;
                *B = smalll;
                return 0;
        }


        if (SPHERICAL_ACC) {
                real r, theta;
                r=exp(X[1])+R0;
                *B = sqrt(8. * M_PI * rho / r) * B_unit; // MODIFIED TO ADD B FIELD TO BONDI ACCRETION
        }

        if(0) {
                printf("\nBcon[0] = %+.15e\n", Bcon[0]);
                printf("Bcon[1] = %+.15e\n", Bcon[1]);
                printf("Bcon[2] = %+.15e\n", Bcon[2]);
                printf("Bcon[3] = %+.15e\n", Bcon[3]);

                printf("Ucon[0] = %+.15e\n", Ucon[0]);
                printf("Ucon[1] = %+.15e\n", Ucon[1]);
                printf("Ucon[2] = %+.15e\n", Ucon[2]);
                printf("Ucon[3] = %+.15e\n", Ucon[3]);

                printf("\nBcov[0] = %+.15e\n", Bcov[0]);
                printf("Bcov[1] = %+.15e\n", Bcov[1]);
                printf("Bcov[2] = %+.15e\n", Bcov[2]);
                printf("Bcov[3] = %+.15e\n", Bcov[3]);

                printf("Ucov[0] = %+.15e\n", Ucov[0]);
                printf("Ucov[1] = %+.15e\n", Ucov[1]);
                printf("Ucov[2] = %+.15e\n", Ucov[2]);
                printf("Ucov[3] = %+.15e\n", Ucov[3]);

                printf("B = %+.15e\n", *B);
                printf("Ne = %+.15e\n", *Ne);
                printf("Thetae = %+.15e\n", *Thetae);
        }

        return 1;

}
