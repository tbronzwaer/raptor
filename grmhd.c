
/*
 * grmhd.c
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


void set_units(real M_unit_)
{
        /** from this, calculate units of length, time, mass,
           and derivative units **/
        L_unit = GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
        T_unit = L_unit / SPEED_OF_LIGHT;

        fprintf(stderr, "\nUNITS\n");
        fprintf(stderr, "L,T,M: %g %g %g\n", L_unit, T_unit, M_unit_);

        RHO_unit = M_unit_ / pow(L_unit, 3);
        U_unit = RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
        B_unit = SPEED_OF_LIGHT * sqrt(4. * M_PI * RHO_unit);

        fprintf(stderr, "rho,u,B: %g %g %g\n", RHO_unit, U_unit, B_unit);

        Ne_unit = RHO_unit / (PROTON_MASS + ELECTRON_MASS);
}

real interp_scalar_2D(real ***var, int i, int j,int k, real coeff[4])
{
        real interp;

        interp =
                var[i][j][k] * coeff[0] +
                var[i][j + 1][k] * coeff[1] +
                var[i + 1][j][k] * coeff[2] + var[i + 1][j + 1][k] * coeff[3];

        return interp;
}

void Xtoij(real X[NDIM], int *i, int *j, real del[NDIM])
{
        *i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 );
        *j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 );

        if (*i < 0) {
                *i = 0;
                del[1] = 0.;
        } else if (*i > N1 - 2) {
                *i = N1 - 2;
                del[1] = 1.;
        } else {
                del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
        }

        if (*j < 0) {
                *j = 0;
                del[2] = 0.;
        } else if (*j > N2 - 2) {
                *j = N2 - 2;
                del[2] = 1.;
        } else {
                del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
        }

        return;
}


// ALLOCATION STUFF BELOW HERE
//////////////////////////////

void init_storage()
{
        int i;

        p = (real****)malloc(NPRIM*sizeof(real***)); //malloc_rank1(NPRIM, sizeof(real *));
        for (i = 0; i < NPRIM; i++) {
                p[i] = (real ***) malloc((N1+1)*sizeof(real**)); //malloc_rank2_cont(N1, N2);
                for(int j =0; j<=N1; j++) {
                        p[i][j]=(real**)malloc((N2+1)*sizeof(real*));
                        for(int k =0; k<=N2; k++) {
                                p[i][j][k]=(real*)malloc((1)*sizeof(real));
                        }
                }
        }

        return;
}
