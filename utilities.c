/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "parameters.h"
#include "constants.h"
#include <stdio.h>

void print_time(int start){
        clock_t diff = clock() - start;
        int msec = diff * 1000 / (CLOCKS_PER_SEC*20);
        printf("Total time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}



void set_constants(real t){
        // Horizon radius for integration cutoff
        real Rh=(1. + sqrt(1. - a * a));
        cutoff_inner = Rh*(1. + horizon_marg); // Cutoff outside or inside BH EH
        R_GRAV = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm
        // Innermost stable circular orbit (ISCO)
        real Z1 = 1. + pow(1. - a * a, (1./3.)) * (pow(1. + a, (1. / 3.)) + pow(1. - a, (1./3.)));
        real Z2 = pow(3. * a * a + Z1 * Z1, 0.5);
        real retro = a < 0. ? 1. : -1.; // 1 for retrograde orbits, -1 for prograde
        R_ISCO = (3. + Z2 + retro * pow((3. - Z1) * (3. + Z1 + 2. * Z2), 0.5));

        // Calibration constant for the spectral irradiance
        // We want Jansky/pixel^2.
        real delta_x, delta_y, d_x,d_y;

#if (NORMCAM)
        d_x     = CAM_SIZE_X * R_GRAV; // Size of image in cm
        delta_x = d_x / source_dist; // it is far away so just a ratio without tan
        d_y     = CAM_SIZE_Y * R_GRAV; // Size of image in cm
        delta_y = d_y / source_dist; // it is far away so just a ratio without tan
#endif


#if (LINEAR_IMPACT_CAM)
        real pixsize = (delta_x / (real) IMG_WIDTH) * (delta_y / (real) IMG_HEIGHT); // Pix size in Sr
#elif (LOG_IMPACT_CAM)
        real pixsize = 1;
#endif

        JANSKY_FACTOR = 1.e23 * pixsize; // 1.e23 is conversion from Jansky to ergs/Sr Hz s cm2

#pragma acc copyin(Rh,cutoff_inner,R_GRAV)
}

void write_image(FILE *imgfile, real **intensityfield,int f, real scalefactor){
        int i, j;

        // Write image to output file
        for(i = 0; i < IMG_WIDTH; i++) {
                for(j = 0; j < IMG_HEIGHT; j++) {
                        real stepx = CAM_SIZE_X / (real) IMG_WIDTH;
                        real stepy = CAM_SIZE_Y / (real) IMG_HEIGHT;
                        real alpha = -CAM_SIZE_X * 0.5 + (i + 0.5) * stepx;
                        real beta  = -CAM_SIZE_Y * 0.5 + (j + 0.5) * stepy;

                        real norm1 = -CAM_SIZE_X * 0.5 + (IMG_WIDTH + 0.5) * stepx;
                        real norm2 = -CAM_SIZE_Y * 0.5 + (IMG_HEIGHT + 0.5) * stepy;

                        fprintf(imgfile, "%+.5e\n",scalefactor * intensityfield[i + j * IMG_WIDTH][f]);
                }
        }
}

void write_VTK_image(FILE *fp, real **intensityfield, int f,real *lambdafield, real scalefactor){
        int i, j;
        real stepx = CAM_SIZE_X / (real) IMG_WIDTH;
        real stepy = CAM_SIZE_Y / (real) IMG_HEIGHT;
        fprintf(fp, "# vtk DataFile Version 2.0\n");
        fprintf(fp, "Image Simulation Result\n");
        fprintf(fp, "ASCII\n");
        fprintf(fp, "DATASET STRUCTURED_GRID\n");
        fprintf(fp, "DIMENSIONS %d %d %d\n", IMG_WIDTH, IMG_HEIGHT, 1);
        fprintf(fp, "POINTS %d float\n", (IMG_WIDTH) *(IMG_HEIGHT));
        for(i = 0; i < IMG_WIDTH; j++)
                for(j = 0; j < IMG_HEIGHT; j++) {
#if (LINEAR_IMPACT_CAM)
                        real xx = -CAM_SIZE_X * 0.5 + (i + 0.5) * stepx;
                        real yy = -CAM_SIZE_Y * 0.5 + (j + 0.5) * stepy;
#elif (LOG_IMPACT_CAM)
                        real r_i = exp(log(20.)*(real)(i+0.5) /(real) IMG_WIDTH) - 1.;
                        real theta_i = 2.*M_PI  * (real)(j+0.5)/ (real)IMG_HEIGHT;

                        real xx = r_i * cos(theta_i);
                        real yy  = r_i * sin(theta_i);
#endif

                        fprintf(fp, "%.20e %.20e %.20e\n", xx, yy, 0.0);
                }
        fprintf(fp, "\nPOINT_DATA %d\n", IMG_WIDTH * IMG_HEIGHT);
        fprintf(fp, "SCALARS Intensity float\n");
        fprintf(fp, "LOOKUP_TABLE default\n");
        real flux=0.0;
        for(i = 0; i < IMG_WIDTH; i++)
                for(j = 0; j < IMG_HEIGHT; j++) {
                        flux += scalefactor * intensityfield[i + j * IMG_WIDTH][f];
                        fprintf(fp, "%+.15e\n", sqrt(scalefactor * intensityfield[i + j * IMG_WIDTH][f]));
                }
        //   fprintf(fp, "%+.15e\n", sqrt(scalefactor * intensityfield[i + j * IMG_WIDTH][f])); // to close the vtk gap
        //   fprintf(fp, "SCALARS lambda float\n");
        //   fprintf(fp, "LOOKUP_TABLE default\n");
        //   for(j = 0; j < IMG_WIDTH; j++)
        //       for(i = 0; i < IMG_WIDTH; i++){
        //           fprintf(fp, "%+.15e\n", lambdafield[i + j * IMG_WIDTH]);
        //       }
        //     fprintf(stderr,"Integrated flux density = %.5e\n", flux);
}
