/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Moscibrodzka, Ziri Younsi
 */

#include "parameters.h"
#include "constants.h"
#include <stdio.h>
#include <math.h>

void set_constants(){
    // Horizon radius for integration cutoff
    double Rh=(1. + sqrt(1. - a * a));
    CUTOFF_INNER = Rh*(1. + horizon_marg); // Cutoff outside or inside BH EH
    R_GRAV = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm

    // Innermost stable circular orbit (ISCO)
    double Z1 = 1. + pow(1. - a * a, (1./3.)) * (pow(1. + a, (1. / 3.)) + pow(1. - a, (1./3.)));
    double Z2 = pow(3. * a * a + Z1 * Z1, 0.5);
    double retro = a < 0. ? 1. : -1.; // 1 for retrograde orbits, -1 for prograde
    R_ISCO = (3. + Z2 + retro * pow((3. - Z1) * (3. + Z1 + 2. * Z2), 0.5));

    // Calibration constant for the spectral irradiance
    // We want Jansky/pixel^2.
    double d_x     = CAM_SIZE_X * R_GRAV; // Size of image in cm
    double delta_x = d_x / source_dist; // it is far away so just a ratio without tan
    double d_y     = CAM_SIZE_Y * R_GRAV; // Size of image in cm
    double delta_y = d_y / source_dist; // it is far away so just a ratio without tan
    double pixsize = (delta_x / (double) IMG_WIDTH) * (delta_y / (double) IMG_HEIGHT); // Pix size in Sr
    JANSKY_FACTOR = 1.e23 * pixsize; // 1.e23 is conversion from Jansky to ergs/Sr Hz s cm2
}

void write_image(FILE *imgfile, double *intensityfield, double scalefactor){
    int i, j;

    // Write image to output file
    for(i = 0; i < IMG_WIDTH; i++){
        for(j = 0; j < IMG_HEIGHT; j++){
//            double stepx = size / (double) width;
//            double stepy = size / (double) height;
//            double alpha = -size * 0.5 + (i + 0.5) * stepx;
//            double beta  = -size * 0.5 + (j + 0.5) * stepy;
            fprintf(imgfile, "%d\t%d\t%+.15e\n", i, j,
                    scalefactor * intensityfield[i + j * IMG_WIDTH]);
        }
    }
}

void write_image_polarized(FILE *imgfile, double *intensityfield, double *f_x_field, double *f_y_field, double *p_field, double scalefactor){
    int i, j;

    // Write image to output file
    for(i = 0; i < IMG_WIDTH; i++){
        for(j = 0; j < IMG_HEIGHT; j++){
            fprintf(imgfile, "%d\t%d\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n", i, j,
                    scalefactor * intensityfield[i + j * IMG_WIDTH], f_x_field[i + j * IMG_WIDTH], f_y_field[i + j * IMG_WIDTH], p_field[i + j * IMG_WIDTH]);
        }
    }
}

void write_image_IQUV(FILE *imgfile, double *Ifield, double *Qfield, double *Ufield, double *Vfield, double scalefactor){
    int i, j;

    // Write image to output file
    for(i = 0; i < IMG_WIDTH; i++){
        for(j = 0; j < IMG_HEIGHT; j++){
            // Note HACK FIX implementation of weird EHT convention, which demands Q -> -Q and U -> -U w.r.t. IEEE convention
            fprintf(imgfile, "%d\t%d\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n", i, j,
                    scalefactor * Ifield[i + j * IMG_WIDTH], scalefactor * Qfield[i + j * IMG_WIDTH], scalefactor * Ufield[i + j * IMG_WIDTH], scalefactor * Vfield[i + j * IMG_WIDTH]);
        }
    }
}


void write_VTK_image(FILE *fp, double *intensityfield, double *lambdafield, double scalefactor){
    int i, j;
    double stepx = CAM_SIZE_X / (double) IMG_WIDTH;
    double stepy = CAM_SIZE_Y / (double) IMG_HEIGHT;
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Image Simulation Result\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", IMG_WIDTH, IMG_HEIGHT, 1);
    fprintf(fp, "POINTS %d float\n", IMG_WIDTH * IMG_HEIGHT);
    for(j = 0; j < IMG_WIDTH; j++)
        for(i = 0; i < IMG_HEIGHT; i++){
            double xx = -CAM_SIZE_X * 0.5 + (i + 0.5) * stepx;
            double yy = -CAM_SIZE_Y * 0.5 + (j + 0.5) * stepy;
            fprintf(fp, "%g %g %g\n", xx, yy, 0.0);
        }
    fprintf(fp, "\nPOINT_DATA %d\n", IMG_WIDTH * IMG_HEIGHT);
    fprintf(fp, "SCALARS Intensity float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    double flux=0.0;
    for(j = 0; j < IMG_WIDTH; j++)
        for(i = 0; i < IMG_HEIGHT; i++){
            flux += scalefactor * intensityfield[i + j * IMG_WIDTH];
            fprintf(fp, "%+.15e\n", scalefactor * intensityfield[i + j * IMG_WIDTH]);
        }
    fprintf(fp, "SCALARS lambda float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(j = 0; j < IMG_WIDTH; j++)
        for(i = 0; i < IMG_WIDTH; i++){
            fprintf(fp, "%+.15e\n", lambdafield[i + j * IMG_WIDTH]);
        }
    fprintf(stdout,"Integrated flux density = %.5e\n", flux);
}
