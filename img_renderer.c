/*
 * Radboud Polarized Integrator v1.0
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *             ___  ___   ___  __________  ___
 *            / _ \/ _ | / _ \/_  __/ __ \/ _ \
 *           / , _/ __ |/ ___/ / / / /_/ / , _/
 *          /_/|_/_/ |_/_/    /_/  \____/_/|_|
 *
 * This program integrates the equations of motion of General Relativity
 * to compute the trajectories of photon bundles (null geodesics); it then
 * performs radiative transfer along these geodesics to compute an image
 * or spectrum. The gravitational field is defined by the metric selected
 * by the user; plasma models can be supplied in the form of GRMHD
 * simulations or analytic models.
 *
 * CONVENTIONS:
 *
 * (Null) geodesics are parametrized by an (affine) parameter called lambda.
 *
 * Metric sign convention: (-,+,+,+)
 *
 * Indices are labeled: "u" (up)   - contravariant index
 *                      "d" (down) - covariant index
 *
 * Examples: U_u[alpha], U_d[alpha], gamma_udd[mu][alpha][beta]
 *
 * A 'ray' (photon bundle position and wave vector) is represented as:
 *
 * photon_u[0] = X_u[0] // Position
 * photon_u[1] = X_u[1]
 * photon_u[2] = X_u[2]
 * photon_u[3] = X_u[3]
 * photon_u[4] = U_u[0] // Wavevector
 * photon_u[5] = U_u[1]
 * photon_u[6] = U_u[2]
 * photon_u[7] = U_u[3]
 *
 * Indices 0, 1, 2, 3 correspond to t, r, theta, phi (Schwarzschild/Kerr).
 */

#include "functions.h"
#include "constants.h"
#include "parameters.h"

int main(int argc, char *argv[]){
        // INPUT FILE
        /////////////
        clock_t start = clock();

        // INITIALIZE MODEL
        ///////////////////
        read_model(argv);

        // LOAD BACKGROUND IMAGE
        ////////////////////////
#if (BACKGROUND_PNG)
        read_png_file("VR_BG_2.png");
#endif

        // read in camera position and velocity
#if (VRCAM)
        camera_trajectory(TIME_INIT);

#endif
        // Initialize HARM2D grmhd model
        // Note: this sets the black hole spin 'a'
        init_model();

        // Set constants such as R_ISCO, JANSKY_FACTOR
        // These depend on the black hole spin
        set_constants(TIME_INIT);
        fprintf(stderr,"\nDone with initialization\n");
        print_time(start);
        // INITIALIZE DATA STRUCTURES
        /////////////////////////////
        real energy_spectrum[num_indices];
        real frequencies[num_indices];
        real **intensityfield;
        intensityfield = (real**)malloc(((IMG_WIDTH) *(IMG_HEIGHT))*sizeof(real*));

        for(int f = 0; f < num_indices; f++) { // For all frequencies...
                frequencies[f] = FREQ_MIN * pow(10., (real) f / (real) FREQS_PER_DEC);
                energy_spectrum[f] = 0.0;
        }
        for(int i = 0; i < IMG_HEIGHT*IMG_WIDTH; i++) {
                intensityfield[i]=(real*)malloc(num_indices*sizeof(real));
                for(int f = 0; f < num_indices; f++) {
                        intensityfield[i][f]=0;
                }
        }

        // MAIN PROGRAM LOOP
        ////////////////////

        //CALCULATE IMAGE
        calculate_image(intensityfield,energy_spectrum, frequencies);
        fprintf(stderr,"\nDone with ray tracing\n");
        print_time(start);
        // WRITE OUTPUT FILES
        /////////////////////

        output_files(intensityfield,energy_spectrum,frequencies);
        fprintf(stderr,"\nDone with output\n");
        print_time(start);
        free(intensityfield);

        // END OF PROGRAM
        /////////////////

        return 0;
}
