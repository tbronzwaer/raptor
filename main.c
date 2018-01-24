/***********************************************************************************
    Copyright: 2014-2020 Black Hole Cam (ERC Synergy Grant 610058)
    Authors:   Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka,
               Ziri Younsi, Heino Falcke, Michael Kramer, and Luciano Rezzolla

               RAPTOR  version 1.0   (released January 28, 2018)

    This file is part of RAPTOR, a program that computes the paths traversed
    by light rays in arbitrary curved spacetimes, and then performs radiative
    transfer calculations along these paths. RAPTOR was developed specifically
    to produce synthetic observational data of accreting supermassive black
    holes for BlackHoleCam/EHT.

    RAPTOR is free to use under the condition that any scientific literature
    resulting from the use of any part of RAPTOR cites the following paper:

    [1] Bronzwaer, T., Davelaar, J., "RAPTOR I: Time-dependent radiative
        transfer in arbitrary spacetimes", A&A

    We strongly encourage you to obtain the latest version of RAPTOR directly
    from our repository:
    www.github.com/tbronzwaer/raptor

    A detailed readme containing operating instructions is also included there.

    RAPTOR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    RAPTOR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RAPTOR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*
 * Radboud Polarized Integrator v1.0
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant 610058)
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
 *
 * HOW TO RUN: QUICKSTART EXAMPLE
 *
 * To compile RAPTOR for execution on a CPU, and to create an image of the included
 * GRMHD dump file, do the following:
 *
 * -Compile using "make harm CPU=1"
 * -Run using ./RAPTOR model.in dump040 1e19 60 1 1 0
 */

#include "functions.h"
#include "constants.h"
#include "parameters.h"

int main(int argc, char *argv[])
{
        fprintf(stderr, "WELCOME TO RAPTOR 1.0\n\n");

        clock_t start = clock();

        // INITIALIZE MODEL
        ///////////////////

        read_model(argv);

        // INITIALIZE GRMHD MODEL
        /////////////////////////

        init_model();

        // Set constants such as R_ISCO, JANSKY_FACTOR
        //////////////////////////////////////////////

        set_constants(TIME_INIT);

        fprintf(stderr,"\nFinished initialization.\n");
        print_time(start);

        // INITIALIZE DATA STRUCTURES
        /////////////////////////////

        real energy_spectrum[num_indices];
        real frequencies[num_indices];
        real **intensityfield;
        intensityfield = (real**)malloc(((IMG_WIDTH) *(IMG_HEIGHT))*sizeof(real*));

        // Frequencies, spectrum
        for(int f = 0; f < num_indices; f++)
        { // For all frequencies...
                frequencies[f] = FREQ_MIN * pow(10., (real) f / (real) FREQS_PER_DEC);
                energy_spectrum[f] = 0.0;
        }

        // "intensityfield", which holds the image data
        for(int i = 0; i < IMG_HEIGHT*IMG_WIDTH; i++)
        {
                intensityfield[i]=(real*)malloc(num_indices*sizeof(real));
                for(int f = 0; f < num_indices; f++)
                        intensityfield[i][f]=0;
        }

        // CALCULATE IMAGE
        //////////////////

        calculate_image(intensityfield,energy_spectrum, frequencies);

        fprintf(stderr,"\nFinished ray tracing calculations.\n");
        print_time(start);

        // WRITE OUTPUT FILES
        /////////////////////

        output_files(intensityfield,energy_spectrum,frequencies);

        fprintf(stderr,"\nFinished writing output files.\n");
        print_time(start);

        // FREE ALLOCATED POINTERS
        //////////////////////////

        free(intensityfield);

        // END OF PROGRAM
        /////////////////

        return 0;
}
