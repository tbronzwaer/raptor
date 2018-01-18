/*
 * This program is based on a simple RNG devised by professor Ronald Kleiss at Radboud University.
 * It was presented in the class notes of his course "Monte Carlo Techniques".
 * Prof. Kleiss has authorized redistribution of this code as long as he is credit as the algorithm's originator.
 * Implementation author: T. Bronzwaer, 2013
 */

#include "functions.h"
#include "parameters.h"
#include "constants.h"

// PARAMETERS
const int s       = 10;
const int r       = 24;
const int B       = 16777216;
const int seed    = 637553;
const int STARTUP = 2000; // Discard this many numbers before proper use begins

// STATIC VARIABLES (only accessible to functions in this cpp file)
int reg[24];
int carrybit = 0;

#pragma acc copyin(s,r,B,seed,STARTUP,reg,carrybit)

// "genrandrcarry" generates a random real number between 0 and 1.
// WARNING: MUST CALL THE "initrcarry" FUNCTION BEFORE USING "genrandrcarry"!

#pragma acc routine(genrandrcarry)
real genrandrcarry(){
        int y = 0;
        int xn = 0;

        // Compute new number
        y = reg[s-1] - reg[r-1] - carrybit;

        // Possibly offset the new number, compute carry bit
        if (y >= 0) {
                xn = y;
                carrybit = 0;
        }
        else{
                xn = y + B;
                carrybit = 1;
        }

        // Shift the register
        for (int j = r-1; j > 0; j--) {
                reg[j] = reg[j-1];
        }
        reg[0] = xn;

        return (real) xn / (real) B;
}

#pragma acc routine(initrcarry)
void initrcarry(int seed_){
        // Initialize the register & variables
        for (int i = 0; i < r; i++) {
                reg[i] = seed_;
        }

        // Generate a few hundred numbers to get through the start-up
        real temp;
        for (int k = 0; k < STARTUP; k++) {
                temp = genrandrcarry();
        }
}
