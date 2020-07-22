
/* 
 * File:   RCARRY.cpp
 * Author: Thomas Bronzwaer
 *
 * Created on February 21, 2014, 3:51 PM
 * 
 * This C++ program is an implementation of the RCARRY pRNG algorithm.
 * WARNING: The init(int seed) function MUST be called before calling genrand().
 * 
 * RCARRY has two functions: 
 * 
 * initialize(int seed) - must be called before generating numbers!
 * genrand()            - generates a double between 0 and 1.
 * 
 * RCARRY was originally created by professor Ronald Kleiss, Radboud University,
 * and my implementation is distributed with his kind permission.
 */

#include <stdlib.h>
#include <stdio.h>

// PARAMETERS
const int s       = 10;
const int r       = 24;
const int B       = 16777216; // Seed must be smaller than B! (WHY??)
//const int seed    = 637553; 
const int STARTUP = 2000; // Discard this many numbers before proper use begins

// STATIC VARIABLES (only accessible to functions in this cpp file)
static int reg[24];
static int carrybit = 0;

// Generate a random number between 0 and 1 (inclusive?)
// WARNING: MUST CALL THE init FUNCTION BEFORE USING genrand!
double genrand_RCARRY(){   
    int y  = 0;
    int xn = 0;   
    
    // Compute new number
    y = reg[s-1] - reg[r-1] - carrybit;

    // Possibly offset the new number, compute carry bit
    if (y >= 0){
        xn = y;
        carrybit = 0;
    }
    else{
        xn = y + B;
        carrybit = 1;
    }

    // Shift the register
    for (int j = r-1; j > 0; j--){
        reg[j] = reg[j-1];
    } 
    reg[0] = xn;
    
    return (double) xn / (double) B;
}

void init_RCARRY(int seed_){
    // Initialize the register & variables    
    for (int i = 0; i < r; i++){
        reg[i] = seed_;
    }  
    // Generate a few hundred numbers to get through the start-up
    double temp;
    for (int k = 0; k < STARTUP; k++){
        temp = genrand_RCARRY();
    }
}


