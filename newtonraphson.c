#include "functions.h"
#include "constants.h"
#include "parameters.h"

// TRANSFORMATION FUNCTIONS
///////////////////////////


// Returns the value of f(Xg2) given some value for Xr2. For the correct Xg2,
// we have f(Xg2) = 0.
#pragma acc routine (f_Xg2)
real f_Xg2(real Xg2, real Xr2,real lr,int init){
        return M_PI * Xg2 + 0.5 * (1. - hslope) * sin(2. * M_PI * Xg2) - Xr2;
}

// Returns the value of f'(Xg2).
#pragma acc routine (f_primed_Xg2)
real f_primed_Xg2(real Xg2,real lr,int init){
        return M_PI + M_PI * (1. - hslope) * cos(2. * M_PI * Xg2);
}

// This function does "one Newton-Raphson step", i.e. it returns the NEW,
// "better" estimate Xg2_1 based on the input estimate Xg2_0.
#pragma acc routine (NR_stepX)
real NR_stepX(real Xg2_0, real Xr2,real lr,int init){
        real fprime = f_primed_Xg2(Xg2_0,lr,init);

        return Xg2_0 - f_Xg2(Xg2_0, Xr2, lr,init) / f_primed_Xg2(Xg2_0,lr,init);
}

// Returns the value of f(Ug2) given some value for Ur2. For the correct Ug2,
// we have f(Ug2) = 0.
#pragma acc routine (f_Ug2)
real f_Ug2(real Ug2, real Ug1, real Ur2, real Xg2, real lr,int init){
        return M_PI * Ug2 * (1. + (1. - hslope) * cos(2. * M_PI * Xg2)) - Ur2;
}

// Returns the value of f'(Ug2).
#pragma acc routine (f_primed_Ug2)
real f_primed_Ug2(real Ug2, real Xg2, real lr,int init){
        return M_PI * (1. + (1. - hslope) * cos(2. * M_PI * Xg2));
}

// This function does "one Newton-Raphson step", i.e. it returns the NEW,
// "better" estimate Ug2_1 based on the input estimate Ug2_0.
#pragma acc routine (NR_stepU)
real NR_stepU(real Ug2_0,real Ug1, real Ur2, real Xg2, real lr,int init){
        real fprime = f_primed_Ug2(Ug2_0, Xg2,lr,init);

        return Ug2_0 - f_Ug2(Ug2_0, Ug1, Ur2, Xg2,lr,init) / f_primed_Ug2(Ug2_0, Xg2, lr,init);
}

// Given the X2 coordinate in RAPTOR's convention, Xr2, we compute and return
// an estimate for the corresponding coordinate in HARM2D's convention, Xg2.
real Xg2_approx_rand(real Xr2, real lr,int init2){
        real Xg2_current = 0.1; // Initial guess; reasonable b/c Xg2 E [0, 1]
        real Xg2_prev = 1.e-15; // Keeps track of previous estimate to converge
        real tolerance = 1.e-9; // Maximum error
        int steps = 0;
        int maxsteps = 1000;
        int count = 0;

        // Main loop
        while (fabs(Xg2_current - Xg2_prev) > tolerance*fabs(Xg2_current) && Xg2_current!=Xg2_prev) {
                Xg2_current =genrandrcarry();
                steps = 0;
                count++;
                while(steps < maxsteps && fabs(Xg2_current - Xg2_prev) > tolerance*fabs(Xg2_current) && Xg2_current!=Xg2_prev) {
                        Xg2_prev = Xg2_current;
                        Xg2_current = NR_stepX(Xg2_current, Xr2, lr,init2);

                        steps++;
                }
        }

        // Clamp output value between 0 and 1
        return fmin(1., fmax(Xg2_current, 0.));
}

// Given the U2 coordinate in RAPTOR's convention, Ur2, we compute and return
// an estimate for the corresponding vector component in HARM2D's convention, Ug2.
real Ug2_approx_rand(real Ur2, real Ug1,real Xg2, real lr,int init2){
        real Ug2_current = 0.1; // Initial guess; reasonable b/c Xg2 E [0, 1]
        real Ug2_prev = 1.e-15; // Keeps track of previous estimate to converge
        real tolerance = 1.e-9; // Maximum error
        int steps = 0;
        int maxsteps = 1000;
        int count = 0;

        // Main loop
        while (fabs(Ug2_current - Ug2_prev) > tolerance*fabs(Ug2_current) && Ug2_current!=Ug2_prev) {
                Ug2_current=genrandrcarry();
                steps = 0;
                count++;
                while(steps < maxsteps && fabs(Ug2_current - Ug2_prev) > tolerance*fabs(Ug2_current)) {
                        Ug2_prev = Ug2_current;
                        Ug2_current = NR_stepU(Ug2_current,Ug1, Ur2, Xg2, lr,init2);
                        steps++;
                }
        }

        return Ug2_current;
}
