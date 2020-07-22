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


#include "raptor_harm3d_model.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "functions.h"

/* HDF5 v1.8 API */

#include "parameters.h"

void init_model()
{
    /* find dimensional quantities from black hole
     mass and its accretion rate */
    set_units(M_UNIT);

    fprintf(stderr, "getting simulation data...\n");

    init_harm3d_data(GRMHD_FILE);
}


void init_harm3d_data(char *fname)
{
	int nghost;
	double th_end,th_cutout,two_temp_gam;
 int i,j,k,l,m;
 double X[NDIM],UdotU,ufac;
double gcov[NDIM][NDIM],gcon[NDIM][NDIM], g;
double Ucon[NDIM],Ucov[NDIM];
double dV,Thetae,V,Be;
double Th_unit;
double UdotBp,Bcon[NDIM],Bcov[NDIM],Bp[NDIM];
double bsq,beta,beta_trans,b2;
double trat;

FILE *fp;

    fprintf(stderr, "%s\n", fname);
    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "can't open sim data file\n");
        exit(1);
    } else {
        fprintf(stderr, "successfully opened %s\n", fname);
    }

    /* get standard HARM header */
    fscanf(fp, "%d ", &N1);
    fscanf(fp, "%d ", &N2);
    fscanf(fp, "%d ", &N3);
    fscanf(fp, "%lf ", &gam);
    fscanf(fp, "%lf ", &a);
    fscanf(fp, "%lf ", &dx[1]);
    fscanf(fp, "%lf ", &dx[2]);
    fscanf(fp, "%lf ", &dx[3]);
    fscanf(fp, "%lf ", &startx[1]);
    fscanf(fp, "%lf ", &startx[2]);
    fscanf(fp, "%lf ", &startx[3]);
    fscanf(fp, "%lf ", &hslope);
    fscanf(fp, "%lf ", &Rin);
    fscanf(fp, "%lf ", &Rout);


	stopx[0] = 1.;
	stopx[1] = startx[1]+N1*dx[1];
	stopx[2] = startx[2]+N2*dx[2];
	stopx[3] = startx[3]+N3*dx[3];
	fprintf(stderr, "phi limits is %e %e\n",startx[3],stopx[3]);

	init_storage();

	for(int i=0;i<N1;i++){
		for(int j=0;j<N2;j++){
			for(int k=0;k<N3;k++){
			        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
			                &p[KRHO][i][j][k],
			                &p[UU][i][j][k],
			                &p[U1][i][j][k],
			                &p[U2][i][j][k],
			                &p[U3][i][j][k],
			                &p[B1][i][j][k],
					&p[B2][i][j][k],
					&p[B3][i][j][k]);
			}
		}
	}

	fprintf(stderr,"%e\n",p[KRHO][N1-1][N2-1][N3-1]);



	fprintf(stderr,"done reading!\n");

//	exit(1);
}


// TRANSFORMATION FUNCTIONS
///////////////////////////

// WARNING: these are not yet listed in functions.h and are only meant for use
// by other functions in this file.

// Returns the value of f(Xg2) given some value for Xr2. For the correct Xg2,
// we have f(Xg2) = 0.
double f_Xg2(double Xg2, double Xr2){
    return M_PI * Xg2 + 0.5 * (1. - hslope) * sin(2. * M_PI * Xg2) - Xr2;
}

// Returns the value of f'(Xg2).
double f_primed_Xg2(double Xg2){
    return M_PI + M_PI * (1. - hslope) * cos(2. * M_PI * Xg2);
}

// This function does "one Newton-Raphson step", i.e. it returns the NEW,
// "better" estimate Xg2_1 based on the input estimate Xg2_0.
double NR_stepX(double Xg2_0, double Xr2){
    double fprime = f_primed_Xg2(Xg2_0);

    if(fabs(fprime) < 1.e-9)
        printf("fprime = %+.15e\n", fprime);

    return Xg2_0 - f_Xg2(Xg2_0, Xr2) / f_primed_Xg2(Xg2_0);
}

// Returns the value of f(Ug2) given some value for Ur2. For the correct Ug2,
// we have f(Ug2) = 0.
double f_Ug2(double Ug2, double Ur2, double Xg2){
    return M_PI * Ug2 * (1. + (1. - hslope) * cos(2. * M_PI * Xg2)) - Ur2;
}

// Returns the value of f'(Ug2).
double f_primed_Ug2(double Ug2, double Xg2){
    return M_PI * (1. + (1. - hslope) * cos(2. * M_PI * Xg2));
}

// This function does "one Newton-Raphson step", i.e. it returns the NEW,
// "better" estimate Ug2_1 based on the input estimate Ug2_0.
double NR_stepU(double Ug2_0, double Ur2, double Xg2){
    double fprime = f_primed_Ug2(Ug2_0, Xg2);

    if(fabs(fprime) < 1.e-9)
        printf("fprime = %+.15e\n", fprime);

    return Ug2_0 - f_Ug2(Ug2_0, Ur2, Xg2) / f_primed_Ug2(Ug2_0, Xg2);
}

// Given the X2 coordinate in RAPTOR's convention, Xr2, we compute and return
// an estimate for the corresponding coordinate in HARM2D's convention, Xg2.
double Xg2_approx_rand(double Xr2){
    double Xg2_current = 0.1; // Initial guess; reasonable b/c Xg2 E [0, 1]
    double Xg2_prev = 1.e-15;     // Keeps track of previous estimate to converge
    double tolerance = 1.e-9; // Maximum error
    int steps = 0;
    int maxsteps = 100;

    int count = 0;

    // Main loop
    while (fabs(Xg2_current - Xg2_prev) > tolerance){
        Xg2_current = (double) rand() / (double)RAND_MAX;
        //Xg2_current = 1.e-16;
        steps = 0;
        count++;

        while(steps < maxsteps && fabs(Xg2_current - Xg2_prev) > tolerance){
            Xg2_prev = Xg2_current;
            Xg2_current = NR_stepX(Xg2_current, Xr2);
            steps++;
        }
    }

    // Clamp output value between 0 and 1
    return fmin(1., fmax(Xg2_current, 0.));
}

// Given the U2 coordinate in RAPTOR's convention, Ur2, we compute and return
// an estimate for the corresponding vector component in HARM2D's convention, Ug2.
double Ug2_approx_rand(double Ur2, double Xg2){
    double Ug2_current = 0.1; // Initial guess; reasonable b/c Xg2 E [0, 1]
    double Ug2_prev = 1.e-15;     // Keeps track of previous estimate to converge
    double tolerance = 1.e-9; // Maximum error
    int steps = 0;
    int maxsteps = 100;

    int count = 0;

    // Main loop
    while (fabs(Ug2_current - Ug2_prev) > tolerance){
        Ug2_current = (double) rand() / (double)RAND_MAX;
        steps = 0;
        count++;

        while(steps < maxsteps && fabs(Ug2_current - Ug2_prev) > tolerance){
            Ug2_prev = Ug2_current;
            Ug2_current = NR_stepU(Ug2_current, Ur2, Xg2);
            steps++;
        }
    }

    return Ug2_current;
}


// Current metric: modified Kerr-Schild, squashed in theta
// to give higher resolution at the equator


#define DLOOP  for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)

/* mnemonics for dimensional indices */
#define TT      0
#define RR      1
#define TH      2
#define PH      3

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{

	*r = exp(X[1]) + R0;
	*th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

	return;
}

void gcon_func(double *X, double gcon[][NDIM])
{

	int k, l;
	double sth, cth, irho2;
	double r, th;
	double hfac;
	/* required by broken math.h */
	void sincos(double in, double *sth, double *cth);

	DLOOP gcon[k][l] = 0.;

	bl_coord(X, &r, &th);

	sincos(th, &sth, &cth);
	sth = fabs(sth) + 1.e-9;

	irho2 = 1. / (r * r + a * a * cth * cth);

	// transformation for Kerr-Schild -> modified Kerr-Schild
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);

	gcon[TT][TT] = -1. - 2. * r * irho2;
	gcon[TT][1] = 2. * irho2;

	gcon[1][TT] = gcon[TT][1];
	gcon[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
	gcon[1][3] = a * irho2 / r;

	gcon[2][2] = irho2 / (hfac * hfac);

	gcon[3][1] = gcon[1][3];
	gcon[3][3] = irho2 / (sth * sth);
}


void gcov_func(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;
	/* required by broken math.h */
	void sincos(double th, double *sth, double *cth);

	DLOOP gcov[k][l] = 0.;

	bl_coord(X, &r, &th);

	sincos(th, &sth, &cth);
	sth = fabs(sth) + 1.e-9;
	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
	pfac = 1.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}

#undef TT
#undef RR
#undef TH
#undef PH


void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov)
{

	ucov[0] = Gcov[0][0] * ucon[0]
	    + Gcov[0][1] * ucon[1]
	    + Gcov[0][2] * ucon[2]
	    + Gcov[0][3] * ucon[3];
	ucov[1] = Gcov[1][0] * ucon[0]
	    + Gcov[1][1] * ucon[1]
	    + Gcov[1][2] * ucon[2]
	    + Gcov[1][3] * ucon[3];
	ucov[2] = Gcov[2][0] * ucon[0]
	    + Gcov[2][1] * ucon[1]
	    + Gcov[2][2] * ucon[2]
	    + Gcov[2][3] * ucon[3];
	ucov[3] = Gcov[3][0] * ucon[0]
	    + Gcov[3][1] * ucon[1]
	    + Gcov[3][2] * ucon[2]
	    + Gcov[3][3] * ucon[3];

	return;
}

// Get the fluid parameters - IN THE PLASMA FRAME?
void get_fluid_params(double X[NDIM], double *Ne,
                      double *Thetae, double *B, double Bcon[NDIM], double Ucon[NDIM], int *IN_VOLUME)
{

 int i, j, k;
     double del[NDIM];
     double rho, uu;
     double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
     double gcon[NDIM][NDIM],gcov[NDIM][NDIM],Bcov[NDIM],Ucov[NDIM], coeff[4];
     double bsq,beta,beta_trans,b2,trat,two_temp_gam,Th_unit, Be;
     double Rlow=1,Rhigh=1;
     if (X[1] < startx[1] ||
         X[1] > stopx[1]  ||
         X[2] < startx[2] ||
         X[2] > stopx[2]) {
         *Ne = 0.;
	 *IN_VOLUME=0;
         return;
     }
     *IN_VOLUME=1;

     Xtoijk(X, &i, &j,&k, del);

    metric_uu(X, gcon);
    metric_dd(X, gcov);


     coeff[1]=del[1];
     coeff[2]=del[2];
     coeff[3]=del[3];


//now interpolated to geodesic location 
     // interpolate (cubiclinear interp.) like in mibothros
     rho= interp_scalar(p[KRHO], i, j, k, coeff);
     uu = interp_scalar(p[UU], i, j, k, coeff);
     *Ne = rho * Ne_unit + 1e-40;


     //here unlike in mibothros it was interpolating scalars and
     //reconstructing velocity and magnetic field based on interpolated coefficients
     Bp[1] = interp_scalar(p[B1], i, j, k, coeff);
     Bp[2] = interp_scalar(p[B2], i, j, k, coeff);
     Bp[3] = interp_scalar(p[B3], i, j, k, coeff);

     Vcon[1] = interp_scalar(p[U1], i, j, k, coeff);
     Vcon[2] = interp_scalar(p[U2], i, j, k, coeff);
     Vcon[3] = interp_scalar(p[U3], i, j, k, coeff);

    //reconstrueren van de 4 vectoren
    // Get Ucov
    VdotV = 0.;
    for (i = 1; i < NDIM; i++)
        for (j = 1; j < NDIM; j++)
            VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
    Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
    Ucon[0] = -Vfac * gcon[0][0];
    for (i = 1; i < NDIM; i++)
        Ucon[i] = Vcon[i] - Vfac * gcon[0][i];

    lower_index(X, Ucon, Ucov);
//    lower(Ucon, gcov, Ucov); // Gammie's lowering function

   double Utot=0;
    for(int i =0;i<NDIM;i++)
        Utot+=Ucon[i]*Ucov[i];

    // Get B and Bcov
    UdotBp = 0.;
    for (i = 1; i < NDIM; i++)
        UdotBp += Ucov[i] * Bp[i];
    Bcon[0] = UdotBp;
    for (i = 1; i < NDIM; i++)
        Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];

//    lower_index(X, Bcon, Bcov);
    lower(Bcon, gcov, Bcov); // Gammie's lowering function

     bsq=Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
     Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3];


     *B = sqrt(bsq) * B_unit +1e-40;

     /*electron temperature depending on the plasma magnetization*/
     beta=uu*(gam-1.)/0.5/bsq;
     Be=(-(1.+gam * uu/(rho))*Ucov[0]);
     beta_trans=1.;
     b2=pow(beta/beta_trans,2);
     trat =3.;// Rhigh * b2/(1. + b2) + Rlow /(1. + b2);
     two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
     Th_unit = (1.4444444444 - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + trat);
     *Thetae=(2./15.)*(uu/rho)* (PROTON_MASS / ELECTRON_MASS) +1e-40;
      Be=(-(1.+two_temp_gam * uu/rho)*Ucov[0]);
    //if(bsq/rho>0.15){
     if(uu<0)
                fprintf(stderr,"U %e %e\n",uu,p[UU][i][j][k]);; 

     if(*Thetae<0)
		fprintf(stderr,"Te %e\n",*Thetae); 
     if(*B<0)
		fprintf(stderr,"B %e\n",*B); 
     if(*Ne<0)
		fprintf(stderr,"Ne %e %e\n",*Ne,p[KRHO][i][j][k]); 

     if(bsq/rho > 1.||  exp(X[1])>50.){
        *Ne =0;
	 *IN_VOLUME=0;

     }

}


void set_units(double M_unit_)
{
//	double MBH;

	/* set black hole mass */
	/** could be read in from file here,
	    along with M_unit and other parameters **/
//	MBH = 4.e6;

	/** input parameters appropriate to Sgr A* **/
	//double BH_MASS = MBH * MSUN;

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


double interp_scalar(double ***var, int i, int j, int k, double coeff[4])
{

	double interp;
	int ip1,jp1,kp1;
	double b1,b2,b3,del[NDIM];

	del[1]=coeff[1];
	del[2]=coeff[2];
	del[3]=coeff[3];
	if(del[1] >1 || del[2]>1||del[3]>1 || del[1] < 0 || del[2] < 0 || del[3] < 0)
		fprintf(stderr,"del[1] %e \n del[2] %e\n del[3] %e\n",del[1],del[2],del[3]);

    ip1 = i+1;
    jp1 = j+1;
    kp1 = k+1;

    b1 = 1.-del[1];
    b2 = 1.-del[2];
    b3 = 1.-del[3];

	interp = var[i][j][k]*b1*b2 +
	  var[i][jp1][k]*b1*del[2] +
	  var[ip1][j][k]*del[1]*b2 +
	  var[ip1][jp1][k]*del[1]*del[2];

	  /* Now interpolate above in x3 */
	interp = b3*interp +
	  del[3]*(var[i][j][kp1]*b1*b2+
		  var[i][jp1][kp1]*b1*del[2]+
		  var[ip1][j][kp1]*del[1]*b2+
		  var[ip1][jp1][kp1]*del[1]*del[2]);

	return interp;
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
	double phi;
	/* Map X[3] into sim range, assume startx[3] = 0 */
	phi = fmod(X[3], stopx[3]);
	if(phi < 0.) phi = stopx[3]+phi;

	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
	*k = (int) (phi/dx[3] + 1000 -0.5) - 1000;

	if(*i < 0) {
		*i = 0 ;
		del[1] = 0. ;
	}
	else if(*i > N1-2) {
		*i = N1-2 ;
		del[1] = 1. ;
	}
	else {
		del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}

	if(*j < 0) {
		*j = 0 ;
		del[2] = 0. ;
	}
	else if(*j > N2-2) {
		*j = N2-2 ;
		del[2] = 1. ;
	}
	else {
		del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
	}

	if(*k < 0) {
		*k = 0 ;
		del[3] = 0. ;
	}
	else if(*k > N3-2) {
		*k = N3-2 ;
		del[3] = 1. ;
	}
	else {
		del[3] = (phi - ((*k + 0.5) * dx[3])) / dx[3];
	}
	if(del[3]<0)
		fprintf(stderr,"%e %e %e %d\n",del[3],phi,(*k+0.5)*dx[3],*k);

	return;
}


// ALLOCATION STUFF BELOW HERE
//////////////////////////////

static void *malloc_rank1(int n1, int alloc_size)
{
    void *A;

    if ((A = (void *) malloc(n1 * alloc_size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

double ***malloc_rank3(int n1, int n2, int n3)
{
	double ***A;
	double *space;
	int i,j;

	space = malloc_rank1(n1*n2*n3, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));

	for(i = 0; i < n1; i++){
		A[i] = malloc_rank1(n2,sizeof(double *));
		for(j = 0; j < n2; j++){
			A[i][j] = &(space[n3*(j + n2*i)]);
		}
	}

	return A;
}

void init_storage(void)
{
  int i,j,k;


  p = (double ****)malloc(NPRIM*sizeof(double***));
  for(i = 0; i < NPRIM; i++){
   p[i]=(double ***)malloc(N1*sizeof(double **));
   for(j= 0; j < N1; j++){
    p[i][j]=(double **)malloc(N2*sizeof(double *));
    for(k=0;k<N2;k++){
     p[i][j][k]=(double *)malloc(N3*sizeof(double));
    }
   }
  }


fprintf(stderr,"done here with memory, %d %d %d %d\n",N1,N2,N3,NPRIM);
  return;

}
