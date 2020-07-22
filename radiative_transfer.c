/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 */

#include <math.h>
#include "functions.h"
#include "constants.h"
#include "parameters.h"
#include <stdio.h>

//kappa distribution function

int N_theta, N_theta_e, N_B, N_n_e,  N_nuratio, N_nu;
double theta_min, theta_max, d_theta;
double theta_e_min, theta_e_max, d_theta_e;
double B_min, B_max, d_B;
double n_e_min, n_e_max, d_n_e;
double nuratio_min, nuratio_max, d_nuratio;
double nu_min, nu_max, d_nu;

double kappa,gamma_min,gamma_max,gamma_cutoff;
double **j_nu_data, **alpha_nu_data;

#define ME  (9.10956e-28)
#define mc2 (8.18726e-07)
#define kb  (1.38e-16)
#define hpl (6.6262e-27)
#define CL  (2.99792458e10)
#define keV (1.602e-9)
#define alphaf  (7.29735e-3)
#define h__mc2  (8.09e-21)
#define SIGMATH (0.665245873e-24)


/*Approximation from Maxon 1972*/
double Ei(double xx){
  double x,x2,Eii;
  x=fabs(xx);
  x2=x*x;

  if (x<=1.){
    Eii=log(x)+0.577-x+0.25*x2-0.055*x2*x;
  }else if(x>1.){
    Eii=-exp(-x)/x*(x2+2.33*x+0.25)/(x2+3.33*x+1.68);
  }
  return Eii;
}


double Bl(double x){

  return 0.85+1.35*sqrt(x)+0.38*x;

}

/*calculates e-e spectrum using tables and interpolations*/
double Gfun(double x,double Thetae){

    int i,ii;
    double en,step_en;
    double f1A,f1B,f1C,f1D;
    double f2A,f2B,f2C,f2D;
    double x1,x2;
    double AA,BB,CC,DD;
    double alphaa,betab,gammag,deltad;
    double Gfun1;
    double energy[]={50,75,100,150,200,300,400,500,600,700,800,900,1000}; //keV
  
  /* low energy data for different temperatures*/

  double A[] = { 1.584, 1.357, 1.197, 1.023, 0.883, 
		 0.700, 0.572, 0.484, 0.417, 0.361, 
		 0.322, 0.286, 0.259 };
  
  double B[] = { 0.578, 0.437, 0.291, 0.204, 0.0835, 
		 -0.0494, -0.139, -0.181, -0.209, -0.204, 
		 -0.244, -0.257, -0.258 };
  
  double C[] = {4.565, 3.842, 3.506, 3.036, 2.831, 
		2.545, 2.352, 2.175, 2.028, 1.914, 
		1.795, 1.705, 1.617 };
  
  double D[] = {2.091, 1.855, 1.672, 1.593, 1.487, 
		1.364, 1.254, 1.179, 1.108, 1.030, 
		0.982, 0.923, 0.879 };
  

    /* high energy data for different temperatures*/
  double alpha[] = {0.0387,0.0633, 0.0862, 0.128, 0.159, 
		    0.208, 0.234, 0.245, 0.248, 0.247, 
		    0.243, 0.239, 0.235};
  
  double beta[] = {0.523, 0.540, 0.569, 0.596, 0.658, 
		   0.633, 0.643, 0.695, 0.729, 0.756, 
		   0.763, 0.755, 0.735};
    
  double gamma[] = {5.319, 4.412, 3.897, 3.383, 2.974, 
		    2.738, 2.424, 2.025, 1.716, 1.457, 
		    1.271, 1.140, 1.060};

  double delta[] = {0.782, 0.689, 0.633, 0.523, 0.532, 
		    0.326, 0.302, 0.394, 0.453, 0.500, 
		    0.515, 0.508, 0.478};

  en=Thetae*ME*CL*CL/keV;

  /* find ABCD or alpha,beta,gamma,detla for a given temperature*/

  
  if(x > 0.05 && x < 1.1){
      
      if(en <  50){
        x1 = energy[0];
        x2 = energy[1];
	f1A = A[0];
	f2A = A[1];
	f1B = B[0];
	f2B = B[1];
	f1C = C[0];
	f2C = C[1];
	f1D = D[0];
	f2D = D[1];
      }else if(en > 1000){
	x1 = energy[11];
	x2 = energy[12];
	f1A = A[11];
	f2A = A[12];
	f1B = B[11];
	f2B = B[12];
	f1C = C[11];
	f2C = C[12];
	f1D = D[11];
	f2D = D[12];
      }else{

	  for(ii=0;ii<12;ii++){
	      if( en >= energy[ii] && en <= energy[ii+1]){
		  i=ii;
		  x1 = energy[i] ;
		  x2 = energy[i+1];
		  //fprintf(stdout,"1:i=%d ii=%d en=%g x1=%g x2=%g\n",i,ii,en,x1,x2);
		  break;
	      }
	  }
	  
	  
	  f1A = A[i];
	  f2A = A[i + 1];
	  f1B = B[i];
	  f2B = B[i + 1];
	  f1C = C[i];
	  f2C = C[i + 1];
	  f1D = D[i];
	  f2D = D[i + 1];
      }
      AA = (f2A - f1A) / (x2 - x1) * en + (f1A * x2 - x1 * f2A) / (x2 - x1);
      BB = (f2B - f1B) / (x2 - x1) * en + (f1B * x2 - x1 * f2B) / (x2 - x1);
      CC = (f2C - f1C) / (x2 - x1) * en + (f1C * x2 - x1 * f2C) / (x2 - x1);
      DD = (f2D - f1D) / (x2 - x1) * en + (f1D * x2 - x1 * f2D) / (x2 - x1);

      /*interpolation or extrapolation*/
      Gfun1=(AA+BB*x)*log(1./x)+CC+DD*x;
      
    }else if(x > 1.1 && x < 10.){
	    
      if(en < 50){
	x1 = energy[0];
	x2 = energy[1];
	f1A = alpha[0];
	f2A = alpha[1];
	f1B = beta[0];
	f2B = beta[1];
	f1C = gamma[0];
	f2C = gamma[1];
	f1D = delta[0];
	f2D = delta[1];
      }else if(en > 1000){
	x1 = energy[11];
	x2 = energy[12];
	f1A = alpha[11];
	f2A = alpha[12];
	f1B = beta[11];
	f2B = beta[12];
	f1C = gamma[11];
	f2C = gamma[12];
	f1D = delta[11];
	f2D = delta[12];
      }else{
	  for(ii=0;ii<12;ii++){
	      if(en >= energy[ii] && en <= energy[ii+1]){
		  i=ii;
		  x1 = energy[i] ;
		  x2 = energy[i+1];
		  //fprintf(stdout,"2:i=%d ii=%d en=%g x1=%g x2=%g\n",i,ii,en,x1,x2);
		  break;
	      }
	  }
	f1A = alpha[i];
	f2A = alpha[i + 1];
	f1B = beta[i];
	f2B = beta[i + 1];
	f1C = gamma[i];
	f2C = gamma[i + 1];
	f1D = delta[i];
	f2D = delta[i + 1];
      }
      alphaa =(f2A - f1A) / (x2 - x1) * en + (f1A * x2 - x1 * f2A) / (x2 - x1);
      betab  =(f2B - f1B) / (x2 - x1) * en + (f1B * x2 - x1 * f2B) / (x2 - x1);
      gammag =(f2C - f1C) / (x2 - x1) * en + (f1C * x2 - x1 * f2C) / (x2 - x1);
      deltad =(f2D - f1D) / (x2 - x1) * en + (f1D * x2 - x1 * f2D) / (x2 - x1);

      Gfun1=alphaa*x*x+betab*x+gammag+deltad/x;
    }else{

      Gfun1=0.0;       /*for other x'es it is zero???*/

    }
    return Gfun1;
}

double Xmax =1e-25, Xmin=1e25;

double emission_coeff_kappa(double nu,double Ne, double Thetae, double B, double theta){

    //emissivity for the kappa distribution function, see Pandya et al. 2016
    double nus,nuc, sth, x,  w, X_c, factor;
    double  J_s;

    double Rhigh = 100.;
    double Rlow = 100.;
    double b2 =  pow((B/B_unit)/(Ne/Ne_unit),2.);
   // printf("%g\n", b2);
   // kappa = Rhigh * b2/(1+b2) + Rlow / (1 + b2);

    w = Thetae;//sqrt(  2./9./kappa *Thetae * Thetae);
    nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    sth = sin(theta);

    factor = nuc*sth*Ne;// (Ne * pow(ELECTRON_CHARGE, 2.) * nuc * sth)/SPEED_OF_LIGHT;

    //      fprintf(stderr,"sinth %g\n", sth);
    nus = nuc * sth * pow( w * kappa,2);
   // if (nu > 1.e12 * nus || Thetae > 400. || Thetae < 1.)
    //    return (0.);
    X_c = nu / (nuc);
   // Thetae = 2.757851e+00;
    //      fprintf(stderr, "X_kappa %g\n", X_kappa);
    J_s=interpolate_scalar_2d(j_nu_data,X_c,Thetae);
   // 	printf("X_c %e Thetae %e J_s %e\n",X_c, Thetae, J_s); 
   if(J_s>0)   
   	fprintf(stderr,"T %g X %g J_s %g %g\n",Thetae, X_c,  J_s*factor,emission_coeff_THSYNCH(B, theta, Thetae, 1*nuc, Ne) );
    return ( J_s*factor);

}

double absorption_coeff_kappa(double nu,double Ne, double Thetae, double B, double theta){

   //emissivity for the kappa distribution function, see Pandya et al. 2016
   double nus,nuc, sth, x,  w, X_c, factor;
   double alpha_s;

   double Rhigh = 100.;
   double Rlow = 100.;
   double b2 =  pow((B/B_unit)/(Ne/Ne_unit),2.);
   // printf("%g\n", b2);
   kappa = Rhigh * b2/(1+b2) + Rlow / (1 + b2);

   w = Thetae;//sqrt(  2./9./kappa *Thetae * Thetae);
   nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
   sth = sin(theta);

   factor = Ne/(B*sth);// (Ne * pow(ELECTRON_CHARGE, 2.) * nuc * sth)/SPEED_OF_LIGHT;

   //      fprintf(stderr,"sinth %g\n", sth);
   nus = nuc * sth * pow( w * kappa,2);
   // if (nu > 1.e12 * nus || Thetae > 400. || Thetae < 1.)
   //    return (0.);
   X_c = nu / (nuc);
   alpha_s=interpolate_scalar_2d(alpha_nu_data,X_c,Thetae);
   return ( alpha_s*factor);

}

//non thermal emission
double emission_coeff_kappa_FIT(double nu,double Ne, double Thetae, double B, double theta){
    //emissivity for the kappa distribution function, see Pandya et al. 2016
    double nuc, sth, nus, x,  w, X_kappa, factor;
    double J_low, J_high, J_s;

    double Rhigh = 4.5;
    double Rlow = 4.5;
    double b2 =  pow((B/B_unit)/(Ne/Ne_unit),2.);
   // printf("%g\n", b2);
    kappa = Rhigh * b2/(1+b2) + Rlow / (1 + b2);

    w = Thetae;//sqrt(  2./9./kappa *Thetae * Thetae);
    nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    sth = sin(theta);

    factor = (Ne * pow(ELECTRON_CHARGE, 2.) * nuc * sth)/SPEED_OF_LIGHT;
    
    //      fprintf(stderr,"sinth %g\n", sth);
    nus = nuc * sth * pow( w * kappa,2);
    if (nu > 1.e12 * nus || Thetae > 400. || Thetae < .1)
        return (0.);
    X_kappa = nu / nus;
    //      fprintf(stderr, "X_kappa %g\n", X_kappa);
    J_low = pow(X_kappa,1./3.) * sth * 4 * M_PI * tgamma(kappa - 4./3.) /(pow(3,7./3.) * tgamma(kappa - 2.));
    //      fprintf(stderr, "J_low %g\n", J_low);
    J_high = pow(X_kappa,-(kappa-2)/2.) * sth * pow(3,(kappa-1)/2.) * (kappa - 2.)*(kappa - 1.)/4. * tgamma(kappa/4. - 1./3.) * tgamma(kappa/4. + 4./3.);
    //      fprintf(stderr, "J_high %g\n", J_high );
    x = 3 * pow(kappa,-3./2.);
    
    J_s = pow( ( pow(J_low,-x) + pow(J_high,-x) ) , -1./x );
    //      fprintf(stderr,"J_s %g\n", J_s * factor);
    return ( J_s*factor);
    
}


void init_memory(void){
        printf("Init memory\n");
        printf("%i %i\n", N_theta_e, N_nuratio);
        j_nu_data=(double **) malloc((N_theta_e+1) * sizeof(double *));
        alpha_nu_data= (double **) malloc((N_theta_e+1) * sizeof(double *));
        for(int i = 0; i<= N_theta_e; i++){
                j_nu_data[i]=(double *) malloc((N_nuratio+1) * sizeof(double ));
                alpha_nu_data[i]=(double *) malloc((N_nuratio+1) * sizeof(double ));
        }
        printf("done!\n");
}

void read_in_table(char* filename){

        FILE *fp;
        fp = fopen(filename, "r");
        printf("reading file %s\n", filename);
//read in header; min and max values of the 5 variables, sizes of these variables

        fscanf(fp, "%lf ", &nuratio_min);
        fscanf(fp, "%lf ", &nuratio_max);
        fscanf(fp, "%d ", &N_nuratio);
        fscanf(fp, "%lf ", &theta_e_min);
        fscanf(fp, "%lf ", &theta_e_max);
        fscanf(fp, "%d ", &N_theta_e);
        fscanf(fp, "%lf ", &theta_min);
        fscanf(fp, "%lf ", &theta_max);
        fscanf(fp, "%d ", &N_theta);
        fscanf(fp, "%lf ", &kappa);
        fscanf(fp, "%lf ", &gamma_min);
        fscanf(fp, "%lf ", &gamma_max);
        fscanf(fp, "%lf ", &gamma_cutoff);
        printf("Reading header completed\n");
        printf("N_nu= %d N_theta_e=%d\n", N_nuratio,  N_theta_e);
        printf("kappa=%lf\n", kappa);
        printf("gamma_min = %lf gamma_max=%lf gamma_cutoff=%lf\n", gamma_min, gamma_max, gamma_cutoff);

        d_nuratio = (nuratio_max - nuratio_min)/(N_nuratio-1);
        d_theta_e = (theta_e_max - theta_e_min)/(N_theta_e-1);

        printf("%lf %lf\n", d_nuratio, d_theta_e);

//init memory
        init_memory();
        double dump1,dump2,dump3;
//read in body; j_nu and alpha_nu
        for(int i = 0; i<= N_theta_e; i++){
                for(int j = 0; j<= N_nuratio; j++){
                        fscanf(fp, "%lf %lf %lf", &dump1,&dump2,&dump3);
                        fscanf(fp, "%lf ", &j_nu_data[i][j]);
                        fscanf(fp, "%lf ", &alpha_nu_data[i][j]);
                }
        }

}

double interpolate_scalar_2d(double **var, double nuratio, double Thetae)
{

        double interp, coeff[4],del[2];
        printf("%e %e\n",nuratio,Thetae);
        int j0 = (int) ((log10(nuratio)-log10(nuratio_min))*10 + 1000) - 1000;
        int i0 = (int) ((log10(Thetae)-log10(theta_e_min))*10 + 1000) - 1000;
        printf("i0 %i j0 %i\n",j0,i0);
        if(j0 < 0) {
                j0 = 0 ;
                del[1] = 0. ;
	//	return 0;
        }
        else if(j0 > N_nuratio-2) {
                j0 = N_nuratio-2 ;
                del[1] = 1. ;
        //	return 0;
	}
        else {
		del[1] =(nuratio - nuratio_min*pow(10., (double)j0/10.))/(nuratio_min*pow(10., (double)(j0+1)/10.) - nuratio_min*pow(10., (double)j0/10.));
        }

       // printf("%e\n",del[1]);

        if(i0 < 0) {
                i0 = 0 ;
                del[0] = 0. ;
       	//	return 0;
	 }
        else if(i0 > N_theta_e-2) {
                i0 = N_theta_e-2 ;
                del[0] = 1. ;
        //	return 0;
	}
        else {
                del[0] =(Thetae - theta_e_min*pow(10., (double)i0/10.))/(theta_e_min*pow(10., (double)(i0+1)/10.)-theta_e_min*pow(10., (double)i0/10.));
        }
        //printf("%e\n",del[0]);


        coeff[0] = (1. - del[0]) * (1. - del[1]);
        coeff[1] = (1. - del[0]) * del[1];
        coeff[2] = del[0] * (1. - del[1]);
        coeff[3] = del[0] * del[1];

       // printf("%e %e %e %e\n",coeff[0],coeff[1],coeff[2],coeff[3]);

        interp =
            var[i0][j0] * coeff[0] +
            var[i0][j0 + 1] * coeff[1] +
            var[i0 + 1][j0] * coeff[2] + var[i0 + 1][j0 + 1] * coeff[3];
        printf("value %e\n",var[i0][j0]);
	return interp;
}

double interpolate_scalar_5d(double *****A, double nu, double Ne, double Thetae, double B, double theta){
    
    //interpolate
    
    double value;
    int i0, j0, k0, l0, m0;
    int i1, j1, k1, l1, m1;
    double del_i, del_j, del_k, del_l, del_m;
    
    i0 = (int) ((nu-nu_min)/(d_nu) + 1000) - 1000;
    j0 = (int) ((Ne - n_e_min)/(d_n_e) + 1000) - 1000;
    k0 = (int) ((Thetae - theta_e_min)/d_theta_e + 1000) - 1000;;
    l0 = (int) ((B - B_min) /(d_B)  + 1000) - 1000;
    m0 = (int) ((theta-theta_min)/d_theta + 1000) - 1000;
    
    if(i0 < 0) {
        i0 = 0 ;
        del_i = 0. ;
    }
    else if(i0 > N_nu-1) {
        i0 = N_nu-1 ;
        del_i = 1. ;
    }
    else {
        del_i =( nu - ((i0 ) * d_nu + nu_min)) / d_nu;
    }
    
    if(j0 < 0) {
        j0 = 0 ;
        del_j = 0. ;
    }
    else if(j0 > N_n_e-1) {
        j0 = N_n_e-1 ;
        del_j = 1. ;
    }
    else {
        del_j = (Ne - ((j0 ) * d_n_e + n_e_min)) / d_n_e;
    }
    
    if(k0 < 0) {
        k0 = 0 ;
        del_k = 0. ;
    }
    else if(k0 > N_theta_e-1) {
        k0 = N_theta_e-1 ;
        del_k = 1. ;
    }
    else {
        del_k = (Thetae - ((k0 ) * d_theta_e + theta_e_min)) / d_theta_e;
    }
    
    if(l0 < 0) {
        l0 = 0 ;
        del_l = 0. ;
    }
    else if(l0 > N_B-1) {
        l0 = N_B-1 ;
        del_l = 1. ;
    }
    else {
        del_l = (B - ((l0 ) * d_B + B_min)) / d_B;
    }
    
    if(m0 < 0) {
        m0 = 0 ;
        del_m = 0. ;
    }
    else if(m0 > N_theta-1) {
        m0 = N_theta-1 ;
        del_m = 1. ;
    }
    else {
        del_m = (theta - ((m0 ) * d_theta + theta_min)) / d_theta;
    }
    
    printf("%d %d %d %d %d\n", i0,j0,k0,l0,m0);
    printf("%lf %lf %lf %lf %lf\n", del_i,del_j,del_k,del_l,del_m);
    i1=i0+1;
    j1=j0+1;
    k1=k0+1;
    l1=l0+1;
    m1=m0+1;
    
    value = A[i0][j0][k0][l0][m0]*(1-del_i)*(1-del_j)*(1-del_k)*(1-del_l)*(1-del_m) +
    A[i1][j0][k0][l0][m0]*(del_i)*(1-del_j)*(1-del_k)*(1-del_l)*(1-del_m) +
    A[i0][j1][k0][l0][m0]*(1-del_i)*(del_j)*(1-del_k)*(1-del_l)*(1-del_m) +
    A[i0][j0][k1][l0][m0]*(1-del_i)*(1-del_j)*(del_k)*(1-del_l)*(1-del_m) +
    A[i0][j0][k0][l1][m0]*(1-del_i)*(1-del_j)*(1-del_k)*(del_l)*(1-del_m) +
    A[i0][j0][k0][l0][m1]*(1-del_i)*(1-del_j)*(1-del_k)*(1-del_l)*(del_m) +
    A[i1][j1][k0][l0][m0]*(del_i)*(del_j)*(1-del_k)*(1-del_l)*(1-del_m) +
    A[i1][j0][k1][l0][m0]*(del_i)*(1-del_j)*(del_k)*(1-del_l)*(1-del_m) +
    A[i1][j0][k0][l1][m0]*(del_i)*(1-del_j)*(1-del_k)*(del_l)*(1-del_m) +
    A[i1][j0][k0][l0][m1]*(del_i)*(1-del_j)*(1-del_k)*(1-del_l)*(del_m) +
    A[i0][j1][k1][l0][m0]*(1-del_i)*(del_j)*(del_k)*(1-del_l)*(1-del_m) +
    A[i0][j1][k0][l1][m0]*(1-del_i)*(del_j)*(1-del_k)*(del_l)*(1-del_m) +
    A[i0][j1][k0][l0][m1]*(1-del_i)*(del_j)*(1-del_k)*(1-del_l)*(del_m) +
    A[i0][j0][k1][l1][m0]*(1-del_i)*(1-del_j)*(del_k)*(del_l)*(1-del_m) +
    A[i0][j0][k1][l0][m1]*(1-del_i)*(1-del_j)*(del_k)*(1-del_l)*(del_m) +
    A[i0][j0][k0][l1][m1]*(1-del_i)*(1-del_j)*(1-del_k)*(del_l)*(del_m) +
    A[i1][j1][k1][l0][m0]*(del_i)*(del_j)*(del_k)*(1-del_l)*(1-del_m) +
    A[i1][j1][k0][l1][m0]*(del_i)*(del_j)*(1-del_k)*(del_l)*(1-del_m) +
    A[i1][j1][k0][l0][m1]*(del_i)*(del_j)*(1-del_k)*(1-del_l)*(del_m) +
    A[i1][j0][k1][l1][m0]*(del_i)*(1-del_j)*(del_k)*(del_l)*(1-del_m) +
    A[i1][j0][k1][l0][m1]*(del_i)*(1-del_j)*(del_k)*(1-del_l)*(del_m) +
    A[i1][j0][k0][l1][m1]*(del_i)*(1-del_j)*(1-del_k)*(del_l)*(del_m) +
    A[i0][j1][k1][l1][m0]*(1-del_i)*(del_j)*(del_k)*(del_l)*(1-del_m) +
    A[i0][j1][k1][l0][m1]*(1-del_i)*(del_j)*(del_k)*(1-del_l)*(del_m) +
    A[i0][j1][k0][l1][m1]*(1-del_i)*(del_j)*(1-del_k)*(del_l)*(del_m) +
    A[i0][j0][k1][l1][m1]*(1-del_i)*(1-del_j)*(del_k)*(del_l)*(del_m) +
    A[i1][j1][k1][l1][m0]*(del_i)*(del_j)*(del_k)*(del_l)*(1-del_m) +
    A[i1][j1][k1][l0][m1]*(del_i)*(del_j)*(del_k)*(1-del_l)*(del_m) +
    A[i1][j1][k0][l1][m1]*(del_i)*(del_j)*(1-del_k)*(del_l)*(del_m) +
    A[i1][j0][k1][l1][m1]*(del_i)*(1-del_j)*(del_k)*(del_l)*(del_m) +
    A[i0][j1][k1][l1][m1]*(1-del_i)*(del_j)*(del_k)*(del_l)*(del_m) +
    A[i1][j1][k1][l1][m1]*(del_i)*(del_j)*(del_k)*(del_l)*(del_m);
    
    
    return value;
}



// Return emissivity j_nu which depends on local plasma parameters
// Ref. Dolence & Moscibrodzka 2009
double emission_coeff_THSYNCH(double B, double theta, double THETA_e, double nu_plasma, double n_e){
    double nu_c = ELECTRON_CHARGE * B /
                  (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    
    double nu_s = 2. / 9. * nu_c * THETA_e * THETA_e * sin(theta);

    double X    =nu_plasma /( nu_s);
    double f    = pow(pow(X , 0.5) + pow(2., 11. / 12.) * pow(X, 1. / 6.), 2.);
    double j_nu = n_e * sqrt(2.) * M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE *
                  nu_s / (6. * THETA_e * THETA_e * SPEED_OF_LIGHT) * f *
                  exp(-pow(X,1./3.));

    return j_nu;
}

// Return emission constant j_nu as described in Dexter (2009) (the geokerr paper)
double emission_coeff_THSYNCHAV(double B, double THETA_e, double nu_p, double n_e){
    double nu_c = 3. * ELECTRON_CHARGE * B * THETA_e * THETA_e /
                  (4. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double x_M = nu_p / nu_c;
    double I = 4.0505 / pow(x_M, 1. / 6.) * (1. + 0.4 / pow(x_M, 1. / 4.) +
               0.5316 / sqrt(x_M)) * exp(-1.8899 * pow(x_M, 1. / 3.));

    double j_nu = nu_p * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                  (2. * sqrt(3.) * SPEED_OF_LIGHT) * 1. / (THETA_e * THETA_e) * I;

//    j_nu = THETA_e * (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / BOLTZMANN_CONSTANT;
//    j_nu = n_e;

    return  j_nu;
}

// Return emission coefficient for thermal free-free radiation
double emission_coeff_FFTHERMAL(double nu, double n_e, double T){
    double n_i = n_e; // Assume neutral hydrogen plasma
    double Z = 1.; // For H, Z = 1
    double g_ff = 1.; // Gaunt factor

    double j_nu = 5.44e-39 * (Z * Z / sqrt(T)) * n_i * n_e * g_ff *
                  exp(-PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T));

    return j_nu;
}

// Return emissivity for the simple Gaussian hot spot model discussed in Dexter 2009.
double emissivity_hotspot(double *X_u){
    double xspot[4];

    double r  = logscale ? exp(X_u[1]) : X_u[1];

    double Rspot = 0.5; // Spot size (radius)

    double r_spot  = 6.0; // Spot orbit radius
    double th_spot = 0.5 * M_PI;
    double r32     = pow(r_spot, 1.5);
    double omega   = 1. / (r32+a);
    double P       = 2. * M_PI/omega; // Period of the spot on a Keplerian orbit[M]

    //spot currrent position
    xspot[0] = X_u[0]; //current coordinate time
    xspot[1] = logscale ? log(r_spot) : r_spot;
    xspot[2] = th_spot;                          //equator 0.5*pi
    xspot[3] = fmod(X_u[0] / P, 1.) * 2. * M_PI + M_PI; //spot current phi at t=X[0]

    // Pseudo-Cartesian coordinates
    double xc = sqrt(r * r + a * a) * cos(X_u[3]);
    double yc = sqrt(r * r + a * a) * sin(X_u[3]);
    double zc = exp(X_u[1]) * cos(X_u[2]);

    double xs = sqrt(r_spot * r_spot + a * a) * cos(xspot[3]);
    double ys = sqrt(r_spot * r_spot + a * a) * sin(xspot[3]);
    double zs = r_spot * cos(xspot[2]);

    //distance^2 between photon position and spot center
    double xx = fabs(pow(xc - xs, 2) + pow(yc - ys, 2) + pow(zc - zs, 2));

    if (xx <= 4.)
        return exp(-(xx) / 2. / Rspot / Rspot);
    return 0.;
}

// Return emissivity for the thin disk model described in Dexter 2009
double emissivity_thindisk(double *X_u){
    double r = logscale ? exp(X_u[1]) : X_u[1];
    return 1 / r / r;
}

// Planck function
double planck_function(double nu, double THETA_e){
    double T = THETA_e * ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT / BOLTZMANN_CONSTANT;
    return 2. * PLANCK_CONSTANT * nu * nu * nu / (SPEED_OF_LIGHT * SPEED_OF_LIGHT) * 1. / (exp(PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T)) - 1.);
}

double planck_function2(double nu, double Thetae)
{
    double X = PLANCK_CONSTANT * nu / (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT * Thetae);

    if (X < 2.e-3)
        return ((2. * PLANCK_CONSTANT / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)) /
        (X / 24. * (24. + X * (12. + X * (4. + X)))));

    return ((2. * PLANCK_CONSTANT / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)) / (exp(X) - 1.));
}

// Return absorption coefficient - assume local thermodynamical equilibrium so that Kirchoff's Law applies
double absorption_coeff_TH(double j_nu, double nu, double THETA_e){
    double B_nu = planck_function(nu, THETA_e); // Planck function
    return j_nu / B_nu;
}

// Compute the photon frequency in the plasma frame:
double freq_in_plasma_frame(double Uplasma_u[4], double k_d[4]){
    double nu_plasmaframe = 0.;
    int i;
    LOOP_i nu_plasmaframe += Uplasma_u[i] * k_d[i];
    nu_plasmaframe *= -(ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
                       PLANCK_CONSTANT;

    if(isnan(nu_plasmaframe))
	fprintf(stderr,"NAN in plasma frame %e %e %e %e %e\n",nu_plasmaframe,Uplasma_u[0],Uplasma_u[1],Uplasma_u[2],Uplasma_u[3]);  
    return nu_plasmaframe;
}

// See eqn 73 in Dexter 2016
double pitch_angle(double *X_u, double *k_u, double *B_u, double *Uplasma_u){
    double b_dot_k = inner_product(X_u, B_u, k_u);
    double b_dot_b = inner_product(X_u, B_u, B_u);
    double k_dot_u = inner_product(X_u, k_u, Uplasma_u);

    // Compute and clamp result (result can slightly exceed domain of acos due to numerics)
    double result = acos(b_dot_k / (-k_dot_u * sqrt(fabs(b_dot_b) + 1.e-15)));
    result = fmax(fmin(result, M_PI), 0.);

//    return result;

//NEW VERSION
    double B, k, mu;


    B = sqrt(fabs(inner_product(X_u, B_u, B_u)));

    if (B == 0.)
        return (M_PI / 2.);

    k = fabs(inner_product(X_u, k_u, Uplasma_u));

    mu = inner_product(X_u, k_u, B_u) / (k * B);

    if (fabs(mu) > 1.)
    mu /= fabs(mu);

    if (isnan(mu))
    fprintf(stderr, "isnan get_bk_angle\n");

    return (acos(mu));
}

/*
double radiative_transfer(double *lightpath, int steps, double frequency){
    int IN_VOLUME, path_counter;
    double I_current = 0.;
    double dI        = 0.;
    double j_nu      = 0.;
    double B, THETA_e, pitch_ang, nu_p, n_e, nu_p2, dl_current;
    int i;
    double X_u[4], k_u[4], k_d[4], B_u[4], Uplasma_u[4];
    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm

    double tau = 0.;

    double a_nu = 0.;

    double K_inv_old = 0, j_inv_old=0, dtau_old=0;

    // Move backward along constructed lightpath
    for (path_counter = steps - 1; path_counter > 0; path_counter--){
        // Current position, wave vector, and dlambda
        LOOP_i{
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter-1) * 9 + 8]);

        // Obtain the parameters n_e, THETA_e, B, and Uplasma_u at X_u
        //get_plasma_parameters(X_u, &n_e, &THETA_e, &B, Uplasma_u);
        get_fluid_params(X_u, &n_e, &THETA_e, &B, B_u, Uplasma_u, &IN_VOLUME);

        // Check whether the ray is currently in the GRMHD simulation volume
        if(IN_VOLUME){
            // Obtain pitch angle: still no units (geometric)
            pitch_ang = pitch_angle(X_u, k_u, B_u, Uplasma_u);

            // CGS UNITS USED FROM HERE ON OUT
            //////////////////////////////////

            // Scale the wave vector to correct energy
            LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                             (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // Convert distance dlambda accordingly
            dl_current *= (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / (PLANCK_CONSTANT * frequency);

            //lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(Uplasma_u, k_d);
            nu_p2 = nu_p * nu_p;

            // Obtain emission coefficient in current plasma conditions
           j_nu = emission_coeff_THSYNCHAV(B, THETA_e, nu_p, n_e);

	    // Obtain absorption coefficient
            if (ABSORPTION){
    	        a_nu = absorption_coeff_TH(j_nu, nu_p, THETA_e);
	    }

            // Constant used in integration (to produce correct units)
            double C = Rg * PLANCK_CONSTANT / (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            double redshift = frequency / nu_p;

            double dtau  = (nu_p * a_nu * dl_current * C + dtau_old);
            double K_inv = (nu_p * a_nu);
            double j_inv = (j_nu / nu_p2);

            // Only add I_current if it is not NaN
            if(j_nu == j_nu && exp(X_u[1]) < RT_OUTER_CUTOFF){ // I_current += exp(-tau) * j_nu / nu_p / nu_p * dl_current * C;
       //         I_current += dI; // Old way of integrating
       		double Ii=I_current;
		double S = j_inv/K_inv; 
		if(K_inv == 0 )
			I_current = Ii;
		else if(dtau < 1.e-5)
                	I_current = Ii - (Ii - S) * ( 0.166666667*dtau * (6. - dtau * (3. - dtau)));
		else{
                	double efac = exp(-dtau);
                	I_current = Ii*efac + S*(1. - efac);
        	}
	     }

        }
    }

    // Store integrated intensity in the image
    return I_current * pow(frequency, 3.);
}
*/
