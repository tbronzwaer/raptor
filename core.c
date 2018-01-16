#include "parameters.h"
#include "functions.h"

void read_model( char *argv[]){
        char inputfile[100];
        char inputfile2[100];

        // model to read
        sscanf(argv[1], "%s", inputfile);
        printf("Reading %s...",inputfile);
        FILE *input;
        input = fopen(inputfile, "r");
        if (input == NULL) {
                printf ("Cannot read input file %s \n",inputfile);
                exit(1);
        }

        char temp[100], temp2[100];


        // Model parameters
        fscanf(input, "%s %s %lf", temp, temp2, &MBH);
        fscanf(input, "%s %s %lf", temp, temp2, &M_UNIT);
        fscanf(input, "%s %s %d",  temp, temp2, &ABSORPTION);

        // Observer parameters
        fscanf(input, "%s %s %d",  temp, temp2, &IMG_WIDTH);
        fscanf(input, "%s %s %d",  temp, temp2, &IMG_HEIGHT);
        fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_X);
        fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_Y);
        fscanf(input, "%s %s %d",  temp, temp2, &FREQS_PER_DEC);
        fscanf(input, "%s %s %lf", temp, temp2, &FREQ_MIN);
        fscanf(input, "%s %s %lf", temp, temp2, &FREQ_MAX);
        fscanf(input, "%s %s %lf", temp, temp2, &STEPSIZE);

        // INPUT AND OUTPUT FILES
        sscanf(argv[2], "%s", GRMHD_FILE);
        sscanf(argv[3], "%lf", &M_UNIT);
        sscanf(argv[4], "%lf", &INCLINATION);
        sscanf(argv[5], "%lf", &R_HIGH);
        sscanf(argv[6], "%lf", &R_LOW);
        sscanf(argv[7], "%lf", &TIME_INIT);

        fprintf(stderr,"Model parameters:\n");
        fprintf(stderr,"MBH \t\t= %g \n", MBH);
        fprintf(stderr,"M_UNIT \t\t= %g \n", M_UNIT);
        fprintf(stderr,"R_HIGH \t\t= %g \n", R_HIGH);
        fprintf(stderr,"R_LOW \t\t= %g \n", R_LOW);
        fprintf(stderr,"INCLINATION \t= %g \n", INCLINATION);

        fprintf(stderr,"ABSORPTION \t= %d \n", ABSORPTION);

        fprintf(stderr,"Outer R \t= %lf \n\n", RT_OUTER_CUTOFF);

        fprintf(stderr,"Observer parameters:\n");
        fprintf(stderr,"IMG_WIDTH \t= %d \n", IMG_WIDTH);
        fprintf(stderr,"IMG_HEIGHT \t= %d \n", IMG_HEIGHT);
        fprintf(stderr,"CAM_SIZE_X \t= %g \n", CAM_SIZE_X);
        fprintf(stderr,"CAM_SIZE_Y \t= %g \n", CAM_SIZE_Y);
        fprintf(stderr,"FREQS_PER_DEC \t= %d \n", FREQS_PER_DEC);
        fprintf(stderr,"FREQ_MIN \t= %.03g \n", FREQ_MIN);
        fprintf(stderr,"FREQ_MAX \t= %.03g \n", FREQ_MIN * pow(10,(num_indices-1.)/FREQS_PER_DEC));
        fprintf(stderr,"NUM OF FREQS \t= %d \n",num_indices);
        fprintf(stderr,"STEPSIZE \t= %g \n", STEPSIZE);
        fclose (input);

        printf("Done!\n");

        initrcarry(982451653);
        srand(982451653);
}

void calculate_image( real ** intensityfield, real energy_spectrum[num_indices],real frequencies[num_indices]){

        static real intensityfield2[maxsize][num_indices];

        int lmax = (int)((real)IMG_HEIGHT*IMG_WIDTH/(real)maxsize + 0.5);
        if(lmax==0)
                lmax=1.;

        for(int i = 0; i < maxsize; i++) {
                for(int f = 0; f < num_indices; f++) {
                        intensityfield2[i][f]=0;
                }
        }

        real start=clock();
        real diff = clock() - start;
        clock_t startgpu=clock();

        int msec;
        int l1,l2;
        for(int l=0; l<lmax; l++) {
                l1 =(int)l*maxsize;
                l2 =(int)(l+1)*maxsize;

                if(l2 >(IMG_WIDTH)*(IMG_HEIGHT))
                        l2 =(IMG_WIDTH)*(IMG_HEIGHT);

                //#pragma acc copyin(intensityfield2[0:1000][0:(num_indices)])
#pragma omp parallel for shared(energy_spectrum,frequencies,intensityfield,p) schedule(static,1)
                //#pragma acc kernels vector_length(1) copyin(Xcam[0:4],Ucam[0:4],IMG_WIDTH,IMG_HEIGHT,intensityfield2[0:maxsize][0:num_indices],p[0:NPRIM][0:N1][0:N2][0:N3],frequencies[0:num_indices],l1,l2) copyout(intensityfield2[0:maxsize][0:(num_indices)])
                for(int i=l1; i < l2; i++) { // For all pixel rows (distributed over threads)...
                        int y=(int)(i/IMG_WIDTH);
                        int x=(int)(i%IMG_WIDTH);

                        // INTEGRATE THIS PIXEL'S GEODESIC AND PERFORM RADIATIVE TRANSFER AT DESIRED FREQUENCIES, STORE RESULTS
                        int icur=(int)(i-l1);
                        integrate_geodesic(icur,x,y,intensityfield2,frequencies,p,TIME_INIT,Xcam,Ucam);
                }

                diff=clock()-startgpu;
                int msec = diff *1000/ (CLOCKS_PER_SEC*20);
                printf("Done: %.02g %%, speed: %.02g [geodesics/sec]\n", 100.*(real)l2/((real)(IMG_WIDTH)*(IMG_HEIGHT)),(real)l2 /((real)msec/1000.));

#pragma omp parallel for shared(energy_spectrum,frequencies,intensityfield,p) schedule(static,1)
                for(int k  = l1; k < l2; k++) { // For all pixel rows (distributed over threads)...
                        for(int fr=0; fr<num_indices; fr++) {
#if (LOG_IMPACT_CAM)
                                int y=(int)(k/IMG_WIDTH);
                                int x=(int)(k%IMG_WIDTH);
                                real r_i = exp(log(20.)*(real)(x+0.5) /(real) IMG_WIDTH) - 1.;
                                real theta_i = 2.*M_PI  * (real)(y+0.5)/ (real)IMG_HEIGHT;
                                real alpha = r_i * cos(theta_i);
                                real beta  = r_i * sin(theta_i);
                                real d_r = R_GRAV * (r_i +1.) * log(20.)/(real)IMG_WIDTH;
                                real d_theta = R_GRAV*2.* M_PI / (real)IMG_HEIGHT;

                                intensityfield[k][fr]=intensityfield2[k-l1][fr] * pow(frequencies[fr], 3.)* r_i* d_r * d_theta/(source_dist*source_dist); // * e2_c;
                                intensityfield2[k-l1][fr]=0;
#elif (LINEAR_IMPACT_CAM)
                                intensityfield[k][fr]=intensityfield2[k-l1][fr] * pow(frequencies[fr], 3.); //* e2_c;
                                intensityfield2[k-l1][fr]=0;

#endif
                        }
                }
        }
        free(p);
#pragma acc wait

        for(int i=0; i<IMG_WIDTH*IMG_HEIGHT; i++) {
                for(int f=0; f<num_indices; f++) {
                        energy_spectrum[f]+=intensityfield[i][f];
                }
        }
}

void output_files(real ** intensityfield,real energy_spectrum[num_indices],real frequencies[num_indices]){
        struct stat st = {0};
        char spec_folder[256]="";
        sprintf(spec_folder, "output");

        if (stat(spec_folder, &st) == -1) {
                mkdir(spec_folder, 0700);
        }

#if (SPECFILE)
        char spec_filename[256] = "";
        sprintf(spec_filename, "%s/spectrum_%d_%.02lf.dat",spec_folder,(int)TIME_INIT,INCLINATION);
        FILE *specfile    = fopen(spec_filename, "w");
#endif

        for(int f = 0; f < num_indices; f++) { // For all frequencies...
                // Create filenames, open files

#if (IMGFILE)
                char dat_filename[256] = "";
                sprintf(dat_filename, "%s/img_data_%d_%e_%.02lf.dat",spec_folder,(int)TIME_INIT,frequencies[f],INCLINATION);
                FILE *imgfile     = fopen(dat_filename, "w");
                printf("%s\n",GRMHD_FILE);
                fprintf(imgfile,"%d %d %e %e %s %e %e %e %e\n",IMG_WIDTH, IMG_HEIGHT, JANSKY_FACTOR *energy_spectrum[f],frequencies[f],GRMHD_FILE,M_UNIT,INCLINATION,R_LOW,R_HIGH);
                write_image(imgfile, intensityfield,f, JANSKY_FACTOR);
                fclose(imgfile);

#endif
#if (VTKFILE)
                char vtk_filename[256] = "";

                sprintf(vtk_filename, "%s/img_data_%d_%e_%.02lf.vtk",spec_folder,(int)TIME_INIT,frequencies[f],INCLINATION);
                FILE *vtkfile          = fopen(vtk_filename, "w");
                write_VTK_image(vtkfile, intensityfield,f,JANSKY_FACTOR);
                fclose(vtkfile);
#endif


                fprintf(stderr,"Frequency %.5e Integrated flux density = %.5e\n", frequencies[f],JANSKY_FACTOR * energy_spectrum[f]);

#if (SPECFILE)
                fprintf(specfile, "%+.15e\t%+.15e\n", frequencies[f], JANSKY_FACTOR * energy_spectrum[f]);
#endif
        }
        fclose(specfile);

}
