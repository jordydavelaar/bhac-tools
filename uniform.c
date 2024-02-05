#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include <fftw3.h>
#include <stdlib.h>
#include <math.h>

#define UNI1 2000
#define UNI2 2000

void perform_fft(fftw_complex *quant, fftw_complex* out);
void write_plasma_file(fftw_complex *quant,double TIME_INIT);
void fftfreq(double* freq, int N, double size);

void uniform_data_array(double TIME_INIT){

	fftw_complex *quant;
	fftw_complex *out;

	struct GRMHD modvar;

    quant = (fftw_complex *) fftw_malloc(UNI1*UNI2*sizeof(fftw_complex));
    out = (fftw_complex *) fftw_malloc(UNI1*UNI2*sizeof(fftw_complex));

    modvar.igrid_c = -1;

	double X[4] = {0,0,0,0};

	double XMAX[2];
	double XMIN[2];
	double dx[2];

	XMAX[0] = 1.;
	XMAX[1] = 1.;

	XMIN[0] = -1;
	XMIN[1] = 0; 

	dx[0] = (XMAX[0] - XMIN[0])/((double) UNI1);
	dx[1] = (XMAX[1] - XMIN[1])/((double) UNI2);

	for(int i=0;i<UNI1; i++){
		for(int j=0;j<UNI2; j++){
			X[1] =  dx[0] * i + XMIN[0] + 0.5 * dx[0];
			X[2] =  dx[1] * j + XMIN[1] + 0.5 * dx[1];

			get_fluid_params(X, &modvar);
			quant[i*UNI2+j][0]=modvar.B*modvar.B;
			quant[i*UNI2+j][1]=0;

		}
	}

	perform_fft(quant,out);
	fprintf(stderr,"done!\n");
	write_plasma_file(out,TIME_INIT);

	fftw_free(out);
	fftw_free(quant);

}


void write_plasma_file(fftw_complex *quant,double TIME_INIT){

    struct stat st = {0};
    char folder[64] = "output";

    if (stat(folder, &st) == -1) {
        mkdir(folder, 0700);
    }

    char filename[256] = "";
    fprintf(stderr,"%d\n",(int)TIME_INIT);
    sprintf(filename, "%s/FFT_%d.dat", folder,(int)TIME_INIT);

    double intg[UNI1];

    for(int i=0;i<UNI1;i++)
		intg[i]=0;

    FILE *specfile = fopen(filename, "w");

        for(int i=0;i<UNI1; i++){
                for(int j=0;j<UNI2; j++){
//						fprintf(stderr,"%e\n",quant[i][j][0]);
                        intg[i]+=(quant[i*UNI2+j][0])/((double) UNI2);
                }
        }
	double freq[UNI1];
	fftfreq(freq,UNI1,8.);
	for(int i=0;i<UNI1/2;i++){
	        fprintf(specfile, "%+.15e %+.15e\n",freq[i], intg[i]);
	}

}

void perform_fft(fftw_complex* quant, fftw_complex* out){
  fftw_plan p, q;
  int i,j;

  /* forward Fourier transform, save the result in 'out' */
  p = fftw_plan_dft_2d(UNI1, UNI2, quant, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  fftw_destroy_plan(p);

}


void fftfreq(double* freq, int N, double size){
	if(N%2==0){
		for(int i=0;i<N/2;i++){
			freq[i] = i/size;
		}

		for(int i=N/2;i<N;i++){
			freq[i] = (-N/2 + (i-N/2))/size;
		}
	}
	else{
		for(int i=0;i<=(N-1)/2;i++){
			freq[i] = i/size;
		}

		for(int i=(N-1)/2;i<N;i++){
			freq[i] = (-(N-1)/2 + (i-(N-1)/2))/size;
		}

	}
}