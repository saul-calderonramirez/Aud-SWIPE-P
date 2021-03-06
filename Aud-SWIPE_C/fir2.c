/*
    Copyright (c) 2008, Alexander Iakovlev

    VIORD MicroSCADA is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    See file '../Copyright' for full details.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>

//   FIR filter design using the window method - arbitrary filter shape.
//   FIR2(N,F,M,NF,B) designs an N'th order FIR digital filter with the
//   frequency response specified by vectors F and M (size NF) and returns the
//   filter coefficients in length N+1 vector B.  Vectors F and M specify
//   the frequency and magnitude breakpoints for the filter such that
//   PLOT(F,M) would show a plot of the desired frequency response.
//   The frequencies in F must be between 0.0 < F < 1.0, with 1.0
//   corresponding to half the sample rate. They must be in increasing
//   order and start with 0.0 and end with 1.0. 
//
//   The filter B is real, and has linear phase, i.e., even symmetric 
//   coefficients obeying B(k) =  B(N+2-k), k = 1,2,...,N+1.
//
//   FIR2 uses a Hamming window.

#define pi 3.1415926535897932384626433832795

void printfir(double* arr, int size, char* texto){
	
	int i;
	for(i = 0; i < size; i++){
		printf("%f \n", arr[i]);
	}
	printf("\n%s\n",texto);
	

}

void Hanning(double *y,int np){
	int i;
	int a = 0;
     for(i=1; i<=np; i++, a++){ 
     	y[a] = .5*(1 - cos((2*pi*i)/(np+1)));//Hanning en lugar de Hamming
	}
}

void Hanning2(double *y,int np){
		int i;
        for(i=0; i<np; i++) 
                y[i] = .54 - .46*cos((2*pi*i)/np);
}

void diff(double *y,double *f,int sz)
{
    --sz;
    int i;
    for(i=0;i<sz;i++)
        f[i]=y[i+1]-y[i];
}

int isneg(double *y,int sz)
{
	int i;
    for(i=0;i<sz;i++)
            if (y[i] < 0) return 1;
    return 0;
}

double fix(double n){return (int)n;}
 
void fir2(int nn, double *ff,double *aa, int ffsz, double *b)
{ 
		//printf("size ffirC\n");
		//out(ff, ffsz, "ffirC.txt");


int npt,lap;
        nn = nn + 1;
        if (nn < 1024)
                npt = 512;
        else
                npt = pow(2,ceil(log(nn)/log(2)));
        lap = fix(npt/25);
   
        int nbrk=ffsz;
        if (abs(ff[0]) > 0 || abs(ff[nbrk - 1] - 1) > 1){
           perror("The first frequency must be 0 and the last 1");
           abort();
        }
// interpolate breakpoints onto large grid
  		
        int nint=nbrk-1;
        double *df = (double *)malloc(sizeof(double)*(ffsz-1));//new double[ffsz];
        diff(ff,df,ffsz);
        //printfir(ff, ffsz, "ff");
        //printfir(df, ffsz-1, "df");
        if (isneg(df,ffsz-1)){
           perror("Frequencies must be non-decreasing");
           abort();
        }
        npt = npt + 1;   // Length of [dc 1 2 ... nyquist] frequencies.
		
        int  nb = 0;
        int  ne=0;
        double  inc;
        double *H = (double*)malloc(sizeof(double)*npt);//new double[npt];// last index overflaw !!!!!
        H[0]=aa[0];
        int i;
        for(i=0;i<nint;i++)
        {
            if (df[i] == 0){
               nb = nb - lap/2;
               ne = nb + lap;
            }else
               ne = fix(ff[i+1]*npt)-1;
            if (nb < 0 || ne > npt){
                   perror("Too abrupt an amplitude change near end of frequency interval");
                   abort();
            }
			int j;
           for(j=nb;j<=ne;j++){
                if (nb == ne)
                        inc = 0;
                else
                        inc = (double)(j-nb)/(double)(ne-nb);
                H[j] = inc*aa[i+1] + (1.0 - inc)*aa[i];
           }
           nb = ne + 1;
        }
	
// Fourier time-shift.
        
        double dt = 0.5 * (nn - 1);
        double complex *Hz1 = (double complex*)malloc(sizeof(double complex)*npt);//new complex[npt];
        double complex *Hz2 = (double complex*)malloc(sizeof(double complex)*npt);//new complex[npt];
        //double complex *Hz = (double complex*)malloc(sizeof(double complex)*2*npt);//new complex[2*npt];
        fftw_complex* Hz = fftw_malloc(sizeof(fftw_complex) * 2*npt); 
        int sc;
        for(i=0;i<npt;i++){
                double rad = -dt * pi * (double)i / ((double)npt-1);
                double complex z = (H[i]*cos(rad)) +(H[i]*sin(rad)*I);
                Hz1[i]=z;
        }
                        

        for(i=0;i<npt;i++)
                Hz2[npt-1-i]=conj(Hz1[i]);


        for(i=0;i<npt;i++)
                Hz[i]=Hz1[i];
        for(i=npt;i<npt*2;i++)
                Hz[i]=Hz2[i-npt];
		
	   int nfft=npt*2-2;
        double* fo = fftw_malloc(sizeof(double) * 2*npt);
        fftw_plan plan = fftw_plan_dft_c2r_1d(nfft, Hz, fo, FFTW_ESTIMATE); 
        fftw_execute(plan);
        double *wind= (double*)malloc(sizeof(double)*nn);//new double[nn];
        Hanning2(wind,(2*ceil(nn/2))+1);//En lugar de llamar a Hamming, lo modificamos para que sea Hanning, mandamos otra cantidad de puntos ( 2*ceil(n/2)+1 )
  
        double kfft=1./nfft;
        for(i=0;i<nn;i++)
                b[i]=wind[i]*fo[i]*kfft;
        
        fftw_destroy_plan(plan); 
    	fftw_free(Hz); 
    	free(wind);
    	fftw_free(fo); 
        free(df);
        free(H);
        free(Hz1);
        free(Hz2);

}

