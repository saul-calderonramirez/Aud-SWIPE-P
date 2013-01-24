	/*
	 * Name: AudSWIPE-P
	 * Authors: Prof. Arturo Camacho
	 * 			Bsc. Saul Calderon
	 * 			Bsc. Gabriel Alvarado
	 * General Description: Parallel C implementation of the AudSWIPE algorithm, originally
	 * implemented in MATLAB
	 * Archive Description: Defines the ERB filters to apply in the first stage of AudSWIPE algorithm
	 */

#include <math.h>
#include <stdio.h>
#include "ERBFilters.h"
#include "complex.h"

/*
 * Half wave rectification
 * Receives the array
 * */
void MaxRCerosVector(double* x, int size_x){
	int  c;
	for(c = 0; c < size_x; c++){
		if(x[c] < 0)
			x[c] = 0;
	}
}

/*
*Writes the given matrix in a file
*@param m, the matrix to write
*@param x, number of rows
*@param y, number of columns
*@param file, file name
*/
void outBinaryM(double** m, int x, int y, char file[]){
	FILE * out;
	out = fopen (file, "wb");
	int i;
	int j;
	for(i = 0; i < x; ++i){
		for(j = 0; j < y; ++j)
		fprintf(out," %f", m[i][j]);
	}
	fprintf(out," \n");
	fclose(out);
}

/*
*Generates the ERB filters to be used, to obtain the inner ear response and the missing harmonics
*@param fs, sample frequency
*@param cf,
*@param fcoefs, array that will contain ERBs coeficients
*@param texto, the name of the array
*/
void ERBFilters(double fs, vector cf, matrix fcoefs){
	double T = 1/fs;
	double EarQ = 9.26449;
	double minBW = 24.7;
	double order = 1;
	vector ERB = zerov(cf.x);
	vector B = zerov(cf.x);
	int i;
	for(i = 0; i < cf.x; i++){
		ERB.v[i] = pow(pow(cf.v[i]/EarQ, order) + pow(minBW, order), 1/order);
		B.v[i] = 1.019*2*pi*ERB.v[i];
		fcoefs.m[0][i] = T;//A0
		fcoefs.m[1][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					+ 2*sqrt(3+pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A11
		fcoefs.m[2][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					- 2*sqrt(3+pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A12
		fcoefs.m[3][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					+ 2*sqrt(3-pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A13
		fcoefs.m[4][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					- 2*sqrt(3-pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A14
		fcoefs.m[5][i] = 0;//A2
		fcoefs.m[6][i] = 1;//B0
		fcoefs.m[7][i] = (-2*cos(2*cf.v[i]*pi*T))/exp(B.v[i]*T);//B1
		fcoefs.m[8][i] = exp(-2*B.v[i]*T);//B2
		fcoefs.m[9][i] = cabs((-2*cexp(4*I*cf.v[i]*pi*T)*T +
								2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
									(ccos(2*cf.v[i]*pi*T)- csqrt(3 - pow(2, (double)3/2))*
									csin(2*cf.v[i]*pi*T)))*
							(-2*cexp(4*I*cf.v[i]*pi*T)*T +
								2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
								(ccos(2*cf.v[i]*pi*T) + csqrt(3 - pow(2, (double)3/2))*
								csin(2*cf.v[i]*pi*T)))*
							(-2*cexp(4*I*cf.v[i]*pi*T)*T +
								2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
								(ccos(2*cf.v[i]*pi*T) -
								csqrt(3 + pow(2, (double)3/2))*csin(2*cf.v[i]*pi*T)))*
								(-2*cexp(4*I*cf.v[i]*pi*T)*T + 2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
           						(ccos(2*cf.v[i]*pi*T) + csqrt(3 + pow(2,(double)3/2))*csin(2*cf.v[i]*pi*T))) /
          						cpow((-2 / cexp(2*B.v[i]*T) - 2*cexp(4*I*cf.v[i]*pi*T) +
           						2*(1 + cexp(4*I*cf.v[i]*pi*T))/cexp(B.v[i]*T)),(double)4));//gain
		
	}
	freev(ERB);
	freev(B);
}

/*
*Applies the ERB filters, and also Upsamples, makes the half wave rectification and finally downsamples the signal, to obtain the missing harmonics
*@param x, input signal
*@param size_x, x's size
*@param fcoefs, array that contains ERBs coeficients
*@param X, the matrix that will contain a signal per row, wich will be the cochlea's response to a specific band 
*/
void ERBFilterBank(vector x, matrix fcoefs, matrix X){
	int i, j;
	double a[3];// = double[3];
	double b[3];// = double[3];
	matrix Y = zerom( fcoefs.y, x.x*2);
	//Filter-bank parallelization, One filter - One Thread
	#pragma omp parallel for private(a,b)
	for(i = 0; i < fcoefs.y; i++){
		a[0] = fcoefs.m[0][i]/fcoefs.m[9][i];
		a[1] = fcoefs.m[1][i]/fcoefs.m[9][i];
		a[2] = fcoefs.m[5][i]/fcoefs.m[9][i];
		b[0] = fcoefs.m[6][i];
		b[1] = fcoefs.m[7][i];
		b[2] = fcoefs.m[8][i];
		Filter_AB(a, b, 3,  x.v, x.x, X.m[i]);
		a[0] = fcoefs.m[0][i];
		a[1] = fcoefs.m[2][i];
		a[2] = fcoefs.m[5][i];
		Filter_AB(a, b, 3,  X.m[i], x.x, X.m[i]);
		a[1] = fcoefs.m[3][i];
		Filter_AB(a, b, 3, X.m[i], x.x, X.m[i]);
		a[1] = fcoefs.m[4][i];
		Filter_AB(a, b, 3, X.m[i], x.x, X.m[i]);
		//Half wave rectification procedure, to obtain the missing harmonics
		Upsample_Vector(X.m[i], Y.m[i], x.x*2);
		MaxRCerosVector(Y.m[i], x.x*2 );
		Downsample_Vector(Y.m[i], X.m[i], x.x );
	}
}






