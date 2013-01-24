/*
 * Name: AudSWIPE-P
 * Authors: Prof. Arturo Camacho
 * 			Bsc. Saul Calderon
 * 			Bsc. Gabriel Alvarado
 * General Description: Parallel C implementation of the AudSWIPE algorithm, originally
 * implemented in MATLAB
 * Archive Description: Resampling is necessary to make the half wave rectificaction
 * and recover the missing harmonics.  Contains Resample filters
 */

#include "resample.h"
/*
*Upsample the signals given in the X matrix.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Upsampled matrix
*@param size_x,
*@param size_c,
*/
void Upsample(double** X, double** Y, int size_x, int size_c){

	int tam_h = 42;
	double  h[] = {0,
		-1.4311e-18,-0.0021018,3.7067e-18,0.0050153,
		-6.9847e-18,-0.0097846,1.1207e-17,0.017104,
		-1.6174e-17,-0.027965,2.1551e-17,0.044023,
		-2.6905e-17,-0.068645,3.1753e-17,0.11052,
		-3.5622e-17,-0.20176,3.8119e-17,0.63307,1,
	0.63307,3.8119e-17,-0.20176,-3.5622e-17,0.11052,
	3.1753e-17,-0.068645,-2.6905e-17,0.044023,2.1551e-17,
	-0.027965,-1.6174e-17,0.017104,1.1207e-17,-0.0097846,
	-6.9847e-18,0.0050153,3.7067e-18,-0.0021018,-1.4311e-18 };
	int f, c, i;
	//Fills with 0's	
	for(f = 0; f < size_c; f++){
	i = 0;
		for(c = 0; c < size_x; c = c+2){
			Y[f][c] = X[f][i]; 
			i++;
		}
	}
	//Applies the filter	
	for(i = 0; i < size_c; i++){
		Filter_B(h, tam_h, Y[i], size_x, Y[i]);
	}
}

/*
*Upsample the signals given in the X vector.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Upsampled matrix
*@param size_x
*/
void Upsample_Vector(double* X, double* Y, int size_x){
	int tam_h = 42;
	double  h[] = {0,
		-1.4311e-18,-0.0021018,3.7067e-18,0.0050153,
		-6.9847e-18,-0.0097846,1.1207e-17,0.017104,
		-1.6174e-17,-0.027965,2.1551e-17,0.044023,
		-2.6905e-17,-0.068645,3.1753e-17,0.11052,
		-3.5622e-17,-0.20176,3.8119e-17,0.63307,1,
	0.63307,3.8119e-17,-0.20176,-3.5622e-17,0.11052,
	3.1753e-17,-0.068645,-2.6905e-17,0.044023,2.1551e-17,
	-0.027965,-1.6174e-17,0.017104,1.1207e-17,-0.0097846,
	-6.9847e-18,0.0050153,3.7067e-18,-0.0021018,-1.4311e-18 };
	int f, c, i;
	//Fills with 0's
	i = 0;
	for(c = 0; c < size_x; c = c+2){
		Y[c] = X[i];
		i++;
	}
	//Applies the filter
	Filter_B(h, tam_h, Y, size_x, Y);
}


/*
*Reduces the sample rate of the signals given in the X Vector.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Downsampled matrix
*@param size_x
*/
void Downsample_Vector(double* X, double* Y, int size_x){
	int f, c, i;
	int tam_h = 43;
	double  h[] = {0,0,-7.1553e-19,-0.0010509,1.8533e-18,0.0025077,
	-3.4923e-18,-0.0048923,5.6036e-18,0.0085521,-8.0868e-18,
	-0.013983,1.0776e-17,0.022012,-1.3453e-17,-0.034322,1.5876e-17,
	0.05526,-1.7811e-17,-0.10088,1.9059e-17,0.31654,0.5,0.31654,1.9059e-17,
	-0.10088,-1.7811e-17,0.05526,1.5876e-17,-0.034322,-1.3453e-17,0.022012,
	1.0776e-17,-0.013983,-8.0868e-18,0.0085521,5.6036e-18,-0.0048923,-3.4923e-18,
	0.0025077,1.8533e-18,-0.0010509,-7.1553e-19 };
	//Applies the filter
	Filter_B(h, tam_h, X, size_x*2, X);
	//Removes samples
	i = 0;
	for(c = 0; c < size_x*2; c = c+2){
		Y[i] = X[c];
		i++;
	}
}

/*
*Reduces the sample rate of the signals given in the X matrix.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Downsampled matrix
*@param size_x,
*@param size_c,
*/
void Downsample(double** X, double** Y, int size_x, int size_c){
	int f, c, i;
	int tam_h = 43;
	double  h[] = {0,0,-7.1553e-19,-0.0010509,1.8533e-18,0.0025077,
	-3.4923e-18,-0.0048923,5.6036e-18,0.0085521,-8.0868e-18,
	-0.013983,1.0776e-17,0.022012,-1.3453e-17,-0.034322,1.5876e-17,
	0.05526,-1.7811e-17,-0.10088,1.9059e-17,0.31654,0.5,0.31654,1.9059e-17,
	-0.10088,-1.7811e-17,0.05526,1.5876e-17,-0.034322,-1.3453e-17,0.022012,
	1.0776e-17,-0.013983,-8.0868e-18,0.0085521,5.6036e-18,-0.0048923,-3.4923e-18,
	0.0025077,1.8533e-18,-0.0010509,-7.1553e-19 };
	//Applies the filter
	for(i = 0; i < size_c; i++){
		Filter_B(h, tam_h, X[i], size_x*2, X[i]);
	}
	//Removes samples
	for(f = 0; f < size_c; f++){
	i = 0;
		for(c = 0; c < size_x*2; c = c+2){
			Y[f][i] = X[f][c]; 
			i++;
		}
	}
}
