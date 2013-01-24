/*
 * Name: AudSWIPE
 * Authors: Prof. Arturo Camacho
 * 			Bsc. Saul Calderon
 * 			Bsc. Gabriel Alvarado
 * General Description: C implementation of the AudSWIPE algorithm, originally
 * implemented in MATLAB
 * Archive Description: Resampling filters, necessary to apply half wave rectification, and recover missing harmonics
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>

/*
*Reduces the sample rate of the signals given in the X matrix.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Downsampled matrix
*@param size_x,
*@param size_c,
*/
void Downsample(double** X, double** Y, int size_x, int size_c);

/*
*Upsample the signals given in the X matrix.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Upsampled matrix
*@param size_x,
*@param size_c,
*/
void Upsample(double** X, double** Y, int size_x, int size_c);

/*
*Reduces the sample rate of the signals given in the X vector.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Downsampled matrix
*@param size_x
*/
void Downsample_Vector(double* X, double* Y, int size_x);

/*
*Upsample the signals given in the X vector.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Upsampled matrix
*@param size_x
*/
void Upsample_Vector(double* X, double* Y, int size_x);
