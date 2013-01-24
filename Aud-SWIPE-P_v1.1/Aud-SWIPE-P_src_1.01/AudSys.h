/*
 * Name: AudSWIPE-P
 * Authors: Prof. Arturo Camacho
 * 			Bsc. Saul Calderon
 * 			Bsc. Gabriel Alvarado
 * General Description: Parallel C implementation of the AudSWIPE algorithm, originally
 * implemented in MATLAB
 * Archive Description: Applies the filter of the auditive system and produces differents signals
 * in response to the different cochlea's segments.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include "resample.h"
#include <fftw3.h>   // http://www.fftw.org/
#include <sndfile.h> // http://www.mega-nerd.com/libsndfile/
#include <time.h>
#include "fir2.h"
#include "ERBFilters.h"


/*
*Prints the received array
*@param arr, the array to print
*@param size, the size of the array
*@param texto, the name of the array
*/
void print(double* arr, int size, char* texto);

/*
*Prints the received array in an archive
*@param arr, the array to print
*@param size, the size of the array
*@param texto, the name of the array
*/
void out(double* arr, int size, char* arch);

/*
*Calculates the maximum element in a given array
*@param r, the array
*@param l, the maximum element to take into account in the search
*/
double Max_v(vector r);

/*
*Aligns the channels of the signal, necessary due to the different response times of each band of frequencies on the cochlea
*@param X, the array
*@param f,
*@param fs, the sampling frequency
*@return double**, the matrix with the channels aligned
*/
matrix Align_Channels(matrix X, vector f, double fs);

/*
*Calculates the maximum element in a given array
*@param r, the array
*@param l, the maximum element to take into account in the search
*/
void Max(double** x, int size_x, int size_c, double lim);

/*
*Converts the hz scale value in to a ERB scale value
*@param hz, the value in hz
*@return double, the value in ERB scale
*/
double hz2erbs(double hz);

/*
*Converts the ERB scaled array in to a HZ scale array
*@param erbs, the array in ERBS
*@param size, the array's size
*/
void v_erbs2hz(vector erbs);


/*
*Separates groups of consecutive harmonics into channels
*It uses the ERB filters to produce each channel as a response of a different cochlea's segment
*@param x, the matlab array x
*@param fs, the sampling frequency
*/
matrix Harmonics_into_channels(vector x, double fs, vector f);


/*
*Filter function that receives the coeficients of  the transfer function
*MATLAB like filter function
*@param b, filter coeficients
*@param L,
*@param x, input signal
*@param y, output signal
*@param tam, the size of the vector x
*/
void Filter_B(double * b, int L, double * x  ,int size_x, double* y);

/*
*Filter function that receives the coeficients of A and B's polynoms from the transfer function
*MATLAB like filter function
*@param b, b's filter coeficients
*@param a, a's filter coeficients
*@param L,
*@param x, input signal
*@param y, output signal
*@param tam, the size of the vector x
*/
void Filter_AB(double * b, double * a, int L, double * x  ,int size_x, double* y);

/*
*Normalizes an array, with the maximum at f[size-1]
*@param f, array of frecuencies
*@param size, f's size
*/
void normalize(vector f);

/*
*Uses lineal interpolation to interpolate the function given in f
*@param f, function to interpolate
*@param g,
*@param x,
*@param end,
*@return double, the interpolated function
*/
double interp(double* f, double* g, double x, int end);

/*
*Calculate the filter of the middle and outter ear
*Flattens the spectral envelope
*@param n,
*@param fs, the sampling frequency
*@return vector, the resulting filter
*/
vector outmidear( int n, double fs);

/*
*Models the auditive system, archive's main function
*@param x, input signal
*@param samplerate, signal samplerate
*@param time1,
*@return matrix, A matrix with a signal per row, representing the cochlea's response to a different frequency band.
*/
matrix audsys(vector x, double samplerate, clocksArray *time1);

