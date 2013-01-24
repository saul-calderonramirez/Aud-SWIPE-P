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
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "MPI_Matlab_Communicator.h"
#define pi 3.1415926535897932384626433832795

/*
*Generates the ERB filters to be used 
*@param fs, sample frequency
*@param cf, 
*@param size_cf, cf's size
*@param fcoefs, array that will contain ERBs coeficients
*@param texto, the name of the array
*/
void ERBFilters(double fs, vector cf,  matrix fcoefs);

/*
*Applies the ERB filters, and also Upsamples, makes the half wave rectification and finally downsamples the signal
*@param x, input signal
*@param size_x, x's size
*@param fcoefs, array that contains ERBs coeficients
*@param X, the matrix that will contain a signal per row, wich will be the cochlea's response to a specific band 
*/
void ERBFilterBank(vector x, matrix fcoefs, matrix X);
