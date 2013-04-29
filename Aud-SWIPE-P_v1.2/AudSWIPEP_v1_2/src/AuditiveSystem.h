/*
 * AuditiveSystem.h
 * Models the auditive system response
 *  Created on: Apr 23, 2013
 *      Author: saul
 */
#include "includes.h"
#include "Clocks.h"
#include "SPUtilities.h"
#ifndef AUDITIVESYSTEM_H_
#define AUDITIVESYSTEM_H_
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
/*
*Separates groups of consecutive harmonics into channels
*It uses the ERB filters to produce each channel as a response of a different cochlea's segment
*@param x, the matlab array x
*@param fs, the sampling frequency
*/
matrix Harmonics_into_channels(vector x, double fs, vector f);
/*
*Aligns the channels of the signal, necessary due to the different response times of each band of frequencies on the cochlea
*@param X, the array
*@param f,
*@param fs, the sampling frequency
*@return double**, the matrix with the channels aligned
*/
matrix alignChannels(matrix X, vector f, double fs);
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

#endif
