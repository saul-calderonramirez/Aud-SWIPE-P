/*
 * Clocks.h
 *
 *  Created on: Apr 23, 2013
 *      Author: saul
 */

#ifndef CLOCKS_H_
#define CLOCKS_H_

#include "includes.h"


typedef struct {
	struct timespec start, end;
	vector values;
	int currentPosition;
	char** tags;
} clocksArray;

/*
 * Gets the local clock started
 * @param arr1, array of clocks to start
 * */
void startLocalClock(clocksArray* clocks);

/*
 * Inits a clocks array
 * @param numClocks, number of clocks to init in the array
 * @return array of clocks
 * */
clocksArray initClocks(int numClocks);
/*
 * Stops a local clock
 * @param isFinal, 1 if is a final time to print, or 0 if a sum of times is going to be performed
 * @param clocksArray clocks, clocks to finish
 * @param tags, array of tags
 * */
double endLocalClock(clocksArray* clocks, int isFinal, char* tag);
/*
 * Generates an array of a union of the two arrays received
 *@param arr1, first array
 *@param arr2, second array
 *@param clocksArray an array with the union of the two arrays
 * */
clocksArray unionArrays(clocksArray* arr1, clocksArray* arr2);
/*
 * Writes the array of clocks to a file
 * @param testName, name of the test file
 * @param clocks, array of clocks to write
 * @param numProcs, number of Processes that executed the piece of code
 * @param soundLength
 * @param wav, corresponding wav file
 * */
void writeTimesToFile(char* testName, clocksArray* clocks, int numProcs, int winSizes, int fs, double soundLength, char* wav);
/*
 * Calculates the time difference
 * @param timeA_p, start
 * @param timeB_p, end
 * @return time difference
 * */
double timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p);
#endif
