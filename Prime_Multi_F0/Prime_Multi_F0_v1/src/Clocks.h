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
void startLocalClock(clocksArray* clocks);

clocksArray initClocks(int numClocks);

double endLocalClock(clocksArray* clocks, int isFinal, char* tag);
clocksArray unionArrays(clocksArray* arr1, clocksArray* arr2);
void writeTimesToFile(char* testName, clocksArray* clocks, int numProcs, double soundLength, char* wav);
double timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p);
#endif
