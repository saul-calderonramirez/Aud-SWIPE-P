/*
 * Clocks.c
 *
 *  Created on: Apr 23, 2013
 *      Author: saul
 */

#include "Clocks.h"


double totalTime;
/*
 * Calculates the time difference
 * @param timeA_p, start
 * @param timeB_p, end
 * */
double timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p){
  return ((double)(timeA_p->tv_sec ) + (double)timeA_p->tv_nsec/(double)1000000000) -
           (((double)timeB_p->tv_sec ) + (double)timeB_p->tv_nsec/(double)1000000000);
}

/*
 * Gets the local clock started
 * */
void startLocalClock(clocksArray* arr1){
	clock_gettime(CLOCK_MONOTONIC, &(arr1 -> start));
}

/*
 * Inits a clocks array
 * @param numClocks, number of clocks to init in the array
 * */
clocksArray initClocks(int numClocks){
	clocksArray clocks;
	clocks.values = zerov(numClocks);
	clocks.currentPosition = 0;
	clocks.tags = malloc(sizeof(char*) * numClocks);

	return clocks;
}

/*
 * Generates an array of a union of the two arrays received
 *@param arr1, first array
 *@param arr2, second array
 * */
clocksArray unionArrays(clocksArray* arr1, clocksArray* arr2){
	clocksArray unionArr;
	unionArr = initClocks(arr1 -> values.x + arr2 -> values.x);
	int i = 0;
	int u = 0;
	for(i = 0; i < arr1->values.x; ++i){
		unionArr.values.v[u] = arr1 ->values.v[i];
		unionArr.tags[u++] = arr1 ->tags[i];
	}
	for(i = 0; i < arr2->values.x; ++i){
		unionArr.values.v[u] = arr2 ->values.v[i];
		unionArr.tags[u++] = arr2 ->tags[i];
	}
	return unionArr;
}

/*
 * Writes the array of clocks to a file
 * */
void writeTimesToFile(char* testName, clocksArray* clocks, int numProcs, double soundLength, char* wav){
	if(testName != NULL){
		FILE * time_test = NULL;
		time_test = fopen (testName, "r");

		//must write header
		//if file did not exist
		int i = 0;
		if(time_test == NULL){
			time_test = fopen (testName, "w+");
			fprintf(time_test, "File;");
			for(i = 0; i <  clocks->values.x; ++i){
				fprintf(time_test, "%s;", clocks->tags[i]);

			}

			fprintf(time_test, "NumProcs;");
			fprintf(time_test, "soundLength;");
			fprintf(time_test, "Diff;");
			fprintf(time_test, "SR\n");

		}
		else{
			time_test = fopen (testName, "a+");
		}

		fprintf(time_test, "%s ;", wav);
		for(i = 0; i < clocks->values.x; ++i){
			fprintf(time_test, " %lf ; ", clocks->values.v[i]);
		}
		fprintf(time_test, " %d ; ", numProcs);
		fprintf(time_test, " %lf ; ", soundLength);
		double diff = soundLength - clocks->values.v[clocks->values.x-1];

		int pass = diff > 0;

		fprintf(time_test, " %lf ; ", diff);
		fprintf(time_test, "%d\n", pass);

		fclose(time_test);
	}
}


/*
 * Stops a local clock
 * @param isFinal, 1 if is a final time to print, or 0 if a sum of times is going to be performed
 * */

double endLocalClock(clocksArray* clocks, int isFinal, char* tag ){
	clock_gettime(CLOCK_MONOTONIC, &(clocks -> end));
	double timeDiff = timespecDiff(&(clocks -> end), &(clocks -> start));
	clocks->tags[clocks->currentPosition] = tag;

	if(isFinal == 1 && clocks->currentPosition < clocks->values.x ){
		//printf("Time: %f %s pos: %d\n", timeDiff, tag, clocks->currentPosition );
		clocks->values.v[clocks->currentPosition++ ] = timeDiff;
	}
	else{
		//printf("Time: %f %s pos: %d\n", timeDiff, tag, clocks->currentPosition );
		clocks->values.v[clocks->currentPosition] += timeDiff;
		if(clocks->currentPosition != clocks->values.x - 1){
			clocks->currentPosition++;
		}
		else{
			clocks->currentPosition = 0;
		}
	}
	//by default it starts the clock again
	startLocalClock(clocks);
	return timeDiff;
}

