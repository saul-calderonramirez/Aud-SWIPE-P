

#include <fftw3.h>   // http://www.fftw.org/
#include <sndfile.h> // http://www.mega-nerd.com/libsndfile/
#include "AuditiveSystem.h"
#include "includes.h"
//#include "MPIMatlabCommunicator.h"




char* ptrTestFile;
double soundLength;

//MPI variables needed
int myid;

double timeCochlea, timeSendData, timeMisc, timeWsCycles, timePitchRec, timeTotal;
double timesAudsys[7];
/*
 * Executes the PrimeMultiF0
 * @param argc, number of arguments
 * @param argv, list of arguments
 * */
void executePrimeMultiF0(int argc, char* argv[]);
/*
 * Calculates the spectogram of the zero padded signal in parallel, with one thread per window
 * @param zpSignal, zeropadded signal
 * @param ws, window size
 * @param w, window, could be a hanning window
 * @param woverlap, window overlap
 * @param f, frequencies vector, filled by the specgram function
 * */
matrix specgram(vector zpSignal, int ws, double fs, vector w, int woverlap, double dn, vector f, vector ti);

/*
* Primary utility function for multi pitch extraction
* @param wav, name of the input signal
* @param paralelism, Paralelism level
* @param min, Start of frequency range
* @param max, End of frequency range
* @param dt, time delta
* @return vector, pitch array
*/
vector primeMultiF0(char wav[], double min, double max, double st, double dt, char* argv[], double paralelism, char* testName);

//matrix primeMulti_F0(vector x, double fs, vector pc, double dt);
/*
void getWsScoreMat(int i, vector ws, vector pc, vector d, matrix S, double fs, vector p0, vector t);*/
/*
 * Gets the pitch candidates associated to the current window size
 * @param i, current window size in the window size array
 * @param ws, window sizes array
 * @param d,
 * @param ptrJ
 * @param ptrK
 * */
void getPitchCandidatesOfWs(int i, vector ws, vector d, vector* ptrJ, vector* ptrK);

/*
*Returns the pitch strength of several pitch candidate, applying the kernel to the matrix of signals
*@param f, the erb scale
*@param L, the loudness matrix of signals
*@param pc, the pitch candidate
*@param j,
*@return matrix, computed pitch strenghts per pitch candidate
*/
matrix pitchStrengthAllCandidates(vector f, matrix L, vector pc, vector j);

/*
 * Receives the score matrix S and applies post processing required by prime multi F0 algorithm
 * @param S, score matrix, rows are pitch candidates, columns are the values in time
 */
matrix postprocessS(matrix S, vector pc);
