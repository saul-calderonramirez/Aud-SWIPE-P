

#include <fftw3.h>   // http://www.fftw.org/
#include <sndfile.h> // http://www.mega-nerd.com/libsndfile/
#include "AuditiveSystem.h"
#include "includes.h"
#include "MPIMatlabCommunicator.h"
#define EXE_NAME "aswipep"

struct timespec start, end, startTotal, endTotal;
clocksArray clocks, cycleClocks, totalTime;

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
 * */
matrix specgram(vector zpSignal, int ws, double fs, vector w, int woverlap, double dn);

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
 * Gets the Score matrix S for the current window size
 * @param i, current window size in the window size array
 * @param pc, Pitch candidates array
 * @param X, matrix with the cochleas response
 * @param fERBs, frequency array in the ERB scale
 * @param d,
 * @param S, Score matrix
 * @param fs, sampling frequency
 * @param p0,
 * @param t, time vector
 * */
void getWsScoreMat(int i, vector ws, vector pc, matrix X, vector fERBs, vector d, matrix S, double fs, vector p0, vector t);
/*
 * Gets the Loudness Matrix associated to a window size
 * @param i, current number of window size, in the window sizes array
 * @param ws, window sizes array
 * @param n,
 * @param Xz,
 * @param w,
 * @param fERBs, ERB frequency array
 * @param fi,
 * @param W,
 * @param tamX, size in x of the cochleas Matrix
 * */
matrix getWsLoudnessMat(int i, vector ws, vector n, matrix Xz, vector w, vector fERBs, vector fi, matrix W, int tamX);

/*
 * Creates the matrix of channels X
 * Must be done by the parent only, after this, it is sent to the slaves
 * @param wav, wav file name
 * @param min, minimum freq. to look for
 * @param max, max. freq. to look for
 * @param dt, pitch detection time interval
 * */
matrix calculateX(char wav[], double min, double max, double dt,  char* argv[], double paralelism, matrix* X, vector* fsInfo );
/*
*Peforms polynomial tuning on the strength matrix to determine the pitch
*@param f, the erb scale
*@param NL, the normalized loudness matrix of signals
*@param pc, the pitch candidate
*@return vector, computed pitch strenghts
*/
vector pitch(matrix S, vector pc, double st);

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
matrix processMatrixSprimeMultiF0(matrix S, vector pc);
