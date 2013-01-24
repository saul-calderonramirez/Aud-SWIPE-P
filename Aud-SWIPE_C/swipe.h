#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <fftw3.h>   // http://www.fftw.org/
#include <sndfile.h> // http://www.mega-nerd.com/libsndfile/
//#include "vector.h"  // comes with release
#include "AudSys.h"
	  

#define NOK                          0

#define TRUE                         1
#define FALSE                        0

#define DERBS                        .1 
#define POLYV                        .0013028 //  1 / 12 / 64 = 1 / 768
#define DLOG2P                       .0104167 // 1/96

#define ST                           .0001  // Feel free to change these
#define DT                           .001
#define MIN                          100.
#define MAX                          600.

#define VNUM                         1.0 // Current version

#ifndef NAN                          
    #define NAN                      sqrt(-1.)
#endif



double hz2mel(double hz);

double hz2erb(double hz);

double erb2hz(double erb);

double fixnan(double x);

double maxim(double x, double y);

matrix repmat(vector A, int x, int y);

void interp1(vector f, vector g, vector x, double* y);
vector pitchStrengthOneCandidate( vector f, matrix NL, double pc );
matrix pitchStrengthAllCandidates(vector f, matrix L, vector pc, vector j);
vector pitch(matrix S, vector pc, double st);
vector swipe(char wav[], double min, double max, double st, double dt, char filename[], FILE * time_test);
void printp(vector p, char out[], double dt, int mel, int vlo);
void outBinaryM(double** m,int x, int y, char file[]);
void outBinaryV(double* v, int x, char file[]);


