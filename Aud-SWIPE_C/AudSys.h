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
//#include "vector.h"  // comes with release
#include "fir2.h"
//#include "ERBFilters.h"


void print(vector vec, char* texto);

void out(vector vec, char* arch);

double Max_v(vector r);

//double** Align_Channels(double ** X, int size_x, int size_c, double *f, double fs, int *tam_y);
matrix Align_Channels(matrix X, vector f, double fs);
void Max(matrix X);

double hz2erbs(double hz);

void v_erbs2hz(vector erbs/*double* erbs, int size*/);

matrix Harmonics_into_channels(vector x, double fs,/* int tam,*/ vector f);

void Filter_B(vector b, double* x  ,int size_x, double* y);

//Funcion que imita Filter de Matlab con coeficientes A y B
void Filter_AB(vector b, vector a, /*int L,*/ double* x  ,int size_x, double* y);

//Funcion para normalizar un vector, suponemos que la frecuencia minima está en f[0] y la maxima en f[size-1]; 
void normalizar(vector f);

//Interpolación lineal
double interp(double* f, double* g, double x, int end);

//Funcion que calcula el filtro del oido medio y externo
//Flatten the spectral envelope
vector outmidear( int n, double fs);


//funcion que modela el sistema auditivo, desde aqui se llaman a todas las demas funciones de este archivo
matrix audsys(vector x, double samplerate, double *time1);
double timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p);
