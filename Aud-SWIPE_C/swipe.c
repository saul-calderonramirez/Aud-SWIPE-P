/* Copyright (c) 2009-2011 Kyle Gorman
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*
*  SWIPE' pitch estimator
*  Camacho, Arturo. A sawtooth waveform inspired pitch estimator for speech   
*  and music. Doctoral dissertation, University of Florida. 2007.             
*  Implemented in C by Kyle Gorman <kgorman@ling.upenn.edu>
*
*  HOW TO CITE:
*
*  Please cite this dissertation, and if possible include a footnote link to the
*  source of this program, the most-recent version of which will always be at:
*
*  http://ling.upenn.edu/~kgorman/c/swipe/
*
*  This program depends on several free ("libre", not "gratuit") libraries. To 
*  obtain them,  follow the instructions below for your platform.
*
*  LINUX:
*
*  All the large libraries should be available as packages if you're using a 
*  "modern" distro. For instance, on a current Debian/Ubuntu system (Ubuntu 
*  9.04, "Jaunty Jackalope", kernel 2.6.28-13-generic), run (as superuser):
*
*  apt-get install libblas-dev liblapack-dev libfftw3-dev libsndfile1-dev
* 
*  This installs the BLAS, (C)LAPACK, fftw3, and sndfile libraries. Installing 
*  the most recent packages on a Fedora, Slackware, etc. should have a similar
*  effect, assuming dependencies are satisfied in the process.
*
*  MAC OS X:
*
*  The linear algebra libraries ([C]LAPACK, a BLAS implementation) ship with Mac
*  OS X. You will need to install the newest versions of fftw3 and libsndfile, 
*  however. They are available for free online:
*
*  http://www.fftw.org/
*  http://www.mega-nerd.com/libsndfile/
*
*  If you are superuser and wish to install globally the autoconf method
*  should work fine:
*
*  tar -xvzf downloadedPackage.tar.gz
*  cd folderOfPackageCreatedByTAR/
*  ./configure; make; make install;
*
*  If you're not superuser, or don't want to install globally, make sure to 
*  use '--prefix=PATH/TO/LOCATION' as an argument to 'configure'. You may 
*  need to alter the #include statements as well. 
*
*  WINDOWS/CYGWIN:
*
*  Unsupported. Send details of any successes, however.
*
*  THANKS:
*  Arturo Camacho, Stephen Isard, Mark Liberman, Chandan Narayan, Dan Swingley
*/

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

#define ST                           .000000001  // Feel free to change these
#define DT                           .001
#define MIN                          100.
#define MAX                          600.

#define VNUM                         1.0 // Current version

#ifndef NAN                          
    #define NAN                      sqrt(-1.)
#endif

#ifndef isnan
int isnan(double x) { 
    return(x !=x);
}
#endif
struct timespec start, end, startTotal, endTotal;
#ifndef log2
double log2(double x) { // A base-2 log function
    return(log(x) / log(2.));
}
#endif

#ifndef round
double round(double x) { // Rounds a double to the nearest integer value
    return(x >= 0. ? floor(x + .5) : floor(x - .5));
}
#endif

double hz2mel(double hz) { // Converts from hertz to Mel frequency
    return(1127.01048 * log(1. + hz / 700.));
}

double hz2erb(double hz) { // Converts from hertz to ERBs
    return(21.4 * log10(1. + hz / 229.));
}

double erb2hz(double erb) { // Converts from ERBs to hertz 
    return((pow(10, erb / 21.4) - 1.) * 229.);
}

double fixnan(double x) { // A silly function that treats NaNs as 0.
    return(isnan(x) ? 0. : x);
}

double maxim(double x, double y){
	return(x>=y ? x : y);
}
double time1, time2, time3;
//
void outBinaryM(double** m, int x, int y, char file[]){
	
	FILE * out;
	out = fopen (file, "wb");
	int i;
	for(i = 0; i < x; ++i){
		fwrite(m[i], sizeof(double), y, out);
	}
	fclose(out);
	
	
	
}

void outBinaryV(double* v, int x, char file[]){
	FILE * out;
	out = fopen (file, "wb");
	fwrite(v, sizeof(double), x, out);
	
	fclose(out);
}



void interp1(vector f, vector g, vector x, double* y){
	int i = 0;
	int a = 0;
	double x0 ;
	double x1;
	double y0;
	double y1;
	x0 = 0;
	x1 = f.v[i];
	y0 = 0;
	y1 = g.v[i];
	while(x.v[a]<f.v[i] && a < x.x){//Esta seccion calcula cuando los indices de x son menores que los de f
		y[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Formula de interpolacion lineal
		a++;
	}
	while(i < f.x && a < x.x){//Esta seccion calcula cuando los valores de x estan en el rango de f
		if(f.v[i] < x.v[a]){
			i++;
		}else{
			x0 = f.v[i-1];
			x1 = f.v[i];
			y0 = g.v[i-1];
			y1 = g.v[i];
			y[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Formula de interpolacion lineal
			a++;
		}
	}
	x0 = f.v[i-1];
	x1 = 0;
	y0 = g.v[i-1];
	y1 = 0;
	if(a < x.x){
		y[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Formula de interpolacion lineal
		a++;
	}
	while(a < x.x){//Esta seccion calcula cuando los valores de x estan por encima del rango de f
		y[a] = 0;//Formula de interpolacion lineal
		a++;
	}
}

vector pitchStrengthOneCandidate( vector f, matrix NL, double pc ){
	int n,i,j;
	vector S = zerov(NL.x);
	double p = f.v[f.x-1]/pc -0.75;//Number of harmonics
	vector q = copyv(f);
	if(p >= 0){
		n = floor(p);
	}else{
		n = ceil(p);
	}
	if(n == 0){
		return S;
	}
	vector k = zerov(f.x);//Kernel
	for(i = 0; i < q.x; i++){
		q.v[i] = q.v[i]/pc;//Normalize frequency w.r.t. candidate
	}
	intvector primos = primes(n);//This is A-SWIPE'
	vector prim = zerov(primos.x+1);
	for(i = primos.x-1;i >= 0; i--){
		prim.v[i+1] = primos.v[i];
	}
	prim.v[i+1] = 1;
	double val;
	for(i = 0; i < prim.x; i++){///PARALELIZAR AQUI **
		for(j = 0; j < q.x; j++){
			val = fabs(q.v[j]-prim.v[i]);
			if(val < 0.25){//Peaks weights
				k.v[j] = cos(2*M_PI*q.v[j]);
			}else{
				if(val > 0.25 && val < 0.75){//Valleys weights
					k.v[j] = k.v[j] + cos(2*M_PI*q.v[j]) / 2;
				}
			}
		}
	}
	//Apply envelope
	double norm = 0;
	for(i = 0; i < k.x; i++){
		k.v[i] = k.v[i]*sqrt(1./f.v[i]);
		if(k.v[i] > 0)
			norm += k.v[i]*k.v[i];
	}
	norm = sqrt(norm);
	//K+- normalize kernel
	for(i = 0; i < k.x; i++){
		k.v[i] = k.v[i]/norm;
	}
	//Compute pitch strength
	for(i = 0; i < NL.x; i++){
		for(j = 0; j < k.x; j++){
		S.v[i] += NL.m[i][j]*k.v[j];
		}
	}
	freeiv(primos);
	freev(k);
	freev(prim);
	freev(q);
	return S;
}


matrix pitchStrengthAllCandidates(vector f, matrix L, vector pc, vector j){
	//Create pitch salience matrix
	printf("\n		Create pitch salience matrix...\n");
	matrix S = zerom(j.x, L.x);
	vector k = zerov(j.x);
	int q, a, find, i, c;
	int ant = 0;
	vector pc2 = zerov(j.x);
	for(q = 0; q < j.x; q++){
		pc2.v[q] = pc.v[(int)j.v[q]];
	}
	for(q = 0; q < k.x; q++){//Find
		for(a = ant; f.x; a++){
			if(f.v[a] > pc2.v[q]/4.){
				find = a;
				break;
			}	
		}
		k.v[q] = find;
		ant = k.v[q];
	}
	//Loudness normalization factor
	printf("\n		Loudness normalization factor...\n");
	matrix N = zerom(L.x, L.y);
	int fil,col;
	double suma = 0;
	for(fil = 0; fil<N.x; fil++, suma = 0){
		for(col = N.y-1; col >=0 ; col--){
			suma += L.m[fil][col]*L.m[fil][col];
			N.m[fil][col] = sqrt(suma);
		}
	} 
	double val;
	vector n;
	matrix NL;
	vector f2;
	printf("\n		Normalize Loudness...\n");
	printf("\n		Compute each candidate's pitch strength...\n");
	clock_gettime(CLOCK_MONOTONIC, &start);
	for(q = 0; q < pc2.x; q++){
		//printf("hilo %d\n", omp_get_thread_num());
		NL = zerom(L.x, L.y-(int)k.v[q]);
		f2 = zerov(L.y-(int)k.v[q]);
		n = zerov(N.x);
		//Normalize Loudness
		for(i = 0; i < n.x; i++){
			val = N.m[i][(int)k.v[q]];
			if(val != 0){
				n.v[i] = val;
			}else{
				n.v[i] = INFINITY;
			}
		}
		for(i = 0; i < NL.x; i++){
			for(a = 0,c = (int)k.v[q]; a < NL.y; a++, c++){
				NL.m[i][a] = L.m[i][c] / n.v[i];
				f2.v[a] = f.v[c];
			}
		}
		//Compute pitch strength
		
		vector Si = pitchStrengthOneCandidate( f2, NL, pc2.v[q] );
		for(i = 0; i < Si.x; i++){
			S.m[q][i] = Si.v[i]; 
		}
		freev(Si);
		freev(f2);
		freev(n);
		freem(NL);
	}
	clock_gettime(CLOCK_MONOTONIC, &end);//end of execution time
	//printf("tiempo: %f\n",timespecDiff(&end, &start));
	time3 += timespecDiff(&end, &start);
	freev(k);
	freev(pc2);
	freem(N);
	return S;
}

// Peforms polynomial tuning on the strength matrix to determine the pitch
vector pitch(matrix S, vector pc, double st) {
    int i;
    int j;
    int maxi; 
    int search = (int) round((log2(pc.v[2]) - log2(pc.v[0])) / POLYV + 1.);
    double nftc; 
    double maxv;
    double log2pc;
    double tc2 = 1. / pc.v[1];
    vector coefs;
    vector s = makev(3);
    vector ntc = makev(3);
    ntc.v[0] = ((1. / pc.v[0]) / tc2 - 1.) * 2. * M_PI; 
    ntc.v[1] = (tc2 / tc2 - 1.) * 2. * M_PI; 
    ntc.v[2] = ((1. / pc.v[2]) / tc2 - 1.) * 2. * M_PI;
    vector p = makev(S.y);
    for (j = 0; j < S.y; j++) {
        maxv = SHRT_MIN;  
        for (i = 0; i < S.x; i++) { 
            if (S.m[i][j] > maxv) {
                maxv = S.m[i][j];
                maxi = i;
            }
        }
        if (maxv > st) { // Make sure it's big enough
            if (maxi == 0 || maxi == S.x - 1) { // First or last? 
                p.v[j] = pc.v[0];
            }
            else { // Generic case
                tc2 = 1. / pc.v[maxi];
                log2pc = log2(pc.v[maxi - 1]); 
                s.v[0] = S.m[maxi - 1][j];
                s.v[1] = S.m[maxi][j];
                s.v[2] = S.m[maxi + 1][j]; 
                coefs = polyfit(ntc, s, 2); 
                maxv = SHRT_MIN; 
                for (i = 0; i < search; i++) { // Check the nftc space
                    nftc = polyval(coefs, ((1. / pow(2, i * POLYV + log2pc)) / 
                                                          tc2 - 1) * 2 * M_PI);
                    if (nftc > maxv) {
                        maxv = nftc;
                        maxi = i;
                    }
                } // Now we've got the pitch numbers we need
                freev(coefs);
                p.v[j] = pow(2, log2pc + (maxi * POLYV));
            }
        }
        else { // If not voiced during that interval, then...
            p.v[j] = NAN;
        } 
    }
    freev(ntc);
    freev(s);
    return(p);
}

// Primary utility function for each pitch extraction
vector swipe(char wav[], double min, double max, double st, double dt/*, char filename[], FILE* time_test*/) {
    
    
	clock_gettime(CLOCK_MONOTONIC, &startTotal);//Beginning of execution time
	int i, ind, tamX, rango, maxvent, minvent; 
	double a,td,dlog2p, dlog2p_max, nyquist, nyquist2, fs;
	double dERBs = 0.1;
	double total_time = 0;
	time1 = 0, time2 = 0, time3 = 0;
	double woverlap = 0.5;
	FILE* wavf; 
	SF_INFO info;
	SNDFILE* source;
	if (strcmp(wav, "<STDIN>") == 0) { // i.e., is coming from STDIN
	   wavf = stdin;
	}
	else { // is specified
	   wavf = fopen(wav, "r");
	}
	printf("\nReading WAV file: %s...\n", wav);
	source = sf_open_fd(fileno(wavf), SFM_READ, &info, TRUE);
	// Perform checks on the wav header
	if (info.sections < 1) {
	   fprintf(stderr, "File or stream %s not read as audio ... \n", wav);
	   return(makev(0)); // This will be detected as an error
	}
	nyquist = info.samplerate / 2.; 
	nyquist2 = info.samplerate; // Used so g.d. often here...
	fs = info.samplerate;
	if (max > nyquist) { 
	   max = nyquist;
	   fprintf(stderr, "Max pitch > Nyquist ... max set to %.2f Hz.\n", max);
	}
	if (dt > nyquist2) {
	   dt = nyquist2;
	   fprintf(stderr, "Timestep > SR ... timestep set to %f.\n", nyquist2);
	}
	vector x = makev(info.frames); // This reads the signal in
	sf_read_double(source, x.v, x.x);
	sf_close(source); // Takes FILE* wavf with it, too
	fclose(wavf);
	vector t = zerov(ceil((double)((double)info.frames/(double)fs)*(1./dt)+1));
	a = 0;
	for(i = 0; i < t.x; i++, a = a+dt){	
		t.v[i] = a;
	}
	printf("\nFiltering signal into segments of the cochlea...\n");
	
	//clock_gettime(CLOCK_MONOTONIC, &start);
	int validSample = 0;
	for(i = 0; i < x.x && validSample == 0; ++i){
		if(x.v[i] != 0)validSample = 1;
	}
	if(validSample == 0){
		printf("Error: not a valid sample\n ");
		 //exit(EXIT_FAILURE);
	}
	
	out(x, "x_read_c.csv");
	
	
	if(validSample == 1){
	
	matrix X = audsys(x, (double)fs, &time1);//Este metodo retorna la señal 
									//luego de ser filtrada por los segmentos de la coclea
	
	
	//clock_gettime(CLOCK_MONOTONIC, &end);//end of execution time
	//time1 = timespecDiff(&end, &start);
	//fprintf(time_test, "%s\t%.15f\t", filename, time);
	
	//Print matrix for testing
	outm(X.m, X.x, X.y, "testX.txt");								
	
	printf("\nSignal divided in segments...\n");
	printf("\nCreating Aud-SwipeP variables...\n");
	tamX = X.x;
	vector f = zerov(tamX);
	//Define pitch candidates 
	dlog2p = 1./48.;
	rango = ceil((log2(max)-log2(min))*(1./dlog2p)); 
	vector log2pc = zerov(rango);
	vector pc = zerov(rango);
	dlog2p_max = log2(max);
	i = 0;
	for(a = log2(min); a <= dlog2p_max; a = a+dlog2p, i++){
		log2pc.v[i] = a;
		pc.v[i] = pow(2, log2pc.v[i]);
	}
	matrix S = zerom(pc.x, t.x);//Pitch Strength matrix
	
	//Determine power-of-two window sizes
	printf("\n	Determine power-of-two window sizes...\n");
	maxvent = round(log2(8*fs/min));
	minvent = round(log2(8*fs/max));	
	vector logWs = zerov(2);
	logWs.v[0] = maxvent;
	logWs.v[1] = minvent;
	//power-of-two window sizes
	vector ws = zerov(maxvent - minvent + 1);
	for(a = maxvent,i = 0; a >= minvent; a--, i++) ws.v[i] = pow(2, a);
	//Optimal pitches for power-of-two window sizes
	printf("\n	Optimal pitches for power-of-two window sizes...\n");
	vector p0 = zerov(ws.x);
	for(i = 0; i < p0.x; i++) p0.v[i] = 8*fs/ws.v[i];
	//Determine window sizes used by each pitch candidate
	printf("\n	Determine window sizes used by each pitch candidate...\n");
	vector d = zerov(log2pc.x);
	for(i = 0; i < d.x; i++) d.v[i] = 1 + log2pc.v[i] - log2(8.*fs/ws.v[0]);
	//Create ERB-scale uniformly-spaced frequencies (in Hertz)
	printf("\n	Create ERB-scale uniformly-spaced frequencies (in Hertz)...\n");
	vector fERBs = zerov(ceil((hz2erb(fs/2.))*(1./dERBs)));
	for(i = 0, a = 0; i < fERBs.x; i++, a = a+dERBs) fERBs.v[i] = erb2hz(a);
	
	for(i = 0; i < ws.x; i++){
	printf("\nCalculating with window size = %f...\n", ws.v[i]);
		//Determine pitch candidates that use this window size
		//printf("\n	Determine pitch candidates that use this window size...\n");
		vector j = zerov(pc.x);
		vector k = zerov(pc.x);
		int cont;
		if(ws.x == 1){
			resizev(&k, 0);
			for(ind = 0; ind < j.x; ind++){
				j.v[ind] = ind+1;
			}
			
		}else{ 
			if(i == ws.x-1){
				cont = 0;
				for(ind = 0; ind < d.x; ind++){
					if(d.v[ind] - (double)(i+1.) > -1){
						j.v[cont] = ind;
						cont++;
					}
				}
				resizev(&j, cont);
				cont = 0;
				for(ind = 0; ind < j.x; ind++){
					if(d.v[(int)j.v[ind]] - (double)(i+1.) < 0){
						k.v[cont] = ind;
						cont++;
					}
				}
				resizev(&k, cont);
			}else{ 
				if(i == 0){
					cont = 0;
					for(ind = 0; ind < d.x; ind++){
						if(d.v[ind] - (double)(i+1) < 1.){
							j.v[cont] = ind;
							cont++;
						}
					}
					resizev(&j, cont);
					cont = 0;
					for(ind = 0; ind < j.x; ind++){
						if(d.v[(int)j.v[ind]] - (double)(i+1.) > 0.){
							k.v[cont] = ind;
							cont++;
						}
					}
					resizev(&k, cont);
				}else{
					cont = 0;
					for(ind = 0; ind < d.x; ind++){
						if(abs(d.v[ind] - (double)(i+1)) < 1.){
							j.v[cont] = ind;
							cont++;
						}
					}
					resizev(&j, cont);
					cont = 0;
					for(ind = 0; ind < j.x; ind++){
						k.v[ind] = ind;
					}
					resizev(&k,j.x);
				}
			}	
		}
		// Zero pad signal
		printf("\n	Zero pad signal...\n");
		double dn = maxim(1, round(8*(1-woverlap)  * fs/p0.v[i] ));//Hop size
		matrix Xz = zerom(X.x, X.y + (ws.v[i]/2)+ (dn+ws.v[i]/2) );
		int fil, col, col2;
		for(fil = 0; fil < X.x; fil++){
			col2 = 0;
			for(col = (ws.v[i]/2); col < X.y + (ws.v[i]/2); col++,col2++){
				Xz.m[fil][col] = X.m[fil][col2];
			}
		}
		
		
		//Compute specific loudness
		printf("\n	Compute specific loudness...\n");
		vector w = zerov(ws.v[i]);
		Hanning(w.v,w.x);
		vector n = zerov(ceil((Xz.y - ws.v[i] + 1.)*(1./dn)));//Centers of the windows
		for(ind = 0, a = 0; ind < n.x; a = a + dn, ind++){
			n.v[ind] = a;
		}
		//Specific-loudness matrix
		printf("\n	Creating specific-loudness matrix...\n");
		matrix L = zerom(n.x, fERBs.x);
		double df = fs/ws.v[i];
		vector fi = zerov(ws.v[i]);
		for(ind = 0; ind < fi.x; ind++){
			fi.v[ind] = (double)ind*df;
		}
		vector l = zerov(tamX);
		for(ind = 0; ind < f.x; ind++) f.v[ind] = ind + 1.5;
		v_erbs2hz(f/*.v, f.x*/);
		int m;
		int ant = 0;
		//Compute characteristic frequencies indices
		printf("\n	Compute characteristic frequencies indices...\n");
		for(m = 0; m < l.x;m++){
			//find
			int find =0;
			for(ind = fi.x-1; ind >=ant; ind--){
				if(fi.v[ind] < 1.25*f.v[m]){
					find = ind;
					break;
				}
			}
			l.v[m] = find; 
			ant = l.v[m];
		}
		int q;
		int q2;
		matrix W = zerom(f.x, fi.x);
		printf("\n	Compute raised-cosine weights...\n");
		for(q = 0; q < f.x; q++){
			//Compute raised-cosine weights
			for(q2 = 0; q2 <= l.v[q]; q2++){
				W.m[q][q2] = (1. - cos(M_PI * hz2erb(fi.v[q2])/(hz2erb(f.v[q]))))/2;
			}
		}
		int p;
		int segm;
		int copia;
		double* window;
		fftw_complex *fo;
		
		printf("\n	Calculating FFT for each segment of the cochlea at every window...\n");
		printf("\n		Creating and Executing FFTW Plan (PARALLEL)...\n");
		
		
		
		clock_gettime(CLOCK_MONOTONIC, &start);
		
		
		
		fftw_plan plan = fftw_plan_dft_r2c_1d(ws.v[i], window, fo, FFTW_ESTIMATE); 
		for(p = 0; p < n.x; p++){
		
			vector sl = zerov(ws.v[i]);//Specific loudness
			//Compute specific loudness
			
			for(segm = 0; segm < tamX; segm++){/////PARALELIZAR AQUI
				window = fftw_malloc(sizeof(double) * ws.v[i]); 
				fo = fftw_malloc(sizeof(fftw_complex)*ws.v[i]);
				
				
				ind = 0;
				for(col = n.v[p]; col < n.v[p] + ws.v[i]; col++, ind++){
					window[ind] = (w.v[ind]*Xz.m[segm][col]);
				}
				//Se calcula la fft a la ventana calculada anteriormente
				////Magnitudes
				fftw_execute_dft_r2c(plan, window, fo);
				copia = ws.v[i]-1;
				for(ind = 1; ind < ws.v[i]/2; ind++, copia--){//fft solo calcula los primeros n/2 + 1, por nyquist
					fo[copia][0] = fo[ind][0];
					fo[copia][1] = fo[ind][1];
				}
				
				for(ind = 0; ind < ws.v[i]; ind++){//PARALELIZAR AQUI
					
					
					sl.v[ind] += sqrt(
									sqrt((fo[ind][0]*fo[ind][0])	
										+(fo[ind][1]*fo[ind][1]))
									* W.m[segm][ind]);//Specific loudness
					
					fo[ind][0] = 0;
					fo[ind][1] = 0;
				
				}
				fftw_free(window); 
				fftw_free(fo); 
			}
			
			// in ERB scale
			interp1(fi, sl, fERBs, L.m[p]);
			freev(sl);
		}
		
		
		
		fftw_destroy_plan(plan); 
		
		
		clock_gettime(CLOCK_MONOTONIC, &end);//end of execution time
		time2 += timespecDiff(&end, &start);
		//fprintf(time_test, "%.15f\t", time);
		
		
		
		//Se eliminan los negativos de la matriz L
		Max(L);
		//Compute pitch strength
		printf("\n	Compute pitch strength...\n");
		
		
		//clock_gettime(CLOCK_MONOTONIC, &start);
		
		
		matrix Si = pitchStrengthAllCandidates(fERBs, L, pc, j);
		

		
		//clock_gettime(CLOCK_MONOTONIC, &end);//end of execution time
		//time3 += timespecDiff(&end, &start);
		//fprintf(time_test, "%.15f\t", time);
		
		//Interpolate pitch strength at desired times
		printf("\n	Interpolate pitch strength at desired times...\n");
		matrix Si2 = zerom(Si.x, t.x);
		vector g = zerov(Si.y);
		if(Si.y > 1){
			for(ind = 0; ind < n.x; ind++){
				n.v[ind] = (n.v[ind]-1)/fs;
			}	
			for(ind = 0; ind < Si2.x; ind++){
				for(p = 0; p < g.x; p++){
					g.v[p] = Si.m[ind][p];
				}
				interp1(n, g,t,Si2.m[ind]);
			}
		}
		//Compute contribution of this window size to pitch strength
		printf("\n	Compute contribution of this window size to pitch strength...\n");
		vector mu = onesv(j.x);
		for(ind = 0; ind < k.x; ind++){
			mu.v[(int)k.v[ind]] = 1 - fabs(d.v[(int)j.v[(int)k.v[ind]]]-(i+1));
		}
		for(ind = 0; ind < j.x; ind++){
			for(p = 0; p < S.y; p++){
			 S.m[(int)j.v[ind]][p]  += mu.v[ind]*Si2.m[ind][p];
			}
		}
		freev(w);
		freev(mu);
		freev(g);
		freem(Si2);
		freem(Si);
		freem(W);
		freem(L);
		freev(l);
		freev(j);
		freev(k);
		freev(n);
		freev(fi);
		freem(Xz);
		
	}
	
	
	printf("\nCompute pitch...\n");
			
	vector p = pitch(S,pc, st);//por mientras
	printf("\nFree everything...\n");
	printf("\nAud-SwipeP finalized!!!\n");
	freev(f);
	freev(p0);
	freev(fERBs);
	freev(d);
	freev(x);
	freem(X);
	freev(log2pc);
	freev(pc);
	freev(logWs);
	freev(ws);
	freev(t);
	freem(S);
	
	
	clock_gettime(CLOCK_MONOTONIC, &endTotal);//end of execution time
	total_time = timespecDiff(&endTotal, &startTotal);
	printf( "%.10f\t%.10f\t%.10f\t%.10f\n", time1, time2, time3, total_time);
	//printf("\nExecution Time = %.10f\n", ((double)end - (double)start)/(double)(CLOCKS_PER_SEC));
	return(p);
	}
	else{
		vector p = zerov(3);
		freev(x);
		freev(t);
		return(p);
	}
}

// Function for printing the pitch vector returned by swipe()
void printp(vector p, char out[], double dt, int mel, int vlo) {
    int i;
    double t = 0.; 
    FILE* sink; // Handle for printing to file/STDOUT
    if (strcmp(out, "<STDOUT>") == 0) {
        sink = stdout;
    }
    else {
        sink = fopen(out, "w");
        if (sink == NULL) {
            fprintf(stderr, "File or stream %s not writable, aborting.\n", out);
            exit(EXIT_FAILURE);
        }
    }
    if (mel) {
        if (vlo) {
            for (i = 0; i < p.x; i++) {
                fprintf(sink, "%4.7f %5.4f\n", t, hz2mel(p.v[i])); 
                t += dt;
            }
        }
        else { // Default case
            for (i = 0; i < p.x; i++) {
                if (!isnan(p.v[i])) {
                    fprintf(sink, "%4.7f %5.4f\n", t, hz2mel(p.v[i])); 
                }
                t += dt;
            }
        }
    }
    else {
        if (vlo) {
            for (i = 0; i < p.x; i++) {
                fprintf(sink, "%4.7f %5.4f\n", t, p.v[i]); 
                t += dt;
            }
        }
        else { 
            for (i = 0; i < p.x; i++) {
                if (!isnan(p.v[i])) {
                    fprintf(sink, "%4.7f %5.4f\n", t, p.v[i]); 
                }
                t += dt;
            }
        }
    } 
    fclose(sink); 
}


// Main method, interfacing with user arguments
int main(int argc, char* argv[]) {

	//start = clock();//Beginning of execution time
	//Se inicia MPI
    char output[] = "OUTPUT:\npitch_0\ttime_0\npitch_1\ttime_1\n...\t...\
    \npitch_N\ttime_N\n\n"; 
    char header[] = "SWIPE' pitch tracker, implemented in C by Kyle Gorman \
<kgorman@ling.upenn.edu>.\nBased on: Camacho, Arturo (2007). A sawtooth \
waveform inspired pitch estimator\nfor speech and music. Doctoral \
dissertation, University of Florida.\n\n\
\tmore information: <http://ling.upenn.edu/~kgorman/c/swipe/>\n\n";
    char synops[] = "SYNPOSIS:\n\n\
swipe [-i INPUT] [-b LIST] [-o OUTPUT] [-r MIN:MAX] [-s ST] [-t DT] [-mnhv]\n\n\
FLAG:\t\tDESCRIPTION:\t\t\t\t\tDEFAULT:\n\n\
-i FILE\t\tinput file\t\t\t\t\tSTDIN\n\
-o FILE\t\toutput file\t\t\t\t\tSTDOUT\n\
-b LIST\t\tbatch mode: [LIST is a file containing\n\
\t\tone \"INPUT OUTPUT\" pair per line]\n\n\
-r MIN:MAX\tpitch range in Hertz\t\t\t\t100:600\n\
-s THRSHLD\tstrength threshold  [0 <= x <= 1]\t\t0.300\n\
-t SECONDS\ttimestep in seconds [must be < SF / 2]\t\t0.001\n\n\
-m\t\tOutput Mel pitch\t\t\t\tno\n\
-n\t\tDon't output voiceless frames\t\t\tno\n\
-h\t\tDisplay this message, then quit\n\
-v\t\tDisplay version number, then quit\n\n";
    double st = ST; // All set by #defines
    double dt = DT;
    int vlo = TRUE;
    int mel = FALSE; 
    double min = MIN;
    double max = MAX; 
    int ch;
    FILE* batch = NULL; // not going to be read that way,
    char wav[FILENAME_MAX] = "<STDIN>";
    char out[FILENAME_MAX] = "<STDOUT>";
    printf("\nInitializing Aud-SwipeP...\n");
    while ((ch = getopt(argc, argv, "i:o:r:s:t:b:mnhv")) != -1) {
        switch(ch) {
            case 'b':
                batch = fopen(optarg, "rt"); 
                break;
            case 'i':
                strcpy(wav, optarg); 
                break; 
            case 'o':
                strcpy(out, optarg); 
                break;
            case 'r':
                min = atof(strtok(optarg, ":"));
                max = atof(strtok(NULL, ":")); 
                break;
            case 's':
                st = atof(optarg);
                break;
            case 't':
                dt = atof(optarg);
                break;
            case 'm':
                mel = TRUE; 
                break;
            case 'n':
                vlo = FALSE; 
                break;
            case 'h':
                fprintf(stderr, "%s", header); 
                fprintf(stderr, "%s", synops); 
                fprintf(stderr, "%s", output);
                exit(EXIT_SUCCESS);
            case 'v':
                fprintf(stderr, "This is SWIPE', v. %1.1f.\n", VNUM); 
                exit(EXIT_SUCCESS);
            case '?': 
            default:  // Would like to do clever things here, but no ideas yet
                fprintf(stderr, "%s", header);
                fprintf(stderr, "%s", synops);
                exit(EXIT_FAILURE);
            argc -= optind; 
            argv += optind;
        }
    }
    if (min < 1.) { // Santiny-check the args
        fprintf(stderr, "Min pitch < 1 Hz, aborting.\n"); 
        exit(EXIT_FAILURE);
    }
    if (max - min < 1.) {
        fprintf(stderr, "Max pitch <= min pitch, aborting.\n"); 
        exit(EXIT_FAILURE);
    } 
    if (st < 0. || st > 1.) { 
        fprintf(stderr, "Strength must be 0 <= x <= 1, set to %.3f.\n", ST); 
        st = ST;
    }
    if (dt < .001) {
        fprintf(stderr, "Timestep must be >= 0.001 (1 ms), set to %.3f.\n", DT);
        dt = DT;
    }
    if (batch != NULL) { // Iterate through batch pairs
        while (fscanf(batch, "%s %s", wav, out) != EOF) {
            fprintf(stderr, "%s -> %s ... ", wav, out);
            vector p = swipe(wav, min, max, st, dt);
            if (p.x == NOK) {
                fprintf(stderr, "File or stream %s failed.\n", wav);
                fclose(batch); 
                exit(EXIT_FAILURE);
            }
            else {
                printp(p, out, dt, mel, vlo);
                printf("done.\n");
            }
            freev(p);
        }
        fclose(batch);
        exit(EXIT_SUCCESS);
    }
    else {
        vector p = swipe(wav, min, max, st, dt);
        if (p.x == NOK) {
            fprintf(stderr, "File or stream %s failed.\n", wav);
            exit(EXIT_FAILURE);
        }
        else {
            printp(p, out, dt, mel, vlo); 
        } 
        freev(p);
    }
    
    
    /*double f[] = {	0, 0.6, 0.6, 1 };
    double m[] = {	1, 1, 0, 0 };
    
    
    vector b = makev(31);
    fir2(30, f, m, 4, b.v);
    out(b.v, b.x, "nuevoB.txt");*/
    
    exit(EXIT_SUCCESS);
}

