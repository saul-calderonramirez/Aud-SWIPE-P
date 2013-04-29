/*
 * AudSWIPEP.c
 *
 *  Created on: Apr 24, 2013
 *      Author: saul
 */

#include "AudSWIPEP.h"
#include <float.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <fftw3.h>   // http://www.fftw.org/
#include <sndfile.h> // http://www.mega-nerd.com/libsndfile/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <fftw3.h>   // http://www.fftw.org/
#include <sndfile.h> // http://www.mega-nerd.com/libsndfile/
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>

#define MI 100
#define MA 500
#define ST .0001
#define DT .001

#define LEADER						0

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


/*
*Returns the pitch strength of one pitch candidate, applying the kernel to the matrix of signals
*@param f, the erb scale
*@param NL, the normalized loudness matrix of signals
*@param pc, the pitch candidate
*@return vector, computed pitch strenghts
*/
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

/*
*Returns the pitch strength of several pitch candidate, applying the kernel to the matrix of signals
*@param f, the erb scale
*@param L, the loudness matrix of signals
*@param pc, the pitch candidate
*@param j,
*@return matrix, computed pitch strenghts per pitch candidate
*/
matrix pitchStrengthAllCandidates(vector f, matrix L, vector pc, vector j){
	//Create pitch salience matrix
	if(DEBUG == 1)printf("\n		Create pitch salience matrix...\n");
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
	if(DEBUG == 1)printf("\n		Loudness normalization factor...\n");
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
	if(DEBUG == 1)printf("\n		Normalize Loudness...\n");
	if(DEBUG == 1)printf("\n		Compute each candidate's pitch strength...\n");
	/*
	 * Thread parallelization 3, each thread calculates the score of a pitch candidate
	 * */
	#pragma omp parallel for private (NL, f2, n, val, i, a, c)
	for(q = 0; q < pc2.x; q++){

		NL = zerom(L.x, L.y-(int)k.v[q]);
		f2 = zerov(L.y-(int)k.v[q]);
		n = zerov(N.x);
		//Normalize Loudness
		for(i = 0; i < n.x; i++){
			val = N.m[i][(int)k.v[q]];
			if(val != 0){
				n.v[i] = val;
			}else{
				n.v[i] = 88888;
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
	freev(k);
	freev(pc2);
	freem(N);
	return S;
}

/*
*Peforms polynomial tuning on the strength matrix to determine the pitch
*@param f, the erb scale
*@param NL, the normalized loudness matrix of signals
*@param pc, the pitch candidate
*@return vector, computed pitch strenghts
*/
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

/*
 * Creates the matrix of channels X
 * Must be done by the parent only, after this, it is sent to the slaves
 * @param wav, wav file name
 * @param min, minimum freq. to look for
 * @param max, max. freq. to look for
 * @param dt, pitch detection time interval
 * */
matrix calculateWindowSizesAndCochleasResponses(char wav[], double min, double max, double dt,  char* argv[], double paralelism, matrix* X, vector* fsInfo ){
		clock_gettime(CLOCK_MONOTONIC, &startTotal);
		int maxvent, minvent;
		double  nyquist, nyquist2, fs;
		*fsInfo = zerov(2);
		matrix ranges;
		FILE* wavf;
		SF_INFO info;
		SNDFILE* source;
		//opening file time
		// Clock 1: Wav read
		startLocalClock(&clocks);
		if (strcmp(wav, "<STDIN>") == 0) { // i.e., is coming from STDIN
		   wavf = stdin;
		}
		else { // is specified
		   wavf = fopen(wav, "r");
		}
		if(DEBUG == 1){printf("\nReading WAV file: %s\n", wav);}
		source = sf_open_fd(fileno(wavf), SFM_READ, &info, TRUE);

		// Perform checks on the wav header
		if (info.sections < 1) {
		   fprintf(stderr, "File or stream %s not read as audio ... \n", wav);
		  // return(makev(0)); // This will be detected as an error
		}
		nyquist = info.samplerate / 2.;
		nyquist2 = info.samplerate; // Used so g.d. often here...
		fs = info.samplerate;
		int frames = (int)info.frames;
		soundLength = (double)((double)frames/ (double)info.samplerate);

		if(DEBUG == 1){printf("\nInfo READ %s\n", wav);}
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
		endLocalClock(&clocks, 1, "FileReading");
		if(DEBUG == 1){printf("\nFiltering signal into segments of the cochlea...\n");}
		//AudSWIPE pre-processing before
		//Calculates the response of the different cochleas segments
		//For the total time we will ignore the read file time.
		*X = audsys(x, (double)fs, &clocks);
		if(DEBUG == 1){printf("\nSignal divided in segments...\n");}
		if(DEBUG == 1){printf("\nCreating Aud-SwipeP variables...\n");}
		fsInfo->v[0] = fs;
		fsInfo->v[1] = (double)info.frames;
		//In here the work amount for every process is defined
		maxvent = round(log2(8*fs/min));
		minvent = round(log2(8*fs/max));
		int wSizes = maxvent - minvent + 1;
		ranges = generateWindowSizesRanges( paralelism,  wSizes);
		if(DEBUG == 1){printf("\n	Cochleas calculation finished...\n");}
		freev(x);
		return ranges;
}



/*
 * Sends processed  data back to the leader or master process
 * @param mySi2
 * @param myJ
 * @param myMu
 * @param ds
 * @param nX
 * */
void sendDataToLeader(matrix mySi2, vector myJ, vector myMu, double ds, int nX){
	//Sends the data to the leader
	int numSlaves = getNumberProcesses() - 1;
	if(numSlaves > 0 && !isParent()){
		sendMatrix(LEADER, mySi2);
		sendVector(LEADER, myMu);
		sendVector(LEADER, myJ);
		sendDouble(LEADER, ds);
		sendDouble(LEADER, (double)nX);
	}
}

/*
 * Receives all the data from the children or slaves processes and builds the scores Matrix S
 * @param S, score matrix
 * @param t, time vector of the sound
 * @param fs, sampling frequency
 * @param ranges, ranges of the matrix S from each slave, to rebuild the final matrix S, from
 * */
matrix receiveDataBuildS(matrix S, vector t, double fs, matrix ranges){
	int ind;
	int p;
	int counter;
	int slave = 0;
	int numSlaves = getNumberProcesses() - 1;
	int a;
	double dn;
	double sizeN;
	//if it is a slave, sends the data
	if(isParent() && numSlaves > 0){//if its the master, receives the data
		vector mu = zerov(2);
		matrix Si = zerom(2,2);
		vector j = zerov(2);
		//receives and uses the part calculated by the slaves
		for(slave = 1; slave < (numSlaves + 1); ++slave){
			//i = (int)ranges.m[ myid ][ 0 ]; i < (int)ranges.m[ myid ][ 1 ]; i++
			for(counter = (int)ranges.m[ slave ][ 0 ]; counter < (int)ranges.m[ slave ][ 1 ]; ++counter ){
				receiveMatrix(slave, &Si);
				receiveVector(slave, &mu);
				receiveVector(slave, &j);
				receiveDouble(slave, &dn);
				receiveDouble(slave, &sizeN);
				vector n = zerov(((int) sizeN));
				for(ind = 0, a = 0; ind < n.x; a = a + dn, ind++){
						n.v[ind] = a;
				}
				for(ind = 0; ind < n.x; ind++){
					n.v[ind] = (n.v[ind]-1)/fs;
				}

				matrix Si2 = zerom(Si.x, t.x);
				vector g = zerov(Si.y);
				if(Si.y > 1){
					for(ind = 0; ind < Si2.x; ind++){
						for(p = 0; p < g.x; p++){
							g.v[p] = Si.m[ind][p];
						}
						interp1(n, g,t,Si2.m[ind]);
					}
				}
				for(ind = 0; ind < j.x; ind++){
					for(p = 0; p < S.y; p++){
						 S.m[(int)j.v[ind]][p]  += mu.v[ind]*Si2.m[ind][p];
					}
				}
			}
		}
	}
	return S;
}


/*
*Function for printing the pitch vector returned by AudSwipe
*@param vector p
*@param out
*@param dt
*@param mel
*@param vlo
*/
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


matrix getWsLoudnessMat(int i, vector ws, vector n, matrix Xz, vector w, vector fERBs, vector fi, matrix W, int tamX){
	//End clock 8
	int p;
	int ind;
	int segm, col;
	int copia;
	double* window;

	matrix L = zerom(n.x, fERBs.x);
	//OJO! SE CAMBIO DE UN POINTER A DOBLE POINTER!
	fftw_complex** fo;
	if(DEBUG==1)printf("\n	Calculating FFT for each segment of the cochlea at every window...\n");
	if(DEBUG==1)printf("\n		Creating and Executing FFTW Plan (PARALLEL)...\n");

	 /* Thread parallelization 2, each thread calculates the Enhanced spectrum or ESRAS for a window
	 * */
	fftw_plan plan = fftw_plan_dft_r2c_1d(ws.v[i], window, fo, FFTW_ESTIMATE);
	#pragma omp parallel for private(segm,ind,col, window, copia, fo)
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
			//Calculates the fft of the window calculated previously
			////Magnitudes
			fftw_execute_dft_r2c(plan, window, fo);
			copia = ws.v[i]-1;
			for(ind = 1; ind < ws.v[i]/2; ind++, copia--){//fft Calculates only the first n/2 + 1, because of nyquist
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
	endLocalClock(&cycleClocks, 0, "FourierTransforms");
	//End Clock 9
	//Eliminates the negatives from L
	Max(L.m, L.x, L.y, 0);
	return L;
}

void getWsScoreMat(int i, vector ws, vector pc, matrix X, vector fERBs, vector d, matrix S, double fs, vector p0, vector t){
	int ind, p;
	double a = 0;
	double woverlap = 0.5;
	int tamX = X.x;
	vector f = zerov(tamX);
	if(DEBUG==1)printf("\nCalculating with window size = %f...\n", ws.v[i]);
	//Determine pitch candidates that use this window size
	if(DEBUG==1)printf("\n	Determine pitch candidates that use this window size...\n");
	startLocalClock(&cycleClocks);
	//starts the clock
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
	endLocalClock(&cycleClocks, 0, "DefPitchCandidates");
	// Zero pad signal
	if(DEBUG==1)printf("\n	Zero pad signal...\n");
	double dn = maxim(1, round(8*(1-woverlap)  * fs/p0.v[i] ));//Hop size
	matrix Xz = zerom(X.x, X.y + (ws.v[i]/2)+ (dn+ws.v[i]/2) );
	int fil, col, col2;
	for(fil = 0; fil < X.x; fil++){
		col2 = 0;
		for(col = (ws.v[i]/2); col < X.y + (ws.v[i]/2); col++,col2++){
			Xz.m[fil][col] = X.m[fil][col2];
		}
	}
	endLocalClock(&cycleClocks, 0, "ZeroPadSignal");
	//End clock 5
	//Compute specific loudness
	if(DEBUG==1)printf("\n	Compute specific loudness...\n");
	vector w = zerov(ws.v[i]);
	Hanning(w.v,w.x);
	vector n = zerov(ceil((Xz.y - ws.v[i] + 1.)*(1./dn)));//Centers of the windows
	for(ind = 0, a = 0; ind < n.x; a = a + dn, ind++){
		n.v[ind] = a;
	}
	//Specific-loudness matrix
	if(DEBUG==1)printf("\n	Creating specific-loudness matrix...\n");
	matrix L = zerom(n.x, fERBs.x);
	double df = fs/ws.v[i];
	vector fi = zerov(ws.v[i]);
	for(ind = 0; ind < fi.x; ind++){
		fi.v[ind] = (double)ind*df;
	}
	vector l = zerov(tamX);
	for(ind = 0; ind < f.x; ind++) f.v[ind] = ind + 1.5;
	v_erbs2hz( f );
	endLocalClock(&cycleClocks, 0, "LoudnessMatrix");
	//End clock 6
	int m;
	int ant = 0;
	//Compute characteristic frequencies indices
	if(DEBUG==1)printf("\n	Compute characteristic frequencies indices...\n");
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

	// End clock 7
	endLocalClock(&cycleClocks, 0, "ComputeCharacteristicFrequencies");
	int q;
	int q2;
	matrix W = zerom(f.x, fi.x);
	if(DEBUG==1)printf("\n	Compute raised-cosine weights...\n");
	for(q = 0; q < f.x; q++){
		//Compute raised-cosine weights
		for(q2 = 0; q2 <= l.v[q]; q2++){
			W.m[q][q2] = (1. - cos(M_PI * hz2erb(fi.v[q2])/(hz2erb(f.v[q]))))/2;
		}
	}
	endLocalClock(&cycleClocks, 0, "ComputeRaised-cosineWeights");

	//Calculates the Loudness Matrix

	L = getWsLoudnessMat( i,  ws,  n,  Xz, w, fERBs, fi, W, tamX);


	//Compute pitch strength
	if(DEBUG==1)printf("\n	Compute pitch strength...\n");
	matrix Si = pitchStrengthAllCandidates(fERBs, L, pc, j);
	endLocalClock(&cycleClocks, 0, "ComputePitchStrength");
	//end clock 10
	//Interpolate pitch strength at desired times
	if(DEBUG==1)printf("\n	Interpolate pitch strength at desired times...\n");
	matrix Si2 = zerom(Si.x, t.x);
	vector g = zerov(Si.y);
	if(Si.y > 1 && isParent()){
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
	//End clock 11
	endLocalClock(&cycleClocks, 0,"InterpolatePitchStrength");
	//Compute contribution of this window size to pitch strength
	if(DEBUG==1)printf("\n	Compute contribution of this window size to pitch strength...\n");
	vector mu = onesv(j.x);
	for(ind = 0; ind < k.x; ind++){
		mu.v[(int)k.v[ind]] = 1 - fabs(d.v[(int)j.v[(int)k.v[ind]]]-(i+1));
	}
	//The parent process makes its own job, in case that there is only one process, runs normally
	if(isParent()){
		for(ind = 0; ind < j.x; ind++){
			for(p = 0; p < S.y; p++){
			 S.m[(int)j.v[ind]][p]  += mu.v[ind]*Si2.m[ind][p];
			}
		}
	}
	endLocalClock(&cycleClocks, 0, "ComputeContributionWindowSizePitchStrength");
	//sends data to parent
	sendDataToLeader(Si, j, mu, dn, n.x);
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

