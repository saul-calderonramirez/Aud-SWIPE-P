/*
 * SPUtilities.cpp
 *
 *  Created on: Apr 23, 2013
 *      Author: saul
 */

#include "SPUtilities.h"


#define VNUM                         1.0 // Current version

#ifndef NAN
    #define NAN                      sqrt(-1.)
#endif

#ifndef isnan
int isnan(double x) {
    return(x !=x);
}
#endif



/*
 * Rounds to the nearest integer value
 * @param x, double value to round up
 * */

double round2(double x) { // Rounds a double to the nearest integer value
    return(x >= 0. ? floor(x + .5) : floor(x - .5));
}


/*
*Returns the interpolated function
*@param f, the function to interpolate
*@param g
*@param x, x values corresponding to the function f
*@param y, interpolated function
*/
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
	y1 = g.v[i];//This section calculates when x's indexes are inferior than the ones in f
	while(x.v[a]<f.v[i] && a < x.x){
		y[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Lineal interpolation formula
		a++;
	}
	//This section calculates when the values from x are in f's the range.
	while(i < f.x && a < x.x){
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
		y[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Lineal interpolation formula
		a++;
	}
	//This section calculates when the values of x are above the ones from f
	while(a < x.x){
		y[a] = 0;
		a++;
	}
}

/*
*Calculates the logarithm in base 2
*@param x, number to calculate
*@return double, the logatithm in base 2
*/
double log2(double x) { // A base-2 log function
    return(log(x) / log(2.));
}


/*
*Rounds a double to the nearest integer value
*@param x, number to round
*@return double, rounded number
*/
#ifndef round
double round(double x) { // Rounds a double to the nearest integer value
    return(x >= 0. ? floor(x + .5) : floor(x - .5));
}
#endif

/*
*Converts hz value to the mel scale
*@param hz, value in hz
*@return double, mel value
*/
double hz2mel(double hz) { // Converts from hertz to Mel frequency
    return(1127.01048 * log(1. + hz / 700.));
}

/*
*Converts from hertz to ERBs
*@param hz, value in hz
*@return double, ERBs value
*/
double hz2erb(double hz) { // Converts from hertz to ERBs
    return(21.4 * log10(1. + hz / 229.));
}

/*
*Converts from ERBs to hertz
*@param hz, value in hz
*@return double, erb value
*/
double erb2hz(double erb) { // Converts from ERBs to hertz
    return((pow(10, erb / 21.4) - 1.) * 229.);
}

/*
*Converts the hz scale value in to a ERB scale value
*@param hz, the value in hz
*@return double, the value in ERB scale
*/
double hz2erbs(double hz){
	return (21.4 * log10( 1 + hz/229 ));
}

/*
*Converts the ERB scaled array in to a HZ scale array
*@param erbs, the array in ERBS
*@param size, the array's size
*/
void v_erbs2hz(vector erbs){
	int i;
	for(i = 0; i < erbs.x; i++){
		erbs.v[i] = (pow(10, (erbs.v[i]/21.4)) - 1 ) * 229;
	}
}

/*
*A silly function that treats NaNs as 0.
*@param hz, value in hz
*@return double, erb value
*/
double fixnan(double x) { // A silly function that treats NaNs as 0.
    return(isnan(x) ? 0. : x);
}

/*
*Returns the maximum between the received numbers
*@param x, number to compare
*@param y, number to compare
*@return double,result
*/
double maxim(double x, double y){
	return(x>=y ? x : y);
}


/*
*Filter function that receives the coeficients of the characteristic function of the transfer function
*MATLAB like filter function
*@param b, filter coeficients
*@param L,
*@param x, input signal
*@param y, output signal
*@param tam, the size of the vector x
*/
void Filter_B(double * b, int L, double * x  ,int size_x, double* y){
	int i;
	int j = 0;
	double * xp = (double*)malloc(sizeof(double)*(size_x+L+L));
	for(i = 0; i < L; i++){ xp[i] = 0; xp[size_x+L+L-i-1] = 0;}
	for(i = L; i < size_x+L; i++){
		xp[i] = x[i-L];
	}
	int k = 0;
//	int chunkSize = (int)((size_x+L+((L+1)/2)-1)- (L+((L+1)/2)-1));
	for(i = L+((L+1)/2)-1; i < size_x+L+((L+1)/2)-1; i++){//Applies the filter
		y[k] = 0;
		for(j = 0; j < L ; j++){
			y[k] = y[k] + (b[j]*xp[i-j]);
		}
		k++;
	}
	free(xp);
}
/*
*Filter function that receives the coeficients of A and B's polynoms from the transfer function
*MATLAB like filter function
*@param b, b's filter coeficients
*@param a, a's filter coeficients
*@param L,
*@param x, input signal
*@param y, output signal
*@param tam, the size of the vector x
*/
void Filter_AB(double * b, double * a, int L, double * x  ,int size_x, double* y){
	int i;
	int j = 0;
	double * xp = (double*)malloc(sizeof(double)*(size_x+L));
	double * yp = (double*)malloc(sizeof(double)*(size_x+L));
	for(i = 0; i < L; i++){ xp[i] = 0; yp[i] = 0;}
	for(i = L; i < size_x+L; i++){
		xp[i] = x[i-L];
		y[i-L] = 0;
		yp[i] = 0;
	}
	int k = 0;
	for(i = L; i < size_x+L; i++, k++){
		for(j = 0; j < L ; j++){
			yp[i] = yp[i] + (b[j]*xp[i-j]);
		}
		y[k] = yp[i];
		for(j = 1; j < L ; j++){
			y[k] = y[k] -(a[j]*yp[i-j]);
			yp[i] = y[k];
		}
	}
	free(xp);
	free(yp);
}



/*
*Normalizes an array, with the minimum frequency at f[0] and the maximum at f[size-1]
*@param f, array of frecuencies
*@param size, f's size
*/
void normalize(vector f){
	int i;
	//double min = f.v[0];
	double max = f.v[f.x-1];
	for(i = 0; i < f.x; i++){
		//f[i] = (f[i]-min) / (max - min);, bug corrected
		f.v[i] = f.v[i] / max;
	}
}

/*
*Uses lineal interpolation to interpolate the function given in f
*@param f, function to interpolate
*@param g,
*@param x,
*@param end,
*@return double, the interpolated function
*/
double interp(double* f, double* g, double x, int end){
	double y;
	int i = 0;
	while(f[i] < x){
		i++;
	}
	double x0 = f[i-1];
	double x1 = f[i];
	double y0 = g[i-1];
	double y1 = g[i];
	y = y0 + (((x- x0)*y1 -(x -x0)*y0)/(x1-x0));//Formula de interpolacion lineal
	return y;
}


/*
*Upsample the signals given in the X vector.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Upsampled matrix
*@param size_x
*/
void Upsample_Vector(double* X, double* Y, int size_x){
	int tam_h = 42;
	double  h[] = {0,
		-1.4311e-18,-0.0021018,3.7067e-18,0.0050153,
		-6.9847e-18,-0.0097846,1.1207e-17,0.017104,
		-1.6174e-17,-0.027965,2.1551e-17,0.044023,
		-2.6905e-17,-0.068645,3.1753e-17,0.11052,
		-3.5622e-17,-0.20176,3.8119e-17,0.63307,1,
	0.63307,3.8119e-17,-0.20176,-3.5622e-17,0.11052,
	3.1753e-17,-0.068645,-2.6905e-17,0.044023,2.1551e-17,
	-0.027965,-1.6174e-17,0.017104,1.1207e-17,-0.0097846,
	-6.9847e-18,0.0050153,3.7067e-18,-0.0021018,-1.4311e-18 };
	int c, i;
	//Fills with 0's
	i = 0;
	for(c = 0; c < size_x; c = c+2){
		Y[c] = X[i];
		i++;
	}
	//Applies the filter
	Filter_B(h, tam_h, Y, size_x, Y);
}


/*
*Reduces the sample rate of the signals given in the X Vector.
*@param X, Matrix with a signal per row, representing the cochlea's response in frequency band
*@param Y, the Downsampled matrix
*@param size_x
*/
void Downsample_Vector(double* X, double* Y, int size_x){
	int c, i;
	int tam_h = 43;
	double  h[] = {0,0,-7.1553e-19,-0.0010509,1.8533e-18,0.0025077,
	-3.4923e-18,-0.0048923,5.6036e-18,0.0085521,-8.0868e-18,
	-0.013983,1.0776e-17,0.022012,-1.3453e-17,-0.034322,1.5876e-17,
	0.05526,-1.7811e-17,-0.10088,1.9059e-17,0.31654,0.5,0.31654,1.9059e-17,
	-0.10088,-1.7811e-17,0.05526,1.5876e-17,-0.034322,-1.3453e-17,0.022012,
	1.0776e-17,-0.013983,-8.0868e-18,0.0085521,5.6036e-18,-0.0048923,-3.4923e-18,
	0.0025077,1.8533e-18,-0.0010509,-7.1553e-19 };
	//Applies the filter
	Filter_B(h, tam_h, X, size_x*2, X);
	//Removes samples
	i = 0;
	for(c = 0; c < size_x*2; c = c+2){
		Y[i] = X[c];
		i++;
	}
}


/*
 *Returns a Hanning window in y
 *@param y, array to store the hanning window
 *@param np,
 */
void Hanning(double *y, int np) {
	int i;
	int a = 0;
	for (i = 1; i <= np; i++, a++) {
		y[a] = .5 * (1 - cos((2 * pi * i) / (np + 1))); //Hanning en lugar de Hamming
	}
}

/*
 *Returns the difference between the consecutive elements of an array
 *@param y, array y
 *@param f, array f
 *@param sz, final index to look for
 */
void diff(double *y, double *f, int sz) {
	--sz;
	int i;
	for (i = 0; i < sz; i++)
		f[i] = y[i + 1] - y[i];
}

/*
 * Returns 1 if the array has a negative value
 * @param y, array to check
 * @param sz,
 * */
int isneg(double *y, int sz) {
	int i;
	for (i = 0; i < sz; i++)
		if (y[i] < 0)
			return 1;
	return 0;
}
/*
 *Rounds the double number
 @param n, number to round
 @return int
 */
double fix(double n) {
	return (int) n;
}
/*
 *   FIR filter design using the window method - arbitrary filter shape.
 *   FIR2(N,F,M,NF,B) designs an N'th order FIR digital filter with the
 *   frequency response specified by vectors F and M (size NF) and returns the
 *   filter coefficients in length N+1 vector B.  Vectors F and M specify
 *   the frequency and magnitude breakpoints for the filter such that
 *   PLOT(F,M) would show a plot of the desired frequency response.
 *   The frequencies in F must be between 0.0 < F < 1.0, with 1.0
 *   corresponding to half the sample rate. They must be in increasing
 *   order and start with 0.0 and end with 1.0.
 *
 *   The filter B is real, and has linear phase, i.e., even symmetric
 *   coefficients obeying B(k) =  B(N+2-k), k = 1,2,...,N+1.
 *
 *   FIR2 uses a Hamming window.
 *@param n, filter order
 *@param f, frequency breakpoints of the filter
 *@param m,
 *@param nf,
 *@param b,
 */
void fir2(int nn, double *ff, double *aa, int ffsz, double *b) {
	if (nn % 2 > 0) {
		nn--;
	}
	int npt, lap;
	nn = nn + 1;
	if (nn < 1024)
		npt = 512;
	else
		npt = pow(2, ceil(log(nn) / log(2)));
	lap = fix(npt / 25);

	int nbrk = ffsz;
	if (abs(ff[0]) > 0 || abs(ff[nbrk - 1] - 1) > 1) {
		perror("The first frequency must be 0 and the last 1");
		abort();
	}
	// interpolate breakpoints onto large grid
	int nint = nbrk - 1;
	double *df = (double *) malloc(sizeof(double) * (ffsz - 1)); //new double[ffsz];
	diff(ff, df, ffsz);
	//printfir(ff, ffsz, "ff");
	//printfir(df, ffsz-1, "df");
	if (isneg(df, ffsz - 1)) {
		perror("Frequencies must be non-decreasing");
		abort();
	}
	npt = npt + 1;   // Length of [dc 1 2 ... nyquist] frequencies.
	int nb = 0;
	int ne = 0;
	double inc;
	double *H = (double*) malloc(sizeof(double) * npt); //new double[npt];// last index overflaw !!!!!
	H[0] = aa[0];
	int i;
	for (i = 0; i < nint; i++) {
		if (df[i] == 0) {
			nb = nb - lap / 2;
			ne = nb + lap;
		} else
			ne = fix(ff[i + 1] * npt) - 1;
		if (nb < 0 || ne > npt) {
			perror(
					"Too abrupt an amplitude change near end of frequency interval");
			abort();
		}
		int j;
		for (j = nb; j <= ne; j++) {
			if (nb == ne)
				inc = 0;
			else
				inc = (double) (j - nb) / (double) (ne - nb);
			H[j] = inc * aa[i + 1] + (1.0 - inc) * aa[i];
		}
		nb = ne + 1;
	}
	// Fourier time-shift.
	double dt = 0.5 * (nn - 1);
	double complex
	*Hz1 = (double complex*)malloc(sizeof(double complex)*npt); //new complex[npt];
	double complex
	*Hz2 = (double complex*)malloc(sizeof(double complex)*npt); //new complex[npt];
	//double complex *Hz = (double complex*)malloc(sizeof(double complex)*2*npt);//new complex[2*npt];
	fftw_complex* Hz = fftw_malloc(sizeof(fftw_complex) * 2 * npt);
	//int sc;
	for (i = 0; i < npt; i++) {
		double rad = -dt * pi * (double) i / ((double) npt - 1);
		double complex
		z = (H[i] * cos(rad)) + (H[i] * sin(rad) * I);
		Hz1[i] = z;
	}
	for (i = 0; i < npt; i++)
		Hz2[npt - 1 - i] = conj(Hz1[i]);
	for (i = 0; i < npt; i++)
		Hz[i] = Hz1[i];
	for (i = npt; i < npt * 2; i++)
		Hz[i] = Hz2[i - npt];
	int nfft = npt * 2 - 2;
	double* fo = fftw_malloc(sizeof(double) * 2 * npt);
	fftw_plan plan = fftw_plan_dft_c2r_1d(nfft, Hz, fo, FFTW_ESTIMATE);
	fftw_execute(plan);
	double *wind = (double*) malloc(sizeof(double) * nn);      //new double[nn];
	Hanning(wind, (2 * ceil(nn / 2)) + 1);

	double kfft = 1. / nfft;
	for (i = 0; i < nn; i++)
		b[i] = wind[i] * fo[i] * kfft;

	fftw_destroy_plan(plan);
	fftw_free(Hz);
	free(wind);
	fftw_free(fo);
	free(df);
	free(H);
	free(Hz1);
	free(Hz2);
}

/*
 * Half wave rectification
 * Receives the array
 * */
void MaxRCerosVector(double* x, int size_x){
	int  c;
	for(c = 0; c < size_x; c++){
		if(x[c] < 0)
			x[c] = 0;
	}
}
/*
*Calculates the maximum element in a given array
*@param r, the array
*@param l, the maximum element to take into account in the search
*/
void Max(double** x, int size_c, int size_x, double lim){
	int f, c;
	for(f = 0; f < size_c; f++){
		for(c = 0; c < size_x; c++){
			if(x[f][c] < 0)
				x[f][c] = 0;
		}
	}
}

vector readSoundFile(char ptrName[], double* ptrNyquist, double* ptrFs, double* ptrSoundLength, int* ptrFrames){
	FILE* wavf;
	SF_INFO info;
	SNDFILE* source;
	vector x;
	x.v = NULL;
	x.x = 0;
	//opening file time
	if (strcmp(ptrName, "<STDIN>") == 0) { // i.e., is coming from STDIN
	   wavf = stdin;
	}
	else { // is specified
	   wavf = fopen(ptrName, "r");
	}
	if(wavf != NULL){
		//Stores the wav file info on the
		source = sf_open_fd(fileno(wavf), SFM_READ, &info, TRUE);
		// Perform checks on the wav header
		if (info.sections < 1) {
		   fprintf(stderr, "File or stream %s not read as audio ... \n", ptrName);
		}
		else{
			*ptrNyquist = info.samplerate / 2.;
			*ptrFs = info.samplerate;
			int frames = (int)info.frames;
			*ptrSoundLength = (double)((double)frames/ (double)info.samplerate);
			x = makev(info.frames); // This reads the signal in
			sf_read_double(source, x.v, x.x);
			sf_close(source); // Takes FILE* wavf with it, too
			fclose(wavf);
			*ptrFrames = frames;
		}
	}
	return x;
}



