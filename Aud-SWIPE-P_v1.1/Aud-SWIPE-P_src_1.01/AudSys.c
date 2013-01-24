	/*
	 * Name: AudSWIPE-P
	 * Authors: Prof. Arturo Camacho
	 * 			Bsc. Saul Calderon
	 * 			Bsc. Gabriel Alvarado
	 * General Description: Parallel C implementation of the AudSWIPE algorithm, originally
	 * implemented in MATLAB
	 * Archive Description: Applies the filter of the auditive system and produces differents signals
	 * in response to the different cochlea's segments.
	 */
#include "AudSys.h"
#define DEBUG 0


/*
 * Rounds to the nearest integer value
 * @param x, double value to round up
 * */

double round2(double x) { // Rounds a double to the nearest integer value
    return(x >= 0. ? floor(x + .5) : floor(x - .5));
}


/*
*Prints the received array
*@param arr, the array to print
*@param size, the size of the array
*@param text, the name of the array
*/
void print(double* arr, int size, char* text){
	int i;
	for(i = 0; i < size; i++){
		printf("%f \n", arr[i]);
	}
	printf("\n%s\n",text);
}


/*
*Prints the received array in an archive
*@param arr, the array to print
*@param size, the size of the array
*@param texto, the name of the array
*/
void out(double* arr, int size, char* arch){
	FILE * f;
	int i;
	if(f = fopen(arch, "w")){
		for(i = 0; i < size; i++){
			if(arr[i]>=0)
				fprintf(f, " %.15f\n", arr[i]);
			else
				fprintf(f, "%.15f \n", arr[i]);
	}
	fclose(f);
	}else{
		printf("unable to create file %s", arch);
	}
	
}

/*
*Prints the received matrix in an archive
*@param arr, the matrix to print
*@param sizex, number of arr's rows
*@param sizey, numer of arr's columns
*@param arch, archive's name
*/
void outm(double** arr, int sizex, int sizey, char* arch){
	FILE * f;
	int i,j;
	if(f = fopen(arch, "w")){
		for(i = 0; i < sizex; i++){
			for(j = 0; j < sizey; j++){
				if(arr[i][j]>=0)
					fprintf(f, "%f;", arr[i][j]);
				else
					fprintf(f, "%f;", arr[i][j]);
			}
			fprintf(f, " \n ", arr[i][j]);
		}
		fclose(f);
	}else{
		printf("unable to create file %s", arch);
	}
	
}

/*
*Calculates the maximum element in a given array
*@param r, the array
*@param l, the maximum element to take into account in the search
*/
double Max_v(vector x){
	int f;
	double max = 0;
	for(f = 0; f < x.x; f++){
		if(x.v[f] > max)
			max = x.v[f];
	}
	return max;
}

/*
*Aligns the channels of the signal, necessary due to the different response times of each band of frequencies on the cochlea
*@param X, the array
*@param f, the frequency scale
*@param fs, the sampling frequency
*@return X, the matrix with the channels aligned
*/
matrix Align_Channels(matrix X, vector f, double fs){
	vector r = zerov(X.x);
	int l = 0;
	int i;
	for(i = 0; i < r.x; i++){
		r.v[i] = round( fs * 4 / ( 2*pi*1.019*( 24.7+0.108*f.v[i] ) ) );
	}
	l = round(X.y - Max_v(r)+1);
	matrix Y = zerom(X.x, l);
	int c;
	int k;
	for(i = 0; i < Y.x; i++){
		k = r.v[i];
		for(c = 0; c < Y.y; c++, k++){
			Y.m[i][c] = X.m[i][k];
		}
	}
	freev(r);
	return Y;
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
*Separates groups of consecutive harmonics into channels
*It uses the ERB filters to produce each channel as a response of a different cochlea's segment
*@param x, the matlab array x
*@param fs, the sampling frequency
*/
matrix Harmonics_into_channels(vector x, double fs, vector f){
	int i;
	matrix fcoefs = zerom(10, f.x);
	//matrix X, channels per time
	matrix X = makem(f.x, x.x);
	for(i = 0; i < f.x; i++) f.v[i] = i + 1.5;
	v_erbs2hz(f);
	if(DEBUG == 1)printf("\n		Creating ERBFilters...\n");
	ERBFilters(fs, f, fcoefs);
	if(DEBUG == 1)printf("\n		Filtering with signal (Paralelized into Segments of the cochlea)...\n");
	ERBFilterBank(x, fcoefs, X);
	freem(fcoefs);
	return X;
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
	int chunkSize = (int)((size_x+L+((L+1)/2)-1)- (L+((L+1)/2)-1));
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
	double min = f.v[0];
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
*Calculate the filter of the middle and outter ear
*Flattens the spectral envelope
*@param n,
*@param fs, the sampling frequency
*@return vector, the resulting filter
*/
vector outmidear( int n, double fs){
	vector b = zerov(n+1);
	int end = 20;
	double f[] = {	0 *1000,.02 *1000 ,.05   *1000, 
					.1 *1000 ,.2 *1000 , .5  *1000, 
					.6 *1000 , .7  *1000, 1  *1000, 
					2  *1000, 3 *1000,4.0 *1000 , 
					5 *1000,  8 *1000  , 9  *1000, 
					10 *1000 ,12 *1000 , 13 *1000,  
					14 *1000 ,15 *1000, 0 };
	
	double g[] = { 	0, -39 ,-19 ,-12.5 ,-8 , -2 , 
					-1  , 0 , 0 , 4,  8 ,  8 , 5, 
					-9 ,-11 ,-11, -9 ,-13 ,-19 ,-15 };
					
	double m[] = { 	0  ,  0.0112  ,  0.1122  ,  0.2371 ,
					0.3981 ,   0.7943 ,   0.8913 ,   1.0000 ,   
					1.0000  ,  1.5849  ,  2.5119  ,  2.5119 ,   
					1.7783 , 0.3548  ,  0.2818 ,   0.2818  ,  
					0.3548  ,  0.2239 ,   0.1122 ,   0.1778 , 0 };

	vector f2 = zerov(end);
	vector m2 = zerov(end);
	vector m_min;
	vector f_min;
	if(fs/2 > f[end-1]){
		f[end] = fs/2;
		m[end] = 0;
		int i;
		//We needed to normalize f previously
		for(i = 0; i <= end; i++) f[i] = f[i]/(fs/2);
		fir2(n, f, m, end+1, b.v);
		
	}else{
		double mend = interp(f, m, fs/2, end);
		int cont = 0;
		int i;
		for(i = 0; i < end; i++	){
			if(f[i] < fs/2){
				f2.v[cont] = f[i];
				m2.v[cont] = m[i];
				cont++;
			}
		}
		f_min = zerov(cont+1);
		m_min = zerov(cont+1);
		for(i = 0; i < cont; i++){
			f_min.v[i] = f2.v[i]/(fs/2);
			m_min.v[i] = m2.v[i];
		}
		f_min.v[i] = (fs/2)/(fs/2);
		m_min.v[i] = 	mend;
		//We needed to normalize f previously
		fir2(n, f_min.v, m_min.v, cont+1, b.v);
		freev(f_min);
		freev(m_min);

	}
	freev(f2);
	freev(m2);
	return b;
}



/*
*Models the auditive system response in different channel
*@param x, input signal
*@param samplerate, signal samplerate
*@param clocks, the clocks array to measure execution time
*@return matrix, A matrix with a signal per row, representing the cochlea's response to a different frequency band.
*/
matrix audsys(vector x, double samplerate, clocksArray *clocks){
	int i;
	double fs = samplerate;
	double time2, time3, time4, time5, time6;
	//Flatten the spectral envelope
	if(DEBUG == 1)printf("\n	Creating outmidear filter...\n");
	vector b = outmidear(round2(fs/100), fs);
	endLocalClock(clocks, 1, "Creating outmidear filter");
	// End Clock 14
	vector y = zerov(x.x);
	if(DEBUG == 1)printf("\n	Filtering signal with outmidear filter...\n");
	Filter_B(b.v, b.x, x.v, x.x, y.v);
	endLocalClock(clocks, 1, "Filtering signal with outmidear filter");
	// End Clock 15
	//Separate groups of consecutive harmonics into channels
	int tam = 0;
	tam = floor(hz2erbs(fs/2)-.5);//Tamaño de un vector tipo [1.5: hz2erbs(fs/2)];
	vector f = zerov(tam);
	if(DEBUG == 1)printf("\n	Harmonics into channels...\n");
	//matrix X = Harmonics_into_channels(y, fs, tam, f);OJO
	matrix X = Harmonics_into_channels(y, fs, f);
	endLocalClock(clocks, 1, "Harmonics into channels");
	// End Clock 16
	//Align the channels
	int tam_y;
	if(DEBUG == 1)printf("\n	Align Channels...\n");
	matrix Y2;
	Y2 = Align_Channels(X, f, fs);
	//Frees everything
	freev(f);
	for(i = 0; i < tam; i++) free(X.m[i]);
	free(X.m);
	freev(b);
	freev(y);
	endLocalClock(clocks, 1, "AlignChannels");
	// End Clock 19
	return Y2;
}
