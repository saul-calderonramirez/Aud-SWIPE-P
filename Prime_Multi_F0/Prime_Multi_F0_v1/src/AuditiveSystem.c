/*
 * AuditiveSystem.c
 *
 *  Created on: Apr 23, 2013
 *      Author: saul
 */

#include "AuditiveSystem.h"

/*
*Calculate the filter of the middle and outer ear
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

	/*double g[] = { 	0, -39 ,-19 ,-12.5 ,-8 , -2 ,
					-1  , 0 , 0 , 4,  8 ,  8 , 5,
					-9 ,-11 ,-11, -9 ,-13 ,-19 ,-15 };*/

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
	//erbs2hz is not needed in Prime multi F0
	//no deberia ser hz2erbs
	v_erbs2hz(f);
	if(DEBUG == 1)printf("\n		Creating ERBFilters...\n");
	ERBFilters(fs, f, fcoefs);
	if(DEBUG == 1)printf("\n		Filtering with signal (Paralelized into Segments of the cochlea)...\n");
	ERBFilterBank(x, fcoefs, X);
	freem(fcoefs);
	return X;
}

/*
*Aligns the channels of the signal, necessary due to the different response times of each band of frequencies on the cochlea
*@param X, the array
*@param f, the frequency scale
*@param fs, the sampling frequency
*@return X, the matrix with the channels aligned
*/
matrix alignChannels(matrix X, vector f, double fs){
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
*Generates the ERB filters to be used, to obtain the inner ear response and the missing harmonics
*@param fs, sample frequency
*@param cf,
*@param fcoefs, array that will contain ERBs coeficients
*@param texto, the name of the array
*/
void ERBFilters(double fs, vector cf, matrix fcoefs){
	double T = 1/fs;
	double EarQ = 9.26449;
	double minBW = 24.7;
	double order = 1;
	vector ERB = zerov(cf.x);
	vector B = zerov(cf.x);
	int i;
	for(i = 0; i < cf.x; i++){
		ERB.v[i] = pow(pow(cf.v[i]/EarQ, order) + pow(minBW, order), 1/order);
		B.v[i] = 1.019*2*pi*ERB.v[i];
		fcoefs.m[0][i] = T;//A0
		fcoefs.m[1][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					+ 2*sqrt(3+pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A11
		fcoefs.m[2][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					- 2*sqrt(3+pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A12
		fcoefs.m[3][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					+ 2*sqrt(3-pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A13
		fcoefs.m[4][i] = -(2*T*cos(2*cf.v[i]*pi*T)/exp(B.v[i]*T)
					- 2*sqrt(3-pow(2, 1.5))*T*sin(2*cf.v[i]*pi*T)/exp(B.v[i]*T))/2;//A14
		fcoefs.m[5][i] = 0;//A2
		fcoefs.m[6][i] = 1;//B0
		fcoefs.m[7][i] = (-2*cos(2*cf.v[i]*pi*T))/exp(B.v[i]*T);//B1
		fcoefs.m[8][i] = exp(-2*B.v[i]*T);//B2
		fcoefs.m[9][i] = cabs((-2*cexp(4*I*cf.v[i]*pi*T)*T +
								2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
									(ccos(2*cf.v[i]*pi*T)- csqrt(3 - pow(2, (double)3/2))*
									csin(2*cf.v[i]*pi*T)))*
							(-2*cexp(4*I*cf.v[i]*pi*T)*T +
								2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
								(ccos(2*cf.v[i]*pi*T) + csqrt(3 - pow(2, (double)3/2))*
								csin(2*cf.v[i]*pi*T)))*
							(-2*cexp(4*I*cf.v[i]*pi*T)*T +
								2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
								(ccos(2*cf.v[i]*pi*T) -
								csqrt(3 + pow(2, (double)3/2))*csin(2*cf.v[i]*pi*T)))*
								(-2*cexp(4*I*cf.v[i]*pi*T)*T + 2*cexp(-(B.v[i]*T) + 2*I*cf.v[i]*pi*T)*T*
           						(ccos(2*cf.v[i]*pi*T) + csqrt(3 + pow(2,(double)3/2))*csin(2*cf.v[i]*pi*T))) /
          						cpow((-2 / cexp(2*B.v[i]*T) - 2*cexp(4*I*cf.v[i]*pi*T) +
           						2*(1 + cexp(4*I*cf.v[i]*pi*T))/cexp(B.v[i]*T)),(double)4));//gain

	}
	freev(ERB);
	freev(B);
}

/*
*Applies the ERB filters, and also Upsamples, makes the half wave rectification and finally downsamples the signal, to obtain the missing harmonics
*@param x, input signal
*@param size_x, x's size
*@param fcoefs, array that contains ERBs coeficients
*@param X, the matrix that will contain a signal per row, wich will be the cochlea's response to a specific band
*/
void ERBFilterBank(vector x, matrix fcoefs, matrix X){
	int i;
	double a[3];// = double[3];
	double b[3];// = double[3];
	matrix Y = zerom( fcoefs.y, x.x*2);
	/*
	 * Thread Parallelization 1, One filter per Thread
	 * */
	//Filter-bank parallelization, One filter - One Thread
	#pragma omp parallel for private(a,b)
	for(i = 0; i < fcoefs.y; i++){
		a[0] = fcoefs.m[0][i]/fcoefs.m[9][i];
		a[1] = fcoefs.m[1][i]/fcoefs.m[9][i];
		a[2] = fcoefs.m[5][i]/fcoefs.m[9][i];
		b[0] = fcoefs.m[6][i];
		b[1] = fcoefs.m[7][i];
		b[2] = fcoefs.m[8][i];
		Filter_AB(a, b, 3,  x.v, x.x, X.m[i]);
		a[0] = fcoefs.m[0][i];
		a[1] = fcoefs.m[2][i];
		a[2] = fcoefs.m[5][i];
		Filter_AB(a, b, 3,  X.m[i], x.x, X.m[i]);
		a[1] = fcoefs.m[3][i];
		Filter_AB(a, b, 3, X.m[i], x.x, X.m[i]);
		a[1] = fcoefs.m[4][i];
		Filter_AB(a, b, 3, X.m[i], x.x, X.m[i]);
		//Half wave rectification procedure, to obtain the missing harmonics
		Upsample_Vector(X.m[i], Y.m[i], x.x*2);
		MaxRCerosVector(Y.m[i], x.x*2 );
		Downsample_Vector(Y.m[i], X.m[i], x.x );
	}
}

/*
*Models the auditive system response in different channel
*@param x, input signal
*@param samplerate, signal samplerate
*@param clocks, the clocks array to measure execution time
*@return matrix, A matrix with a signal per row, representing the cochlea's response to a different frequency band.
*/
matrix audsys(vector x, double samplerate, clocksArray *clocks){
	double fs = samplerate;
	vector b, y, f;
	matrix X, Y2;
	//Flatten the spectral envelope
	if(DEBUG == 1)printf("\n	Creating outmidear filter...\n");
	b = outmidear(round2(fs/100), fs);
	endLocalClock(clocks, 1, "Creating outmidear filter");
	// End Clock 14
	//y has the sound signal filtered by the out-midear filter
	y = zerov(x.x);
	if(DEBUG == 1)printf("\n	Filtering signal with outmidear filter...\n");
	Filter_B(b.v, b.x, x.v, x.x, y.v);
	endLocalClock(clocks, 1, "Filtering signal with outmidear filter");
	// End Clock 15
	//Separate groups of consecutive harmonics into channels
	int tam = 0;
	tam = floor(hz2erbs(fs / 2)-.5);//Tamaño de un vector tipo [1.5: hz2erbs(fs/2)];
	f = zerov(tam);
	if(DEBUG == 1)printf("\n	Harmonics into channels...\n");
	//matrix X = Harmonics_into_channels(y, fs, tam, f);OJO
	X = Harmonics_into_channels(y, fs, f);
	endLocalClock(clocks, 1, "Harmonics into channels");
	// End Clock 16
	//Align the channels
	if(DEBUG == 1)printf("\n	Align Channels...\n");
	Y2 = alignChannels(X, f, fs);
	//Frees everything
	freev(f);
	//for(i = 0; i < tam; i++) free(X.m[i]);
	freem(X);
	freev(b);
	freev(y);
	endLocalClock(clocks, 1, "AlignChannels");
	// End Clock 19
	return Y2;
}
