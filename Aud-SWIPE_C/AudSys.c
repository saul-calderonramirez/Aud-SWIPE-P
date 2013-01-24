//Aplica el filtro de sistema auditivo y divide la señal en segmentos de la coclea para
//luego ser procesadas por swipe

#include "AudSys.h"

struct timespec start, end;

double tiempo;
//Imprime un vector en la consola
void print(vector vec, char* texto){
	
	int i;
	for(i = 0; i < vec.x; i++){
		printf("%f \n", vec.v[i]);
	}
	printf("\n%s\n",texto);
	

}
double timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p){
  return ((double)(timeA_p->tv_sec ) + (double)timeA_p->tv_nsec/(double)1000000000) -
           (((double)timeB_p->tv_sec ) + (double)timeB_p->tv_nsec/(double)1000000000);
}
//Imprime un vector en un archivo, para facilitar la importacion a matlab
void out(vector vec, char* arch){

	FILE * f;
	int i;
	if(f = fopen(arch, "w")){
		for(i = 0; i < vec.x; i++){
			if(vec.v[i]>=0)
				fprintf(f, " %.15f \n", vec.v[i]);
			else
				fprintf(f, "%.15f \n", vec.v[i]);
	}
	fclose(f);
	}else{
		printf("unable to create file %s", arch);
	}
	
}
void outm(double** arr, int sizex, int sizey, char* arch){

	FILE * f;
	int i,j;
	if(f = fopen(arch, "w")){
		for(i = 0; i < sizex; i++){
			for(j = 0; j < sizey; j++){
				if(arr[i][j]>=0)
					fprintf(f, "%f ", arr[i][j]);
				else
					fprintf(f, "%f ", arr[i][j]);
			}
			fprintf(f, " \n ", arr[i][j]);
		}
		fclose(f);
	}else{
		printf("unable to create file %s", arch);
	}
	
}

//Calcula el maximo en un vector, devuelve el valor del maximo
double Max_v(vector x){

	int f;
	double max = 0;
	for(f = 0; f < x.x; f++){
		if(x.v[f] > max)
			max = x.v[f];
	}
	return max;
}

//Alinea los canales de X luego de dividirla en los segmentos de la coclea
//double** Align_Channels(double ** X, int size_x, int size_c, double *f, double fs, int *tam_y){
matrix Align_Channels(matrix X, vector f, double fs){
	//double* r = (double*)malloc(sizeof(double)*size_c);
	vector r = zerov(X.x);
	int l = 0;
	int i;
	for(i = 0; i < r.x; i++){
		r.v[i] = round( fs * 4 / ( 2*pi*1.019*( 24.7+0.108*f.v[i] ) ) );
	}
	l = round(X.y - Max_v(r)+1);
	
	
	matrix Y = zerom(X.x, l);//(double**)malloc(sizeof(double*)*size_c);
	//for(i = 0; i < size_c; i++) Y[i] = (double*)malloc(sizeof(double)*l);
	
	int c;
	int k;
	for(i = 0; i < Y.x; i++){
		k = r.v[i];
		for(c = 0; c < Y.y/*-1*/; c++, k++){
			Y.m[i][c] = X.m[i][k];
		}
	}
	freev(r);
	return Y;
}

//Imita la funcion max(x,0) de matlab, x es una matriz
void Max(matrix X){
	int f, c;
	for(f = 0; f < X.x; f++){
		for(c = 0; c < X.y; c++){
			if(X.m[f][c] < 0)
				X.m[f][c] = 0;
		}
	}
}

//Hz a erbs
double hz2erbs(double hz){
	return (21.4 * log10( 1 + hz/229 ));
}

//Transforma un vector de erb a un vector de Hz
void v_erbs2hz(vector erbs/*double* erbs, int size*/){
	int i;
	for(i = 0; i < erbs.x/*size*/; i++){
		erbs.v[i] = (pow(10, (erbs.v[i]/21.4)) - 1 ) * 229;
	}
}
//Separate groups of consecutive harmonics into channels
matrix Harmonics_into_channels(vector x, double fs, /*int tam,*/ vector f){
	
	out(x, "vectorXC.csv");
	
	int i;
	//double ** fcoefs = (double**)malloc(sizeof(double*)*10);
	matrix fcoefs = zerom(10,f.x/*tam*/);
	matrix X = makem(f.x/*tam*/, x.x);
	//for(i = 0; i < 10; i++) fcoefs[i] = (double*)malloc(sizeof(double)*tam);
	
	for(i = 0; i < f.x/*tam*/; i++) f.v[i] = i + 1.5;
	v_erbs2hz(f/*, tam*/);
	printf("\n		Creating ERBFilters...\n");
	ERBFilters(fs, f/*,tam*/, fcoefs);
	
	
	//outm(fcoefs, 10, tam, "fcoefsC.csv");
	//Aqui obtiene la matriz gracias a los procesos hijos
	printf("\n		Filtering with signal (Paralelized into Segments of the cochlea)...\n");
	tiempo = 0.0;
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	ERBFilterBank(x, fcoefs, /*tam,*/ X);
	clock_gettime(CLOCK_MONOTONIC, &end);//end of execution time
	tiempo = timespecDiff(&end, &start);
	
	
	//for(i = 0; i < 10; i++) free(fcoefs[i]);
	freem(fcoefs);
	return X;
}

//Funcion que imita Filter de Matlab solo con coeficientes B, agregando zf al final de la entrada
void Filter_B(vector b, double* x  ,int size_x, double* y){
	
	//out(x, "xFilterB.csv");

	int i;
	int j = 0;
	vector xp = zerov(size_x+b.x+b.x);
	//double * xp = (double*)malloc(sizeof(double)*(size_x+L+L));
	//for(i = 0; i < L; i++){ xp[i] = 0; xp[size_x+L+L-i-1] = 0;}
	for(i = b.x; i < size_x+b.x; i++){
		xp.v[i] = x[i-b.x];
	}
	int k = 0;
	for(i = b.x+((b.x+1)/2)-1; i < size_x+b.x+((b.x+1)/2)-1; i++, k++){//Aplica el filtro
		y[k] = 0;
		for(j = 0; j < b.x ; j++){
			y[k] = y[k] + (b.v[j]*xp.v[i-j]);
		}
	}
	freev(xp);
}

//Funcion que imita Filter de Matlab con coeficientes A y B
void Filter_AB(vector b, vector a, /*int L,*/ double* x  ,int size_x, double* y){

	int i;
	int j = 0;
	//double * xp = (double*)malloc(sizeof(double)*(size_x+L));
	vector xp = zerov(size_x+b.x);
	//double * yp = (double*)malloc(sizeof(double)*(size_x+L));
	vector yp = zerov(size_x+b.x);
	
	for(i = 0; i < b.x; i++){ xp.v[i] = 0; yp.v[i] = 0;}
	
	for(i = b.x; i < size_x+b.x/*size_x+L*/; i++){
		xp.v[i] = x[i-b.x];
		y[i-b.x] = 0;
		yp.v[i] = 0;
	}
	
		
	int k = 0;
	for(i = b.x; i < size_x+b.x/*size_x+L*/; i++, k++){
		for(j = 0; j < b.x ; j++){
			yp.v[i] = yp.v[i] + (b.v[j]*xp.v[i-j]);
		}
		
		y[k] = yp.v[i];
		
		for(j = 1; j < b.x ; j++){
			y[k] = y[k] -(a.v[j]*yp.v[i-j]);
			yp.v[i] = y[k]; 
		}
	}
	freev(xp);
	freev(yp);
}

//Funcion para normalizar un vector, suponemos que la frecuencia minima está en f[0] y la maxima en f[size-1]; 
void normalizar(vector f){
	int i;
	double min = f.v[0];
	double max = f.v[f.x-1];
	for(i = 0; i < f.x; i++){
		//f[i] = (f[i]-min) / (max - min);
		f.v[i] = f.v[i] / max;

	}	
}

//Interpolación lineal
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
//Funcion que calcula el filtro del oido medio y externo
//Flatten the spectral envelope
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
					
	//double* f2 = (double*) malloc(sizeof(double)*(end));
	vector f2 = zerov(end);
	//double* m2 = (double*) malloc(sizeof(double)*(end));
	vector m2 = zerov(end);
	//double * m_min;
	vector m_min;
	vector f_min;
	//double * f_min;
	if(fs/2 > f[end-1]){
		f[end] = fs/2;
		m[end] = 0;
		int i;
		//Antes de calcular el filtro hay que normalizar f, ya que a fir2 solo le gustan las frecuencias de 0 a 1
		
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
		f_min = zerov(cont+1);//(double*) malloc(sizeof(double)*(cont+1));
		m_min = zerov(cont+1);//(double*) malloc(sizeof(double)*(cont+1));
		for(i = 0; i < cont; i++){
			f_min.v[i] = f2.v[i]/(fs/2);
			m_min.v[i] = m2.v[i];
			
		}
		
		f_min.v[i] = (fs/2)/(fs/2);
		m_min.v[i] = 	mend;
		
		
		//Antes de calcular el filtro hay que normalizar f, ya que a fir2 solo le gustan las frecuencias de 0 a 1
		//normalizar(f_min, cont+1);
		fir2(n, f_min.v, m_min.v, cont+1, b.v);
		
		freev(f_min);
		freev(m_min);
		
	}
	freev(f2);
	freev(m2);
	
	
	return b;
}



matrix audsys(vector x, double samplerate, double *time1){

	int i;
	double fs = samplerate;
	//Flatten the spectral envelope
	printf("\n	Creating outmidear filter...\n");
	vector b = outmidear(round(fs/10), fs);
	
	//print test b
	out(b, "testB.txt");
	
	vector y = zerov(x.x);
	printf("\n	Filtering signal with outmidear filter...\n");
	Filter_B(b, x.v, x.x, y.v);
	
	out(y, "vectorY_C.csv");
	
	//Separate groups of consecutive harmonics into channels
	int tam = 0;
	tam = floor(hz2erbs(fs/2)-.5);//Tamaño de un vector tipo [1.5: hz2erbs(fs/2)];
	//double * f = (double*)malloc(sizeof(double)*tam);
	vector f = zerov(tam);
	printf("\n	Harmonics into channels...\n");
	matrix X = Harmonics_into_channels(y, fs, f);
	outm(X.m, X.x, X.y, "xc.txt");
	
	*time1 = tiempo;
	//Create new harmonics by applying half-wave rectification
	//double ** Y = (double**)malloc(sizeof(double*)*tam);
	matrix Y = zerom(tam, ceil(y.x*2));
	printf("\n	Resample...\n");
	/*
	int j;
	
	for(i = 0; i < tam; i++){
		Y[i] = (double*)malloc(sizeof(double)*ceil(y.x*2));
		for(j = 0; j < ceil(y.x*2); j++){
		Y[i][j] = 0;//Llenamos con 0's de una vez para el upsampling
		}
	}
	*/
	printf("\n		Upsample...\n");
	Upsample(X, Y/*, y.x*2, tam*/);
	printf("\n		Max...\n");
	Max(Y);
	printf("\n		Downsample...\n");
	Downsample(Y, X/*, y.x, tam*/);
	
	//Align the channels
	int tam_y;
	printf("\n	Align Channels...\n");
	/*
	
	for(i = 0; i < tam; i++) free(Y[i]);
	free(Y);
	*/
	freem(Y);
	matrix Y2;
	//Y2.m = Align_Channels(X.m,x.x,tam, f.v, fs, &tam_y);
	Y2 = Align_Channels(X, f, fs);
	//Liberar todo lo que se usó y guardar la nueva señal en X
	freev(f);
	freem(X);
	//for(i = 0; i < tam; i++) free(X.m[i]);
	
	//free(X.m);
	freev(b);
	freev(y);
	//X.m = Y2.m;
	//X.y = tam_y;
	
	return Y2;
}













