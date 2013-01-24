#include <math.h>
#include <stdio.h>
#include "ERBFilters.h"
#include "complex.h"

void ERBFilters(double fs, vector cf, /*int size_cf,*/ matrix fcoefs){

	double T = 1/fs;
	
	double EarQ = 9.26449;
	double minBW = 24.7;
	double order = 1;
	//double * ERB = (double*)malloc(sizeof(double)*size_cf);
	vector ERB = zerov(cf.x);
	//double * B = (double*)malloc(sizeof(double)*size_cf);
	vector B = zerov(cf.x);
	int i;
	
	for(i = 0; i < cf.x/*size_cf*/; i++){
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


void ERBFilterBank(vector x, /*int size_x,*/ matrix fcoefs/*, int size_cf*/, matrix X){
	int i;
	
	//double * a = (double*)malloc(sizeof(double)*3);
	//double * b = (double*)malloc(sizeof(double)*3);
	//double a[3];// = double[3];
	
	vector a = zerov(3);
	vector b = zerov(3);
	//double b[3];// = double[3];
	//#pragma omp parallel for private(a,b)
	for(i = 0; i < fcoefs.y; i++){
		//printf("numero de hilos: %d\n", omp_get_num_threads());
		a.v[0] = fcoefs.m[0][i]/fcoefs.m[9][i];
		a.v[1] = fcoefs.m[1][i]/fcoefs.m[9][i];
		a.v[2] = fcoefs.m[5][i]/fcoefs.m[9][i];
		b.v[0] = fcoefs.m[6][i];
		b.v[1] = fcoefs.m[7][i];
		b.v[2] = fcoefs.m[8][i];
		//printf("paso aqui!!!%d\n", x.x);
		Filter_AB(a, b, /*3,*/ x.v, x.x, X.m[i]);
		a.v[0] = fcoefs.m[0][i];
		a.v[1] = fcoefs.m[2][i];
		a.v[2] = fcoefs.m[5][i];
		Filter_AB(a, b, /*3,*/ X.m[i], X.x, X.m[i]);
		
		a.v[1] = fcoefs.m[3][i];
		
		Filter_AB(a, b, /*3,*/ X.m[i], X.x, X.m[i]);
		
		a.v[1] = fcoefs.m[4][i];
		
		Filter_AB(a, b, /*3,*/ X.m[i], X.x, X.m[i]);
		
	}
	
	freev(a);
	freev(b);
	
}






