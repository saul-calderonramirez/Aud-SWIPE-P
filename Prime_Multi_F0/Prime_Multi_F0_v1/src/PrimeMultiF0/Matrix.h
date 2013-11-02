/*
 * Arrays.h
 *
 *  Created on: Apr 23, 2013
 *      Author: saul
 *      / vector stuff
 Copyright (c) 2009-2011 Kyle Gorman
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
*  vector: some data structures for swipe
*  Kyle Gorman <kgorman@ling.upenn.edu>
*
*  This is version 1.0., i.e. I think I got all obvious stuff working ideally.
*/


#ifndef ARRAYS_H_
#define ARRAYS_H_
#pragma once
// prime sieve
#define P                          1
#define NP                         0
#define PRIME(x)                   (x == 1)

// cubic spline
#define YP1                        2.
#define YPN                        2.

#include "includes.h"

#define DEBUG 0

// vector stuff
typedef struct{
		int x;
		double* v;
} vector;

// intvector stuff
typedef struct{
	int x;
	int* v;
} intvector;

typedef struct{
	int x;
	int y;
	double** m;
} matrix;

// intmatrix stuff
typedef struct{
	int x;
	int y;
	int** m;
} intmatrix;

void outBinaryM(double** m, int x, int y, char file[]);
vector makev(int);
vector zerov(int);
vector onesv(int);
vector nansv(int);
vector copyv(vector);
void resizev(vector *, int);
int  maxv(vector);
int  minv(vector);
int  bisectv(vector, double);
int  bilookv(vector, double, int);
void freev(vector);
void printv(vector);
intvector makeiv(int);
intvector zeroiv(int);
intvector onesiv(int);
intvector copyiv(intvector);
vector  iv2v(intvector);
int maxiv(intvector);
int miniv(intvector);
int bisectiv(intvector, int);
int bilookiv(intvector, int, int);
void freeiv(intvector);
void printiv(intvector);
// matrix stuff
matrix makem(int, int);
matrix zerom(int, int);
matrix onesm(int, int);
matrix nansm(int, int);
vector* mat2vect(matrix);
matrix  copym(matrix);
void freem(matrix);
void printm(matrix);
matrix add(matrix, matrix);
intmatrix makeim(int, int);
intmatrix zeroim(int, int);
intmatrix onesim(int, int);
intmatrix copyim(intmatrix);
void im2m(intmatrix); // cast
void freeim(intmatrix);
void printim(intmatrix);
int sieve(intvector);
intvector primes(int);
vector spline(vector, vector);
double splinv(vector, vector, vector, double, int);
// polynomial fitting
vector polyfit(vector, vector, int);
double polyval(vector, double);
/*
*Writes the given vector in a file
*@param m, the matrix to write
*@param x, number of rows
*@param file, file name
*/
void outBinaryV(double* v, int x, char file[]);
/*
*Calculates the maximum element in a given array
*@param r, the array
*@param l, the maximum element to take into account in the search
*/
double Max_v(vector x);
#endif
