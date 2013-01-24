#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "vector.h"
#define pi 3.1415926535897932384626433832795

void ERBFilters(double fs, vector cf, /*int size_cf,*/ matrix fcoefs);

void ERBFilterBank(vector x, /*int size_x,*/matrix fcoefs, /*int size_cf,*/ matrix X);
