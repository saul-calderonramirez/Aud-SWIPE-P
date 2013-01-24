#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include "ERBFilters.h"

//X es la señal dividida en segmentos de la coclea, Y es la salida, size_x es el tamaño de la señal, size_c es la cantidad de segmentos de la coclea
void Downsample(matrix X, matrix Y/*, int size_x, int size_c*/);
void Upsample(matrix X, matrix Y/*, int size_x, int size_c*/);
