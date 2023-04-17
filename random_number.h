#ifndef RANDOM_NUMBER_H
#define RANDOM_NUMBER_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define NormRANu (2.3283063671E-10F)
#define PI 3.141592653589
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


void ini_ran(int SEMILLA);/**Funcion para cambiar la semilla**/
float Parisi_Rapuano(long *idum);
double rand_gaussiano(long* idum);



#endif