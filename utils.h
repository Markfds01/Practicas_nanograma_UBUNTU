#ifndef UTILS_HEADER
#define UTILS_HEADER

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "stdbool.h"


void lee_param(int ac,char **av,int* V,double* k_o,double* k_m,double* gamma,double* beta,double* mu,double* prob_travel,int* time_max,double* fraction_people_class_one,double* social_activity,double* mixing_probability);
int ****allocate_memory_4D(int a, int b, int c, int d);
void free_memory_4D(int ****array, int a, int b, int c);
int ***allocate_memory_3D(int a, int b, int c);
void free_memory_3D(int ***array, int a, int b);
void printf_3D_array(int*** pop,int b,int c);
double* allocate_memory_1D(int N);
void total_states(int*** pop,int time,int V,int* total_susceptible,int* total_infected,int* total_recovered);
double** allocate_contact_matrix();
void free_memory_contact_matrix(double** contact_matrix);
double binomial(int n, int k);
bool* initialite_bool_array(int V);
int **allocate_memory_2D(int rows, int cols);
void free_memory_2D(int **array, int rows);
void matriz_a_cero(int a,int b,int** matrix);















#endif