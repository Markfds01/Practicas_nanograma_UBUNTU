#ifndef NANOGRAM_UTILS_HEADER
#define NANOGRAM_UTILS_HEADER

#include "network_structs.h"
#include "utils.h"

bool check_conditions(struct Graph* rows_restrictions,struct Graph* columns_restricions,int* number_of_conditions);
int** initialite_nanogram_beggining(int rows,int columns,int number_of_conditions);
void print_nanogram(int** nanogram,int n_rows,int n_columns);
//int* cluster_the_row(int* row,int n_columns,int* n_clusters);
struct ClusterList* cluster_the_row(int* row,int n_columns);
int* copy_vector(int* vector_tobe_copied,int n);
void free_ClusterList(struct ClusterList* cluster_list);
double calcula_energia_una_fila(int *row,int n_columns,struct AdjList* adjlist);
double total_energy_rows(int** nanogram, struct Graph* rows_restrictions,int n_columns);
int** copy_matrix(int** matrix_tobe_copied,int n_rows,int n_columns);
void propose_change(int*** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,int n_changes,struct Graph* rows_restrictions,double beta,int* aceptancia,struct Graph* columns_restrictions);
bool find_one_near(int**nanogram, int row, int column);
void propose_change_near_one(int*** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,int n_changes,struct Graph* rows_restrictions,double beta, int* aceptancia);
void propose_change_per_row(int** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,struct Graph* rows_restrictions,double beta, int* aceptancia);
int number_of_ones(int** nanogram,int n_rows,int n_columns);
bool differences_in_columns(int* column1,int* column2,int n);
void fprint_nanogram(int** nanogram,int n_rows,int n_columns,FILE* f);
double calcula_energia_una_columna(int **nanogram, int column_index, int n_rows, struct AdjList *adjlist);
double calcula_energia_una_columna_with_a_change(int **nanogram, int column_index, int n_rows, struct AdjList *adjlist,int row_index,int value);
void propose_change_per_position(int** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,struct Graph* rows_restrictions,double beta, int* aceptancia,struct Graph* columns_restrictions);
double total_energy_columns(int** nanogram, int n_rows,struct Graph* columns_restrictions,int n_columns);
double total_energy_nanogram(int** nanogram, struct Graph* rows_restrictions,int n_columns,int n_rows,struct Graph* columns_restrictions);
double calcula_energia_una_columna_with_two_change(int **nanogram, int column_index, int n_rows, struct AdjList *adjlist,int row_index,int value,int row_to_go,int value2);
struct ones_positions*  assign_ones_positions(int** nanogram,int number_of_conditions,int n_rows,int n_columns);
void propose_change_per_position_optimized(int** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,struct Graph* rows_restrictions,double beta, int* aceptancia,struct Graph* columns_restrictions,int number_of_conditions,struct ones_positions* positions_of_ones);
bool check_ones_positions(int** nanogram, struct ones_positions* positions_of_ones,int number_of_conditions);
int** initialite_random_nanogram(int n_rows,int n_columns,int number_of_conditions,long* idum);




























#endif