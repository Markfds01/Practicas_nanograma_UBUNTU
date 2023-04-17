#include "utils.h"

//Function to read the main parameters from a external txt file
void lee_param(int ac,char **av,int* V,double* k_o,double* k_m,double* gamma,double* beta,double* mu,double* prob_travel,int* time_max,double* fraction_people_class_one,double* social_activity,double* mixing_probability)
{
    char Nombrfich[256],dummy[256];
    FILE *input_file;
    switch(ac)
    {
    case 1:
        printf("Abrimos fichero por defecto: datos_iniciales.txt");
        sprintf(Nombrfich,"datos_iniciales.txt");
        break;
    case 2:
        scanf(av[1],"%s",Nombrfich);
        break;
     default:
        printf("uso fr_new2 [parameters_file:Def:datos_iniciales.txt] \n");
        exit(1);
    }
    if((input_file=fopen(Nombrfich,"rt"))==NULL)
    {
        printf("Error opening file %s for reading \n",Nombrfich);
        exit(1);
    }
    if(fscanf(input_file,"%d%s",V,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de N_par");
        exit(2);
    }
    if(fscanf(input_file,"%lf%s",k_o,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de k_o");
        exit(3);
    }
    if(fscanf(input_file,"%lf%s",k_m,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de k_m");
        exit(4);
    }
    if(fscanf(input_file,"%lf%s",gamma,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de gamma");
        exit(5);
    }
    if(fscanf(input_file,"%lf%s",beta,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de beta");
        exit(6);
    }
    if(fscanf(input_file,"%lf%s",mu,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de mu");
        exit(7);
    }
    if(fscanf(input_file,"%lf%s",prob_travel,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de prob_travel");
        exit(8);
    }
     if(fscanf(input_file,"%d%s",time_max,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de time_max");
        exit(9);
    }
    if(fscanf(input_file,"%lf%s",fraction_people_class_one,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de fraction_people_class_one");
        exit(10);
    }
    if(fscanf(input_file,"%lf%s",social_activity,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de social_activity");
        exit(11);
    }
     if(fscanf(input_file,"%lf%s",mixing_probability,dummy)!=2)
    {
        printf("no se ha leido correctamente el valor de mixing_probability");
        exit(11);
    }
}

//Function to allocate memory for the main variable which is pop[][][][] 
//¡¡¡4 DIMENSIONS!!! (Include the types of agents)
int ****allocate_memory_4D(int a, int b, int c, int d) {
  int ****array = (int ****)malloc(a * sizeof(int ***));

  for (int i = 0; i < a; i++) {
    array[i] = (int ***)malloc(b * sizeof(int **));

    for (int j = 0; j < b; j++) {
      array[i][j] = (int **)malloc(c * sizeof(int *));

      for (int k = 0; k < c; k++) {
        array[i][j][k] = (int *)malloc(d * sizeof(int));
      }
    }
  }

  return array;
}

//Function to free memory of a 4D array
void free_memory_4D(int ****array, int a, int b, int c) {
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < b; j++) {
      for (int k = 0; k < c; k++) {
        free(array[i][j][k]);
      }
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
}

//Function to allocate memory for the main variable which is pop[][][] 
//¡¡¡3 DIMENSIONS!!! (Include the types of agents)
int ***allocate_memory_3D(int a, int b, int c) {
  int ***array = (int ***)malloc(a * sizeof(int **));

  for (int i = 0; i < a; i++) {
    array[i] = (int **)malloc(b * sizeof(int *));

    for (int j = 0; j < b; j++) {
      array[i][j] = (int *)malloc(c * sizeof(int));
    }
  }

  return array;
}

//Function to free memory of a 3D array
void free_memory_3D(int ***array, int a, int b) {
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < b; j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
}

int **allocate_memory_2D(int rows, int cols) {
    int **array = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++) {
        array[i] = (int *)malloc(cols * sizeof(int));
    }
    return array;
}

// Function to free memory allocated for a 2D integer array
void free_memory_2D(int **array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

void matriz_a_cero(int a,int b,int** matrix)
{
  for(int i=0;i<a;i++)
  {
    for(int j=0;j<b;j++)
    {
      matrix[i][j] = 0;
    }
  }
}

//Allocate dinamical memory to the contact matrix(2x2)
double** allocate_contact_matrix()
{
  //Initialite the first dimension
  double** array = (double**)malloc(2 * sizeof(double*));

  //Loop for the other dimension
  for(int i=0;i<2;i++)
  {
    array[i] = (double*)malloc(2 * sizeof(double));
  }
  return array;
}

//Free the memory of the contact matrix
void free_memory_contact_matrix(double** contact_matrix)
{
  for(int i=0;i<2;i++)
  {
    free(contact_matrix[i]);
  }
  free(contact_matrix);
}

void printf_3D_array(int*** pop,int b,int c)
{
    FILE* output_file;
    output_file = fopen("pop_file.dat","w");
        for(int j=0;j<b;j++)
        {
            for(int k=0;k<c;k++)
            {
                fprintf(output_file,"%d\t%d\t%d\n",pop[k][j][0],pop[k][j][1],pop[k][j][2]);
            }
        }
        fclose(output_file);
}

//Function to allocate memory in a 1D array
double* allocate_memory_1D(int N)
{
    double* array = (double*) malloc(N * sizeof(double));

    return array;
}

double binomial(int n, int k) {
  return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1));
}

bool* initialite_bool_array(int V)
{
  bool* array = (bool*) malloc(V * sizeof(bool));

  return array;
}

//void printf_4D_array()

