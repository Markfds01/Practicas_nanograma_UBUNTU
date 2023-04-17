#include "nanogram_utils.h"

bool check_conditions(struct Graph* rows_restrictions,struct Graph* columns_restricions,int* number_of_conditions)
{
    struct AdjListNode* pCrawl = NULL;
    int rows_conditions = 0;
    int columns_conditions = 0;

    //We look for the rows restrictions
    for(int i=0;i<rows_restrictions->V;i++)
    {
        pCrawl = rows_restrictions->array[i].head;
        while(pCrawl)
        {
            rows_conditions += pCrawl->dest;
            pCrawl = pCrawl->next;
        }
    }
    //We look for the columns restricitons
    for(int i=0;i<columns_restricions->V;i++)
    {
        pCrawl = columns_restricions->array[i].head;
        while(pCrawl)
        {
            columns_conditions += pCrawl->dest;
            pCrawl = pCrawl->next;
        }
    }
    if(rows_conditions == columns_conditions)
    {
        printf("nanograma bien definido \n");
        *number_of_conditions = rows_conditions;
        return true;
    }
    return false;
}
//Initialite the nanogram putting all the possible squares in a continious row
int** initialite_nanogram_beggining(int rows,int columns,int number_of_conditions)
{
    int** nanogram = allocate_memory_2D(rows,columns);

    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<columns;j++)
        {
            if((i*columns+j)<number_of_conditions)
            {
                nanogram[i][j] = 1;
            }
            else{
                nanogram[i][j] = 0;
            }
        }
    }
    return nanogram;
}

void print_nanogram(int** nanogram,int n_rows,int n_columns)
{
    printf("\n");
    for(int i=0;i<n_rows;i++)
    {
        for(int j=0;j<n_columns;j++)
        {
            printf("%d\t",nanogram[i][j]);
        }
        printf("\n");
    }
}

void fprint_nanogram(int** nanogram,int n_rows,int n_columns,FILE* f)
{
    for(int i=0;i<n_rows;i++)
    {
        for(int j=0;j<n_columns;j++)
        {
            fprintf(f,"%d\t",nanogram[i][j]);
        }
        fprintf(f,"\n");
    }
}


int* copy_vector(int* vector_tobe_copied,int n)
{
    int* vector_to_return = (int*) malloc(n * sizeof(int));
    memcpy(vector_to_return, vector_tobe_copied, n * sizeof(int));
    //for(int i=0;i<n;i++)
    //{
    //    vector_to_return[i] = vector_tobe_copied[i];
    //}
    return vector_to_return;
}
struct ClusterList* cluster_the_row(int* row,int n_columns)
{
    int label= 1;
    int tamagno_cluster = 0;
    //Equalize the vector to the row
    int* clustered_row = copy_vector(row,n_columns);
    struct ClusterList* Cluster_row = newAdjClusterList();
    

    //Going through the line
    for(int i=0;i<n_columns;i++)
    {
        //If it is drawed
        if(row[i] == 1)
        {   
            tamagno_cluster++;
            clustered_row[i] = label;
            if(row[i+1] == 0 || i == (n_columns-1))
            {
                label++;
                Cluster_row->number_of_clusters++;
                addCluster(Cluster_row,tamagno_cluster);
                tamagno_cluster = 0;
            }
        }
    }
    free(clustered_row);
    return Cluster_row;
}

void free_ClusterList(struct ClusterList* cluster_list)
{
      struct AdjListNode *ptr, *ptr2;
  
    ptr = cluster_list->head;
    while(ptr != NULL) {
      ptr2 = ptr->next;
      free(ptr);

      ptr = ptr2;
    }
	free(cluster_list);
}

double calcula_energia_una_fila(int *row,int n_columns,struct AdjList* adjlist)
{
    //We clustered the row
    struct ClusterList* clustered_row = cluster_the_row(row,n_columns);

    double energy_per_row = 0;
    //We analyze the case
    if(clustered_row->number_of_clusters < adjlist->number_restrictions)
    {
        //We have that there are less clusters than coditions
        //So we compare the clusters with the conditions we have and add
        //the others conditions

        //We define two pointers
        struct AdjListNode* ptr_cluster;
        struct AdjListNode* ptr_restrictions;

        //We assign the corresponding sites
        ptr_cluster = clustered_row->head;
        ptr_restrictions = adjlist->head;

        while(ptr_cluster)
        {
            //Compute the enegry of the first condition
            energy_per_row = energy_per_row + fabs((ptr_cluster->dest - ptr_restrictions->dest));

            //Go to the next
            ptr_cluster = ptr_cluster->next;
            ptr_restrictions = ptr_restrictions->next;
        }
        
        //Add the final restrictions
        while(ptr_restrictions)
        {
            //Add the energy
            energy_per_row = energy_per_row + ptr_restrictions->dest;

             ptr_restrictions = ptr_restrictions->next;
        }
 
    }
    else{
        if(clustered_row->number_of_clusters > adjlist->number_restrictions)
        {
            //We have that there are more clusters than coditions
            //So we compare the clusters with the conditions we have and add
            //the others clusters

            //We define two pointers
            struct AdjListNode* ptr_cluster;
            struct AdjListNode* ptr_restrictions;

            //We assign the corresponding sites
            ptr_cluster = clustered_row->head;
            ptr_restrictions = adjlist->head;

            while(ptr_restrictions)
            {
                //Compute the enegry of the first condition
                energy_per_row = energy_per_row + fabs((ptr_cluster->dest - ptr_restrictions->dest));

                //Go to the next
                ptr_cluster = ptr_cluster->next;
                ptr_restrictions = ptr_restrictions->next;
            }
        
            //Add the final clusters
            while(ptr_cluster)
            {
                //Add the energy
                energy_per_row = energy_per_row + ptr_cluster->dest;

                ptr_cluster = ptr_cluster->next;
            }
        }
        else{
            //We have that there are the same amount of clusters than conditions

            //We define two pointers
            struct AdjListNode* ptr_cluster;
            struct AdjListNode* ptr_restrictions;

            //We assign the corresponding sites
            ptr_cluster = clustered_row->head;
            ptr_restrictions = adjlist->head;

            while(ptr_restrictions)
            {
                //Compute the enegry of the first condition
                energy_per_row = energy_per_row + fabs((ptr_cluster->dest - ptr_restrictions->dest));

                //Go to the next
                ptr_cluster = ptr_cluster->next;
                ptr_restrictions = ptr_restrictions->next;
            }
        }
    }
    free_ClusterList(clustered_row);
    return energy_per_row;
}

double total_energy_rows(int** nanogram, struct Graph* rows_restrictions,int n_columns)
{
    double total_energy = 0.0;

    for(int i=0;i<rows_restrictions->V;i++)
    {
        total_energy = total_energy + calcula_energia_una_fila(nanogram[i],n_columns,rows_restrictions->array+i);
    }

    return total_energy;
    
}

int** copy_matrix(int** matrix_tobe_copied,int n_rows,int n_columns)
{
    int** matrix =allocate_memory_2D(n_rows,n_columns);
    for(int i=0;i<n_rows;i++)
    {
        for(int j=0;j<n_columns;j++)
        {
            matrix[i][j] = matrix_tobe_copied[i][j];
        }
    }
    return matrix;
}

void propose_change(int*** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,int n_changes,struct Graph* rows_restrictions,double beta, int* aceptancia,struct Graph* columns_restrictions)
{
    int row_to_change;
    int row_to_go;
    int column_to_change;
    int column_to_go;


    int** temp_nanogram = copy_matrix(*nanogram,n_rows,n_columns);


    for(int i=0;i<n_changes;i++)
    {
        //Declare the points to change
        do{
            row_to_change = Parisi_Rapuano(idum) * n_rows;
            column_to_change = Parisi_Rapuano(idum) * n_columns;
        }while(temp_nanogram[row_to_change][column_to_change]!=1);
        
        do{
             row_to_go = Parisi_Rapuano(idum) * n_rows;
            column_to_go = Parisi_Rapuano(idum) * n_columns;
        }while(temp_nanogram[row_to_go][column_to_go] != 0);
       

        //Change the points
        int aux = temp_nanogram[row_to_change][column_to_change];
        temp_nanogram[row_to_change][column_to_change] = temp_nanogram[row_to_go][column_to_go];
        temp_nanogram[row_to_go][column_to_go] = aux;
    }
   

    double new_energy = total_energy_nanogram(temp_nanogram,rows_restrictions,n_columns,n_rows,columns_restrictions);

    double energy_diff = beta*(new_energy - *old_energy);

    //printf("\n new_energy es %lf \n ",new_energy);
    /*if(exp(-energy_diff)<1)
    {
        printf("ojo que se pasa de 100 y es %lf",exp(-energy_diff));
    }*/
    if(exp(-energy_diff)>Parisi_Rapuano(idum))
    {
        //free_memory_2D(nanogram,n_rows);
        *nanogram = copy_matrix(temp_nanogram,n_rows,n_columns);
        *aceptancia = *aceptancia + 1 ;
       // print_nanogram(temp_nanogram,n_rows,n_columns);
        //print_nanogram(*nanogram,n_rows,n_columns);
        free_memory_2D(temp_nanogram,n_rows);
        *old_energy = new_energy;
        return;
    }

    else {
        free_memory_2D(temp_nanogram,n_rows);
        return;
    }
}


bool find_one_near(int**nanogram, int row, int column)
{
    if(nanogram[row][column + 1] == 1 || nanogram[row][column - 1])
    {
        return true;
    }
    else{
        return false;
    }
}
void propose_change_near_one(int*** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,int n_changes,struct Graph* rows_restrictions,double beta, int* aceptancia)
{
    int row_to_change;
    int row_to_go;
    int column_to_change;
    int column_to_go;


    int** temp_nanogram = copy_matrix(*nanogram,n_rows,n_columns);


    for(int i=0;i<n_changes;i++)
    {
        //Declare the points to change
        do{
            row_to_change = Parisi_Rapuano(idum) * n_rows;
            column_to_change = Parisi_Rapuano(idum) * n_columns;
        }while(temp_nanogram[row_to_change][column_to_change]!=1);
        
        do{
             row_to_go = Parisi_Rapuano(idum) * n_rows;
            column_to_go = Parisi_Rapuano(idum) * n_columns;
        }while(temp_nanogram[row_to_go][column_to_go] != 0 && find_one_near(*nanogram, row_to_go, column_to_go));
       

        //Change the points
        int aux = temp_nanogram[row_to_change][column_to_change];
        temp_nanogram[row_to_change][column_to_change] = temp_nanogram[row_to_go][column_to_go];
        temp_nanogram[row_to_go][column_to_go] = aux;
    }
   

    double new_energy = total_energy_rows(temp_nanogram,rows_restrictions,n_columns);

    double energy_diff = beta*(new_energy - *old_energy);
    

    //printf("\n new_energy es %lf \n ",new_energy);

    if(exp(-energy_diff)>Parisi_Rapuano(idum))
    {
        //free_memory_2D(nanogram,n_rows);
        *nanogram = copy_matrix(temp_nanogram,n_rows,n_columns);
        *aceptancia = *aceptancia + 1 ;
       // print_nanogram(temp_nanogram,n_rows,n_columns);
        //print_nanogram(*nanogram,n_rows,n_columns);
        free_memory_2D(temp_nanogram,n_rows);
        *old_energy = new_energy;
        return;
    }

    else {
        free_memory_2D(temp_nanogram,n_rows);
        return;
    }
}

void propose_change_per_row(int** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,struct Graph* rows_restrictions,double beta, int* aceptancia)
{
    int row_to_change;
    int row_to_go;
    int column_to_change;
    int column_to_go;

    do{
            row_to_change = Parisi_Rapuano(idum) * n_rows;
            column_to_change = Parisi_Rapuano(idum) * n_columns;
    }while(nanogram[row_to_change][column_to_change]!=1);

    do{
             row_to_go = Parisi_Rapuano(idum) * n_rows;
            column_to_go = Parisi_Rapuano(idum) * n_columns;
    }while(nanogram[row_to_go][column_to_go] != 0);

    //Compute the energy of the rows before the change
    double old_energy_row_to_change = calcula_energia_una_fila(nanogram[row_to_change],n_columns,rows_restrictions->array+row_to_change);
    double old_energy_row_to_go = calcula_energia_una_fila(nanogram[row_to_go],n_columns,rows_restrictions->array+row_to_go);

    //Create the new vectors to see the new energy
    int* new_row_to_change = copy_vector(nanogram[row_to_change],n_columns);
    int* new_row_to_go = copy_vector(nanogram[row_to_go],n_columns);

    //Introduce the new energy variables
    double new_energy_row_to_change ;
    double new_energy_row_to_go;

    //Introduce the energy variable
    double energy_diff;

    if(row_to_change == row_to_go)
    {
        //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;
        aux  = new_row_to_change[column_to_go];
        new_row_to_change[column_to_go] = new_row_to_go[column_to_change];
        new_row_to_go[column_to_change] = aux;

        //Compute the new energies, only once because as it is the same row only have one contribution
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);

        energy_diff = (new_energy_row_to_change ) - (old_energy_row_to_change );

        energy_diff *= beta;
    }
    else{
         

        //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;

        //Compute the new energies
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);

        energy_diff = (new_energy_row_to_change + new_energy_row_to_go) - (old_energy_row_to_change + old_energy_row_to_go);

        energy_diff *= beta;
    }
   
    if(exp(-energy_diff)>Parisi_Rapuano(idum))
    {
        //free(nanogram[row_to_change]);
        /*double energy_probe = total_energy_rows(nanogram,rows_restrictions,n_columns);
        if(energy_probe!=*old_energy)
        {
            printf("aqui es");
        }*/
        nanogram[row_to_change][column_to_change] = 0;
        *old_energy = *old_energy + energy_diff/beta;
        
        //free(nanogram[row_to_go]);
        nanogram[row_to_go][column_to_go] = 1;
        *aceptancia = *aceptancia + 1 ;
         /*energy_probe = total_energy_rows(nanogram,rows_restrictions,n_columns);
        if(energy_probe!=*old_energy)
        {
            printf("aqui es");
        }*/
        free(new_row_to_change);
        free(new_row_to_go);
        return;
    }
    else{
        free(new_row_to_change);
        free(new_row_to_go);
    }
}

void propose_change_per_position(int** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,struct Graph* rows_restrictions,double beta, int* aceptancia,struct Graph* columns_restrictions)
{
    int row_to_change;
    int row_to_go;
    int column_to_change;
    int column_to_go;

    do{
            row_to_change = Parisi_Rapuano(idum) * n_rows;
            column_to_change = Parisi_Rapuano(idum) * n_columns;
    }while(nanogram[row_to_change][column_to_change]!=1);

    do{
             row_to_go = Parisi_Rapuano(idum) * n_rows;
            column_to_go = Parisi_Rapuano(idum) * n_columns;
    }while(nanogram[row_to_go][column_to_go] != 0);

    //Compute the energy of the rows before the change
    double old_energy_row_to_change = calcula_energia_una_fila(nanogram[row_to_change],n_columns,rows_restrictions->array+row_to_change);
    double old_energy_row_to_go = calcula_energia_una_fila(nanogram[row_to_go],n_columns,rows_restrictions->array+row_to_go);
    double old_energy_column_change = calcula_energia_una_columna(nanogram,column_to_change,n_rows,columns_restrictions->array + column_to_change);
    double old_energy_column_go = calcula_energia_una_columna(nanogram,column_to_go,n_rows,columns_restrictions->array + column_to_go);

    //Create the new vectors to see the new energy
    int* new_row_to_change = copy_vector(nanogram[row_to_change],n_columns);
    int* new_row_to_go = copy_vector(nanogram[row_to_go],n_columns);

    //Introduce the new energy variables
    double new_energy_row_to_change ;
    double new_energy_row_to_go;
    double new_energy_column_to_go;
    double new_energy_column_to_change;

    //Introduce the energy variable
    double energy_diff;

    if(row_to_change == row_to_go )
    {
        //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;
        aux  = new_row_to_change[column_to_go];
        new_row_to_change[column_to_go] = new_row_to_go[column_to_change];
        new_row_to_go[column_to_change] = aux;

        //Compute the new energies, only once because as it is the same row only have one contribution
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
        new_energy_column_to_change = calcula_energia_una_columna_with_a_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0);
        new_energy_column_to_go = calcula_energia_una_columna_with_a_change(nanogram,column_to_go,n_rows,columns_restrictions->array+column_to_go,row_to_go,1);
        

        energy_diff = (new_energy_row_to_change + new_energy_column_to_go + new_energy_column_to_change) - (old_energy_row_to_change + old_energy_column_go + old_energy_column_change);

        energy_diff *= beta;
    }
    else{
        if(column_to_change == column_to_go)
        {
             //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;

        //Compute the new energies
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
        new_energy_column_to_change = calcula_energia_una_columna_with_two_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0,row_to_go,1);

        energy_diff = (new_energy_row_to_change + new_energy_row_to_go  + new_energy_column_to_change) - (old_energy_row_to_change + old_energy_row_to_go + old_energy_column_change  );

        energy_diff *= beta;
        }
        else{
         

        //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;

        //Compute the new energies
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
        new_energy_column_to_change = calcula_energia_una_columna_with_a_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0);
        new_energy_column_to_go = calcula_energia_una_columna_with_a_change(nanogram,column_to_go,n_rows,columns_restrictions->array+column_to_go,row_to_go,1);

        energy_diff = (new_energy_row_to_change + new_energy_row_to_go + new_energy_column_to_change + new_energy_column_to_go) - (old_energy_row_to_change + old_energy_row_to_go + old_energy_column_change + old_energy_column_go );

        energy_diff *= beta;
    }
    }
    
   
    if(exp(-energy_diff)>Parisi_Rapuano(idum))
    {
        //free(nanogram[row_to_change]);
        /*double energy_probe = total_energy_nanogram(nanogram,rows_restrictions,n_columns,n_rows,columns_restrictions);
        if(energy_probe!=*old_energy)
        {
            printf("aqui es");
        }*/
        nanogram[row_to_change][column_to_change] = 0;
        *old_energy = *old_energy + energy_diff/beta;
        
        //free(nanogram[row_to_go]);
        nanogram[row_to_go][column_to_go] = 1;
        *aceptancia = *aceptancia + 1 ;
        //double energy_columns = total_energy_columns(nanogram,n_rows,columns_restrictions);
         //energy_probe = total_energy_nanogram(nanogram,rows_restrictions,n_columns,n_rows,columns_restrictions);
       /* if(energy_probe!=*old_energy)
        {
            printf("aqui es");
        }*/
       // new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
       // new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
       // new_energy_column_to_change = calcula_energia_una_columna_with_a_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0);
       // new_energy_column_to_go = calcula_energia_una_columna_with_a_change(nanogram,column_to_go,n_rows,columns_restrictions->array+column_to_go,row_to_go,1);
        free(new_row_to_change);
        free(new_row_to_go);
        return;
    }
    else{
        free(new_row_to_change);
        free(new_row_to_go);
    }
}
int number_of_ones(int** nanogram,int n_rows,int n_columns)
{
    int contador = 0;
    for(int i=0;i<n_rows;i++)
    {
        for(int j=0;j<n_columns;j++)
        {
            if(nanogram[i][j] == 1)
            {
                contador++;
            }
        }
    }
    return contador;
}

double calcula_energia_una_columna(int **nanogram, int column_index, int n_rows, struct AdjList *adjlist) {
    int *column = malloc(n_rows * sizeof(int));  // create a new array to store the column
    for (int i = 0; i < n_rows; i++) {
        column[i] = nanogram[i][column_index];   // fill the array with the values from the column
    }
    double energy = calcula_energia_una_fila(column, n_rows, adjlist);  // calculate the energy of the column using the existing function
    free(column);  // free the memory allocated for the column array
    return energy;
}

double calcula_energia_una_columna_with_a_change(int **nanogram, int column_index, int n_rows, struct AdjList *adjlist,int row_index,int value) {
    int *column = malloc(n_rows * sizeof(int));  // create a new array to store the column
    for (int i = 0; i < n_rows; i++) {
        column[i] = nanogram[i][column_index];   // fill the array with the values from the column
    }
    column[row_index] = value;
    double energy = calcula_energia_una_fila(column, n_rows, adjlist);  // calculate the energy of the column using the existing function
    free(column);  // free the memory allocated for the column array
    return energy;
}

double calcula_energia_una_columna_with_two_change(int **nanogram, int column_index, int n_rows, struct AdjList *adjlist,int row_index,int value,int row_to_go,int value2) {
    int *column = malloc(n_rows * sizeof(int));  // create a new array to store the column
    for (int i = 0; i < n_rows; i++) {
        column[i] = nanogram[i][column_index];   // fill the array with the values from the column
    }
    column[row_index] = value;
    column[row_to_go] = value2;
    double energy = calcula_energia_una_fila(column, n_rows, adjlist);  // calculate the energy of the column using the existing function
    free(column);  // free the memory allocated for the column array
    return energy;
}
double total_energy_columns(int** nanogram, int n_rows,struct Graph* columns_restrictions,int n_columns)
{
    double energy= 0;
    for(int i=0;i<n_columns;i++)
    {
        energy = energy + calcula_energia_una_columna(nanogram,i,n_rows,columns_restrictions->array + i);
    }
    return energy ;
}
double total_energy_nanogram(int** nanogram, struct Graph* rows_restrictions,int n_columns,int n_rows,struct Graph* columns_restrictions)
{
    double rows_energy = total_energy_rows(nanogram,rows_restrictions,n_columns);
    double columns_energy = total_energy_columns(nanogram,n_rows,columns_restrictions,n_columns);

    return columns_energy + rows_energy;
}

struct ones_positions*  assign_ones_positions(int** nanogram,int number_of_conditions,int n_rows,int n_columns)
{
    //Define the variable
    struct ones_positions* positions = (struct ones_positions*) malloc(number_of_conditions * sizeof(struct ones_positions));

    //Looking for ones
    int k = 0;
    for(int i=0;i<n_rows;i++)
    {
        for(int j=0;j<n_columns;j++)
        {
            if(nanogram[i][j] == 1)
            {
                positions[k].row_pos = i;
                positions[k].column_pos = j;
                k++;
            }
        }
    }

    return positions;
}

void propose_change_per_position_optimized(int** nanogram, int n_rows, int n_columns, long* idum, double* old_energy,struct Graph* rows_restrictions,double beta, int* aceptancia,struct Graph* columns_restrictions,int number_of_conditions,struct ones_positions* positions_of_ones)
{
    int row_to_change;
    int row_to_go;
    int column_to_change;
    int column_to_go;

    int number_of_one = Parisi_Rapuano(idum) * number_of_conditions;

    row_to_change = positions_of_ones[number_of_one].row_pos;
    column_to_change = positions_of_ones[number_of_one].column_pos;

    do{
             row_to_go = Parisi_Rapuano(idum) * n_rows;
            column_to_go = Parisi_Rapuano(idum) * n_columns;
    }while(nanogram[row_to_go][column_to_go] != 0);

    //Compute the energy of the rows before the change
    double old_energy_row_to_change = calcula_energia_una_fila(nanogram[row_to_change],n_columns,rows_restrictions->array+row_to_change);
    double old_energy_row_to_go = calcula_energia_una_fila(nanogram[row_to_go],n_columns,rows_restrictions->array+row_to_go);
    double old_energy_column_change = calcula_energia_una_columna(nanogram,column_to_change,n_rows,columns_restrictions->array + column_to_change);
    double old_energy_column_go = calcula_energia_una_columna(nanogram,column_to_go,n_rows,columns_restrictions->array + column_to_go);

    //Create the new vectors to see the new energy
    int* new_row_to_change = copy_vector(nanogram[row_to_change],n_columns);
    int* new_row_to_go = copy_vector(nanogram[row_to_go],n_columns);

    //Introduce the new energy variables
    double new_energy_row_to_change ;
    double new_energy_row_to_go;
    double new_energy_column_to_go;
    double new_energy_column_to_change;

    //Introduce the energy variable
    double energy_diff;

    if(row_to_change == row_to_go )
    {
        //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;
        aux  = new_row_to_change[column_to_go];
        new_row_to_change[column_to_go] = new_row_to_go[column_to_change];
        new_row_to_go[column_to_change] = aux;

        //Compute the new energies, only once because as it is the same row only have one contribution
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
        new_energy_column_to_change = calcula_energia_una_columna_with_a_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0);
        new_energy_column_to_go = calcula_energia_una_columna_with_a_change(nanogram,column_to_go,n_rows,columns_restrictions->array+column_to_go,row_to_go,1);
        

        energy_diff = (new_energy_row_to_change + new_energy_column_to_go + new_energy_column_to_change) - (old_energy_row_to_change + old_energy_column_go + old_energy_column_change);

        energy_diff *= beta;
    }
    else{
        if(column_to_change == column_to_go)
        {
             //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;

        //Compute the new energies
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
        new_energy_column_to_change = calcula_energia_una_columna_with_two_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0,row_to_go,1);

        energy_diff = (new_energy_row_to_change + new_energy_row_to_go  + new_energy_column_to_change) - (old_energy_row_to_change + old_energy_row_to_go + old_energy_column_change  );

        energy_diff *= beta;
        }
        else{
         

        //Change the values 
        int aux  = new_row_to_change[column_to_change];
        new_row_to_change[column_to_change] = new_row_to_go[column_to_go];
        new_row_to_go[column_to_go] = aux;

        //Compute the new energies
        new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
        new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
        new_energy_column_to_change = calcula_energia_una_columna_with_a_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0);
        new_energy_column_to_go = calcula_energia_una_columna_with_a_change(nanogram,column_to_go,n_rows,columns_restrictions->array+column_to_go,row_to_go,1);

        energy_diff = (new_energy_row_to_change + new_energy_row_to_go + new_energy_column_to_change + new_energy_column_to_go) - (old_energy_row_to_change + old_energy_row_to_go + old_energy_column_change + old_energy_column_go );

        energy_diff *= beta;
    }
    }
    //energy_diff /= (n_columns * n_rows);
    
   
    if(exp(-energy_diff)>Parisi_Rapuano(idum))
    {
        //free(nanogram[row_to_change]);
        /*double energy_probe = total_energy_nanogram(nanogram,rows_restrictions,n_columns,n_rows,columns_restrictions);
        if(energy_probe!=*old_energy)
        {
            printf("aqui es");
        }*/
        nanogram[row_to_change][column_to_change] = 0;
        *old_energy = *old_energy + energy_diff/beta;
        
        //free(nanogram[row_to_go]);
        nanogram[row_to_go][column_to_go] = 1;
        *aceptancia = *aceptancia + 1 ;
        
        positions_of_ones[number_of_one].row_pos = row_to_go;
        positions_of_ones[number_of_one].column_pos = column_to_go;

        //double energy_columns = total_energy_columns(nanogram,n_rows,columns_restrictions);
         //energy_probe = total_energy_nanogram(nanogram,rows_restrictions,n_columns,n_rows,columns_restrictions);
       /* if(energy_probe!=*old_energy)
        {
            printf("aqui es");
        }*/
       // new_energy_row_to_change = calcula_energia_una_fila(new_row_to_change,n_columns,rows_restrictions->array+row_to_change);
       // new_energy_row_to_go = calcula_energia_una_fila(new_row_to_go,n_columns,rows_restrictions->array+row_to_go);
       // new_energy_column_to_change = calcula_energia_una_columna_with_a_change(nanogram,column_to_change,n_rows,columns_restrictions->array+column_to_change,row_to_change,0);
       // new_energy_column_to_go = calcula_energia_una_columna_with_a_change(nanogram,column_to_go,n_rows,columns_restrictions->array+column_to_go,row_to_go,1);
        free(new_row_to_change);
        free(new_row_to_go);
        return;
    }
    else{
        free(new_row_to_change);
        free(new_row_to_go);
    }
}

bool check_ones_positions(int** nanogram, struct ones_positions* positions_of_ones,int number_of_conditions)
{
    int guard = 0;
    for(int i=0;i<number_of_conditions;i++)
    {
        if(nanogram[positions_of_ones[i].row_pos][positions_of_ones[i].column_pos] == 0)
        {
            guard = 1;
        }
    }
    return guard == 1;
}

int** initialite_random_nanogram(int n_rows,int n_columns,int number_of_conditions,long* idum)
{
    int** nanogram = allocate_memory_2D(n_rows,n_columns);
    matriz_a_cero(n_rows,n_columns,nanogram);
    int i=0;
    for(i=0;i<number_of_conditions;i++)
    {
        int row;
        int column;
        do{
            row = Parisi_Rapuano(idum) * n_rows;
            column = Parisi_Rapuano(idum) * n_columns;
        }while(nanogram[row][column] == 1);
        
        nanogram[row][column] = 1;
    }
    return nanogram;
}
