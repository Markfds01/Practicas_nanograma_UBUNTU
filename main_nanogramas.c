#include "main_nanogramas.h"
#define BIG
#define ALL

int main()
{
    int number_of_conditions;
    long idum=0;//-time(NULL);
    int n_changes = 0;
    int n_iters_termalization = 5000;
    int n_iters_optimization = 1000000;
    FILE* output_nanogram;
    FILE* energy_file;
    energy_file = fopen("energy_file.csv","w");
    fprintf(energy_file,"time_step,energy\n");
    output_nanogram = fopen("output_nanogram.txt","w");
    //Chose between the big and the small nanogram
    #ifdef BIG
    int n_rows = 21;
    int n_columns = 21;
    char rows_file[256] = "big_rows.txt";
    char columns_file[256] = "big_columns.txt";
    struct Graph* rows_graph= readGraphfile(rows_file,n_rows);
    struct Graph* columns_graph= readGraphfile(columns_file,n_columns);
    #endif
    #ifdef SMALL
    int n_rows = 12;
    int n_columns = 10;
    char rows_file[256] = "small_rows.txt";
    char columns_file[256] = "small_columns.txt";
    struct Graph* rows_graph= readGraphfile(rows_file,n_rows);
    struct Graph* columns_graph= readGraphfile(columns_file,n_columns);
    #endif
    printGraph(columns_graph);
    if(check_conditions(rows_graph,columns_graph,&number_of_conditions))
    {
        printf("Condiciones comprobadas y guardadas, empieza la optimizacion\n");
        int n_clusters = 0;
        //Initialite the nanogram as you like, all in a row for example, then we will inspect other initial configurations
       // int** nanogram = initialite_nanogram_beggining(n_rows,n_columns,number_of_conditions);
        int **nanogram = initialite_random_nanogram(n_rows,n_columns,number_of_conditions,&idum);
       int contador2 = number_of_ones(nanogram,n_rows,n_columns);
        //print_nanogram(nanogram,n_rows,n_columns);
        //print_nanogram(nanogram,n_rows,n_columns);
        double enrgy_trial = calcula_energia_una_columna(nanogram,7,n_rows,columns_graph->array+7);
        //Compute the energy of the initial configuration
        #ifdef ONLY_ROWS
            double total_energy = total_energy_rows(nanogram,rows_graph,n_columns);
        #endif
        #ifdef ALL
        double total_energy = total_energy_nanogram(nanogram,rows_graph,n_columns,n_rows,columns_graph);
        //total_energy /= (n_columns * n_rows);
        #endif

        struct ones_positions* positions_of_ones = assign_ones_positions(nanogram,number_of_conditions,n_rows,n_columns);
        if(check_ones_positions(nanogram,positions_of_ones,number_of_conditions))
        {
            printf("debug");
        }

        //Set the temperature and its loop parameters, because the temperature should decrease
        double beta = 0.2;
        double delta_beta = 0.3;
        double beta_final = 0.7;
        int beta_iteration = (beta_final - beta) / delta_beta;

        int aceptancia = 0;
        int contador;
        /*for(int i=0;i<10;i++)
        {
            propose_change_per_row(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia);
        total_energy = total_energy_rows(nanogram,rows_graph,n_columns);
        contador = number_of_ones(nanogram,n_rows,n_columns);
        }*/
        
        //Set the number of changes per iteration, at the beggining is better to apply more changes in order to arrive the minimum energy in fewer steps
        /*n_changes = sqrt(sqrt(n_rows * n_columns));
        for(int k=0;k<beta_iteration;k++)
        {
            int aceptancia = 0;
            for(int j=0;j<n_iters_termalization;j++)
            {   
                propose_change(&nanogram,n_rows,n_columns,&idum,&total_energy,n_changes,rows_graph,beta,&aceptancia);
                //propose_change_near_one(&nanogram,n_rows,n_columns,&idum,&total_energy,n_changes,rows_graph,beta,&aceptancia);
                //print_nanogram(nanogram,n_rows,n_columns);
            }

            printf("\n total energy es %lf \t aceptancia es %lf \t n_cambios es %d\n ",total_energy,(double)aceptancia/100000.,n_changes);
            //print_nanogram(nanogram,n_rows,n_columns);
            beta = beta + delta_beta;
            n_changes-=5;
        }*/

        //Change the configuration and go to little changes in order to find the optimus way
        #ifdef ONLY_ROWS
        beta = 5.0 ;
        n_changes = 1;
        aceptancia = 0;
            for(int j=0;j<n_iters_optimization;j++)
            {   
                //propose_change(&nanogram,n_rows,n_columns,&idum,&total_energy,n_changes,rows_graph,beta,&aceptancia);
                //propose_change_near_one(&nanogram,n_rows,n_columns,&idum,&total_energy,n_changes,rows_graph,beta,&aceptancia);
                propose_change_per_row(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia);
                //propose_change_per_position(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia,columns_graph);
                contador = number_of_ones(nanogram,n_rows,n_columns);
                if(total_energy == 0)
                {
                    double energy = total_energy_rows(nanogram,rows_graph,n_columns);
                    printf("\n total energy es %lf total energy true es: %lf\t aceptancia es %lf \t n_cambios es %d\n ",total_energy,energy, (double)aceptancia/(double)j,n_changes);
                    print_nanogram(nanogram,n_rows,n_columns);
                    fprint_nanogram(nanogram,n_rows,n_columns,output_nanogram);
                    fclose(output_nanogram);
                    return 0;
                }
            }

         printf("\n total energy es %lf total energy true es: %lf\t aceptancia es %lf \t n_cambios es %d\n ",total_energy,total_energy_nanogram(nanogram,rows_graph,n_columns,n_rows,columns_graph), (double)aceptancia/n_iters_optimization,n_changes);

        #endif
        #ifdef ALL
           
            n_changes = 3;
            aceptancia = 0;
            for(int i=0;i<100;i++)
            {
                 beta = 4 ;
                for(int j=0;j<n_iters_optimization;j++)
                {   
                    int minimum_energy = total_energy;
                propose_change_per_position_optimized(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia,columns_graph,number_of_conditions,positions_of_ones);
               
               // positions_of_ones = assign_ones_positions(nanogram,number_of_conditions,n_rows,n_columns);
                //propose_change_per_position(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia,columns_graph);
                //contador = number_of_ones(nanogram,n_rows,n_columns);
                int minimum_energy_contador = 0;
                
                if(j%5==0)
                {
                   // fprintf(energy_file,"%d,%d\n",j+i*n_iters_optimization,(int)total_energy);
                }
                if(total_energy == 0)
                {
                    double energy = total_energy_nanogram(nanogram,rows_graph,n_columns,n_rows,columns_graph);
                    printf("\n total energy es %lf total energy true es: %lf\t aceptancia es %lf \t n_cambios es %d\n ",total_energy,energy, (double)aceptancia/(double)j,n_changes);
                    print_nanogram(nanogram,n_rows,n_columns);
                    fprint_nanogram(nanogram,n_rows,n_columns,output_nanogram);
                    fclose(output_nanogram);
		    free_memory_2D(nanogram,n_rows);
                    free_graph(rows_graph);
                    free_graph(columns_graph);
                    free(positions_of_ones);
                    return 0;
                }
                }
                printf("\n OPTIMIZATION total energy es %lf total energy true es: %lf\t aceptancia es %lf \t n_cambios es %d\n ",total_energy,total_energy_nanogram(nanogram,rows_graph,n_columns,n_rows,columns_graph), (double)aceptancia/n_iters_optimization,n_changes);
                beta = 1.2;
                aceptancia = 0;
                for(int j=0;j<n_iters_termalization;j++)
                {
                    //propose_change(&nanogram,n_rows,n_columns,&idum,&total_energy,n_changes,rows_graph,beta,&aceptancia,columns_graph);
                    //propose_change_per_position(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia,columns_graph);
                   propose_change_per_position_optimized(nanogram,n_rows,n_columns,&idum,&total_energy,rows_graph,beta,&aceptancia,columns_graph,number_of_conditions,positions_of_ones);
                   if(j%5==0)
                    {
                   // fprintf(energy_file,"%d,%d\n",j+(i+1)*n_iters_optimization,(int)total_energy);
                    }
                     if(total_energy == 0)
                    {
                    double energy = total_energy_nanogram(nanogram,rows_graph,n_columns,n_rows,columns_graph);
                    printf("\n total energy es %lf total energy true es: %lf\t aceptancia es %lf \t n_cambios es %d\n ",total_energy,energy, (double)aceptancia/(double)j,n_changes);
                    print_nanogram(nanogram,n_rows,n_columns);
                    fprint_nanogram(nanogram,n_rows,n_columns,output_nanogram);
                    fclose(output_nanogram);
                    return 0;
                    }
                }

            printf("\n TERMALIZATION total energy es %lf total energy true es: %lf\t aceptancia es %lf \t n_cambios es %d\n ",total_energy,total_energy_nanogram(nanogram,rows_graph,n_columns,n_rows,columns_graph), (double)aceptancia/n_iters_termalization,n_changes);
            printf("\n \n La energia es : %lf \n \n ",total_energy);
            }
            
            print_nanogram(nanogram,n_rows,n_columns);
        #endif
        //Print the nanogram
        
        free_memory_2D(nanogram,n_rows);
        free_graph(rows_graph);
        free_graph(columns_graph);
	free(positions_of_ones);
    }
}
