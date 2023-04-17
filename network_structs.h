#ifndef NETWORKS_STRUCTS_HEADER
#define NETWORKS_STRUCTS_HEADER
#include "stdlib.h"
#include "stdio.h"
#include "random_number.h"
#include <stdbool.h>
#include "string.h"
//#include "DFS.h"
//#include "network_types.h"

// A structure to represent an adjacency list node

struct ones_positions{
    int row_pos;
    int column_pos;
};

struct AdjListNode {
    int dest;
    struct AdjListNode* next;
};
 
// A structure to represent an adjacency list
struct AdjList {
    struct AdjListNode* head;
    int number_restrictions;
};

struct ClusterList{
    struct AdjListNode* head;
    int number_of_clusters;
};
 
// A structure to represent a graph. A graph
// is an array of adjacency lists.
// Size of array will be V (number of vertices
// in graph)
struct Graph {
    int V;
    struct AdjList* array;
};

struct Queue {
    int front, rear, size;
    unsigned capacity;
    int* array;
};

//Funciones
struct AdjListNode* newAdjListNode(int dest);
struct Graph* createGraph(int V);
void addEdge(struct Graph* graph, int src, int dest);
void printGraph(struct Graph* graph);
int *create_vector_grades(int V);//Crea un vector con los grados
int *create_vector_to_randomize(int* vector_grades,struct Graph* graph);
int tamagno_vector_randomize(int *vector_grades,int V);
int* calculate_new_vector_grades(struct Graph *graph);
void free_graph(struct Graph* graph);
void free_vector_grades(int* vector_grades);
bool busca_si_esta(int* vector,int N,int sospechoso);
void vector_a_cero(int* vector,int N);
int how_many_visited_nodes(int* vector,int N);
int deep_first_search(struct Graph *graph, int* vector_grades,long* idum);
int visit_order(int* vector,int N,int index);
bool check_loops(struct Graph* graph);
bool check_repeated_links(struct Graph* graph,int k_max);
struct Queue* createQueue(unsigned capacity);
int isFull(struct Queue* queue);
int isEmpty(struct Queue* queue);
void enqueue(struct Queue* queue, int item);
int dequeue(struct Queue* queue);
int* nNeighbours(int node, struct Graph* graph, int* arr, int n, int V);
void fprintGraph(struct Graph* graph,FILE *output_file);
void addRestriction(struct Graph* graph, int src, int dest);
struct Graph* readGraphfile(char* filename,int n);
void addCluster(struct ClusterList* Cluster_Row,  int dest);
struct ClusterList* newAdjClusterList();





























#endif
