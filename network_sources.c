#include "network_structs.h"

// A utility function to create a new adjacency list node
struct AdjListNode* newAdjListNode(int dest)
{
    struct AdjListNode* newNode
        = (struct AdjListNode*)malloc(
            sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

struct ClusterList* newAdjClusterList()
{
    struct ClusterList* newClusterList = (struct ClusterList*) malloc(sizeof(struct ClusterList));

    newClusterList->head = NULL;
    newClusterList->number_of_clusters = 0;
    return newClusterList; 
}
 
// A utility function that creates a graph of V vertices
struct Graph* createGraph(int V)
{
    struct Graph* graph
        = (struct Graph*)malloc(sizeof(struct Graph));
    graph->V = V;
 
    // Create an array of adjacency lists.  Size of
    // array will be V
    graph->array = (struct AdjList*)malloc(
        V * sizeof(struct AdjList));
 
    // Initialize each adjacency list as empty by
    // making head as NULL
    int i;
    for (i = 0; i < V; ++i){
        graph->array[i].head = NULL;
        graph->array[i].number_restrictions = 0;
    }

    return graph;
}

//Free alocate memory of the graph
void free_graph(struct Graph* graph)
{
      struct AdjListNode *ptr, *ptr2;
  for(int i = 0; i < graph->V; i++) {
    ptr = graph->array[i].head;
    while(ptr != NULL) {
      ptr2 = ptr->next;
      free(ptr);

      ptr = ptr2;
    }
    //free(graph->array[i].head);
  }
  free(graph->array);
    free(graph);
}

//Function to initialite the grades vector which contains the grade of each node
int* create_vector_grades(int V)
{
    int *vector_grades = (int*) malloc(V * sizeof(int));
    return vector_grades;
}

//Free the memory of the vector_grades
void free_vector_grades(int* vector_grades)
{
    free(vector_grades);
}

//Function to create the vector where we´re going to put each node head as many times as its grade
int* create_vector_to_randomize(int* vector_grades,struct Graph* graph)
{
    int size = tamagno_vector_randomize(vector_grades,graph->V);
    int *vector_randomize = (int*) malloc(size * sizeof(int));
    int k=0;
    for(int i=0;i<graph->V;i++)
    {
        //Cuantas veces se tiene que repetir el nodo
        int repetition = vector_grades[i];
        for(int j=0;j<repetition;j++)
        {
            //Se asigna el valor del nodo las veces necesarias
            vector_randomize[k] = i;
            k++;
            //printf("\n %d",vector_randomize[k-1]);
        }
    }
    return vector_randomize;
}

//Function to calculate the size of the vector of the function above
int tamagno_vector_randomize(int *vector_grades,int V)
{
    int i;
    int res = 0;
    for(i=0;i<V;i++)
    {
        res+=vector_grades[i];
    }
    return res;
}

int* calculate_new_vector_grades(struct Graph *graph)
{
    int *vector_grades = create_vector_grades(graph->V);//Creo el vector que guarda los grados de cada nodo
    for(int i=0;i<graph->V;i++)
    {
        double grado_indice = 0;//Grado de cada indice
        struct AdjListNode* pCrawl = graph->array[i].head;//Igualo variable a la cabeza de cada nodo
        while(pCrawl)
        {
            grado_indice++;
            pCrawl = pCrawl->next;
        }
        vector_grades[i] = grado_indice;
    }
    return vector_grades;
}
// Adds an edge to an undirected graph
void addEdge(struct Graph* graph, int src, int dest)
{
    // Add an edge from src to dest.  A new node is
    // added to the adjacency list of src.  The node
    // is added at the beginning
    struct AdjListNode* check = NULL;
    struct AdjListNode* newNode = newAdjListNode(dest);
 

   if (src == dest) { // check if is a loop 
        //printf("No se permiten loops, no se agrega el enlace\n");
        return;
    }
    
    check = graph->array[src].head;
    while (check != NULL) { // check if already exists
        if (check->dest == dest) {
            //printf("Enlace ya existe, no se agrega nuevamente\n");
            return;
        }
        check = check->next;
    }

    if (graph->array[src].head == NULL) {
        newNode->next = graph->array[src].head;
        graph->array[src].head = newNode;
    }
    else {
 
        check = graph->array[src].head;
        while (check->next != NULL) {
            check = check->next;
        }
        // graph->array[src].head = newNode;
        check->next = newNode;
    }
 
    // Since graph is undirected, add an edge from
    // dest to src also
    newNode = newAdjListNode(src);
    if (graph->array[dest].head == NULL) {
        newNode->next = graph->array[dest].head;
        graph->array[dest].head = newNode;
    }
    else {
        check = graph->array[dest].head;
        while (check->next != NULL) {
            check = check->next;
        }
        check->next = newNode;
    }
 
    // newNode = newAdjListNode(src);
    // newNode->next = graph->array[dest].head;
    // graph->array[dest].head = newNode;
}
 
// A utility function to print the adjacency list
// representation of graph
void printGraph(struct Graph* graph)
{
    int v;
    for (v = 0; v < graph->V; ++v) {
        struct AdjListNode* pCrawl = graph->array[v].head;
        printf("\n Adjacency list of vertex %d\n head ", v);
        while (pCrawl) {
            printf("-> %d", pCrawl->dest);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

void fprintGraph(struct Graph* graph,FILE *output_file)
{
    int v;
    for (v = 0; v < graph->V; ++v) {
        struct AdjListNode* pCrawl = graph->array[v].head;
        fprintf(output_file,"%d", v);
        while (pCrawl) {
            fprintf(output_file,"\t %d", pCrawl->dest);
            pCrawl = pCrawl->next;
        }
        fprintf(output_file,"\n");
    }
}

int deep_first_search(struct Graph *graph, int* vector_grades,long* idum)
{
    //Start node to do the searching
    int start = 0;
    while(vector_grades[start]==0)
    {
        int start = Parisi_Rapuano(idum) * graph->V;
    }
    //Create the vector of visited nodes
    int* vector_nodos_visitados = (int*) malloc(graph->V * sizeof(int));
    int* orden_de_visita = (int*) malloc(graph->V *1000* sizeof(int));
    orden_de_visita[0] = start;
    int k=1;
    vector_a_cero(vector_nodos_visitados,graph->V);
    //Stablish the first node as visited
    vector_nodos_visitados[start] = 1;
    //Defining destino variable
    int destino = 0;
    //Create the auxiliar variable to look
    struct AdjListNode* actual_node = NULL;
    //Stablish the variable as the start node
    actual_node = graph->array[start].head;
    //Variable del vecino anterior
    int old_one = start;
    for(long int i=0;i<10000000;i++)
    {
        for(int j=0;j<vector_grades[start];j++)
        {
            destino = actual_node->dest;
            //Si no está visitado, se avanza hacia ese nodo
            if(vector_nodos_visitados[destino] == 0)
            {
                old_one = start;
                //Ahora se sale desde ese nodo
                start = destino;
                //Se añade a visitado
                vector_nodos_visitados[destino] = 1;
                actual_node = graph->array[start].head;
                j=0;
                orden_de_visita[k] = start;
                k++;
            }
            else{
                actual_node = actual_node->next;
            }
            if(how_many_visited_nodes(vector_nodos_visitados,graph->V)==1)
            {
                printf("\n Hemos recorrido 90 de la red\n");
                return 0;
            }
        }
    start = orden_de_visita[visit_order(orden_de_visita,graph->V,start) - 1];
    //orden_de_visita[k] = start;
    //k++;
    actual_node = graph->array[start].head;
    }
}
//Comprobar si un elemento esta en un vector
bool busca_si_esta(int* vector,int N,int sospechoso)
{
    for(int i=0;i<N;i++)
    {
        if(vector[i]==sospechoso)
        {
            return true;
        }
    }
    return false;
}
//Establecer el vector a 0
void vector_a_cero(int* vector,int N)
{
    for(int i=0;i<N;i++)
    {
        vector[i] = 0;
    }
}

//Function to compute how many neighbours have been visited
int how_many_visited_nodes(int* vector,int N)
{
    int contador = 0;
    for(int i=0;i<N;i++)
    {
        contador += vector[i];
    }
    if(contador>=0.9*N)
    {
        return 1;
    }
    else{
        return 0;
    }
}

//Function to know which is the previous node
int visit_order(int* vector,int N,int index)
{
    for(int i=0;i<N;i++)
    {
        if(vector[i]==index)
        {
            return i;
        }
    }
}

//Function to look if there are any loops
bool check_loops(struct Graph* graph)
{
    //Iniciamos el nodo auxiliar
    struct AdjListNode* pCrawl = NULL;
    //Recorremos todos los nodos
    for(int i=0;i<graph->V;i++)
    {
        //Asignamos el puntero al nodo
        pCrawl = graph->array[i].head;
        while(pCrawl)
        {
            //Comprobamos que no haya loops
            if(pCrawl->dest == i)
            {
                return false;
            }
            pCrawl = pCrawl->next;
        }
    }
    return true;
}

//Function to check if there are any repeated links()
bool check_repeated_links(struct Graph* graph,int k_max)
{
    //Lista con los posibles vecinos
    int* vecinos = (int*) malloc(k_max * sizeof(int));
    for(int i=0;i<k_max;i++)
    {
        vecinos[i] = -1;
    }
    struct AdjListNode* pCrawl = NULL;
    int k=0;
    for(int i=0;i<graph->V;i++)
    {
        pCrawl = graph->array[i].head;
        while(pCrawl)
        {
            //Comprobamos si esta ya en vecinos
            if(busca_si_esta(vecinos,k_max,pCrawl->dest))
            {
                return false;
            }
            //Lo añadimos a la lista de vecinos checkeada
            vecinos[k] = pCrawl->dest;
            k++;
            //Pasamos al siguiente vecino
            pCrawl = pCrawl->next;
        }
        k=0;
        for(int i=0;i<k_max;i++)
        {
            //Restablecemos la lista de vecinos
            vecinos[i] = -1;
        }
    }
}

// Function to create a queue with given capacity.
// Initializes size of queue as 0.
struct Queue* createQueue(unsigned capacity)
{
    struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));
    queue->capacity = capacity;
    queue->front = queue->size = 0;
    queue->rear = capacity - 1;  // This is important, see the enqueue
    queue->array = (int*)malloc(queue->capacity * sizeof(int));
    return queue;
}

// Queue is full when size becomes equal to the capacity
int isFull(struct Queue* queue)
{
    return (queue->size == queue->capacity);
}
 
// Queue is empty when size is 0
int isEmpty(struct Queue* queue)
{
    return (queue->size == 0);
}
 
// Function to add an item to the queue.
// It changes rear and size
void enqueue(struct Queue* queue, int item)
{
    if (isFull(queue))
        return;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->array[queue->rear] = item;
    queue->size = queue->size + 1;
}

int dequeue(struct Queue* queue) {
    if (isEmpty(queue)) {
        return -1;
    }
    int item = queue->array[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}

int* nNeighbours(int node, struct Graph* graph, int* arr, int n, int V) {
    // Initialize all values in arr to -1 (impossible to reach)
    arr = (int*) malloc(V * sizeof(int));
    for(int i=0;i<V;i++)
    {
        arr[i] = -1;
    }
    // Initialize all values in arr to -1 (impossible to reach)
    //memset(arr, -1, V * sizeof(int));
    
    // Set the seed node's distance to 0
    arr[node] = 0;
    
    // Create a queue for BFS and enqueue the seed node
    struct Queue* queue = createQueue(V);
    enqueue(queue, node);
 
    // Traverse the graph using BFS
    while (!isEmpty(queue)) {
        // Dequeue a vertex from queue
        int currNode = dequeue(queue);
 
        // Find all adjacent vertices of the dequeued vertex currNode
        struct AdjListNode* adjList = graph->array[currNode].head;
        while (adjList != NULL) {
            int adjNode = adjList->dest;
            // If adjNode is not visited before, enqueue it and mark it visited
            if (arr[adjNode] == -1) {
                arr[adjNode] = arr[currNode] + 1;
                enqueue(queue, adjNode);
            }
            adjList = adjList->next;
        }
    }
    
    // Set all values greater than n to -1 (impossible to reach)
    for (int i = 0; i < V; ++i) {
        if (arr[i] > n) {
            arr[i] = -1;
        }
    }
    return arr;
}

void addRestriction(struct Graph* graph, int src, int dest)
{
    // Add an edge from src to dest.  A new node is
    // added to the adjacency list of src.  The node
    // is added at the beginning
    struct AdjListNode* check = NULL;
    struct AdjListNode* newNode = newAdjListNode(dest);
    
    check = graph->array[src].head;

    if (graph->array[src].head == NULL) {
        newNode->next = graph->array[src].head;
        graph->array[src].head = newNode;
    }
    else {
 
        check = graph->array[src].head;
        while (check->next != NULL) {
            check = check->next;
        }
        // graph->array[src].head = newNode;
        check->next = newNode;
    }
}

void addCluster(struct ClusterList* Cluster_Row,  int dest)
{
    // Add an edge from src to dest.  A new node is
    // added to the adjacency list of src.  The node
    // is added at the beginning
    struct AdjListNode* check = NULL;
    struct AdjListNode* newNode = newAdjListNode(dest);
    
    check = Cluster_Row->head;

    if (Cluster_Row->head == NULL) {
        newNode->next = Cluster_Row->head;
        Cluster_Row->head = newNode;
    }
    else {
 
        check = Cluster_Row->head;
        while (check->next != NULL) {
            check = check->next;
        }
        // graph->array[src].head = newNode;
        check->next = newNode;
    }
}

struct Graph* readGraphfile(char* filename, int n) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Unable to open file %s\n", filename);
        exit(1);
    }

    // Create the graph
    struct Graph* graph = createGraph(n);

    // Read each line in the file and add the numbers to the corresponding line list
    int line_number = 0;
    int dest;
    while (fscanf(file, "%d", &dest) != EOF) {
        addRestriction(graph, line_number, dest);
        graph->array[line_number].number_restrictions++;
        char c = fgetc(file);
        if (c == '\n') {
            line_number++;
        }
    }

    fclose(file);
    return graph;
}
