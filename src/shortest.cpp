/**
 * 15618 Project shortest.cpp
 *
 * Xinqi Wang, Yuou Lei
 */

#include "shortest.h"
using namespace std::chrono;

/**
 * Dijkstra's Algorithms with Priority Queue OpenMP
 */
void Dijkstra_pq(int** graph, int source, int n_nodes, int *dist, int *prev, int num_threads)
{
    // printf("Computing using Dijkstra's Algorithm with Priority Queue...\n");

    std::priority_queue<distPair, std::vector<distPair>,compare> pq;
    for(int i=0; i<n_nodes;i++){

        if(i==source){
            dist[i] = 0;
            pq.push({i, 0});
            prev[i] = source;
        }
        else{
            dist[i] = 1000*n_nodes;
            prev[i] = -1;
        }
    }

    std::vector<distPair>::iterator it;
    distPair p;
    int u, v, alt, id, num;
    std::vector<distPair> new_pairs;

    # pragma omp parallel private (num, p, v, alt, new_pairs, id) shared (pq, u, graph, dist)
    {
        id = omp_get_thread_num();
        num = omp_get_num_threads();
        for(int n=0; n<n_nodes; n++){
            // printf("----%d----%d----\n", id, n);
            #pragma omp single
            {
                p = pq.top(); 
                u = p.first;
                pq.pop();
                // size = graph[u].size();
            }
            #pragma omp barrier
            
            // # pragma omp parallel
            for(int i=id; i<n_nodes ; i+=num){
                v = i;
                alt = dist[u] + graph[u][v];
                // printf("%d u=%d v=%d alt=%d+%d=%d dist[v]=%d\n", 
                //         id, u, v, dist[u], graph[u][i].second, alt, dist[v]);
                if (alt<dist[v]){
                    dist[v] = alt;
                    prev[v] = u;
                    new_pairs.push_back({v,alt});
                }
            }
            #pragma omp critical
            {
                it = new_pairs.begin();
                while (it != new_pairs.end())
                {
                    // printf("%d added (%d, %d)\n", id, it->first, it->second);
                    pq.push(*it);
                    it = new_pairs.erase(it);
                }            
            }
            #pragma omp barrier
        }
    }
}

/**
 * Dijkstra's Algorithm OpenMP
 */ 
void Dijkstra(int** graph, int source, int n_nodes, int *dist, int *prev, int num_of_threads)
{
    // printf("Computing using Dijkstra's Algorithm...\n");
    // omp_set_num_threads(num_of_threads);
    // printf("asdasdasdasfasdasd\n");
    // std::priority_queue<distPair, std::vector<distPair>, compare> pq;
    int* visited = (int*)calloc(n_nodes, sizeof(int));
    int* mindist = (int*)malloc(sizeof(int)*num_of_threads);
    int* minInd = (int*)malloc(sizeof(int)*num_of_threads);

    int finished = 0;
    // int id, numThreads;
    int global_min = 1000*n_nodes;
    int u = 0;

    #pragma omp parallel shared(finished, source, dist, visited, global_min, u, mindist, minInd, graph, n_nodes)
    {
        int id = omp_get_thread_num();
        int numThreads = omp_get_num_threads();
    
        for(int i=id; i<n_nodes; i+=numThreads){
            dist[i] = 1000*n_nodes;
            prev[i] = -1;
        }

        #pragma omp single
        {
            dist[source] = 0;
            prev[source] = source;
        }

        for(; finished<n_nodes;)
        {
            int minDis = 1000*n_nodes;
            int minNode = 0;
            
            // start = (id * n_nodes) / numThreads;
            // end = ((id + 1) * n_nodes) / numThreads;

            for (int i = id; i < n_nodes; i+=numThreads)
            {
                if (visited[i] == 0)
                {
                    if (dist[i] < minDis)
                    {
                        minDis = dist[i];
                        minNode = i;
                    }
                }
            }
            // printf("%d %d minDis=%d minInd=%d\n", id, finished, minDis, minNode);
            mindist[id] = minDis;
            minInd[id] = minNode;
            
            #pragma omp single
            {
                global_min = 1000*n_nodes;
                u = 0;
                for(int t=0; t<numThreads; t++){
                    if(mindist[t]<global_min){
                        global_min = mindist[t];
                        u = minInd[t];
                    }
                }
                visited[u] = 1;
            }
            #pragma omp barrier

            for (int v = id; v<n_nodes; v+=num_of_threads)
            {
                // int node = it->first;
                if (visited[v] == 0)
                {
                    int alt = dist[u] + graph[u][v];
                    // printf("%d u=%d v=%d alt=%d+%d=%d dist[v]=%d\n",
                    //        id, u, v, dist[u], graph[u][v], alt, dist[v]);
                    if (alt < dist[v])
                    {
                        dist[v] = dist[u] + graph[u][v];
                        prev[v] = u;
                    }
                }
            }
            #pragma omp barrier
            #pragma omp single
            finished++;
        }
    }

    free(visited);
}

/**
 * Dijkstra's Algorithm MPI
 * 
 */
void Dijkstra_MPI_core(int** graph, int source, int n_nodes, int *dist,
 int *prev, int procID, int nproc, int *global_min, int startIndex, int endIndex) {


    std::vector<distPair>::iterator it;
    std::vector<distPair> new_pairs;
    int visited[n_nodes];
    for (int i = 0; i < n_nodes; i++) {
        visited[i] = 0;
    }
    
    int local_min[2];
    
    for (int n = 0; n < n_nodes; n++) {
        int local_min_dist = 1000*n_nodes+1;
        int local_min_node = 0;
        for(int i=startIndex; i<endIndex ; i++){
            if (visited[i] == 1) continue;

            if (dist[i] < local_min_dist)
            {
                local_min_dist = dist[i];
                local_min_node = i;
            }
        }
        local_min[0] = local_min_dist;
        local_min[1] = local_min_node;
        
        MPI_Allreduce(local_min, global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        visited[global_min[1]] = 1;
        dist[global_min[1]] = global_min[0];

        for(int i=startIndex; i<endIndex ; i++){
            if (visited[i] == 1) continue;
            int alt = dist[global_min[1]] + graph[global_min[1]][i];
            // printf("%d u=%d v=%d alt=%d+%d=%d dist[v]=%d\n",
            //        id, u, v, dist[u], graph[u][v], alt, dist[v]);
            if (alt < dist[i])
            {
                dist[i] = alt;
                prev[i - startIndex] = global_min[1];   
            }
        }
    }
    
}

void Dijkstra_MPI(int source, int nproc, char *inputPath, int end)
{

    double startTime;
    double endTime;
    int procID;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int n_nodes, n_edges;
    int *dist = (int *)malloc(sizeof(int) * n_nodes);
    int *prev = (int *)malloc(sizeof(int) * n_nodes);
    int **graph = readInput(inputPath, &n_nodes, &n_edges);
    int global_min[2];
    if (procID == 0) {
        for(int i=0; i<n_nodes;i++){
            if(i==source){
                dist[i] = 0;
            }
            else{
                dist[i] = 1000*n_nodes;
            }
        }
    }

    // Run computation

    MPI_Bcast(dist, n_nodes, MPI_INT, 0, MPI_COMM_WORLD);
    int span = (n_nodes + nproc - 1) / nproc;
    int startIndex = procID * span;
    int endIndex = min(n_nodes, startIndex + span);
    int local_prev[span];
    for (int i = 0; i < span; i++) {
        local_prev[i] = 0;
    }
    startTime = MPI_Wtime();
    Dijkstra_MPI_core(graph, source, n_nodes, dist, local_prev, procID, nproc, global_min, startIndex, endIndex);
    MPI_Gather(local_prev, span, MPI_INT, prev, span, MPI_INT, 0, MPI_COMM_WORLD);
    endTime = MPI_Wtime();
    if (procID == 0) {
        writeOutput(inputPath, n_nodes, n_edges, dist, prev, nproc, 'd', source, end, true);       
    }
    // MPI_Reduce(local_prev, prev, n_nodes, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Cleanup
    MPI_Finalize();
    printf("Elapsed time for proc %d: %f\n", procID, endTime - startTime);
    free(dist);
    free(prev);
    free(graph);

}

/**
 * Sequential Dijkstra's Algorithm
 * 
 */
void Dijkstra_seq(int **graph, int source, int n_nodes, int *dist, int *prev)
{
    // printf("Computing using Sequential Dijkstra's Algorithm...\n");

    int *visited = (int *)calloc(n_nodes, sizeof(int));

    for (int i = 0; i < n_nodes; i++){
        dist[i] = 1000 * n_nodes;
        prev[i] = -1;
    }

    dist[source] = 0;
    prev[source] = source;

    for (int finished = 0; finished < n_nodes; finished++){
        int minDis = 1000 * n_nodes;
        int u = 0;

        for (int i = 0; i < n_nodes; i ++){
            if (visited[i] == 0 && dist[i] < minDis)
            {
                minDis = dist[i];
                u = i;
            }
        }
        visited[u] = 1;

        for (int v = 0; v < n_nodes; v++){
            if (visited[v] == 0)
            {
                int alt = dist[u] + graph[u][v];
                if (alt < dist[v])
                {
                    dist[v] = dist[u] + graph[u][v];
                    prev[v] = u;
                }
            }
        }
    }

    free(visited);
}

/**
 * Sequential Bellman Ford's Algorithm
 * 
 */
void BellmanFord_seq(Graph *graph, int source, int *dist, int *prev)
{

    Edge* edges = graph->edges;
    int *weights = graph->weights;
    int n_nodes = graph->n_nodes;
    int n_edges = graph->n_edges;

    for (int i = 0; i < n_nodes; i++)
    {
        dist[i] = 1000 * n_nodes;
        prev[i] = -1;
    }
    dist[source] = 0;
    prev[source] = source;
    
    for (int i=0; i<n_nodes; i++){
        for (int j=0; j<n_edges; j++){
            int u = edges[j].u;
            int v = edges[j].v;
            int alt = dist[u] + weights[j];
            if(dist[u]!=1000*n_nodes && alt<dist[v]){
                dist[v] = alt;
                prev[v] = u;
            }
        }
    }
}

/**
 * Bellman Ford Algorithm
 *
 */
void BellmanFord(Graph* graph, int source, int *dist, int *prev)
{
    printf("Computing using Bellman Ford's Algorithm ...\n");

    Edge* edges = graph->edges;
    int *weights = graph->weights;
    int n_nodes = graph->n_nodes;
    int n_edges = graph->n_edges;

    #pragma omp parallel for
    for (int i = 0; i < n_nodes; i++)
    {
        dist[i] = 1000 * n_nodes;
        prev[i] = -1;
    }
    dist[source] = 0;
    prev[source] = source;

    for (int i=0; i<n_nodes; i++){
        #pragma omp parallel for
        for (int j=0; j<n_edges; j++){
            int u = edges[j].u;
            int v = edges[j].v;
            int alt = dist[u] + weights[j];
            if(dist[u]!=1000*n_nodes && alt<dist[v]){
                dist[v] = alt;
                prev[v] = u;
            }
        }
    }
}

/**
 * Bellman Ford Algorithm MPI
 */
void minPair(void *inB, void *inoutB, int *len, MPI_Datatype *dptr)
{
    prevPair *in = (prevPair *)inB;
    prevPair *inout = (prevPair *)inoutB;
    int i;
    for (i = 0; i < *len; i++)
    {
        if (in->dist < inout->dist)
        {
            inout->dist = in->dist;
            inout->prev = in->prev;
        }
        in++;
        inout++;
    }
}

void BellmanFord_MPI_core(int id, int num, Graph *graph, int source, prevPair *dist_prev)
{
    Edge* edges = graph->edges;
    int *weights = graph->weights;
    int n_nodes = graph->n_nodes;
    int n_edges = graph->n_edges;

    MPI_Op minRed;
    MPI_Datatype MPI_PAIR;
    MPI_Type_contiguous(2, MPI_INT, &MPI_PAIR);
    MPI_Type_commit(&MPI_PAIR);
    MPI_Op_create(&minPair, 1, &minRed);

    // assign work
    int n = n_edges / num;
    int start = n * id;
    if (id < n_edges % num)
    {
        n += 1;
        start += id;
    }
    else
    {
        start += n_edges % num;
    }

    // initialize distance
    for (int i = 0; i < n_nodes; i++)
    {
        dist_prev[i].dist = 1000 * n_nodes;
        dist_prev[i].prev = -1;
    }
    dist_prev[source].dist = 0;
    dist_prev[source].prev = source;

    MPI_Barrier(MPI_COMM_WORLD);

    // bool done;
    for (int i = 0; i < n_nodes; i++)
    {
        // done = true;

        for (int e = start; e < start + n; e++)
        {
            // printf("%d (%d, %d, %d) dist[%d]=%d prev[%d]=%d \n", id, u, v, graph[u * n_nodes + v],
            //                 v, dist_prev[v].dist, v, dist_prev[v].prev);
            int u = edges[e].u;
            int v = edges[e].v;
            int alt = dist_prev[u].dist + weights[e];

            if (alt < dist_prev[v].dist)
            {
                dist_prev[v].dist = alt;
                dist_prev[v].prev = u;
                // done = false;
            }
        }
        // MPI_Allreduce(MPI_IN_PLACE, &done, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        // if (done)
        //     break;
        // printf("%d reduce minimum\n", id);
        MPI_Allreduce(MPI_IN_PLACE, dist_prev, n_nodes, MPI_PAIR, minRed, MPI_COMM_WORLD);
    }
}

void BellmanFord_MPI(Graph *graph, int source, char* inputPath, char algorithm, int end)
{
    MPI_Init(NULL, NULL);

    int procID;
    int nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if (procID == 0)
    {
        printf("Computing using MPI Bellman Ford's Algorithm...\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int n_nodes = graph->n_nodes;
    int n_edges = graph->n_edges;
    prevPair *dist_prev = (prevPair*)malloc(sizeof(prevPair) * n_nodes);
    /* DONE Compute */
    auto compute_start = Clock::now();
    double compute_time = 0;

    BellmanFord_MPI_core(procID, nproc, graph, source, dist_prev);

    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("***Computation Time: %lf.\n", compute_time);

    /* DONE Write output */
    if (procID == 0)
    {
        int *dist = (int *)malloc(sizeof(int) * n_nodes);
        int *prev = (int *)malloc(sizeof(int) * n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            dist[i] = dist_prev[i].dist;
            prev[i] = dist_prev[i].prev;
        }
        free(dist_prev);
        writeOutput(inputPath, n_nodes, n_edges, dist, prev, nproc, algorithm, source, end, true);
        free(dist);
        free(prev);
    }

    /* DONE Cleanup */
    MPI_Finalize();
}

/**
 * Read Input
 *
 */
int **readInput(char *inputPath, int *n_nodes, int *n_edges)
{
    // printf("Reading input %s...\n", inputPath);
    FILE *input = fopen(inputPath, "r");

    if (!input)
    {
        printf("Unable to open file: %s.\n", inputPath);
        return NULL;
    }

    fscanf(input, "%d\t%d\n", n_nodes, n_edges);
    // printf("n_nodes = %d\tn_edges = %d\n", *n_nodes, *n_edges);

    int **graph = (int **)malloc(sizeof(int *) * (*n_nodes));
    for (int i = 0; i < *n_nodes; i++)
    {
        graph[i] = (int *)malloc(sizeof(int) * (*n_nodes));
        if (graph[i] == NULL)
        {
            printf("Unable to maloc graph[%d]\n", i);
            return NULL;
        }
    }

    for (int i = 0; i < *n_nodes; i++)
    {
        for (int j = 0; j < *n_nodes; j++)
        {
            if (i == j)
                graph[i][j] = 0;
            else
                graph[i][j] = 1000 * (*n_nodes);
        }
    }

    int u, v, w;
    for (int i = 0; i < *n_edges; i++)
    {
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w);
        // printf("(%d, %d) = %d\n", u, v, w);
        graph[u][v] = w;
    }
    return graph;
}

/**
 * Read Input for Bellman Ford
 *
 */
Graph *readInput_BellmanFord(char *inputPath)
{
    // printf("Reading input Bellman Ford%s...\n", inputPath);
    FILE *input = fopen(inputPath, "r");

    if (!input)
    {
        printf("Unable to open file: %s.\n", inputPath);
        return NULL;
    }

    int n_nodes, n_edges;
    fscanf(input, "%d\t%d\n", &n_nodes, &n_edges);
    // printf("n_nodes = %d\tn_edges = %d\n", *n_nodes, *n_edges);

    Edge *edges = (Edge*)malloc(sizeof(Edge)*n_edges);
    int *weights = (int*)malloc(sizeof(int)*n_edges);

    int u, v, w;
    for (int i = 0; i < n_edges; i++)
    {
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w);
        edges[i].u = u;
        edges[i].v = v;
        weights[i] = w;
    }


    Graph *graph = (Graph*)malloc(sizeof(Graph));
    graph->edges = edges;
    graph->weights = weights;
    graph->n_nodes = n_nodes;
    graph->n_edges = n_edges;

    return graph;
}

/**
 * Write output
 *
 */
void writeOutput(char *inputPath, int n_nodes, int n_edges, int *dist, int *prev, 
                int num_threads, char algorithm, int source, int end, bool mpi)
{
    // printf("Writing output...\n");
    char *inputname = (char*)malloc(sizeof(char)*strlen(inputPath));
    strcpy(inputname, inputPath);
    char *pt;
    pt = strtok(inputname, "../");
    for(int i=0; i<1; i++){
        pt = strtok(NULL, "/.");
    }

    char outputname[64];
    if(num_threads==0)
        sprintf(outputname, "../outputs/seq_%s_%d_%c.txt", pt, num_threads, algorithm);
    else if(mpi)
        sprintf(outputname, "../outputs/mpi_%s_%d_%c.txt", pt, num_threads, algorithm);
    else
        sprintf(outputname, "../outputs/openmp_%s_%d_%c.txt", pt, num_threads, algorithm);

    FILE *output = fopen(outputname, "w");
    if (!output) {
        printf("Unable to open file: %s.\n", outputname);
        return;
    }

    fprintf(output, "%d\t%d\n", n_nodes, n_edges);
    fprintf(output, "Node\tDis\tPrev\n");
    for(int n=0; n<n_nodes; n++){
        // printf("%d\t%d\t%d\n", n, dist[n], prev[n]);
        fprintf(output, "%d\t%d\t%d\n", n, dist[n], prev[n]);  
    }
    fclose(output);

    /* TODO outpur shortest path if end node is specified */
    if(end>=0){
        int p = end;
        printf("%d",p);
        while(p!=source){
            p = prev[p];
            printf("->%d", p);
        }
        printf("\n");
    }
}
