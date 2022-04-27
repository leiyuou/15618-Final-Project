/**
 * 15618 Project Main.cpp
 * 
 * Xinqi Wang, Yuou Lei
 */

#include "shortest.h"
#include "mpi.h"
#include <chrono>
#include <bits/stdc++.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <omp.h>

int **readInput(char *inputPath, int* n_nodes, int *n_edges)
{
    printf("Reading input %s...\n", inputPath);
    FILE *input = fopen(inputPath, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", inputPath);
        return NULL;
    }

    fscanf(input, "%d\t%d\n", n_nodes, n_edges);
    int **graph = (int**)malloc(sizeof(int*)*(*n_nodes));
    for (int i=0; i<*n_nodes; i++){
        graph[i]=(int*)malloc(sizeof(int)*(*n_nodes));
        if(graph[i]==NULL){
            printf("Unable to maloc graph[%d]\n", i);
            return NULL;
        }
    }
    // printf("1done reading...\n"); 
    
    for(int i=0; i<*n_nodes; i++){
        for (int j=0; j<*n_nodes; j++){
            if(i==j)
                graph[i][j] = 0;
            else
                graph[i][j] = 1000*(*n_nodes);
        }
    }
    // printf("2done reading...\n"); 

    int u, v, w;
    for(int i=0; i<*n_edges; i++){
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w); 
        // printf("graph[%d][%d] = %d\n", u, v, w);  
        graph[u][v] = w;
    }
    // printf("3done reading...\n"); 
    return graph;
}

int *readInput_MPI(char *inputPath, int* n_nodes, int *n_edges)
{
    printf("Reading input %s...\n", inputPath);
    FILE *input = fopen(inputPath, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", inputPath);
        return NULL;
    }

    fscanf(input, "%d\t%d\n", n_nodes, n_edges);
    int *graph = (int*)malloc(sizeof(int*)*(*n_nodes)*(*n_nodes));
    
    for(int i=0; i<*n_nodes; i++){
        for (int j=0; j<*n_nodes; j++){
            if(i==j)
                graph[i*(*n_nodes) + j] = 0;
            else
                graph[i*(*n_nodes) + j] = 1000*(*n_nodes);
        }
    }
    int u, v, w;
    for(int i=0; i<*n_edges; i++){
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w);  
        graph[u*(*n_nodes) + v] = w;
    }

    return graph;
}

void writeOutput(char *inputPath, int n_nodes, int n_edges, int *dist, int *prev, 
                int num_threads, char algorithm, int source, int end, bool mpi)
{
    printf("Writing output...\n");
    char *inputname = (char*)malloc(sizeof(char)*strlen(inputPath));
    strcpy(inputname, inputPath);
    char *pt;
    pt = strtok(inputname, "../");
    for(int i=0; i<1; i++){
        pt = strtok(NULL, "/.");
    }

    char outputname[64];
    if(mpi)
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


void minPair (void *inB, void *inoutB, int *len, MPI_Datatype *dptr)
{
    prevPair *in = (prevPair*)inB;
    prevPair *inout = (prevPair*)inoutB;
    int i;
    for(i=0; i<*len; i++){
        if(in->dist<inout->dist){
            inout->dist = in->dist;
            inout->prev = in->prev;
        }
        in++;
        inout++;
    }
}

/**
 * Bellman Ford Algorithm MPI
 */
void BellmanFord_MPI(int id, int num, int* graph, int source, int n_nodes, prevPair *dist_prev){
    printf("Computing using Bellman Ford's Algorithm ...\n");

    MPI_Op minRed; 
    MPI_Datatype MPI_PAIR;
    MPI_Type_contiguous( 2, MPI_INT, &MPI_PAIR); 
    MPI_Type_commit(&MPI_PAIR);
    MPI_Op_create(&minPair, 1, &minRed);

    // broadcast number of nodes and graph
    MPI_Bcast(&n_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(id!=0){
        dist_prev = (prevPair*)malloc(sizeof(prevPair) * n_nodes);
        graph = (int*)malloc(sizeof(int*)*(n_nodes*n_nodes));
    }
    MPI_Bcast(&graph, n_nodes*n_nodes, MPI_INT, 0, MPI_COMM_WORLD);

    // assign work
    int n = n_nodes/num;
    int start = n*id;
    if(id<n_nodes%num){
        n += 1;
        start += id;
    }else{
        start += n_nodes%num;
    }

    // initialize distance
    for (int i = 0; i < n_nodes; i++){
        dist_prev[i].dist = 1000*n_nodes;
        dist_prev[i].prev = -1;
    }
    dist_prev[source].dist = 0;
    dist_prev[source].prev = source;

    MPI_Barrier(MPI_COMM_WORLD);

    bool done;
    for (int i = 0; i < n_nodes; i++)
    {
        done = true;
        for (int u = 0; u < n_nodes; u++)
        {
            for (int v = start; v < start + n; v++)
            {

                int alt = dist_prev[u].dist + graph[u*n_nodes+v];

                if (dist_prev[u].dist != 1000 * n_nodes && alt < dist_prev[v].dist)
                {
                    dist_prev[v].dist = alt;
                    dist_prev[v].prev = u;
                    done = false;
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &done, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        if (done)
            break;

        MPI_Allreduce(MPI_IN_PLACE, &dist_prev, n, MPI_PAIR, minRed, MPI_COMM_WORLD);

    }

    if(id!=0){
        free(dist_prev);
        free(graph);
    }
}

/**
 * Main 
 */ 
int main(int argc, char *argv[])
{
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;

    auto init_start = Clock::now();
    double init_time = 0;

    char *inputPath = NULL;
    char algorithm = 'd';
    int num_threads=1;
    int source=0;
    int end = -1;
    int opt=0;
    bool with_pq = false;
    bool mpi = false;
    
    do {
        opt = getopt(argc, argv, "f:n:s:e:a:q");
        switch (opt) {
        case 'f':
            inputPath = optarg;
            break;
        case 'n':
            num_threads = atoi(optarg);
            break;
        case 's':
            source = atoi(optarg);
            break;
        case 'e':
            end = atoi(optarg);
            break;
        case 'a':
            algorithm = optarg[0];
            break;
        case 'q':
            with_pq = true;
            break;
        case 'i':
            mpi = true;
            break;
        case -1:
            break;
        default:
            break;
        }
    } while (opt != -1);

    if (inputPath == NULL) {
        printf("Usage: %s -f <input filename> [-n <No_threads>] -s <source node>] -e <end node> -a <algorithm>\n", argv[0]);
        return 1;
    }

    if (mpi){
        int procID;
        int nproc;
        MPI_Comm_rank(MPI_COMM_WORLD, &procID);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);

        MPI_Init(NULL,NULL);

        /* DONE Read input */
        int n_nodes, n_edges;
        int *graph;
        prevPair *dist_prev; 
        
        if(procID==0){
            graph = readInput_MPI(inputPath, &n_nodes, &n_edges);
            init_time += duration_cast<dsec>(Clock::now() - init_start).count();
            printf("Initialization Time: %lf.\n", init_time);
            dist_prev = (prevPair*)malloc(sizeof(prevPair) * n_nodes);
        }

        /* DONE Compute */
        // printf("Computing using %c...\n", algorithm);
        auto compute_start = Clock::now();
        double compute_time = 0;

        BellmanFord_MPI(procID, nproc, graph, source, n_nodes, dist_prev);

        compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
        printf("***Computation Time: %lf.\n", compute_time);

        /* DONE Write output */
        if (procID==0){
            int *dist = (int *)malloc(sizeof(int) * n_nodes);
            int *prev = (int *)malloc(sizeof(int) * n_nodes);
            for(int i=0; i<n_nodes; i++){
                dist[i] = dist_prev[i].dist;
                prev[i] = dist_prev[i].prev;
            }
            free(dist_prev);
            writeOutput(inputPath, n_nodes, n_edges, dist, prev, num_threads, algorithm, source, end, mpi); 
        }

        /* DONE Cleanup */
        MPI_Finalize();
    }
    else {
        printf("Number of threads: %d\n", num_threads);
        printf("Source node: %d\n", source);

        /* DONE Read input */
        int n_nodes, n_edges;
        int **graph = readInput(inputPath, &n_nodes, &n_edges);

        init_time += duration_cast<dsec>(Clock::now() - init_start).count();
        printf("Initialization Time: %lf.\n", init_time);

        /* DONE Compute */
        // printf("Computing using %c...\n", algorithm);

        int *dist = (int *)malloc(sizeof(int) * n_nodes);
        int *prev = (int *)malloc(sizeof(int) * n_nodes);

        auto compute_start = Clock::now();
        double compute_time = 0;

        if (algorithm == 'd')
        {
            if (with_pq)
                Dijkstra_pq(graph, source, n_nodes, dist, prev, num_threads);
            else
                Dijkstra(graph, source, n_nodes, dist, prev, num_threads);
        }
        else
            BellmanFord(graph, source, n_nodes, dist, prev);

        compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
        printf("***Computation Time: %lf.\n", compute_time);

        /* DONE Write output */
        writeOutput(inputPath, n_nodes, n_edges, dist, prev, num_threads, algorithm, source, end, mpi);       
    }

    return 0;
}

