/**
 * 15618 Project Main.cpp
 * 
 * Xinqi Wang, Yuou Lei
 */

 
#include "shortest.h"

#include <chrono>
#include <bits/stdc++.h>
# include <cstdio>
# include <cstdlib>
# include <unistd.h>
# include <omp.h>

/**
 * Dijkstra's Algorithms
 */ 
void Dijkstra_pq(int** graph, int source, int n_nodes, int *dist, int *prev, int num_threads)
{
    printf("Computing using Dijkstra's Algorithm with Priority Queue...\n");

    std::priority_queue<distPair, std::vector<distPair>,compare> pq;
    for(int i=0; i<n_nodes;i++){

        if(i==source){
            dist[i] = 0;
            pq.push({i, 0});
        }
        else
            dist[i] = 1000*n_nodes;
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

void Dijkstra(int** graph, int source, int n_nodes, int *dist, int *prev, int num_of_threads)
{
    printf("Computing using Dijkstra's Algorithm...\n");
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
        }

        #pragma omp single
        {
            dist[source] = 0;
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
 * Bellman Ford Algorithm
 */
void BellmanFord(int** graph, int source, int n_nodes, int *dist, int *prev){
    printf("Computing using Bellman Ford's Algorithm ...\n");


    # pragma omp parallel shared (prev, graph, dist)
    {
        int id = omp_get_thread_num();
        int num = omp_get_num_threads();
        for(int i=id; i<n_nodes;i+=num){
            if(i==source)
                dist[i] = 0;
            else
                dist[i] = 1000*n_nodes;
        }

        for (int i=0; i<n_nodes; i++){
            for (int u=0; u<n_nodes; u++){
                for(int v = id; v<n_nodes; v+=num){
                    // int v = it->first;
                    int alt = dist[u] + graph[u][v];

                    if (dist[u]!=1000*n_nodes && alt<dist[v]){
                        dist[v] = alt;
                        prev[v] = u;
                    }
                }
            }
        }
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
    
    do {
        opt = getopt(argc, argv, "f:n:s:e:a:");
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

    printf("Number of threads: %d\n", num_threads);
    // printf("Algorithm to be used: %c\n", algorithm);
    printf("Source node: %d\n", source);

    omp_set_num_threads(num_threads);

    /* DONE Read input */
    printf("Reading input %s...\n", inputPath);
    FILE *input = fopen(inputPath, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", inputPath);
        return 1;
    }

    int n_nodes, n_edges;
    fscanf(input, "%d\t%d\n", &n_nodes, &n_edges);
    // std::vector<adjPair> *graph = new std::vector<adjPair>;
    // std::vector<std::vector<adjPair>> graph(n_nodes);
    int **graph = (int**)malloc(sizeof(int*)*n_nodes);
    for (int i=0; i<n_nodes; i++){
        graph[i]=(int*)malloc(sizeof(int)*n_nodes);
        if(graph[i]==NULL){
            printf("Unable to maloc graph[%d]\n", i);
            return 1;
        }
    }
    // printf("1done reading...\n"); 
    
    for(int i=0; i<n_nodes; i++){
        for (int j=0; j<n_nodes; j++){
            if(i==j)
                graph[i][j] = 0;
            else
                graph[i][j] = 1000*n_nodes;
        }
    }
    // printf("2done reading...\n"); 

    int u, v, w;
    for(int i=0; i<n_edges; i++){
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w); 
        // printf("graph[%d][%d] = %d\n", u, v, w);  
        graph[u][v] = w;
    }
    // printf("3done reading...\n"); 

    init_time += duration_cast<dsec>(Clock::now() - init_start).count();
    printf("Initialization Time: %lf.\n", init_time);

    /* DONE Compute */
    // printf("Computing using %c...\n", algorithm);

    int* dist = (int*)malloc(sizeof(int)*n_nodes);
    int* prev = (int*)malloc(sizeof(int)*n_nodes);

    auto compute_start = Clock::now();
    double compute_time = 0;

    if(algorithm=='d')
        Dijkstra(graph, source, n_nodes, dist, prev, num_threads);
    else
        BellmanFord(graph, source, n_nodes, dist, prev);
    
    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("***Computation Time: %lf.\n", compute_time);

    /* DONE Write output */
    // printf("Writing output...\n");
    char *inputname = (char*)malloc(sizeof(char)*strlen(inputPath));
    strcpy(inputname, inputPath);
    char *pt;
    pt = strtok(inputname, "../");
    for(int i=0; i<1; i++){
        
        pt = strtok(NULL, "/.");
    }

    char outputname[64];
    sprintf (outputname, "./outputs/output_%s_%d.txt", pt, num_threads);
    
    FILE *output = fopen(outputname, "w");
    if (!output) {
        printf("Unable to open file: %s.\n", outputname);
        return 1;
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

    return 0;
}

