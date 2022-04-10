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
void Dijkstra(std::vector<adjPair> *graph, int source, int n_nodes, int *dist, int *prev)
{
    printf("Computing using Dijkstra's Algorithm ...\n");

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
    int u, v, alt, size, id, num;
    std::vector<distPair> new_pairs;

    # pragma omp parallel private (num, p, v, alt, new_pairs, id) shared (size, pq, u, graph, dist)
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
                size = graph[u].size();
            }
            #pragma omp barrier
            
            // # pragma omp parallel
            for(int i=id; i<size ; i+=num){
                // int i = ind + id;
                v = graph[u][i].first;
                alt = dist[u] + graph[u][i].second;
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
 * Bellman Ford Algorithm
 */
void BellmanFord(std::vector<adjPair> *graph, int source, int n_nodes, int *dist, int *prev){
    printf("Computing using Bellman Ford's Algorithm ...\n");

    // for(int i=0; i<n_nodes;i++){
    //     if(i==source)
    //         dist[i] = 0;
    //     else
    //         dist[i] = 1000*n_nodes;
    // }

    // std::list<distPair>::iterator it;
    // for (int i=0; i<n_nodes; i++){
    //     for (int u=0; u<n_nodes; u++){

    //         for(it=graph[u].begin();it!=graph[u].end();it++){
    //             int v = it->first;
    //             int alt = dist[u] + it->second;

    //             if (dist[u]!=1000*n_nodes && alt<dist[v]){
    //                 dist[v] = alt;
    //                 prev[v] = u;
    //             }
    //         }
    //     }
    // }
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
    std::vector<adjPair> *graph = new std::vector<adjPair>[n_nodes];
    
    int u, v, w;
    for(int i=0; i<n_edges; i++){
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w);        
        // printf("(%d\t%d)\t%d\n", u, v, w);
        graph[u].push_back({v, w});
    }

    init_time += duration_cast<dsec>(Clock::now() - init_start).count();
    printf("Initialization Time: %lf.\n", init_time);

    /* DONE Compute */
    // printf("Computing using %c...\n", algorithm);

    int* dist = (int*)malloc(sizeof(int)*n_nodes);
    int* prev = (int*)malloc(sizeof(int)*n_nodes);

    auto compute_start = Clock::now();
    double compute_time = 0;

    if(algorithm=='d')
        Dijkstra(graph, source, n_nodes, dist, prev);
    else
        BellmanFord(graph, source, n_nodes, dist, prev);
    
    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("Computation Time: %lf.\n", compute_time);

    /* DONE Write output */
    printf("Writing output...\n");
    char *inputname = (char*)malloc(sizeof(char)*strlen(inputPath));
    strcpy(inputname, inputPath);
    char *pt;
    pt = strtok(inputname, "../");
    for(int i=0; i<1; i++){
        
        pt = strtok(NULL, "/.");
    }

    char outputname[64];
    sprintf (outputname, "../outputs/output_%s_%d.txt", pt, num_threads);
    
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

