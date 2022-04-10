/**
 * 15618 Project Main.cpp
 * 
 * Xinqi Wang, Yuou Lei
 */

 
#include "shortest.h"

#include <limits>
#include <bits/stdc++.h>
# include <cstdio>
# include <cstdlib>
# include <unistd.h>
# include <omp.h>


void Dijkstra(std::list<adjPair> *graph, int source, int n_nodes, int *dist, int *prev)
{

    std::priority_queue<distPair, std::vector<distPair>,compare> pq;
    for(int i=0; i<n_nodes;i++){

        if(i==source){
            dist[i] = 0;
            pq.push({i, 0});
        }
            
        else
            dist[i] = 1000*n_nodes;
    }

    std::list<distPair>::iterator it;
    while(!pq.empty()){
        distPair p = pq.top(); 
        int u = p.first;
        pq.pop();

        for(it=graph[u].begin();it!=graph[u].end();it++){
            int v = it->first;
            int alt = dist[u] + it->second;

            // printf("u=%d v=%d alt=%d+%d=%d dist[v]=%d\n", 
                    // u, v, dist[u], it->second, alt, dist[v]);
            if (alt<dist[v]){
                dist[v] = alt;
                prev[v] = u;
                pq.push({v, alt});
            }
        }
    }

}

int main(int argc, char *argv[])
{

    char *inputPath = NULL;
    int num_threads=1;
    int source=0;
    int opt=0;
    
    do {
        opt = getopt(argc, argv, "f:n:s:");
        switch (opt) {
        case 'f':
            inputPath = optarg;
            break;
        case 'n':
            num_threads = atoi(optarg);
            break;
        case 's':
            source = atoi(optarg);
        case -1:
            break;
        default:
            break;
        }
    } while (opt != -1);

    if (inputPath == NULL) {
        printf("Usage: %s -f <input filename> [-n <No_threads>] -s <source node>]\n", argv[0]);
        return 1;
    }

    printf("Number of threads: %d\n", num_threads);
    printf("Input file: %s\n", inputPath);

    /* DONE Read input */
    printf("Reading input %s...\n", inputPath);
    FILE *input = fopen(inputPath, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", inputPath);
        return 1;
    }

    int n_nodes, n_edges;
    fscanf(input, "%d\t%d\n", &n_nodes, &n_edges);
    std::list<adjPair> *graph = new std::list<adjPair>[n_nodes];
    
    int u, v, w;
    for(int i=0; i<n_edges; i++){
        fscanf(input, "%d\t%d\t%d\n", &u, &v, &w);        
        // printf("(%d\t%d)\t%d\n", u, v, w);
        graph[u].push_back({v, w});
    }

    /* DONE Compute */
    printf("Computing..\n");

    int* dist = (int*)malloc(sizeof(int)*n_nodes);
    int* prev = (int*)malloc(sizeof(int)*n_nodes);

    Dijkstra(graph, source, n_nodes, dist, prev);

    /* DONE Write output */
    printf("Writing output...\n");
    char *inputname = (char*)malloc(sizeof(char)*strlen(inputPath));
    strcpy(inputname, inputPath);
    char *pt;
    pt = strtok(inputname, "./");
    for(int i=0; i<1; i++){
        
        pt = strtok(NULL, "/.");
    }

    char outputname[64];
    sprintf (outputname, "outputs/output_%s_%d.txt", pt, num_threads);
    
    FILE *output = fopen(outputname, "w");
    if (!output) {
        printf("Unable to open file: %s.\n", outputname);
        return 1;
    }

    fprintf(output, "%d\t%d\n", n_nodes, n_edges);
    fprintf(output, "Node\tDis\tPrev\n");
    for(int n=0; n<n_nodes; n++){
        printf("%d\t%d\t%d\n", n, dist[n], prev[n]);
        fprintf(output, "%d\t%d\t%d\n", n, dist[n], prev[n]); 
    }
    fclose(output);

    return 0;
}

