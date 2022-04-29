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

/**
 * Main
 */
int main(int argc, char *argv[])
{
    using namespace std::chrono;

    char *inputPath = NULL;
    char algorithm = 'd';
    int num_threads = 1;
    int source = 0;
    int end = -1;
    int opt = 0;
    bool with_pq = false;
    bool mpi = false;

    do
    {
        opt = getopt(argc, argv, "f:n:s:e:a:qi");
        switch (opt)
        {
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

    if (inputPath == NULL)
    {
        printf("Usage: %s -f <input filename> [-n <No_threads>] -s <source node>] -e <end node> -a <algorithm>\n", argv[0]);
        return 1;
    }

    if(num_threads==0){

        int n_nodes, n_edges;
        int *dist;
        int *prev;
    
        if (algorithm == 'd')
        {
            /* DONE Read input */
            int **graph = readInput(inputPath, &n_nodes, &n_edges);

            // init_time += duration_cast<dsec>(Clock::now() - init_start).count();
            // printf("Initialization Time: %lf.\n", init_time);

            /* DONE Compute */
            dist = (int *)malloc(sizeof(int) * n_nodes);
            prev = (int *)malloc(sizeof(int) * n_nodes);

            auto compute_start = Clock::now();
            double compute_time = 0;

            printf("Computing using Sequential Dijkstra's Algorithm...\n");
            Dijkstra_seq(graph, source, n_nodes, dist, prev);

            compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
            printf("***Computation Time: %lf.\n", compute_time);
        }else{
            /* DONE Read input */
            Graph *graph = readInput_BellmanFord(inputPath);

            // init_time += duration_cast<dsec>(Clock::now() - init_start).count();
            // printf("Initialization Time: %lf.\n", init_time);

            /* DONE Compute */
            n_nodes = graph->n_nodes;
            n_edges = graph->n_edges;
            dist = (int *)malloc(sizeof(int) * n_nodes);
            prev = (int *)malloc(sizeof(int) * n_nodes);

            auto compute_start = Clock::now();
            double compute_time = 0;

            printf("Computing using Sequential Bellman Ford's Algorithm...\n");
            BellmanFord_seq(graph, source, dist, prev);

            compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
            printf("***Computation Time: %lf.\n", compute_time);
        }

        /* DONE Write output */
        writeOutput(inputPath, n_nodes, n_edges, dist, prev, num_threads, algorithm, source, end, mpi);
        free(dist);
        free(prev);
        
    }
    else if (mpi)
    {
        if (algorithm == 'd')
        {
            printf("Computing using MPI Dijkstra's Algorithm...\n");
            Dijkstra_MPI(source, num_threads, inputPath, end);
        }
        else
        {
            /* DONE Read input */
            // init_time = 0;
            Graph *graph = readInput_BellmanFord(inputPath);
            // init_time += duration_cast<dsec>(Clock::now() - init_start).count();
            // printf("Initialization Time: %lf.\n", init_time);
            BellmanFord_MPI(graph, source, inputPath, algorithm, end);
        }
    }
    else
    {
        // printf("Number of threads: %d\n", num_threads);
        // printf("Source node: %d\n", source);
        omp_set_num_threads(num_threads);
        int n_nodes, n_edges;
        int *dist;
        int *prev;
    
        if (algorithm == 'd')
        {   
            /* DONE Read input */
            int **graph = readInput(inputPath, &n_nodes, &n_edges);
            // init_time += duration_cast<dsec>(Clock::now() - init_start).count();
            // printf("Initialization Time: %lf.\n", init_time);

            /* DONE Compute */
            dist = (int*)malloc(sizeof(int) * n_nodes);
            prev = (int*)malloc(sizeof(int) * n_nodes);

            auto compute_start = Clock::now();
            double compute_time = 0;

            if (with_pq){
                printf("Computing using OpenMP Dijkstra's Algorithm with PQ...\n");
                Dijkstra_pq(graph, source, n_nodes, dist, prev, num_threads);
            }
            else{
                printf("Computing using OpenMP Dijkstra's Algorithm...\n");
                Dijkstra(graph, source, n_nodes, dist, prev, num_threads);
            }
            
            compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
            printf("***Computation Time: %lf.\n", compute_time);
        }
        else{
            /* DONE Read input */
            Graph *graph = readInput_BellmanFord(inputPath);
            // init_time += duration_cast<dsec>(Clock::now() - init_start).count();
            // printf("Initialization Time: %lf.\n", init_time);
            
            n_nodes = graph->n_nodes;
            n_edges = graph->n_edges;
            dist = (int*)malloc(sizeof(int) * n_nodes);
            prev = (int*)malloc(sizeof(int) * n_nodes);
            
            printf("Computing using OpenMP Bellman Ford's Algorithm...\n");
            auto compute_start = Clock::now();
            double compute_time = 0;
            BellmanFord(graph, source, dist, prev);

            compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
            printf("***Computation Time: %lf.\n", compute_time);
        }

        /* DONE Write output */
        writeOutput(inputPath, n_nodes, n_edges, dist, prev, num_threads, algorithm, source, end, mpi);
    }

    return 0;
}
