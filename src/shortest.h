/**
 * 15618 Project shortest.h
 * 
 * Xinqi Wang, Yuou Lei
 */
#ifndef _SHORTEST_H
#define _SHORTEST_H

#include <bits/stdc++.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "mpi.h"
#include <omp.h>

#include <math.h>



typedef std::pair<int,int> distPair;

typedef struct{
    int dist;
    int prev;
} prevPair;

struct compare  
 {  
   bool operator()(const distPair& d1, const distPair& d2)  
   {  
       return d1.second > d2.second;  
   }  
 }; 
static inline int min(int x, int y) { return (x < y) ? x : y; }

void writeOutput(char *inputPath, int n_nodes, int n_edges, int *dist, int *prev, 
                int num_threads, char algorithm, int source, int end, bool mpi);
/**
 * Dijkstra's Algorithms with Priority Queue
 */ 
 void Dijkstra_pq(int** graph, int source, int n_nodes, int *dist, int *prev, int num_threads);
 /**
 * Dijkstra's Algorithms
 */ 
 void Dijkstra(int** graph, int source, int n_nodes, int *dist, int *prev, int num_of_threads);
 /**
 * Bellman Ford Algorithm
 */
 void BellmanFord(int** graph, int source, int n_nodes, int *dist, int *prev);

//  void minPair(void *inB, void *inoutB, int *len, MPI_Datatype *dptr);

/**
 * Bellman Ford Algorithm MPI
 */
// void BellmanFord_MPI(int id, int num, int* graph, int source, int n_nodes, prevPair *dist_prev);

/**
 * @brief Dijkstra's Algo with MPI
 * 
 */
void Dijkstra_MPI_core(int** graph, int source, int n_nodes, int *dist, int *prev, int procID, int nproc, int *global_min);

void Dijkstra_MPI(int **graph, int source, int n_nodes, int *dist, int *prev, int nproc, char *inputPath, int n_edges, int end);

#endif