/**
 * 15618 Project shortest.h
 * 
 * Xinqi Wang, Yuou Lei
 */
#include <bits/stdc++.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "mpi.h"
#include <omp.h>

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