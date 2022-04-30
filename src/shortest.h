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
#include <chrono>
#include <math.h>


typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::duration<double> dsec;
typedef std::pair<int, int> intPair;

typedef struct 
{
  int n_nodes;
  int n_edges;
  int *n_neighbors;
  int** nodes;
  int** weights;
} AdjList;

typedef struct
{
  int u;
  int v;
} Edge;

typedef struct
{
  Edge *edges;
  int *weights;
  int n_nodes;
  int n_edges;
} Graph;

typedef struct
{
  int dist;
  int prev;
} prevPair;

struct compare
{
  bool operator()(const intPair &d1, const intPair &d2)
  {
    return d1.second > d2.second;
  }
};

static inline int min(int x, int y) { return (x < y) ? x : y; }

/**
 * Dijkstra's Algorithm with Priority Queue OpenMP
 */ 
void Dijkstra_pq(AdjList *graph, int source, int *dist, int *prev, int num_threads);

/**
 * Dijkstra's Algorithm OpenMP
 */ 
void Dijkstra(AdjList *graph, int source, int *dist, int *prev, int num_of_threads);

/**
 * Sequential Dijkstra's Algorithm
 */
void Dijkstra_seq(AdjList *graph, int source, int *dist, int *prev);

/**
 * Dijkstra's Algo with MPI
 * 
 */
void Dijkstra_MPI(int source, int nproc, char *inputPath, int end);

/**
 * Sequential Bellman Ford Algorithm
 *
 */
void BellmanFord_seq(Graph *graph, int source, int *dist, int *prev);

/**
 * Bellman Ford Algorithm OpenMP
 *
 */
void BellmanFord(Graph* graph, int source, int *dist, int *prev);

/**
 * Bellman Ford Algorithm MPI
 *
 */
void BellmanFord_MPI(Graph *graph, int source, char* inputPath, char algorithm, int end);

/**
 * Read Input
 *
 */
AdjList *readInput(char *inputPath);

/**
 * Read Input for Bellman Ford
 *
 */
Graph *readInput_BellmanFord(char *inputPath);

/**
 * Write output
 *
 */
void writeOutput(char *inputPath, int n_nodes, int n_edges, int *dist, int *prev, 
                int num_threads, char algorithm, int source, int end, bool mpi);

#endif