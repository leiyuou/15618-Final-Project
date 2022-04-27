#include "shortest.h"

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

void Dijkstra_MPI_core(int** graph, int source, int n_nodes, int *dist,
 int *prev, int procID, int nproc, int *global_min) {

    int span = (n_nodes + nproc - 1) / nproc;
    int startIndex = procID * span;
    int endIndex = min(n_nodes, startIndex + span);
    std::vector<distPair>::iterator it;
    std::vector<distPair> new_pairs;
    int visited[n_nodes];
    
    int local_min[2];
    
    for (int n = 0; n < n_nodes; n++) {
        int local_min_dist = 1000*n_nodes;
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
                prev[i] = global_min[1];
            }
        }
    }
    
}

void Dijkstra_MPI(int **graph, int source, int n_nodes, int *dist, int *prev, int nproc)
{
    double startTime;
    double endTime;
    int procID;
    std::priority_queue<distPair, std::vector<distPair>,compare> pq;
    for(int i=0; i<n_nodes;i++){

        if(i==source){
            dist[i] = 0;
            pq.push({i, 0});
        }

        else
            dist[i] = 1000*n_nodes;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    MPI_Init(NULL, NULL);

    // Run computation
    distPair minPair = pq.top();
    int global_min[2];
    global_min[0] = minPair.second;
    global_min[1] = minPair.first;
    startTime = MPI_Wtime();
    Dijkstra_MPI_core(graph, source, n_nodes, dist, prev, procID, nproc, global_min);
    endTime = MPI_Wtime();

    // Cleanup
    MPI_Finalize();
    
    printf("Elapsed time for proc %d: %f\n", procID, endTime - startTime);

}

/**
 * Bellman Ford Algorithm
 */
void BellmanFord(int** graph, int source, int n_nodes, int *dist, int *prev){
    printf("Computing using Bellman Ford's Algorithm ...\n");

    #pragma omp parallel for
    for (int i = 0; i < n_nodes; i++) {
        dist[i] = 1000*n_nodes;
    }
    dist[source] = 0;

    bool done[n_nodes];
    bool global_done;
    # pragma omp parallel shared (prev, graph, dist, done)
    {
        int id = omp_get_thread_num();
        int num = omp_get_num_threads();
        int n = n_nodes/num;
        int start = n*id;
        if(id<n_nodes%num){
            n += 1;
            start += id;
        }else{
            start += n_nodes%num;
        }

        for (int i=0; i<n_nodes; i++){
            done[id] = true;
            for (int u=0; u<n_nodes; u++){
                for(int v = start; v<start+n; v++){

                    int alt = dist[u] + graph[u][v];

                    if (dist[u]!=1000*n_nodes && alt<dist[v]){
                        dist[v] = alt;
                        prev[v] = u;
                        done[id] = false;
                    }
                }
            }

            #pragma opm barrier
            #pragma omp single 
            {
                global_done = true;
                for(int d=0; d<n_nodes; d++)
                    global_done &= done[d];
            }
            if(global_done)
                break;
        }
    }
}
