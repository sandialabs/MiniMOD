#include "earlycoll.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void all_to_all_rma_init(MPIX_Request *req, int size, int* recvbuf, int chunk_size, MPI_Comm comm) {
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;
    MPI_Win_create(req->recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win);
}

void all_to_all_rma_init_datacopy(MPIX_Request *req, int size, int* recvbuf, int chunk_size, MPI_Comm comm) {
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;
    req->_recvbuf = malloc(sizeof(int)*size);
    MPI_Win_create(req->_recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win);
    MPI_Win_fence(0, req->win);
}

void all_to_all_rma_direct(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    MPI_Win_fence(0, req->win);

    // Direct RMA communication with proper offsets
    for (int i = 0; i < total_processes; i++) {
        MPI_Put(&sendbuf[i * chunk_size], chunk_size, MPI_INT, i, req->rank * chunk_size, chunk_size, MPI_INT, req->win);
    }
}

void all_to_all_rma_pairwise(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    MPI_Win_fence(0, req->win);

    // Pairwise exchange with proper offsets
    for (int step = 0; step < total_processes; step++) {
        int target = (rank + step) % total_processes;
        MPI_Put(&sendbuf[target * chunk_size], chunk_size, MPI_INT, target, rank * chunk_size, chunk_size, MPI_INT, req->win);
    }
}

void all_to_all_rma_datacopy(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    // Pairwise exchange with proper offsets
    for (int step = 0; step < total_processes; step++) {
        int target = (rank + step) % total_processes;
        MPI_Put(&sendbuf[target * chunk_size], chunk_size, MPI_INT, target, rank * chunk_size, chunk_size, MPI_INT, req->win);
    }
}


// NOT IN USE, CONFIRM ALLGATHER FIRST
void all_to_all_rma_init_hierarchical(MPIX_Request *req, int size, int* recvbuf, int chunk_size, MPI_Comm comm) {
    // In the hierarchical rma, there must be an internal recvbuf per leader (per node)
    // TODO: within init, these are what the request for hierarchical needs:
    // PPN, Leader, Leader Window, Leader recvbufs
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;
    req->_recvbuf = (char *)malloc(size * sizeof(int));
    MPI_Win_create(req->_recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win);
}


// NOT IN USE, CONFIRM ALLGATHER FIRST
// -- logic for data positioning needs reconsideration (where to place)
void all_to_all_rma_hierarchical(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int size, int total_processes, int processes_per_node) {
    // Hierarchical communication, 1 leader per node

    int local_rank = rank % processes_per_node;
    int leader_local_rank = 0; // defaulting to leader of 0 for now TODO: either dynamic or smart location.
    int leader_global_rank = (rank / processes_per_node) * processes_per_node; // Leader rank for node
    
    // send from local proc to node leader. TODO: DO WE NEED FENCE PER NODE? -- FENCE OR SYNCH PER LEADER COMM (no leader comm, we use offsets just think of this as the 'connector')
    MPI_Win_fence(0, req->win);
    MPI_Put(&sendbuf[local_rank * chunk_size], chunk_size, MPI_INT, leader_global_rank, local_rank * chunk_size, chunk_size, MPI_INT, req->win);

    // leader direct communication between leaders
    int num_leaders = total_processes / processes_per_node;

    // TODO: Swivel position of data here or at final junction.
    // TODO: Use node win instead?
    MPI_Win_fence(0, req->win);
    for (int target_leader = 0; target_leader < num_leaders; target_leader++) {
        if (target_leader != leader_global_rank / processes_per_node) { // avoid self
            MPI_Put(req->_recvbuf, total_processes, MPI_INT, target_leader * processes_per_node, leader_global_rank, total_processes, MPI_INT, req->win);
        }
    }

    // TODO: Which win to use? -- make a win for recvbuf then?
    MPI_Win_fence(0, req->win);

    // leader put from _recvbufs to recvbuf, trigger when leader is ready (move on from leader win fence?)
    // TODO: Swivel position of data here or at leader junction.
    if (rank == leader_global_rank){
        for (int local_ranks = 0; local_ranks < processes_per_node; local_ranks++) {
            MPI_Put(req->_recvbuf, size, MPI_INT, local_ranks/*get global rank from local rank*/, 0, size, MPI_INT, req->win);
        }
    }
}

void all_to_all_rma_wait_hierarchical(MPIX_Request *req) {
    MPI_Win_fence(0, req->win);
}

void all_to_all_rma_wait(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);
}

void all_to_all_rma_wait_datacopy(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);
    memcpy(req->recvbuf, req->_recvbuf, req->size*sizeof(int));
    MPI_Win_fence(0, req->win);
}




void all_to_all_rma_free(MPIX_Request *req) {
    MPI_Win_free(&req->win);
}

void MPIX_Alltoall_init(int *sendbuf, int chunk_size, int *recvbuf, MPI_Comm comm, MPI_Info info, MPIX_Request *request) {
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &request->rank);

    char algorithm[256];
    int flag;
    MPI_Info_get(info, "algorithm", 256, algorithm, &flag);

    request->sendbuf = sendbuf;
    request->world_size = size;

    request->operation_func = all_to_all_rma_direct;
    request->wait_func = all_to_all_rma_wait;
    request->free_func = all_to_all_rma_free;


    if (flag) {
        if (strcmp(algorithm, "direct") == 0) {
            all_to_all_rma_init(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = all_to_all_rma_direct;
        } else if (strcmp(algorithm, "pairwise") == 0) {
            all_to_all_rma_init(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = all_to_all_rma_pairwise;
        } else if (strcmp(algorithm, "datacopy") == 0) {
            all_to_all_rma_init_datacopy(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = all_to_all_rma_datacopy;
	        request->wait_func = all_to_all_rma_wait_datacopy;
        /* alltoall hierarchical not implemented yet
        } else if (strcmp(algorithm, "hierarchical") == 0) {
            all_to_all_rma_init_hierarchical(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = all_to_all_rma_hierarchical;
            request->wait_func = all_to_all_rma_wait_hierarchical; */
        } else {
            fprintf(stderr, "Unknown algorithm: %s\n", algorithm);
            MPI_Abort(comm, 1);
        }
    } 
}

