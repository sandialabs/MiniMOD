/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

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
        } else {
            fprintf(stderr, "Unknown algorithm: %s\n", algorithm);
            MPI_Abort(comm, 1);
        }
    } 

    
    
}

