/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "earlycoll.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void broadcast_rma_init(MPIX_Request *rma, int *recvbuf, int size, int chunk_size, MPI_Comm comm) {
    rma->size = size;
    rma->chunk_size = chunk_size;
    rma->recvbuf = recvbuf;
    MPI_Win_create(rma->recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &rma->win);
}

void broadcast_rma_init_datacopy(MPIX_Request *rma, int *recvbuf, int size, int chunk_size, MPI_Comm comm) {
    rma->size = size;
    rma->chunk_size = chunk_size;
    rma->recvbuf = recvbuf;
    rma->_recvbuf = malloc(size * sizeof(int));
    MPI_Win_create(rma->_recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &rma->win);
}

void broadcast_rma_direct(MPIX_Request *rma, int *sendbuf, int rank, int chunk_size, int total_processes) {
    MPI_Win_fence(0, rma->win);

    if (rank == rma->root) {
        for (int i = 0; i < total_processes; i++) {
            MPI_Put(sendbuf, chunk_size, MPI_INT, i, 0, chunk_size, MPI_INT, rma->win);
        }
    }
}

void broadcast_rma_datacopy(MPIX_Request *rma, int *sendbuf, int rank, int chunk_size, int total_processes) {
    if (rank == rma->root) {
        for (int i = 0; i < total_processes; i++) {
            MPI_Put(sendbuf, chunk_size, MPI_INT, i, 0, chunk_size, MPI_INT, rma->win);
        }
    }
}


void broadcast_rma_recursive_doubling(MPIX_Request *rma, int *sendbuf, int rank, int chunk_size, int total_processes) {
    
    if(rank == rma->root) memcpy(rma->recvbuf, sendbuf, chunk_size * sizeof(int));	    

    // Recursive doubling broadcast
    int mask = 1;

    // Get the relative rank for the purposes of the tree computation
    int relative_rank = (rank - rma->root + total_processes) % total_processes; 

    while (mask < total_processes) {
	MPI_Win_fence(0, rma->win);
	if(relative_rank < mask)
	{
            // This process is a sender in this step
            int partner = relative_rank | mask;
            if (partner < total_processes) {
                int partner_rank = (partner + rma->root) % total_processes;
                MPI_Put(rma->recvbuf, chunk_size, MPI_INT, partner_rank, 0, chunk_size, MPI_INT, rma->win);
            }
        } 
        mask <<= 1;
    }
}

void broadcast_rma_init_datacopy_rounds (MPIX_Request *req, int size, int *recvbuf, int chunk_size, MPI_Comm comm, int worldsize) {
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;

    // Create the internal buffer for datacopy. 
    req->_recvbuf = malloc(sizeof(int)*size);
    MPI_Win_create(req->_recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win);

    // Set up the rounds flags.
    req->num_rounds = (int)ceil(log2(worldsize));
    req->rounds = calloc(req->num_rounds, sizeof(int));
    MPI_Win_create(req->rounds, req->num_rounds * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win_rounds);
}

void broadcast_rma_datacopy_rounds(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {

    if(rank == req->root) memcpy(req->recvbuf, sendbuf, chunk_size * sizeof(int));

    // Recursive doubling broadcast
    int mask = 1;
    int round = 0;
    int one = 1;
    // Get the relative rank for the purposes of the tree computation
    int relative_rank = (rank - req->root + total_processes) % total_processes;

    while (mask < total_processes) {
        if(round == 0 || req->flags[round-1] == 1)
        {
            //MPI_Win_fence(0, rma->win);
            if(relative_rank < mask)
            {
                // This process is a sender in this step
                int partner = relative_rank | mask;
                if (partner < total_processes) {
                    int partner_rank = (partner + req->root) % total_processes;
                    MPI_Put(req->recvbuf, chunk_size, MPI_INT, partner_rank, 0, chunk_size, MPI_INT, req->win);
                    MPI_Win_flush_all(req->win);
                    MPI_Put(&one, 1, MPI_INT, partner_rank, round, 1, MPI_INT, req->win_rounds);
                    MPI_Win_flush_all(req->win_rounds);
                }
            }
        }
        mask <<= 1;
        round++;
    }
}



void broadcast_rma_wait(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);
}

void broadcast_rma_wait_datacopy(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);
    memcpy(req->recvbuf, req->_recvbuf, req->size*sizeof(int));
    MPI_Win_fence(0, req->win);
}


void broadcast_rma_free(MPIX_Request *rma) {
    MPI_Win_free(&rma->win);
}

void MPIX_Broadcast_init(int *sendbuf, int *recvbuf, int chunk_size, int root, MPI_Comm comm, MPI_Info info, MPIX_Request *request) {
    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    request->sendbuf = sendbuf;
    request->rank = rank;
    request->world_size = size;

    request->root = root;

    char algorithm[256];
    int flag;
    MPI_Info_get(info, "algorithm", 256, algorithm, &flag);

    request->wait_func = broadcast_rma_wait;
    request->free_func = broadcast_rma_free;

    if (flag) {
        if (strcmp(algorithm, "direct") == 0) {
            broadcast_rma_init(request, recvbuf, chunk_size * size, chunk_size, comm);
	    request->operation_func = broadcast_rma_direct;
        } else if (strcmp(algorithm, "recursive_doubling") == 0) {
            broadcast_rma_init(request, recvbuf, chunk_size * size, chunk_size, comm);
   	    request->operation_func = broadcast_rma_recursive_doubling;
        } else if (strcmp(algorithm, "datacopy") == 0) {
            broadcast_rma_init_datacopy(request, recvbuf, chunk_size * size, chunk_size, comm);
            request->operation_func = broadcast_rma_datacopy;
        } else {
            fprintf(stderr, "Unknown algorithm: %s\n", algorithm);
            MPI_Abort(comm, 1);
        }
    } else {
        // Default to direct algorithm if no algorithm is specified
        request->operation_func = broadcast_rma_direct;
    }

    request->root = root;
}

