/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "earlycoll.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

void allgather_rma_init(MPIX_Request *req, int size, int *recvbuf, int chunk_size, MPI_Comm comm) {
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;
    MPI_Win_create(req->recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win);
}

void allgather_rma_init_datacopy(MPIX_Request *req, int size, int *recvbuf, int chunk_size, MPI_Comm comm) {
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;

    // Create the internal buffer for data copy. 
    req->_recvbuf = malloc(sizeof(int)*size);
    MPI_Win_create(req->_recvbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &req->win);
    MPI_Win_fence(0, req->win);
}

void allgather_rma_init_rounds_datacopy(MPIX_Request *req, int size, int *recvbuf, int chunk_size, MPI_Comm comm, int worldsize) {
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

void allgather_p2p_init_bruck_datacopy(MPIX_Request *req, int size, int rank, int *sendbuf, int *recvbuf, int chunk_size, MPI_Comm comm, int worldsize) {
    req->size = size;
    req->chunk_size = chunk_size;
    req->recvbuf = recvbuf;

    // Create the internal buffer for datacopy.
    req->_recvbuf = malloc(sizeof(int)*size);
    

    // Create the requests. 
    int i = 1;
    int round = 0;
    req->round_count = (int)ceil(log2(worldsize)); 

    req->send_reqs = malloc(sizeof(MPI_Request) * req->round_count);
    req->recv_reqs = malloc(sizeof(MPI_Request) * req->round_count);

    while(i<worldsize) {
        int send_rank = (rank + i) % worldsize;
        int recv_rank = (rank - i + worldsize) % worldsize;

        // Calculate the amount of data to send and recv the offset for receiving
        int send_size = (i * 2 > worldsize) ? (worldsize - i) * chunk_size : i * chunk_size;  // Amount of data to send in this round
        int recv_offset = i * chunk_size;  // Offset in recvbuf for received data


        MPI_Send_init(req->_recvbuf, send_size, MPI_INT, send_rank, round, comm, &req->send_reqs[round]);
        MPI_Recv_init(req->_recvbuf + recv_offset, send_size, MPI_INT, recv_rank, round, comm, &req->recv_reqs[round]);
        i <<= 1;
        round++;
    }
  
}

void allgather_rma_direct(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    MPI_Win_fence(0, req->win);

    // Direct RMA allgather with proper offsets
    for (int i = 0; i < total_processes; i++) {
        MPI_Put(sendbuf, chunk_size, MPI_INT, i, rank * chunk_size, chunk_size, MPI_INT, req->win);
    }
}

void allgather_rma_datacopy(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    // Direct RMA allgather with proper offsets
    for (int i = 0; i < total_processes; i++) {
        MPI_Put(sendbuf, chunk_size, MPI_INT, i, rank * chunk_size, chunk_size, MPI_INT, req->win);
    }
}


void allgather_rma_recursive_doubling(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    MPI_Win_fence(0, req->win);

    int offset_rank = rank;

    // Start with a self put (TODO: change this to a mem copy)
    memcpy(req->recvbuf + (rank * chunk_size), sendbuf, chunk_size * sizeof(int));

    // Recursive doubling allgather
    int mask = 1;
    while (mask < total_processes) {
        int partner = rank ^ mask;
        // Puts here are using 
        MPI_Win_fence(0, req->win);

        MPI_Put(req->recvbuf + (offset_rank*chunk_size), mask * chunk_size, MPI_INT, partner, offset_rank * chunk_size, mask * chunk_size, MPI_INT, req->win);
        if(partner < offset_rank) offset_rank = partner; // Update offest rank to the lowest set of the chunk we have. 
        mask <<= 1;
    }
}

// Global array and its size
MPI_Request request;

// Function to initialize the global array
void initialize_dummy() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int dest = (rank + 1) % size;
    MPI_Irecv(0, 0, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &request);
}

// Function to perform computations on the global array
void perform_dummy() {
    int flag;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
}

void end_dummy()
{
    MPI_Cancel(&request);
}

void allgather_rma_recursive_doubling_rounds_datacopy(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    int offset_rank = rank;
    initialize_dummy();
    memcpy(req->_recvbuf + (rank * chunk_size), sendbuf, chunk_size * sizeof(int));
    // Recursive doubling allgather
    int mask = 1;
    int round = 0;
    int one = 1; // Also Hacky, but we need a value for 1 in a pointer to incriment the rounds. 

    // Enter a pasive target epoch, to allow win flush
    MPI_Win_lock_all(0, req->win);
    MPI_Win_lock_all(0, req->win_rounds);


    while (mask < total_processes) {
        // Lock the window for the local process for win_rounds
        int local_round_value = (round == 0) ? 1 : req->rounds[round-1];

        if(local_round_value == 1)
        {
            int partner = rank ^ mask;
            MPI_Put(req->_recvbuf + (offset_rank*chunk_size), mask * chunk_size, MPI_INT, partner, offset_rank * chunk_size, mask * chunk_size, MPI_INT, req->win);
            MPI_Win_flush_all(req->win);

            // Lock the window for the target process for win_rounds
            MPI_Put(&one, 1, MPI_INT, partner, round, 1, MPI_INT, req->win_rounds);
            MPI_Win_flush_all(req->win_rounds);

            if(partner < offset_rank) offset_rank = partner; // Update offset rank to the lowest set of the chunk we have.
            mask <<= 1;
            round++;
        }
        perform_dummy();
    }
    end_dummy();
    MPI_Win_unlock_all(req->win);
    MPI_Win_unlock_all(req->win_rounds);
}


void allgather_rma_bruck(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {

    memcpy(req->recvbuf, sendbuf, req->chunk_size * sizeof(int));

    for (int i = 1; i < req->world_size; i <<= 1) {
        int put_rank = (rank + i) % total_processes;
        MPI_Win_fence(0, req->win);
        MPI_Put(req->recvbuf, i * chunk_size, MPI_INT, put_rank, i * chunk_size, i * chunk_size, MPI_INT, req->win);
    }
}

void allgather_rma_bruck_rounds_datacopy(MPIX_Request *req, int *sendbuf, int rank, int chunk_size, int total_processes) {
    initialize_dummy();
    //
    memcpy(req->_recvbuf, sendbuf, req->chunk_size * sizeof(int));
    int i = 1;
    int round = 0;
    int one = 1; // Also Hacky, but we need a value for 1 in a pointer to incriment the rounds.

    // Enter a pasive target epoch, to allow win flush
    MPI_Win_lock_all(0, req->win);
    MPI_Win_lock_all(0, req->win_rounds);

    while (i < total_processes) {
        // Lock the window for the local process for win_rounds
        int local_round_value = (round == 0) ? 1 : req->rounds[round-1];

        if(local_round_value == 1)
        {
            int partner = (rank + i) % total_processes;
            MPI_Put(req->_recvbuf, i * chunk_size, MPI_INT, partner, i * chunk_size, i * chunk_size, MPI_INT, req->win);
            MPI_Win_flush_all(req->win);

            MPI_Put(&one, 1, MPI_INT, partner, round, 1, MPI_INT, req->win_rounds);
            MPI_Win_flush_all(req->win_rounds);
            i <<= 1;
            round++;
        }
        perform_dummy();
    }
    end_dummy();
    MPI_Win_unlock_all(req->win);
    MPI_Win_unlock_all(req->win_rounds);
}


void allgather_rma_wait(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);
}

void allgather_rma_wait_datacopy(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);
    memcpy(req->recvbuf, req->_recvbuf, req->size*sizeof(int));
    MPI_Win_fence(0, req->win);
}

void allgather_rma_wait_rounds_datacopy(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    //while(1 != req->rounds[req->num_rounds - 1]) ; // Spin check

    while (1) {
        // Lock the window for the local process for win_rounds
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, req->rank, 0, req->win_rounds);

        // Check the completion status
        if (1 == req->rounds[req->num_rounds - 1]) {
            // Unlock the window for the local process for win_rounds
            MPI_Win_unlock(req->rank, req->win_rounds);
            break;
        }

        MPI_Win_unlock(req->rank, req->win_rounds);
    }

    memcpy(req->recvbuf, req->_recvbuf, req->size*sizeof(int));
    memset(req->rounds, 0, req->num_rounds * sizeof(int)); // Reset rounds.
    MPI_Barrier(req->comm); // Force synch to prevent overtaking and make sure that everyone is on the same allgather. 
}

void allgather_rma_wait_bruck(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    MPI_Win_fence(0, req->win);

    // Reorder the data back to its original order
    int *tempbuf = (int *)malloc(req->chunk_size * req->world_size * sizeof(int));
    for (int i = 0; i < req->world_size; i++) {
        int src_rank = (req->rank - i + req->world_size) % req->world_size;
        memcpy(tempbuf + i * req->chunk_size, req->recvbuf + src_rank * req->chunk_size, req->chunk_size * sizeof(int));
    }
    memcpy(req->recvbuf, tempbuf, req->chunk_size * req->world_size * sizeof(int));
    free(tempbuf);
}


void allgather_rma_wait_bruck_rounds_datacopy(MPIX_Request *req) {
    // Wait for the RMA communication to complete
    //while(1 != req->rounds[req->num_rounds - 1]) ; // Spin check

    while (1) {
        // Lock the window for the local process for win_rounds
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, req->rank, 0, req->win_rounds);

        // Check the completion status
        if (1 == req->rounds[req->num_rounds - 1]) {
            // Unlock the window for the local process for win_rounds
            //MPI_Win_unlock(req->rank, req->win_rounds); // This now happens below.
            break;
        }

        MPI_Win_unlock(req->rank, req->win_rounds);
    }

    
    // Lock the window for the local process for win
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, req->rank, 0, req->win);

    for (int i = 0; i < req->world_size; i++) {
        int src_rank = (req->rank - i + req->world_size) % req->world_size;
        // Prepare to print the values being moved
    //printf("Moving data from rank %d to index %d (size: %d): ", src_rank, i, req->chunk_size);
    
    // Print the values being copied
    //for (int j = 0; j < req->chunk_size; j++) {
    //    int value = req->_recvbuf[src_rank * req->chunk_size + j];
    //    printf("%d ", value);
    //}
    
   // printf("\n"); // New line after printing all values
        memcpy(req->recvbuf + i * req->chunk_size, req->_recvbuf + src_rank * req->chunk_size, req->chunk_size * sizeof(int));
    }
    memset(req->rounds, 0, req->num_rounds * sizeof(int)); // Reset rounds.

    MPI_Win_unlock(req->rank, req->win);
    MPI_Win_unlock(req->rank, req->win_rounds);
    MPI_Barrier(req->comm); // Force synch to prevent overtaking and make sure that everyone is on the same allgather. 
}

void allgather_rma_free(MPIX_Request *req) {
    MPI_Win_free(&req->win);
}

void allgather_rma_free_rounds(MPIX_Request *req) {
    MPI_Win_free(&req->win);
    MPI_Win_free(&req->win_rounds);
}

// TODO: We've cleaned the MPI state, but have some extra memory floating around. Fix that. 

void MPIX_Allgather_init(int *sendbuf, int chunk_size, int *recvbuf, MPI_Comm comm, MPI_Info info, MPIX_Request *request) {
    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    allgather_rma_init(request, chunk_size * size, recvbuf, chunk_size, comm);

    request->sendbuf = sendbuf;
    request->rank = rank;
    request->world_size = size;
    request->comm = comm;
    char algorithm[256];
    int flag;
    MPI_Info_get(info, "algorithm", 256, algorithm, &flag);

    request->wait_func = allgather_rma_wait;
    request->free_func = allgather_rma_free;

    if (flag) {
        if (strcmp(algorithm, "direct") == 0) {
            allgather_rma_init(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = allgather_rma_direct;
        } else if (strcmp(algorithm, "recursive_doubling") == 0) {
            allgather_rma_init(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = allgather_rma_recursive_doubling;
        } else if (strcmp(algorithm, "bruck") == 0) {
            allgather_rma_init(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = allgather_rma_bruck;
            request->wait_func = allgather_rma_wait_bruck;
        } else if (strcmp(algorithm, "datacopy") == 0) {
            allgather_rma_init_datacopy(request, chunk_size * size, recvbuf, chunk_size, comm);
            request->operation_func = allgather_rma_datacopy;
            request->wait_func = allgather_rma_wait_datacopy;
        } else if (strcmp(algorithm, "rd_dc_rnds") == 0) {
            allgather_rma_init_rounds_datacopy(request, chunk_size * size, recvbuf, chunk_size, comm, size);
            request->operation_func = allgather_rma_recursive_doubling_rounds_datacopy;
            request->wait_func = allgather_rma_wait_rounds_datacopy;
            request->free_func = allgather_rma_free_rounds;
        } else if (strcmp(algorithm, "b_dc_rnds") == 0) {
            allgather_rma_init_rounds_datacopy(request, chunk_size * size, recvbuf, chunk_size, comm, size);
            request->operation_func = allgather_rma_bruck_rounds_datacopy;
            request->wait_func = allgather_rma_wait_bruck_rounds_datacopy;
            request->free_func = allgather_rma_free_rounds;
        } else {
            fprintf(stderr, "MPIX_allgather_init Unknown algorithm: %s\n", algorithm);
            MPI_Abort(comm, 1);
        }
    } else {
        // Default to direct algorithm if no algorithm is specified
        request->operation_func = allgather_rma_direct;
    }

}

