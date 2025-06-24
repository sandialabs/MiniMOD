/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */


#ifndef EARLYCOLL_H
#define EARLYCOLL_H

#include <mpi.h>

typedef struct MPIX_Request MPIX_Request;

struct MPIX_Request {
    MPI_Win win;
    int *sendbuf, *recvbuf;
    int size;
    int chunk_size;
    int rank;
    int world_size;
    int root;
    void (*operation_func)(MPIX_Request *, int *, int, int, int);
    void (*wait_func)(MPIX_Request *);
    void (*free_func)(MPIX_Request *);

    // The following are algorithm specific and may not be needed/used.
    int num_rounds;
    int* rounds;
    MPI_Comm comm;
    MPI_Win win_rounds; // this is a round counter, it allows for the operations to continue once they have received data from the previous round (TODO - Evaluate if we might use this for direct as an completion counter)
    int* flags; 
    MPI_Win win_flags; // This is for R2S flags
    int *_recvbuf;  // This is a secondary buffer which allows us to do datacopy approaches. 

    // The following are exclusively for point to point. We should eventually make a union of this and the rma since they're mutually exclusive but for research purposes the extra memory probably doesn't hurt. 
    int round_count; 
    MPI_Request* recv_reqs;
    MPI_Request* send_reqs;
};

void MPIX_Alltoall_init(int *sendbuf, int chunk_size, int *recvbuf, MPI_Comm comm, MPI_Info info, MPIX_Request *request);
void MPIX_Broadcast_init(int *sendbuf, int *recvbuf, int chunk_size, int root, MPI_Comm comm, MPI_Info info, MPIX_Request *request);
void MPIX_Allgather_init(int *sendbuf, int chunk_size, int *recvbuf,  MPI_Comm comm, MPI_Info info, MPIX_Request *request);

// Currently there's some funkyness with these functions, we need to replace it

void MPIX_Start(MPIX_Request *request); //*sendbuf, int rank, int size);

//void MPIX_Wait(MPIX_Request *request) {
//    request->wait_func(request);
//}

//void MPIX_Request_free(MPIX_Request *request) {
//    request->free_func(request);
//}


#endif // EARLYCOLL_H

