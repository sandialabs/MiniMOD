/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "earlycoll.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void verify_allgather(int *recvbuf, int chunk_size, int rank, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < chunk_size; j++) {
            int expected = i * chunk_size + j;
            if (recvbuf[i * chunk_size + j] != expected) {
                printf("Allgather verification failed at rank %d, index %d: expected %d, got %d\n",
                       rank, i * chunk_size + j, expected, recvbuf[i * chunk_size + j]);
                return;
            }
        }
    }
    printf("Allgather verification passed at rank %d\n", rank);
}

void print_buffer(const char *label, int *buffer, int size, int rank) {
    printf("%s at rank %d: ", label, rank);
    for (int i = 0; i < size; i++) {
        printf("%d ", buffer[i]);
    }
    printf("\n");
}

void test_allgather(MPI_Comm comm, const char *algorithm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int chunk_size = 4;
    int total_elements = chunk_size * size;
    int *sendbuf = (int *)calloc(chunk_size, sizeof(int));
    int *recvbuf = (int *)calloc(total_elements, sizeof(int));
    for (int i = 0; i < chunk_size; i++) {
        sendbuf[i] = rank * chunk_size + i;
    }

    MPIX_Request request;
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "algorithm", algorithm);

    MPIX_Allgather_init(sendbuf, chunk_size, recvbuf, comm, info, &request);

    MPIX_Start(&request);

    request.wait_func(&request);

    request.free_func(&request);

    print_buffer("Send buffer", sendbuf, chunk_size, rank);
    print_buffer("Receive buffer", recvbuf, total_elements, rank);

    verify_allgather(recvbuf, chunk_size, rank, size);

    MPI_Info_free(&info);
    free(sendbuf);
    free(recvbuf);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <algorithm>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char *algorithm = argv[1];

    MPI_Init(&argc, &argv);

    test_allgather(MPI_COMM_WORLD, algorithm);

    MPI_Finalize();
    return 0;
}

