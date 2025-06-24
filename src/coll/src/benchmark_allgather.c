/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "earlycoll.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> // For usleep()
#include <mpi.h>

void benchmark_allgather_custom(MPI_Comm comm, const char *algorithm, int chunk_size, int delay_rank, int delay_microseconds, int iterations) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

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

    double total_time = 0.0;

    double start_time = MPI_Wtime();

    for (int iter = 0; iter < iterations; iter++) {
        // Introduce delay for the specified rank
        if (rank == delay_rank) {
            usleep(delay_microseconds); // Delay in microseconds
        }

        MPIX_Start(&request);
        request.wait_func(&request);

    }

    double end_time = MPI_Wtime();
    total_time = end_time - start_time;
    total_time -= iterations * (delay_microseconds / 1e9); // Don't measure the delay.  

    request.free_func(&request);

    MPI_Info_free(&info);
    free(sendbuf);
    free(recvbuf);

    // Print summary statistics (only once, after all iterations)
    double avg_time = total_time / iterations;
    if (rank == 0) { // Only rank 0 prints the summary to avoid excessive output
        printf("Benchmark summary (custom algorithm '%s', chunk size %d, iterations %d):\n", algorithm, chunk_size, iterations);
        printf("  Average time: %.6f seconds\n", avg_time);
    }
}

void benchmark_allgather_default(MPI_Comm comm, int chunk_size, int delay_rank, int delay_microseconds, int iterations) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int total_elements = chunk_size * size;
    int *sendbuf = (int *)calloc(chunk_size, sizeof(int));
    int *recvbuf = (int *)calloc(total_elements, sizeof(int));
    for (int i = 0; i < chunk_size; i++) {
        sendbuf[i] = rank * chunk_size + i;
    }

    double total_time = 0.0;

    double start_time = MPI_Wtime();


    for (int iter = 0; iter < iterations; iter++) {
        // Introduce delay for the specified rank
        if (rank == delay_rank) {
            usleep(delay_microseconds); // Delay in microseconds
        }

        MPI_Allgather(sendbuf, chunk_size, MPI_INT, recvbuf, chunk_size, MPI_INT, comm);
    }

    double end_time = MPI_Wtime();
    total_time = end_time - start_time;
    total_time -= iterations * (delay_microseconds / 1e9); // Don't measure the delay.

    free(sendbuf);
    free(recvbuf);

    // Print summary statistics (only once, after all iterations)
    double avg_time = total_time / iterations;
    if (rank == 0) { // Only rank 0 prints the summary to avoid excessive output
        printf("Benchmark summary (default MPI persistent implementation, chunk size %d, iterations %d):\n", chunk_size, iterations);
        printf("  Average time: %.6f seconds\n", avg_time);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <algorithm|default> <chunk_size> <delay_rank> <delay_microseconds> <iterations>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char *algorithm = argv[1];
    int chunk_size = atoi(argv[2]);
    int delay_rank = atoi(argv[3]);
    int delay_microseconds = atoi(argv[4]);
    int iterations = atoi(argv[5]);

    if (chunk_size <= 0 || iterations <= 0 || delay_microseconds < 0) {
        fprintf(stderr, "Error: chunk_size and iterations must be positive integers, and delay_microseconds must be non-negative.\n");
        exit(EXIT_FAILURE);
    }

    MPI_Init(&argc, &argv);

    if (strcmp(algorithm, "default") == 0) {
        benchmark_allgather_default(MPI_COMM_WORLD, chunk_size, delay_rank, delay_microseconds, iterations);
    } else {
        benchmark_allgather_custom(MPI_COMM_WORLD, algorithm, chunk_size, delay_rank, delay_microseconds, iterations);
    }

    MPI_Finalize();
    return 0;
}

