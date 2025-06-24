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

void benchmark_algorithm(const char *operation, const char *algorithm, int *sendbuf, int *recvbuf, int chunk_size, int iterations, MPI_Comm comm, int delay_rank, double delay_amount) {
    int rank, size;
    MPI_Info info;
    MPIX_Request request;
    double total_time = 0.0;
    double *times = (double *)malloc(iterations * sizeof(double));
    double total_delay = 0.0;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Create MPI_Info object and set the algorithm
    MPI_Info_create(&info);
    MPI_Info_set(info, "algorithm", algorithm);

    // Initialize the collective operation
    if (strcmp(operation, "alltoall") == 0) {
        MPIX_Alltoall_init(sendbuf, chunk_size, recvbuf, comm, info, &request);
    } else if (strcmp(operation, "broadcast") == 0) {
        MPIX_Broadcast_init(sendbuf, recvbuf, chunk_size, 0, comm, info, &request);
    } else if (strcmp(operation, "allgather") == 0) {
        MPIX_Allgather_init(sendbuf, chunk_size, recvbuf, comm, info, &request);
    } else {
        fprintf(stderr, "Unknown operation: %s\n", operation);
        MPI_Abort(comm, 1);
    }

    for (int i = 0; i < iterations; i++) {
        MPI_Barrier(comm);

        // Introduce delay for the selected rank
        if (rank == delay_rank) {
            usleep((useconds_t)(delay_amount * 1e6));
            total_delay += delay_amount;
        }

        double start_time = MPI_Wtime();

        // Start the collective operation
        MPIX_Start(&request);

        // Wait for the operation to complete
        request.wait_func(&request);

        double end_time = MPI_Wtime();
        times[i] = end_time - start_time;
        total_time += times[i];
    }

    // Calculate the average time
    double average_time = total_time / iterations;

    // Calculate the standard deviation
    double variance = 0.0;
    for (int i = 0; i < iterations; i++) {
        variance += (times[i] - average_time) * (times[i] - average_time);
    }
    variance /= iterations;
    double stddev = sqrt(variance);

    // Adjust the average time by subtracting the total delay
    if (rank == 0) {
        average_time -= total_delay / iterations;
    }

    // Print the average time and standard deviation
    if (rank == 0) {
        printf("Operation: %s, Algorithm: %s, Chunk Size: %d, Average Time: %f seconds, Standard Deviation: %f seconds\n",
               operation, algorithm, chunk_size, average_time, stddev);
    }

    // Free resources
    request.free_func(&request);
    MPI_Info_free(&info);
    free(times);
}

void benchmark_operation(const char *operation, int *sendbuf, int *recvbuf, int chunk_size, int iterations, MPI_Comm comm, int delay_rank, double delay_amount) {
    if (strcmp(operation, "alltoall") == 0) {
        benchmark_algorithm(operation, "direct", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
        benchmark_algorithm(operation, "pairwise", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
    } else if (strcmp(operation, "broadcast") == 0) {
        benchmark_algorithm(operation, "direct", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
        benchmark_algorithm(operation, "recursive_doubling", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
    } else if (strcmp(operation, "allgather") == 0) {
        benchmark_algorithm(operation, "direct", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
        benchmark_algorithm(operation, "recursive_doubling", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
        benchmark_algorithm(operation, "bruck", sendbuf, recvbuf, chunk_size, iterations, comm, delay_rank, delay_amount);
    } else {
        fprintf(stderr, "Unknown operation: %s\n", operation);
        MPI_Abort(comm, 1);
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    int *sendbuf, *recvbuf;
    int min_chunk_size, max_chunk_size, iterations;
    int delay_rank;
    double delay_amount;

    if (argc != 7) {
        fprintf(stderr, "Usage: %s <operation> <min_chunk_size> <max_chunk_size> <iterations> <delay_rank> <delay_amount>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char *operation = argv[1];
    min_chunk_size = atoi(argv[2]);
    max_chunk_size = atoi(argv[3]);
    iterations = atoi(argv[4]);
    delay_rank = atoi(argv[5]);
    delay_amount = atof(argv[6]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int chunk_size = min_chunk_size; chunk_size <= max_chunk_size; chunk_size *= 2) {
        int total_elements = chunk_size * size;

        // Allocate and initialize send buffer
        sendbuf = (int *)malloc(total_elements * sizeof(int));
	recvbuf = (int *)malloc(total_elements * sizeof(int));
        for (int i = 0; i < total_elements; i++) {
            sendbuf[i] = rank * total_elements + i; // Example data
        }

        // Benchmark the specified operation with valid algorithms
        benchmark_operation(operation, sendbuf, recvbuf, chunk_size, iterations, MPI_COMM_WORLD, delay_rank, delay_amount);

        free(sendbuf);
    }

    MPI_Finalize();
    return 0;
}
