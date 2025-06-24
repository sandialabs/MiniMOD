/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "earlycoll.h" // Include your custom library for MPIX implementations

// Function to print help message
void print_help() {
    printf("Usage: mpirun -np <num_processes> ./benchmark [options]\n");
    printf("Options:\n");
    printf("  -o <operation>       Specify the operation (alltoall, allgather, broadcast) [REQUIRED]\n");
    printf("  -a <algorithm>       Specify the algorithm (default, direct, pairwise, etc.) [default: default]\n");
    printf("  -m <message_size>    Specify the message size in integers [default: 1024]\n");
    printf("  -i <iterations>      Specify the number of iterations [default: 10]\n");
    printf("  -distr <distribution> Specify the delay distribution (laggard, gaussian, burst) [default: laggard]\n");
    printf("  -p1 <param1>         For laggard: delay in nanoseconds.\n");
    printf("                       For gaussian: standard deviation of computational time in nanoseconds.\n");
    printf("                       For burst: computational spike in nanoseconds.\n");
    printf("  -p2 <param2>         For laggard: rank to introduce serialized work.\n");
    printf("                       For gaussian: mean computational time in nanoseconds.\n");
    printf("                       For burst: number of processes experiencing the spike.\n");
    printf("  -verify              Enable verification of results after each iteration [default: off]\n");
    printf("  -d                   Enable debug messages [default: off]\n");
    printf("  -h                   Show this help message and exit\n");
    printf("\nExample:\n");
    printf("  mpirun -np 4 ./benchmark -o alltoall -a default -m 1024 -i 10 -distr gaussian -p1 500000 -p2 1000000 -verify -d\n");
    exit(EXIT_SUCCESS);
}

// Function to induce delay based on the specified distribution and rank
void induce_delay(const char *distribution, long param1, long param2, int rank, int delayed_rank, int iteration, int debug) {
    struct timespec req, rem;
    req.tv_sec = 0;
    req.tv_nsec = 0;

    if (strcmp(distribution, "laggard") == 0) {
        // Apply delay only to the specified delayed rank
        if (delayed_rank != -1 && rank != delayed_rank) {
            if (debug) {
                printf("Rank %d: No delay applied (not the delayed rank).\n", rank);
            }
            return;
        }
        req.tv_nsec = param1;
    } else if (strcmp(distribution, "gaussian") == 0) {
        // Gaussian delay using Box-Muller transform
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        req.tv_nsec = param2 + (long)(z * param1); // Mean + Gaussian noise
        if (req.tv_nsec < 0) req.tv_nsec = 0; // Ensure non-negative delay
    } else if (strcmp(distribution, "burst") == 0) {
        // Apply burst delay to a subset of processes
        if (rank < param2) { // Apply delay to the first param2 ranks
            req.tv_nsec = param1;
        }
    }

    if (debug) {
        printf("Rank %d: Inducing delay: %ld ns\n", rank, req.tv_nsec);
    }

    nanosleep(&req, &rem);
}

// Function to verify the results of the collective operation
int verify_results(const char *operation, int *sendbuf, int *recvbuf, int message_size, int rank, int size, int debug) {
    int total_elements = message_size * size;

    if (strcmp(operation, "alltoall") == 0) {
        for (int i = 0; i < total_elements; i++) {
            int expected = (i / message_size) * total_elements + (rank * message_size) + (i % message_size);
            if (recvbuf[i] != expected) {
                if (debug) {
                    printf("Rank %d: Verification failed at index %d. Expected %d, got %d.\n", rank, i, expected, recvbuf[i]);
                }
                return 0; // Verification failed
            }
        }
    } else if (strcmp(operation, "allgather") == 0) {
        for (int i = 0; i < total_elements; i++) {
            if (recvbuf[i] != i) {
                if (debug) {
                    printf("Rank %d: Verification failed at index %d. Expected %d, got %d.\n", rank, i, i, recvbuf[i]);
                }
                return 0; // Verification failed
            }
        }
    } else if (strcmp(operation, "broadcast") == 0) {
        for (int i = 0; i < message_size; i++) {
            if (recvbuf[i] != i) {
                if (debug) {
                    printf("Rank %d: Verification failed at index %d. Expected %d, got %d.\n", rank, i, i, recvbuf[i]);
                }
                return 0; // Verification failed
            }
        }
    } else {
        fprintf(stderr, "Unknown operation: %s\n", operation);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return 1; // Verification passed
}

// Benchmark function for a single operation and algorithm
void benchmark(const char *operation, const char *algorithm, int *sendbuf, int *recvbuf, int message_size, int iterations, MPI_Comm comm, const char *distribution, long param1, long param2, int delayed_rank, int debug, int verify) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double start_time, end_time, total_time = 0.0;

    if (debug && rank == 0) {
        printf("Starting benchmark: operation=%s, algorithm=%s, message_size=%d, iterations=%d, delayed_rank=%d\n", operation, algorithm, message_size, iterations, delayed_rank);
    }

    MPI_Info info;
    MPIX_Request request;
    MPI_Info_create(&info);
    MPI_Info_set(info, "algorithm", algorithm); // Set the algorithm in MPI_Info

    // Initialize the collective operation
    if (strcmp(operation, "alltoall") == 0) {
        MPIX_Alltoall_init(sendbuf, message_size, recvbuf, comm, info, &request);
    } else if (strcmp(operation, "allgather") == 0) {
        MPIX_Allgather_init(sendbuf, message_size, recvbuf, comm, info, &request);
    } else if (strcmp(operation, "broadcast") == 0) {
        MPIX_Broadcast_init(sendbuf, recvbuf, message_size, 0, comm, info, &request);
    } else {
        fprintf(stderr, "Unknown operation: %s\n", operation);
        MPI_Abort(comm, 1);
    }

    if (debug && rank == 0) { printf("Start iterations"); }

    for (int i = 0; i < iterations; i++) {
        MPI_Barrier(comm);
        start_time = MPI_Wtime();

        // Introduce delay for the current rank
        induce_delay(distribution, param1, param2, rank, delayed_rank, i, debug);

        // Start the collective operation
        MPIX_Start(&request);

        // Wait for the operation to complete
        request.wait_func(&request);

        end_time = MPI_Wtime();
        total_time += (end_time - start_time);

        if (verify) {
            if (!verify_results(operation, sendbuf, recvbuf, message_size, rank, size, debug)) {
                fprintf(stderr, "Rank %d: Verification failed in iteration %d.\n", rank, i + 1);
                MPI_Abort(comm, 1);
            }
        }

        if (debug && rank == 0) {
            printf("Iteration %d completed, time=%f seconds\n", i + 1, end_time - start_time);
        }
    }

    // Free resources
    request.free_func(&request);
    MPI_Info_free(&info);

    // Print the average time taken for the operation
    if (rank == 0) {
        printf("Operation: %s, Algorithm: %s, Message Size: %d, Average Time: %f seconds\n", operation, algorithm, message_size, total_time / iterations);
    }
}

int main(int argc, char *argv[]) {
    const char *operation = NULL;
    const char *algorithm = "default";
    int message_size = 1024;
    int iterations = 10;
    const char *distribution = "laggard";
    long param1 = 1000000; // Default delay in nanoseconds for laggard, or standard deviation for Gaussian
    long param2 = 0;       // Default mean delay for Gaussian or delayed rank for laggard
    int debug = 0;
    int verify = 0;

    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            print_help();
        } else if (strcmp(argv[i], "-o") == 0) {
            operation = argv[++i];
        } else if (strcmp(argv[i], "-a") == 0) {
            algorithm = argv[++i];
        } else if (strcmp(argv[i], "-m") == 0) {
            message_size = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-i") == 0) {
            iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-distr") == 0) {
            distribution = argv[++i];
        } else if (strcmp(argv[i], "-p1") == 0) {
            param1 = atol(argv[++i]);
        } else if (strcmp(argv[i], "-p2") == 0) {
            param2 = atol(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0) {
            debug = 1;
        } else if (strcmp(argv[i], "--verify") == 0) {
            verify = 1;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if (operation == NULL) {
        fprintf(stderr, "Error: Operation (-o) must be specified.\n");
        exit(EXIT_FAILURE);
    }

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



    int *sendbuf = (int *)malloc(2*message_size * size * sizeof(int));
    if (sendbuf == NULL) {
        fprintf(stderr, "Rank %d: Failed to allocate memory for sendbuf (size: %lu bytes)\n", rank, message_size * size * sizeof(int));
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    } else {
        printf("Rank %d: Successfully allocated memory for sendbuf (size: %lu bytes)\n", rank, message_size * size * sizeof(int));
    }

    int *recvbuf = (int *)malloc(2* message_size * size * sizeof(int));
    if (recvbuf == NULL) {
        fprintf(stderr, "Rank %d: Failed to allocate memory for recvbuf (size: %lu bytes)\n", rank, message_size * size * sizeof(int));
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    } else {
        printf("Rank %d: Successfully allocated memory for recvbuf (size: %lu bytes)\n", rank, message_size * size * sizeof(int));
    }


    //int *sendbuf = (int *)malloc(message_size * size * sizeof(int));
    //int *recvbuf = (int *)malloc(message_size * size * sizeof(int));

    // Initialize send buffer
    if (strcmp(operation, "allgather") == 0) {
        for (int i = 0; i < message_size; i++) {
            sendbuf[i] = rank * message_size + i;
        }
    } else if (strcmp(operation, "alltoall") == 0) {
        for (int i = 0; i < message_size * size; i++) {
            sendbuf[i] = rank * message_size + i;
        }
    } else if (strcmp(operation, "broadcast") == 0) {
        if (rank == 0) {
            for (int i = 0; i < message_size; i++) {
                sendbuf[i] = i;
            }
        }
    }

    benchmark(operation, algorithm, sendbuf, recvbuf, message_size, iterations, MPI_COMM_WORLD, distribution, param1, param2, param2, debug, verify);

    free(sendbuf);
    free(recvbuf);

    MPI_Finalize();
    return 0;
}
