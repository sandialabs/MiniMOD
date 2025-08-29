#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void rma_buffers(int rank, int size) {
    int *firstbuf = (int *)malloc((size * sizeof(int)) / 2);
    int *secondbuf = (int *)malloc((size * sizeof(int)) / 2);
    int *targetbuf = (int *)malloc(size * sizeof(int));

    // we want to split ranks to put into a shared targetbuf -- using 2 windows
    if (rank % 2 == 1) {
        firstbuf[rank] = rank;
    } else {
        secondbuf[rank] = rank;
    }

    // Second communicator
    MPI_Comm comm2;
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &comm2);

    // MPI_Win win1, win2; // this is how it would be done with 2 windows
    MPI_Win targetwin, targetwin2;
    // MPI_Win_create(firstbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_create(targetbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm2, &targetwin);
    MPI_Win_create(targetbuf, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm2, &targetwin2);

    // if (rank % 2 == 1) {
    //     MPI_Win_fence(0, targetwin);
    //     for (int i = 0; i < size; i++){
    //         MPI_Put(&firstbuf[rank], 1, MPI_INT, i, rank, 1, MPI_INT, targetwin);
    //     }
    //     MPI_Win_fence(0, targetwin);
    // } else {
    //     MPI_Win_fence(0, targetwin2);
    //     for (int i = 0; i < size; i++){
    //         MPI_Put(&secondbuf[rank], 1, MPI_INT, i, rank, 1, MPI_INT, targetwin2);
    //     }
    //     MPI_Win_fence(0, targetwin2);
    // }

    MPI_Win_free(&targetwin);
    MPI_Win_free(&targetwin2);

    // if (rank == 0){
    //     printf("targetbuf at %d size:\n", size);
    //     for (int i = 0; i < size; i++) {
    //         printf("%d\n", targetbuf[i]);
    //     }
    // }

    free(firstbuf);
    free(secondbuf);
    free(targetbuf);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    rma_buffers(rank, size);

    MPI_Finalize();
    return 0;
}