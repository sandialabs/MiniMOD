/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_ERROR_H
#define MINIMOD_ERROR_H

#if defined(__cplusplus)
extern "C" {
#endif



#define MINIMOD_MODULE_RETURN_VALUES                                    \
    X(MINIMOD_MODULE_SUCCESS, "Module Success")                         \
        X(MINIMOD_MODULE_ERR_COMM_CRAYSHMEM, "Comm Craysmem error")            \
        X(MINIMOD_MODULE_ERR_COMM_MPI_PART, "Comm MPI Part error")             \
        X(MINIMOD_MODULE_ERR_COMM_MPI_PT2PT, "Comm MPI pt2pt error")           \
        X(MINIMOD_MODULE_ERR_COMM_MPI_RMA, "Comm MPI RMA error")               \
        X(MINIMOD_MODULE_ERR_COMM_OPENSHMEM, "Comm Open Shmem error")          \
        X(MINIMOD_MODULE_ERR_GRAN_BINS, "Gran Bins error")                     \
        X(MINIMOD_MODULE_ERR_GRAN_BULK, "Gran Bulk error")                     \
        X(MINIMOD_MODULE_ERR_GRAN_FINE, "Gran Fine error")                     \
        X(MINIMOD_MODULE_ERR_THREADS_OMP, "Threads OpenMP error")              \
        X(MINIMOD_MODULE_ERR_THREADS_PTHREADS, "Threads pthreads error")




#define MINIMOD_KERNEL_RETURN_VALUES                                    \
    X(MINIMOD_KERNEL_SUCCESS, "Kernel Success")                         \
        X(MINIMOD_KERNEL_ERR_MINIFE, "Kernel Minife error")



#define MINIMOD_ERROR_VALUES                              \
    X(MINIMOD_SUCCESS, "Minimod Success")                 \
        X(MINIMOD_ERR_MALLOC, "Malloc error")             \
        X(MINIMOD_ERR_MPI, "MPI error")                   \
        X(MINIMOD_ERR_MODULE_LOAD, "Module load error")


#define X(__err__, __msg__) __err__,
    typedef enum { MINIMOD_MODULE_RETURN_VALUES } MODULE_ERROR;
#undef X

#define X(__err__, __msg__) __err__,
    typedef enum { MINIMOD_KERNEL_RETURN_VALUES } KERNEL_ERROR;
#undef X

#define X(__err__, __msg__) __err__,
    typedef enum { MINIMOD_ERROR_VALUES } MINIMOD_ERROR;
#undef X


static const char* const minimod_module_error_messages[] = {
#define X(__err__, __msg__) __msg__,
    MINIMOD_MODULE_RETURN_VALUES
#undef X
};

static const char* const minimod_kernel_error_messages[] = {
#define X(__err__, __msg__) __msg__,
    MINIMOD_KERNEL_RETURN_VALUES
#undef X
};

static const char* const minimod_error_messages[] = {
#define X(__err__, __msg__) __msg__,
    MINIMOD_ERROR_VALUES
#undef X
};
    
#endif
