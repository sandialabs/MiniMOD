/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_COMM
#define MINIMOD_COMM

#include <stddef.h>


typedef void comm_request;

typedef int comm_datatype;
#define COMM_DATATYPE_INT 0
#define COMM_DATATYPE_LONGLONG 1
#define COMM_DATATYPE_FLOAT 2 
#define COMM_DATATYPE_DOUBLE 3

typedef int comm_reduce_op;
#define COMM_REDUCE_OP_SUM 0
#define COMM_REDUCE_OP_MAX 1
#define COMM_REDUCE_OP_MIN 2

#define COMM_ANY_ID -2

typedef struct
{
  int (*init)               (int);
  int (*get_id)             (int*);
  int (*job_size)           (int*);

  int (*channel_init)       (comm_request**, void** recv_buf, void* send_buf, int recv_id, int send_id, int num_entries, int entry_size);
  int (*channel_start)      (comm_request*);
  int (*channel_send)       (comm_request*, int index);
  int (*channel_send_range) (comm_request*, int index_start, int index_end);
  int (*channel_end)        (comm_request*);
  int (*channel_finalize)   (comm_request*);

  int (*barrier)            ();
  int (*broadcast)          (void* buffer, int count, comm_datatype datatype, int root_id);
  int (*reduce)             (void* send_buffer, void* recv_buffer, int count, comm_datatype datatype, comm_reduce_op op, int root);
  int (*allreduce)          (void* send_buffer, void* recv_buffer, int count, comm_datatype datatype, comm_reduce_op op);
  int (*allgather)          (void* send_buf, int send_count, comm_datatype send_datatype, void* recv_buf, int recv_count, comm_datatype recv_datatype);
  int (*clean)              ();
} comm_funcs;

int comm_load(char* name);
int comm_init(int threads);
int comm_get_id(int*);
int comm_get_job_size(int*);
int comm_channel_init(comm_request**, void** recv_buf, void* send_buf, int recv_id, int send_id, int num_entries, int entry_size); 
int comm_channel_start(comm_request*);
int comm_channel_send(comm_request*, int index);
int comm_channel_send_range(comm_request*, int index_start, int index_end);
int comm_channel_end(comm_request*);
int comm_channel_finalize (comm_request*);
int comm_barrier();
int comm_broadcast(void*, int, comm_datatype, int);
int comm_reduce(void*, void*, int, comm_datatype, comm_reduce_op, int);
int comm_allreduce(void*, void*, int, comm_datatype, comm_reduce_op);
int comm_allgather(void*, int, comm_datatype, void*, int, comm_datatype);
int comm_clean();

#include"communication/comm_mpi_pt2pt.h" 
#ifdef HAVE_MPI_PART
#include"communication/comm_mpi_part.h" 
#endif
#ifdef HAVE_MPI_PCL
#include"communication/comm_mpi_pcl.h" 
#endif
#ifdef HAVE_OSHMEM
#include"communication/comm_openshmem.h" 
#endif
#ifdef HAVE_CRAY_SHMEM
#include"communication/comm_crayshmem.h" 
#endif
#include"communication/comm_mpi_pers.h" 
#include"communication/comm_mpi_rma.h"

#endif
