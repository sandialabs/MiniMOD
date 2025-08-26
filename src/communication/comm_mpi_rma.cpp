/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include <stdlib.h>
#include "comm_mpi_rma.h"
#include <assert.h>

#ifdef HAVE_MPI


#define DEBUG 0 
#define VERBOSE 0

namespace Comm_MPI_RMA
{
  int comm_mpi_rma_init(int threads)
  {
    int err;
    if(threads)
    {
      int provided;
      err = MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
      if(provided != MPI_THREAD_MULTIPLE)
      {
        printf("MPI_Init_thread did not provide MPI_THREAD_MULTIPLE\n"); 
        return -9;
      }
    }
    else
    {
      err = MPI_Init(0, 0);
    }
    if(err != MPI_SUCCESS)
    {
      printf("MPI_Init_thread returned an error\n"); 
      return -3;
    }
    return 0;
  }

  int comm_mpi_rma_get_id(int* id)
  {
    int err = MPI_Comm_rank(MPI_COMM_WORLD, id);
    if(err != MPI_SUCCESS) return -5;
    return 0;
  }

  int comm_mpi_rma_job_size(int* size)
  {
    int err = MPI_Comm_size(MPI_COMM_WORLD, size);
    if(err != MPI_SUCCESS) return -5;
    return 0;
  }

  int comm_mpi_rma_clean()
  {
    int err = MPI_Finalize();
    if(err != MPI_SUCCESS) return -4;
    return 0;
  }

  int comm_mpi_rma_barrier()
  {
    int err = MPI_Barrier(MPI_COMM_WORLD);
    if(err != MPI_SUCCESS) return -10;
    return 0;
  }

  MPI_Datatype comm_mpi_rma_convert_datatype(comm_datatype datatype)
  {
    switch(datatype)
    {
      case COMM_DATATYPE_INT:
        return MPI_INT;
      case COMM_DATATYPE_FLOAT:
        return MPI_FLOAT;
      case COMM_DATATYPE_LONGLONG:
        return MPI_LONG_LONG;
      case COMM_DATATYPE_DOUBLE:
        return MPI_DOUBLE;
    }
    printf("Error: Datatype not recognized!");
    MPI_Abort(MPI_COMM_WORLD,-1); 
    return MPI_INT; 
  }

  MPI_Op comm_mpi_rma_convert_op(comm_reduce_op op)
  {
    switch(op)
    {
      case COMM_REDUCE_OP_SUM:
        return MPI_SUM;
      case COMM_REDUCE_OP_MAX:
        return MPI_MAX;
      case COMM_REDUCE_OP_MIN:
        return MPI_MIN;
    }
    printf("Error: Datatype not recognized!");
    MPI_Abort(MPI_COMM_WORLD,-1); 
    return MPI_SUM; 
  }

  int comm_mpi_rma_broadcast(void* buffer, int count, comm_datatype datatype, int root)
  {
    return MPI_Bcast(
      buffer,
      count,
      comm_mpi_rma_convert_datatype(datatype),
      root,
      MPI_COMM_WORLD);
  }

  int comm_mpi_rma_reduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op, int root)
  {
    return MPI_Reduce(
      send_buf,
      recv_buf,
      count,
      comm_mpi_rma_convert_datatype(datatype), 
      comm_mpi_rma_convert_op(op),
      root,
      MPI_COMM_WORLD);
  }

  int comm_mpi_rma_allreduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op)
  {
    return MPI_Allreduce(
      send_buf,
      recv_buf,
      count,
      comm_mpi_rma_convert_datatype(datatype),
      comm_mpi_rma_convert_op(op),
      MPI_COMM_WORLD);
  }

  int comm_mpi_rma_allgather(void* send_buf, int send_count, comm_datatype send_datatype, void* recv_buf, int recv_count, comm_datatype recv_datatype)
  {
      return MPI_Allgather(
      send_buf,
      send_count,
      comm_mpi_rma_convert_datatype(send_datatype), 
      recv_buf,
      recv_count,
      comm_mpi_rma_convert_datatype(recv_datatype), 
      MPI_COMM_WORLD);
  }

  int comm_mpi_rma_channel_init(comm_request** comm_req, void** recv_buf, void* send_buf, int recv_id, int send_id, int num_entries, int entry_size)
  { 
    
    comm_mpi_rma_req* req = (comm_mpi_rma_req*)malloc(sizeof(comm_mpi_rma_req)); 
    if(req == NULL) return -4; 
 
    
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);


    int err = MPI_Win_allocate(num_entries * entry_size, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &req->recv_buf, &req->win);
    if(err != MPI_SUCCESS)
    {
      printf("MPI_Win_allocate returned an error\n"); 
      return -1;
    }
    if(req->recv_buf == NULL) return -4; 
 
    
    req->send_buf = send_buf;
    req->recv_id = recv_id;
    req->send_id = send_id;
    comm_mpi_rma_get_id(&req->my_id);
    req->num_entries = num_entries;
    req->entry_size = entry_size;
    req->ready = 0;
    
    *comm_req = req;
    *recv_buf = req->recv_buf;
 

    return 0;
  }

  int comm_mpi_rma_channel_start(comm_request* comm_req)
  {
    comm_mpi_rma_req* req = (comm_mpi_rma_req*)comm_req; 
    return MPI_Win_fence(0, req->win);







  }

  int comm_mpi_rma_channel_send(comm_request* comm_req, int index)
  {
    comm_mpi_rma_req* req = (comm_mpi_rma_req*)comm_req;  
    if(req->send_id != -1)
    {
      MPI_Put(((unsigned char*)req->send_buf) + (index * req->entry_size), req->entry_size, MPI_BYTE, req->send_id, index * req->entry_size, req->entry_size, MPI_BYTE, req->win);
    }
    return 0;
  }

  int comm_mpi_rma_channel_send_range(comm_request* comm_req, int index_start, int index_end)
  {
    comm_mpi_rma_req* req = (comm_mpi_rma_req*)comm_req;  
    if(req->send_id != -1)
    {
      if(index_start - index_end > 1) return -1;
      MPI_Put(((unsigned char*)req->send_buf) + (index_start * req->entry_size), (index_end - index_start + 1) * req->entry_size, MPI_BYTE, req->send_id, index_start * req->entry_size, (index_end - index_start + 1) * req->entry_size, MPI_BYTE, req->win);
    }
    return 0;
  }

  int comm_mpi_rma_channel_end(comm_request* comm_req)
  {
    comm_mpi_rma_req* req = (comm_mpi_rma_req*)comm_req;
    return MPI_Win_fence(0, req->win);






  }

  int comm_mpi_rma_channel_finalize(comm_request* comm_req)
  {
    comm_mpi_rma_req* req = (comm_mpi_rma_req*)comm_req;
    MPI_Win_free(&(req->win));
    free(comm_req);
    return 0;
  }

  int comm_mpi_rma_load(comm_funcs* funcs)
  {
    funcs->init = &comm_mpi_rma_init;
    funcs->get_id = &comm_mpi_rma_get_id;
    funcs->job_size = &comm_mpi_rma_job_size;
    funcs->channel_init = &comm_mpi_rma_channel_init;
    funcs->channel_start = &comm_mpi_rma_channel_start;
    funcs->channel_send = &comm_mpi_rma_channel_send;
    funcs->channel_send_range = &comm_mpi_rma_channel_send_range;
    funcs->channel_end = &comm_mpi_rma_channel_end;
    funcs->channel_finalize = &comm_mpi_rma_channel_finalize;
    funcs->barrier = &comm_mpi_rma_barrier;
    funcs->broadcast = &comm_mpi_rma_broadcast;
    funcs->reduce = &comm_mpi_rma_reduce;
    funcs->allreduce = &comm_mpi_rma_allreduce;
    funcs->allgather = &comm_mpi_rma_allgather;
    funcs->clean = &comm_mpi_rma_clean;
    return 0;
  }
} 
#endif
