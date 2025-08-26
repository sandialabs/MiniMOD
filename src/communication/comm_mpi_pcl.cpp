/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifdef HAVE_MPI

#include <stdlib.h>
#include "comm_mpi_pcl.h"
#include <assert.h>


#define DEBUG 0 
#define VERBOSE 0




namespace Comm_MPI_PCL
{
  int comm_mpi_pcl_init(int threads)
  {
    int err;
    if(threads)
    {
      int provided;
      err = MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
      if(provided != MPI_THREAD_MULTIPLE)
      {
        return -9;
      }
    }
    else
    {
      err = MPI_Init(0, 0);
    }
    if(err != MPI_SUCCESS) return -3;



    return 0;
  }

  int comm_mpi_pcl_get_id(int* id)
  {
    int err = MPI_Comm_rank(MPI_COMM_WORLD, id);
    if(err != MPI_SUCCESS) return -5;
    return 0;
  }

  int comm_mpi_pcl_job_size(int* size)
  {
    int err = MPI_Comm_size(MPI_COMM_WORLD, size);
    if(err != MPI_SUCCESS) return -5;
    return 0;
  }

  int comm_mpi_pcl_clean()
  {
    int err = MPI_Finalize();
    if(err != MPI_SUCCESS) return -4;
    return 0; 
  }

  int comm_mpi_pcl_barrier()
  {
    int err = MPI_Barrier(MPI_COMM_WORLD);
    if(err != MPI_SUCCESS) return -10;
    return 0;
  }

  MPI_Datatype comm_mpi_pcl_convert_datatype(comm_datatype datatype)
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

  MPI_Op comm_mpi_pcl_convert_op(comm_reduce_op op)
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

  int comm_mpi_pcl_broadcast(void* buffer, int count, comm_datatype datatype, int root)
  {
    return MPI_Bcast(
      buffer, 
      count, 
      comm_mpi_pcl_convert_datatype(datatype), 
      root, 
      MPI_COMM_WORLD);
  }

  int comm_mpi_pcl_reduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op, int root)
  {
    return MPI_Reduce(
      send_buf,
      recv_buf,
      count,
      comm_mpi_pcl_convert_datatype(datatype), 
      comm_mpi_pcl_convert_op(op),
      root,
      MPI_COMM_WORLD);
  }

  int comm_mpi_pcl_allreduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op)
  {
    return MPI_Allreduce(
      send_buf,
      recv_buf,
      count,
      comm_mpi_pcl_convert_datatype(datatype), 
      comm_mpi_pcl_convert_op(op),
      MPI_COMM_WORLD);
  }

  int comm_mpi_pcl_allgather(void* send_buf, int send_count, comm_datatype send_datatype, void* recv_buf, int recv_count, comm_datatype recv_datatype)
  {
      return MPI_Allgather(
      send_buf,
      send_count,
      comm_mpi_pcl_convert_datatype(send_datatype), 
      recv_buf,
      recv_count,
      comm_mpi_pcl_convert_datatype(recv_datatype), 
      MPI_COMM_WORLD);
  }

  int comm_mpi_pcl_channel_init(comm_request** comm_req, void** recv_buf, void* send_buf, int recv_id, int send_id, int num_entries, int entry_size)
  { 
    
    comm_mpi_pcl_req* req = (comm_mpi_pcl_req*) malloc(sizeof(comm_mpi_pcl_req)); 
    if(req == NULL) { printf("req didn't malloc"); return -4; } 

    
    req->recv_buf = malloc(num_entries * entry_size);
    if(req->recv_buf == NULL) return -4; 
    
    


    if(recv_id != -1)
    {

      MPIX_Precv_init(req->recv_buf, num_entries, entry_size, MPI_BYTE, recv_id, 1, MPI_COMM_WORLD, &req->recv_req);
    }
    if(send_id != -1)
    {

      MPIX_Psend_init(send_buf, num_entries, entry_size, MPI_BYTE, send_id, 1, MPI_COMM_WORLD, &req->send_req);
    }

    
    req->send_buf = send_buf;
    req->recv_id = recv_id;
    req->send_id = send_id;
    req->num_entries = num_entries;
    req->entry_size = entry_size;

    
    *comm_req = req;
    *recv_buf = req->recv_buf;
    return 0;
  }

  int comm_mpi_pcl_channel_start(comm_request* comm_req)
  {
    comm_mpi_pcl_req* req = (comm_mpi_pcl_req*) comm_req;
    if(req->send_id != -1)
    {
      MPIX_Start(&req->send_req);
    }
    if(req->recv_id != -1)
    {
      MPIX_Start(&req->recv_req);
    }
    return 0;
  }

  int comm_mpi_pcl_channel_send(comm_request* comm_req, int index)
  {
    comm_mpi_pcl_req* req = (comm_mpi_pcl_req*) comm_req;
    if(req->send_id != -1)
    {
      return MPIX_Pready(index, &req->send_req);
    }
    else
    {
      return 0;
    }
  }

  int comm_mpi_pcl_channel_end(comm_request* comm_req)
  {
    comm_mpi_pcl_req* req = (comm_mpi_pcl_req*)comm_req;
    int err; 
    if(req->send_id != -1)
    {
      err = MPIX_Wait(&req->send_req, MPI_STATUSES_IGNORE);
      if(err != MPI_SUCCESS) return -4; 
    }
    if(req->recv_id != -1)
    {
      err = MPIX_Wait(&req->recv_req, MPI_STATUSES_IGNORE);
      if(err != MPI_SUCCESS) return -4; 
    }
    return 0;
  }

  int comm_mpi_pcl_channel_finalize(comm_request* comm_req)
  {
    comm_mpi_pcl_req* req = (comm_mpi_pcl_req*) comm_req;
    if(req->send_id != -1) {
      MPIX_Request_free(&req->send_req);
    }
    if(req->recv_id != -1) {
      MPIX_Request_free(&req->recv_req);
    }
    free(req->recv_buf);
    free(req);
    return 0;
  }

  int comm_mpi_pcl_load(comm_funcs* funcs)
  {
    funcs->init = &comm_mpi_pcl_init;
    funcs->get_id = &comm_mpi_pcl_get_id;
    funcs->job_size = &comm_mpi_pcl_job_size;
    funcs->channel_init = &comm_mpi_pcl_channel_init;
    funcs->channel_start = &comm_mpi_pcl_channel_start;
    funcs->channel_send = &comm_mpi_pcl_channel_send;
    funcs->channel_end = &comm_mpi_pcl_channel_end;
    funcs->channel_finalize = &comm_mpi_pcl_channel_finalize;
    funcs->barrier = &comm_mpi_pcl_barrier;
    funcs->broadcast = &comm_mpi_pcl_broadcast;
    funcs->reduce = &comm_mpi_pcl_reduce;
    funcs->allreduce = &comm_mpi_pcl_allreduce;
    funcs->allgather = &comm_mpi_pcl_allgather;
    funcs->clean = &comm_mpi_pcl_clean;
    return 0;
  }
} 
#endif
