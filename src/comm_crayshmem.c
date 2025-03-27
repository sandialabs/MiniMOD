/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifdef HAVE_CRAY_SHMEM

#include "comm_crayshmem.h"

#include <shmem.h>
#include <shmemx.h>
#include <stdlib.h>
#include <string.h>


#define DEBUG 0 
#define VERBOSE 0
#if 1
namespace Comm_Crayshmem
{
  int comm_crayshmem_init(int threads)
  {
    if(threads)
    {
      int provided;
      int required = SHMEM_THREAD_MULTIPLE;
      provided = shmemx_init_thread(required);

      if(provided == required)
      {
        return 0;
      }
      else
      {
        return -9;
      }
    }
    else
    {
      shmem_init();
      return 0;
    }
  }

  int comm_crayshmem_get_id(int* id)
  {
    *id = shmem_my_pe();
    return 0;
  }

  int comm_crayshmem_job_size(int* size)
  {
    *size = shmem_n_pes();
    return 0;
  }

  
  

  int comm_crayshmem_buffer_source_size = 0;
  int comm_crayshmem_buffer_dest_size = 0;
  void* comm_crayshmem_buffer_source = 0;
  void* comm_crayshmem_buffer_dest   = 0;
  long* comm_crayshmem_broadcast_psynch = 0;
  int comm_crayshmem_pwrk_size = 0;
  long* comm_crayshmem_pwrk = 0;

  int comm_crayshmem_clean()
  {
    
    
    if( comm_crayshmem_pwrk != nullptr ) {
      shmem_free(comm_crayshmem_pwrk);
    }
    if ( comm_crayshmem_broadcast_psynch != nullptr ) {
      shmem_free(comm_crayshmem_broadcast_psynch);
    }
    if ( comm_crayshmem_buffer_dest != nullptr ) {
      shmem_free(comm_crayshmem_buffer_dest); 
    }
    if ( comm_crayshmem_buffer_source != nullptr ) {
      shmem_free(comm_crayshmem_buffer_source);
    }
    shmem_finalize();
    return 0; 
  }

  int comm_crayshmem_barrier()
  {
    shmem_barrier_all();
    return 0;
  }

  int comm_crayshmem_broadcast(void* buffer, int count, comm_datatype datatype, int root)
  {
    if(comm_crayshmem_broadcast_psynch == 0)
    {
      comm_crayshmem_broadcast_psynch = (long*)shmem_malloc(sizeof(long)*SHMEM_BCAST_SYNC_SIZE);
    }

    for(int i = 0 ; i < SHMEM_BCAST_SYNC_SIZE; i++) comm_crayshmem_broadcast_psynch[i] = SHMEM_SYNC_VALUE;

    shmem_barrier_all(); 

    switch(datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:
        











        if(comm_crayshmem_buffer_source_size < count * 4)
        {
          comm_crayshmem_buffer_source_size = count * 4;
          if(comm_crayshmem_buffer_source != 0) shmem_free(comm_crayshmem_buffer_source); 
          comm_crayshmem_buffer_source = (void*)shmem_malloc(comm_crayshmem_buffer_source_size);
        }
        if(comm_crayshmem_buffer_dest_size < count * 4)
        {
          comm_crayshmem_buffer_dest_size = count * 4;
          if(comm_crayshmem_buffer_dest != 0) shmem_free(comm_crayshmem_buffer_dest); 
          comm_crayshmem_buffer_dest = (void*)shmem_malloc(comm_crayshmem_buffer_dest_size);
        }

        
        if(root == shmem_my_pe())
        {
          memcpy(comm_crayshmem_buffer_source, buffer, count * 4);
          shmem_broadcast32(comm_crayshmem_buffer_dest, comm_crayshmem_buffer_source, count, root, 0, 0, shmem_n_pes(), comm_crayshmem_broadcast_psynch);
        }
        else
        {
          shmem_broadcast32(comm_crayshmem_buffer_dest, comm_crayshmem_buffer_source, count, root, 0, 0, shmem_n_pes(), comm_crayshmem_broadcast_psynch);
          memcpy(buffer,  comm_crayshmem_buffer_dest, count * 4);
        }
        return 0;
      case COMM_DATATYPE_DOUBLE:











        if(comm_crayshmem_buffer_source_size < count * 8)
        {
          comm_crayshmem_buffer_source_size = count * 8;
          if(comm_crayshmem_buffer_source != 0) shmem_free(comm_crayshmem_buffer_source); 
          comm_crayshmem_buffer_source = (void*)shmem_malloc(comm_crayshmem_buffer_source_size);
        }
        if(comm_crayshmem_buffer_dest_size < count * 8)
        {
          comm_crayshmem_buffer_dest_size = count * 8;
          if(comm_crayshmem_buffer_dest != 0) shmem_free(comm_crayshmem_buffer_dest); 
          comm_crayshmem_buffer_dest = (void*)shmem_malloc(comm_crayshmem_buffer_dest_size);
        }

       if(root == shmem_my_pe())
        {
          memcpy(comm_crayshmem_buffer_source, buffer, count * 8);
          shmem_broadcast64(comm_crayshmem_buffer_dest, comm_crayshmem_buffer_source, count, root, 0, 0, shmem_n_pes(), comm_crayshmem_broadcast_psynch);
        }
        else
        {
          shmem_broadcast64(comm_crayshmem_buffer_dest, comm_crayshmem_buffer_source, count, root, 0, 0, shmem_n_pes(), comm_crayshmem_broadcast_psynch);
          memcpy(buffer, comm_crayshmem_buffer_dest, count * 8);
        }
        return 0;
    }
    return -1;
  }

  int comm_crayshmem_allreduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op)
  {
    
    if(comm_crayshmem_broadcast_psynch == 0)
    {
      comm_crayshmem_broadcast_psynch = (long*)shmem_malloc(sizeof(long)*SHMEM_BCAST_SYNC_SIZE);
    }

    
    if(comm_crayshmem_pwrk_size < SHMEM_REDUCE_MIN_WRKDATA_SIZE || comm_crayshmem_pwrk_size < (count/2 + 1))
    {
      int comm_crayshmem_pwrk_size = (count/2) + 1;
      if(comm_crayshmem_pwrk_size < SHMEM_REDUCE_MIN_WRKDATA_SIZE)
      {
        comm_crayshmem_pwrk_size = SHMEM_REDUCE_MIN_WRKDATA_SIZE;
      }
      if( comm_crayshmem_pwrk != nullptr ) {
        shmem_free(comm_crayshmem_pwrk);
      }
      comm_crayshmem_pwrk = (long*)shmem_malloc(sizeof(long)*comm_crayshmem_pwrk_size);
    }

    
    for(int i = 0; i < SHMEM_BCAST_SYNC_SIZE; i++) comm_crayshmem_broadcast_psynch[i] = SHMEM_SYNC_VALUE;
 
    shmem_barrier_all(); 

    
    switch(datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:











        if(comm_crayshmem_buffer_source_size < count * 4)
        {
          comm_crayshmem_buffer_source_size = count * 4;
          if(comm_crayshmem_buffer_source != 0) shmem_free(comm_crayshmem_buffer_source); 
          comm_crayshmem_buffer_source = (void*)shmem_malloc(comm_crayshmem_buffer_source_size);
        }
        if(comm_crayshmem_buffer_dest_size < count * 4)
        {
          comm_crayshmem_buffer_dest_size = count * 4;
          if(comm_crayshmem_buffer_dest != 0) shmem_free(comm_crayshmem_buffer_dest); 
          comm_crayshmem_buffer_dest = (void*)shmem_malloc(comm_crayshmem_buffer_dest_size);
        }
        memcpy(comm_crayshmem_buffer_source, send_buf, count * 4);
        break;
      case COMM_DATATYPE_LONGLONG:
      case COMM_DATATYPE_DOUBLE:











        if(comm_crayshmem_buffer_source_size < count * 8)
        {
          comm_crayshmem_buffer_source_size = count * 8;
          if(comm_crayshmem_buffer_source != 0) shmem_free(comm_crayshmem_buffer_source); 
          comm_crayshmem_buffer_source = (void*)shmem_malloc(comm_crayshmem_buffer_source_size);
        }
        if(comm_crayshmem_buffer_dest_size < count * 8)
        {
          comm_crayshmem_buffer_dest_size = count * 8;
          if(comm_crayshmem_buffer_dest != 0) shmem_free(comm_crayshmem_buffer_dest); 
          comm_crayshmem_buffer_dest = (void*)shmem_malloc(comm_crayshmem_buffer_dest_size);
        }
        memcpy(comm_crayshmem_buffer_source, send_buf, count * 8);
        break;
    }

    switch(op)
    {
      case COMM_REDUCE_OP_SUM:
        switch(datatype)
        {
          case COMM_DATATYPE_INT:
            shmem_int_sum_to_all((int*)comm_crayshmem_buffer_dest, (int*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(int*)comm_crayshmem_pwrk,comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_LONGLONG:
            shmem_longlong_sum_to_all((long long*)comm_crayshmem_buffer_dest, (long long*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(long long*)comm_crayshmem_pwrk,comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_FLOAT:
            shmem_float_sum_to_all((float*)comm_crayshmem_buffer_dest, (float*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(), (float*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_DOUBLE:
            shmem_double_sum_to_all((double*)comm_crayshmem_buffer_dest, (double*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(), (double*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
        }
        break;
      case COMM_REDUCE_OP_MAX:
        switch(datatype)
        {
          case COMM_DATATYPE_INT:
            shmem_int_max_to_all((int*)comm_crayshmem_buffer_dest, (int*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(int*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_LONGLONG:
            shmem_longlong_max_to_all((long long*)comm_crayshmem_buffer_dest, (long long*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(long long*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_FLOAT:
            shmem_float_max_to_all((float*)comm_crayshmem_buffer_dest, (float*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(), (float*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_DOUBLE:
            shmem_double_max_to_all((double*)comm_crayshmem_buffer_dest, (double*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(), (double*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
        }
        break;
      case COMM_REDUCE_OP_MIN:
        switch(datatype)
        {
          case COMM_DATATYPE_INT:
            shmem_int_min_to_all((int*)comm_crayshmem_buffer_dest, (int*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(int*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_LONGLONG:
            shmem_longlong_min_to_all((long long*)comm_crayshmem_buffer_dest, (long long*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(long long*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_FLOAT:
            shmem_float_min_to_all((float*)comm_crayshmem_buffer_dest, (float*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(),(float*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break;
          case COMM_DATATYPE_DOUBLE:
            shmem_double_min_to_all((double*)comm_crayshmem_buffer_dest, (double*)comm_crayshmem_buffer_source, count, 0, 0, shmem_n_pes(), (double*)comm_crayshmem_pwrk, comm_crayshmem_broadcast_psynch);
            break; 
        }
        break;
    }

    
    switch(datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:
        memcpy(recv_buf, comm_crayshmem_buffer_dest, count * 4);
        break;
      case COMM_DATATYPE_LONGLONG:
      case COMM_DATATYPE_DOUBLE:
        memcpy(recv_buf, comm_crayshmem_buffer_dest, count * 8);
        break;
    }

    return 0; 
  }

  int comm_crayshmem_reduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op, int root){
    return comm_crayshmem_allreduce(
      send_buf,
      recv_buf, 
      count, 
      datatype, 
      op
    );
  }

  int comm_crayshmem_allgather(void* send_buf, int send_count, comm_datatype send_datatype, void* recv_buf, int recv_count, comm_datatype recv_datatype)
  {
    if(comm_crayshmem_broadcast_psynch == 0)
    {
      comm_crayshmem_broadcast_psynch = (long*)shmem_malloc(sizeof(long)*SHMEM_BCAST_SYNC_SIZE);
    }

    for(int i = 0 ; i < SHMEM_BCAST_SYNC_SIZE; i++) comm_crayshmem_broadcast_psynch[i] = SHMEM_SYNC_VALUE;

    shmem_barrier_all(); 

    
    switch(send_datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:
        if(comm_crayshmem_buffer_source_size < send_count * 4)
        {
          comm_crayshmem_buffer_source_size = send_count * 4;
          if(comm_crayshmem_buffer_source != 0) shmem_free(comm_crayshmem_buffer_source); 
          comm_crayshmem_buffer_source = (void*)shmem_malloc(comm_crayshmem_buffer_source_size);
        }
        break;
      case COMM_DATATYPE_DOUBLE:
        if(comm_crayshmem_buffer_source_size < send_count * 8)
        {
          comm_crayshmem_buffer_source_size = send_count * 8;
          if(comm_crayshmem_buffer_source != 0) shmem_free(comm_crayshmem_buffer_source); 
          comm_crayshmem_buffer_source = (void*)shmem_malloc(comm_crayshmem_buffer_source_size);
        }
        break;
    }
    switch(recv_datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:
        if(comm_crayshmem_buffer_dest_size < shmem_n_pes() * recv_count * 4)
        {
          comm_crayshmem_buffer_dest_size = shmem_n_pes() * recv_count * 4;
          if(comm_crayshmem_buffer_dest != 0) shmem_free(comm_crayshmem_buffer_dest); 
          comm_crayshmem_buffer_dest = (void*)shmem_malloc(comm_crayshmem_buffer_dest_size);
        }
        break;
      case COMM_DATATYPE_DOUBLE:
        if(comm_crayshmem_buffer_dest_size < shmem_n_pes() * recv_count * 8)
        {
          comm_crayshmem_buffer_dest_size = shmem_n_pes() * recv_count * 8;
          if(comm_crayshmem_buffer_dest != 0) shmem_free(comm_crayshmem_buffer_dest); 
          comm_crayshmem_buffer_dest = (void*)shmem_malloc(comm_crayshmem_buffer_dest_size);
        }
        break;
    }

    
    switch(send_datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:
        memcpy(comm_crayshmem_buffer_source, send_buf, send_count * 4);
        shmem_fcollect32(comm_crayshmem_buffer_dest, comm_crayshmem_buffer_source, send_count, 0, 0, shmem_n_pes(), comm_crayshmem_broadcast_psynch);
      break;
      case COMM_DATATYPE_DOUBLE:
        memcpy(comm_crayshmem_buffer_source, send_buf, send_count * 8);
        shmem_fcollect64(comm_crayshmem_buffer_dest, comm_crayshmem_buffer_source, send_count, 0, 0, shmem_n_pes(), comm_crayshmem_broadcast_psynch);
      break;
    }
    switch(recv_datatype)
    {
      case COMM_DATATYPE_INT:
      case COMM_DATATYPE_FLOAT:
        memcpy(recv_buf, comm_crayshmem_buffer_dest, shmem_n_pes() * recv_count * 4);
        return 0;
      case COMM_DATATYPE_DOUBLE:
        memcpy(recv_buf, comm_crayshmem_buffer_dest, shmem_n_pes() * recv_count * 8);
        return 0;
    }
    return -1;
  }

  int comm_crayshmem_channel_init(comm_request** comm_req, void** recv_buf, void* send_buf, int recv_id, int send_id, int num_entries, int entry_size)
  { 
    
    comm_crayshmem_req* req = (comm_crayshmem_req*)malloc(sizeof(comm_crayshmem_req)); 
    if(req == NULL) return -4; 
 
    
    req->recv_buf = shmem_malloc(num_entries * entry_size);
    if(req->recv_buf == NULL) return -4; 
 
    
    req->send_buf = send_buf;
    req->recv_id = recv_id;
    req->send_id = send_id;
    req->num_entries = num_entries;
    req->entry_size = entry_size;
    req->ready = 0;
 
    *comm_req = req;
    *recv_buf = req->recv_buf;
    return 0;
  }

  int comm_crayshmem_channel_start(comm_request* comm_req)
  {
    comm_crayshmem_barrier();
    return 0;
  }

  int comm_crayshmem_channel_send(comm_request* comm_req, int index)
  {
    comm_crayshmem_req* req = (comm_crayshmem_req*)comm_req;  
    if(req->send_id != -1)
    {

      shmem_putmem(((unsigned char*)req->recv_buf) + (index * req->entry_size), ((unsigned char*)req->send_buf) + (index * req->entry_size), req->entry_size, req->send_id);


    }
    return 0;
  }

  int comm_crayshmem_channel_end(comm_request* comm_req)
  {

    comm_crayshmem_barrier();
    return 0;
  }

  int comm_crayshmem_channel_finalize(comm_request* comm_req)
  {
    comm_crayshmem_req* req = (comm_crayshmem_req*)comm_req;
    shmem_free(req->recv_buf);
    free(req);
    return 0;
  }

  
  int comm_crayshmem_channel_start_unknown (comm_request* comm_req)
  {
    comm_crayshmem_barrier();
    return 0;
  }

  int comm_crayshmem_channel_send_unknown (comm_request* comm_req, int index, int entry_size)
  {
    comm_crayshmem_req* req = (comm_crayshmem_req*)comm_req;
    if(req->send_id != -1)
    {
      shmem_putmem(((unsigned char*)req->recv_buf) + (index * entry_size), ((unsigned char*)req->send_buf) + (index * entry_size), entry_size, req->send_id);
    }
    return 0;
  }

  int comm_crayshmem_channel_match_unknown (comm_request*, int, int)
  {
    return 0; 
  }

  int comm_crayshmem_channel_end_unknown (comm_request* comm_req, int, int)
  {
    comm_crayshmem_barrier();
    return 0;
  }

#endif

  int comm_crayshmem_load(comm_funcs* funcs)
  {
    funcs->init = &comm_crayshmem_init;
    funcs->get_id = &comm_crayshmem_get_id;
    funcs->job_size = &comm_crayshmem_job_size;
    funcs->channel_init = &comm_crayshmem_channel_init;
    funcs->channel_start = &comm_crayshmem_channel_start;
    funcs->channel_send = &comm_crayshmem_channel_send;
    funcs->channel_finalize = &comm_crayshmem_channel_finalize;
    funcs->channel_end = &comm_crayshmem_channel_end;
    funcs->channel_start_unknown = &comm_crayshmem_channel_start_unknown;
    funcs->channel_send_unknown = &comm_crayshmem_channel_send_unknown;
    funcs->channel_match_unknown = &comm_crayshmem_channel_match_unknown;
    funcs->channel_end_unknown = &comm_crayshmem_channel_end_unknown;
    funcs->barrier = &comm_crayshmem_barrier;
    funcs->broadcast = &comm_crayshmem_broadcast;
    funcs->reduce = &comm_crayshmem_reduce;
    funcs->allreduce = &comm_crayshmem_allreduce;
    funcs->allgather = &comm_crayshmem_allgather;
    funcs->clean = &comm_crayshmem_clean;
    return 0;
  }
} 
#endif
