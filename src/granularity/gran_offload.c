/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "gran_offload.h"
#include <stdlib.h>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <pthread.h>


#define DEBUG 0
#define VERBOSE 0

namespace Gran_Offload
{
  Util_Ring::util_ring_data* ring;
  pthread_t* thread;

  void* gran_offload_func(void* args)
  {
    gran_offload_info* info;

    int err;
    while(1)
    {
      Util_Ring::util_ring_get_event(ring, (void**)(&info));
      
      if(NULL != info)
      {
        if(VERBOSE) 
        {
          int id;
          comm_get_id(&id);
          printf("(%d)util_ring_get_event<dir:%d, ind:%d, rdy: %d>\n",id,info->send_dir,info->index,info->req->ready_count.load());
        }
 
        err = comm_channel_send(info->req->comm_req[info->send_dir], info->index);
        if (err) return (void*)-1;
        
        
        info->req->ready_count++;


        free(info);
      }
    }

    return 0;
  }

  int gran_offload_init()
  {
    
    int err = Util_Ring::util_ring_init(&ring);
    if(err)
    {
      printf("ring init error\n");
      return err;
    }

    
    thread = (pthread_t*)malloc(sizeof(pthread_t));
    if (thread == NULL) return -1;
    
    err = pthread_create(thread, NULL, &gran_offload_func, NULL);
    if(err)
    {
      printf("pthread error\n");
      return -1;
    }
    return 0;
  }

  int gran_offload_clean()
  {
    return -1; 
  }

  int gran_offload_thread_req()
  {
    return 0;
  }

  int gran_offload_stencil_init(gran_request** gran_req, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int stencil)
  {
    int err; 

    
    gran_offload_req* req = (gran_offload_req*) malloc(sizeof(gran_offload_req));
    if(req == NULL) return -4; 
 
    
    req->comm_req = (comm_request**) malloc(sizeof(comm_request*) * stencil);
    for(int i = 0; i < stencil; i++)
    {
      err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], num_entries[i], entry_size[i]);
      if (err) exit(-1);
    }

    
    req->part_count = 0;
    for(int i = 0; i < stencil; i++)
    {







      req->part_count += num_entries[i];
    }

    req->stencil = stencil;
 
    *gran_req = (gran_request*)req;
    
    return 0;
  }

  int gran_offload_stencil_start(gran_request* gran_req)
  {
    if (gran_req == NULL) {
      printf("%s:%d: ERROR: In function %s, the provided gran_req is NULL\n", __FILE__, __LINE__, __func__);
      exit(-1);
    }

    gran_offload_req* req = (gran_offload_req*) gran_req;

    for(int i = 0; i < req->stencil; i++)
    {
      int err =  comm_channel_start(req->comm_req[i]);
      if(err) return err;
    }

    
    req->ready_count.store(0);

    if(VERBOSE) 
    {
      int id;
      comm_get_id(&id);
      printf("(%d)[starting iteration] <part_count:%d>\n", id, req->part_count);
    }

    return 0;
  }

  int gran_offload_stencil_ready(gran_request* gran_req, int send_dir, int index)
  {
    gran_offload_req* req = (gran_offload_req*) gran_req;

    gran_offload_info* info = (gran_offload_info*)calloc(1, sizeof(gran_offload_info));
    if (info == NULL) return -1;
    info->req = req;
    info->send_dir = send_dir;
    info->index = index;

    if(VERBOSE) 
    {
      int id;
      comm_get_id(&id);
      printf("(%d)util_ring_add_event<dir:%d,ind:%d>\n",id,send_dir,index);
    }

    Util_Ring::util_ring_add_event(ring, info);
    return 0;
  }

  int gran_offload_stencil_end(gran_request* gran_req)
  {
    gran_offload_req* req = (gran_offload_req*) gran_req;

    
    while(req->ready_count.load() < req->part_count);

    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_end(req->comm_req[i]);
      if(err) return err;
    }

    return 0;
  }

  int gran_offload_stencil_finish(gran_request* gran_req)
  {
    gran_offload_req* req = (gran_offload_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_finalize(req->comm_req[i]);
    }
    free(req->comm_req);

#if 0
    gran_offload_info* tmp;
    Util_Ring::util_ring_get_event(ring, (void**)(&tmp));
    while (tmp != NULL)
    {
      free(tmp);
      Util_Ring::util_ring_get_event(ring, (void**)(&tmp));
    }
#endif

    free(req);

    return err;
  }

  int gran_offload_load(gran_funcs* funcs)
  {
    funcs->init = &gran_offload_init;
    funcs->clean = &gran_offload_clean;
    funcs->thread_req = &gran_offload_thread_req;
    funcs->stencil_init = &gran_offload_stencil_init;
    funcs->stencil_start = &gran_offload_stencil_start;
    funcs->stencil_ready = &gran_offload_stencil_ready;
    funcs->stencil_end = &gran_offload_stencil_end;
    funcs->stencil_finish = &gran_offload_stencil_finish;
    return 0;
  }
} 
