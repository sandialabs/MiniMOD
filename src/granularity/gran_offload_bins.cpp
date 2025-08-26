/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "gran_offload_bins.h"
#include <stdlib.h>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <pthread.h>


#define DEBUG 0
#define VERBOSE 0

namespace Gran_Offload_Bins
{
  int num_bins = -1;
  Util_Ring::util_ring_data* ring;
  pthread_t* thread;

  void* gran_offload_bins_func(void* args)
  {
    gran_offload_bins_info* info;

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

        gran_offload_bins_req* req = (gran_offload_bins_req*)info->req; 

        if(req->bin_thresh[info->send_dir] == -1){
          err = comm_channel_send(req->comm_req[info->send_dir], info->index);
          if (err) return (void*)-1;
        }
        else
        {
          int capture = 0;

          capture = req->bin_ready[info->send_dir][info->index/req->bin_thresh[info->send_dir]].fetch_add(1) + 1;
          
          if(capture == req->bin_thresh[info->send_dir]) 
          {
            err = comm_channel_send(req->comm_req[info->send_dir], info->index/req->bin_thresh[info->send_dir]);
            if (err) return (void*)-1;
          }
          else if(capture > req->bin_thresh[info->send_dir])
          {
            printf("Error caputure %d greater than thresh hold %d\n", capture,  req->bin_thresh[info->send_dir]);
          }
        }

        
        req->ready_count++;

        free(info);
      }
    }

    return 0;
  }

  int gran_offload_bins_init()
  {
    
    int err = Util_Ring::util_ring_init(&ring);
    if(err)
    {
      printf("ring init error\n");
      return err;
    }

    
    thread = (pthread_t*)malloc(sizeof(pthread_t));
    if (thread == NULL) return -1;
    
    err = pthread_create(thread, NULL, &gran_offload_bins_func, NULL);
    if(err)
    {
      printf("pthread error\n");
      return -1;
    }

    
    if(const char* str = std::getenv("MINIMOD_BINS_ENV"))
    {
      num_bins = atoi(str);

    }
    else
    {
      num_bins = GRAN_BINS_NUMBER;
    }
    return 0;
  }

  int gran_offload_bins_clean()
  {
    return -1; 
  }

  int gran_offload_bins_thread_req()
  {
    return 0;
  }

  int gran_offload_bins_stencil_init(gran_request** gran_req, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int stencil)
  {
    int err; 

    for(int i = 0; i < stencil; i++)
    {
      if(num_entries[i] % num_bins != 0 && num_entries[i] != 1) 
      {
        printf("Error: GRAN_OFFLOAD_BINS module currently doesn't support number_entries %% num_bins != 0\n");
        return -10;
      }
    }

    
    gran_offload_bins_req* req = (gran_offload_bins_req*) malloc(sizeof(gran_offload_bins_req));
    if(req == NULL) return -4; 
 
    
    req->comm_req = (comm_request**) malloc(sizeof(comm_request*) * stencil);

    req->bin_ready =  (std::atomic<int>**)malloc(sizeof(std::atomic<int>*)*stencil); 
    req->bin_thresh = (int*)malloc(sizeof(int)*stencil);
    for(int i = 0; i < stencil; i++)
    {
      if(num_entries[i] == 1){
        err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], num_entries[i], entry_size[i]);
        if (err) exit(-1);

        req->bin_ready[i]  = (std::atomic<int>*) malloc(sizeof(std::atomic<int>)*num_bins);
        req->bin_thresh[i] = -1;
        continue;
      }

      err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], num_bins, (num_entries[i] * entry_size[i]) / num_bins);
      if (err) exit(-1);

      req->bin_ready[i]  = (std::atomic<int>*) malloc(sizeof(std::atomic<int>)*num_bins);
      req->bin_thresh[i] = num_entries[i]/num_bins;
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

  int gran_offload_bins_stencil_start(gran_request* gran_req)
  {
    gran_offload_bins_req* req = (gran_offload_bins_req*) gran_req;

    for(int i = 0; i < req->stencil; i++)
    {
      int err =  comm_channel_start(req->comm_req[i]);
      if(err) return err;
      for(int j = 0; j < num_bins; j++)
      {

        req->bin_ready[i][j].store(0);
      }
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

  int gran_offload_bins_stencil_ready(gran_request* gran_req, int send_dir, int index)
  {
    gran_offload_bins_req* req = (gran_offload_bins_req*) gran_req;

    gran_offload_bins_info* info = (gran_offload_bins_info*)calloc(1, sizeof(gran_offload_bins_info));
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

  int gran_offload_bins_stencil_end(gran_request* gran_req)
  {
    gran_offload_bins_req* req = (gran_offload_bins_req*) gran_req;

    
    while(req->ready_count.load() < req->part_count);

    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_end(req->comm_req[i]);
      if(err) return err;
    }

    if(VERBOSE) 
    {
      int id;
      comm_get_id(&id);
      printf("(%d)gran_offload_bins_stencil_end: completed busy loop!\n",id);
    }

    return 0;
  }

  int gran_offload_bins_stencil_finish(gran_request* gran_req)
  {
    gran_offload_bins_req* req = (gran_offload_bins_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_finalize(req->comm_req[i]);
    }
    free(req->comm_req);

#if 0
    gran_offload_bins_info* tmp;
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

  int gran_offload_bins_load(gran_funcs* funcs)
  {
    funcs->init = &gran_offload_bins_init;
    funcs->clean = &gran_offload_bins_clean;
    funcs->thread_req = &gran_offload_bins_thread_req;
    funcs->stencil_init = &gran_offload_bins_stencil_init;
    funcs->stencil_start = &gran_offload_bins_stencil_start;
    funcs->stencil_ready = &gran_offload_bins_stencil_ready;
    funcs->stencil_end = &gran_offload_bins_stencil_end;
    funcs->stencil_finish = &gran_offload_bins_stencil_finish;
    return 0;
  }
} 
