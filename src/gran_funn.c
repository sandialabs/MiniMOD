/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "gran_funn.h"
#include <stdlib.h>
#include <cstring>


#define DEBUG 0
#define VERBOSE 0

namespace Gran_Funn
{
  int gran_funn_init()
  {
    return 0;
  }

  int gran_funn_clean()
  {
    return -1; 
  }

  int gran_funn_thread_req()
  {
    return 1;
  }

  int gran_funn_stencil_init(gran_request** gran_req, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int stencil)
  {
    
    gran_funn_req* req = (gran_funn_req*) malloc(sizeof(gran_funn_req));
    if(req == NULL) return -4; 
    req->comm_req = (comm_request**) malloc(sizeof(comm_request*) * stencil);
 
    
    int err; 
    for(int i = 0; i < stencil; i++)
    {
      err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], num_entries[i], entry_size[i]);
      if (err) return err;
    }

    err = Util_Queue::util_queue_init(&(req->queue));
    if (err) return err;

    req->stencil = stencil;
 
    *gran_req = (gran_request*)req;
    return 0;
  }

  int gran_funn_stencil_start(gran_request* gran_req)
  {
    if (gran_req == NULL) {
      printf("%s:%d: ERROR: In function %s, the provided gran_req is NULL\n", __FILE__, __LINE__, __func__);
      exit(-1);
    }
    gran_funn_req* req = (gran_funn_req*) gran_req;
    for(int i = 0; i < req->stencil; i++)
    {
      int err =  comm_channel_start(req->comm_req[i]);
      if(err) return err;
    }
    return 0;
  }

  int gran_funn_stencil_ready(gran_request* gran_req, int send_dir, int index)
  {
    gran_funn_req* req = (gran_funn_req*) gran_req;

    gran_funn_info* info = (gran_funn_info*)malloc(sizeof(gran_funn_info));
    if (info == NULL) return -1;
    info->send_dir = send_dir;
    info->index = index;

    Util_Queue::util_queue_event* event;
    int err = Util_Queue::util_queue_create_event(&event, info);
    if(err) return err;

    if(VERBOSE) 
    {
      int id;
      comm_get_id(&id);
      printf("(%d)util_queue_add_event<dir:%d,ind:%d>\n",id,send_dir,index);
    }

    Util_Queue::util_queue_add_event(req->queue, event);
    return 0;
  }








  int gran_funn_stencil_end(gran_request* gran_req)
  {
    gran_funn_req* req = (gran_funn_req*) gran_req;

    Util_Queue::util_queue_event* event;
    Util_Queue::util_queue_get_event(req->queue, &event);
    int i, err;
    while(event != NULL)
    {
      gran_funn_info* info = (gran_funn_info*)event->info;

      if(VERBOSE) 
      {
        int id;
        comm_get_id(&id);
        printf("(%d)util_queue_get_event<dir:%d, ind:%d, nxt:%p>\n",id,info->send_dir,info->index, event->next);
      }

      err = comm_channel_send(req->comm_req[info->send_dir], info->index);
      if (err) return err;

      free(info);
      free(event);

      Util_Queue::util_queue_get_event(req->queue, &event);
    }

    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_end(req->comm_req[i]);
      if(err) return err;
    }
    return 0;
  }

  int gran_funn_stencil_finish(gran_request* gran_req)
  {
    gran_funn_req* req = (gran_funn_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_finalize(req->comm_req[i]);
    }
    free(req->comm_req);

    Util_Queue::util_queue_event* tmp;
    Util_Queue::util_queue_get_event(req->queue, &tmp);
    while (tmp != NULL)
    {
      free(tmp);
      Util_Queue::util_queue_get_event(req->queue, &tmp);
    }

    free(req);

    return err;
  }

  int gran_funn_load(gran_funcs* funcs)
  {
    funcs->init = &gran_funn_init;
    funcs->clean = &gran_funn_clean;
    funcs->thread_req = &gran_funn_thread_req;
    funcs->stencil_init = &gran_funn_stencil_init;
    funcs->stencil_start = &gran_funn_stencil_start;
    funcs->stencil_ready = &gran_funn_stencil_ready;
    funcs->stencil_end = &gran_funn_stencil_end;
    funcs->stencil_finish = &gran_funn_stencil_finish;
    return 0;
  }
} 
