/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "gran_fine.h"
#include <stdlib.h>


#define DEBUG 0
#define VERBOSE 0

namespace Gran_Fine
{
  int gran_fine_init()
  {
    return 0;
  }

  int gran_fine_clean()
  {
    return -1; 
  }

  int gran_fine_thread_req()
  {
    return 1;
  }

  int gran_fine_stencil_init(gran_request** gran_req, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int stencil)
  {
    
    gran_fine_req* req = (gran_fine_req*) malloc(sizeof(gran_fine_req));
    if(req == NULL) return -4; 
 
    
    req->comm_req = (comm_request**) malloc(sizeof(comm_request*) * stencil);
    int err; 
    for(int i = 0; i < stencil; i++)
    {
      err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], num_entries[i], entry_size[i]);
      if (err) return err;
    }
 
    req->stencil = stencil;
 
    *gran_req = (gran_request*)req;
    return 0;
  }

  int gran_fine_stencil_start(gran_request* gran_req)
  {
    if (gran_req == NULL) {
      printf("%s:%d: ERROR: In function %s, the provided gran_req is NULL\n", __FILE__, __LINE__, __func__);
      exit(-1);
    }
    gran_fine_req* req = (gran_fine_req*) gran_req;
    for(int i = 0; i < req->stencil; i++)
    {
      int err =  comm_channel_start(req->comm_req[i]);
      if(err) return err;
    }
    return 0;
  }

  int gran_fine_stencil_ready(gran_request* gran_req, int send_dir, int entry)
  {
    gran_fine_req* req = (gran_fine_req*) gran_req;
    int err = comm_channel_send(req->comm_req[send_dir], entry);
    return err;
  }

  int gran_fine_stencil_end(gran_request* gran_req)
  {
    gran_fine_req* req = (gran_fine_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_end(req->comm_req[i]);
      if(err) return err;
    }
    return 0;
  }

  int gran_fine_stencil_finish(gran_request* gran_req)
  {
    gran_fine_req* req = (gran_fine_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++)
    {
      err = comm_channel_finalize(req->comm_req[i]);
    }
    free(req->comm_req);
    free(req);

    return err;
  }

  int gran_fine_load(gran_funcs* funcs)
  {
    funcs->init = &gran_fine_init;
    funcs->clean = &gran_fine_clean;
    funcs->thread_req = &gran_fine_thread_req;
    funcs->stencil_init = &gran_fine_stencil_init;
    funcs->stencil_start = &gran_fine_stencil_start;
    funcs->stencil_ready = &gran_fine_stencil_ready;
    funcs->stencil_end = &gran_fine_stencil_end;
    funcs->stencil_finish = &gran_fine_stencil_finish;
    return 0;
  }
} 
