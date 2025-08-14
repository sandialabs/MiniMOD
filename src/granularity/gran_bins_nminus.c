/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "gran_bins_nminus.h"
#include <stdlib.h>


#define DEBUG 0
#define VERBOSE 0

namespace Gran_Bins_Nminus {
  int num_bins = -1;

  int gran_bins_nminus_init() {
    if(const char* str = std::getenv("MINIMOD_BINS_ENV")) {
      num_bins = atoi(str);

    }
    else {
      num_bins = GRAN_BINS_NUMBER;
    }
    return 0;
  }

  int gran_bins_nminus_clean() {
    return -1; 
  }

  int gran_bins_nminus_thread_req() {
    return 1;
  }

  int gran_bins_nminus_stencil_init(gran_request** gran_req, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int stencil) {
    for(int i = 0; i < stencil; i++) {
      if(num_entries[i] % num_bins != 0 && num_entries[i] != 1)  {
        printf("Error: GRAN_BINS module currently doesn't support number_entries %% num_bins != 0\n");
        return -10;
      }
    }
    
    gran_bins_nminus_req* req = (gran_bins_nminus_req*) malloc(sizeof(gran_bins_nminus_req));
    if(req == NULL) return -4; 
    req->comm_req = (comm_request**) malloc(sizeof(comm_request*) * stencil);
 

    req->bin_ready =  (std::atomic<int>**)malloc(sizeof(std::atomic<int>*) * stencil); 
    req->bin_thresh = (int*)malloc(sizeof(int)*stencil);

    req->bin_count = (std::atomic<int>*)malloc(sizeof(std::atomic<int>) * stencil);
    
    int err; 
    for(int i = 0; i < stencil; i++) {
      if(num_entries[i] == 1){
        err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], 1, num_entries[i] * entry_size[i]);
        if (err) return err;

        req->bin_ready[i]  = (std::atomic<int>*) malloc(sizeof(std::atomic<int>)*num_bins);
        req->bin_thresh[i] = -1;
	continue;
      }
      err = comm_channel_init(&(req->comm_req[i]), recv_bufs[i], send_bufs[i], recv_ids[i], send_ids[i], num_bins, (num_entries[i] * entry_size[i])/num_bins);
      if (err) return err;

      req->bin_ready[i]  = (std::atomic<int>*) malloc(sizeof(std::atomic<int>)*num_bins);
      req->bin_thresh[i] = num_entries[i]/num_bins;
    }
    req->stencil = stencil;
 
    *gran_req = (gran_request*)req;
    return 0;
  }

  int gran_bins_nminus_stencil_start(gran_request* gran_req) {
    gran_bins_nminus_req* req = (gran_bins_nminus_req*) gran_req;
    for(int i = 0; i < req->stencil; i++) {
      int err =  comm_channel_start(req->comm_req[i]);
      req->bin_count[i] = 0;
      for(int j = 0; j < num_bins; j++) {

        req->bin_ready[i][j].store(0);
      }
      if(err) return err;
    }
    return 0;
  }

  int gran_bins_nminus_stencil_ready(gran_request* gran_req, int send_dir, int entry) {
    int err = 0;
    gran_bins_nminus_req* req = (gran_bins_nminus_req*) gran_req;
 
    
    if(req->bin_thresh[send_dir] == -1){
      return comm_channel_send(req->comm_req[send_dir], entry);
    }

    int bin = entry/req->bin_thresh[send_dir];
    int capture = 0;
 

    capture = req->bin_ready[send_dir][bin].fetch_add(1) + 1;

    
    if(capture == req->bin_thresh[send_dir])  {


      int bc_capture = req->bin_count[send_dir].fetch_add(1) + 1;

      

      if(bc_capture == num_bins - 1) {
        
        int last_bin = -1;
        for(int j = 0; j < num_bins; j++) {
          if(req->bin_ready[send_dir][j].load() != req->bin_thresh[send_dir]) {
            last_bin = j;
	    break;
          }
        }

        
        if (-1 == last_bin) {
          
          err = comm_channel_send_range(req->comm_req[send_dir], 0, num_bins - 1);
          if(err) {
            int id;
            comm_get_id(&id);
            printf("(%d) Error: complete buffer send in bins_nminus: comm_channel_send_range(%p,%d,%d)\n", id, req->comm_req[send_dir], 0, num_bins - 1);
            return -1;
          }
        }
        else {
          err = comm_channel_send_range(req->comm_req[send_dir], 0, (last_bin - 1));
          if(err) {
            int id;
            comm_get_id(&id);
            printf("(%d) Error: first send in bins_nminus: comm_channel_send_range(%p,%d,%d)\n", id, req->comm_req[send_dir], 0, last_bin - 1);
            return -1;
          }
          if(last_bin != num_bins - 1) {
            err = comm_channel_send_range(req->comm_req[send_dir], (last_bin + 1), (num_bins - 1));
            if(err) {
              int id;
              comm_get_id(&id);
              printf("(%d) Error: second send in bins_nminus: comm_channel_send_range(%p,%d,%d)\n", id, req->comm_req[send_dir], last_bin + 1, num_bins - 1);
              return -1;
            }
          }
        }
      }
      

      else if(bc_capture == num_bins) {
        err = comm_channel_send(req->comm_req[send_dir], bin);
        if(err) {
          int id;
          comm_get_id(&id);
          printf("(%d) Error: final send in bins_nminus: comm_channel_send(%p,%d)\n", id, req->comm_req[send_dir], bin);
          return -1;
        }
      }
    }
    else if(capture > req->bin_thresh[send_dir]) {
      printf("Error caputure %d greater than thresh hold %d\n", capture,  req->bin_thresh[send_dir]);
    }
    return 0;
  }

  int gran_bins_nminus_stencil_end(gran_request* gran_req) {
    gran_bins_nminus_req* req = (gran_bins_nminus_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++) {
      err = comm_channel_end(req->comm_req[i]);
      if(err) return err;
    }
    return 0;
  }

  int gran_bins_nminus_stencil_finish(gran_request* gran_req) {
    gran_bins_nminus_req* req = (gran_bins_nminus_req*) gran_req;
    int i, err;
    for (i = 0; i < req->stencil; i++) {
      err = comm_channel_finalize(req->comm_req[i]);
      free(req->bin_ready[i]);
    }
    free(req->comm_req);
    free(req->bin_ready);
    free(req->bin_thresh);
    free(req->bin_count);
    free(req);

    return err;
  }

  int gran_bins_nminus_load(gran_funcs* funcs) {
    funcs->init = &gran_bins_nminus_init;
    funcs->clean = &gran_bins_nminus_clean;
    funcs->thread_req = &gran_bins_nminus_thread_req;
    funcs->stencil_init = &gran_bins_nminus_stencil_init;
    funcs->stencil_start = &gran_bins_nminus_stencil_start;
    funcs->stencil_ready = &gran_bins_nminus_stencil_ready;
    funcs->stencil_end = &gran_bins_nminus_stencil_end;
    funcs->stencil_finish = &gran_bins_nminus_stencil_finish;
    return 0;
  }
} 

