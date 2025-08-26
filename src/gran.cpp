/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "gran.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


#define DEBUG 0
#define VERBOSE 0

gran_funcs* active_gran;

int gran_load(char* name) {
  active_gran = (gran_funcs*) malloc(sizeof(gran_funcs));
  
  
  
  if(strcmp(name,"bulk") == 0) {
    using namespace Gran_Bulk;
    return gran_bulk_load(active_gran);
  }

  if(strcmp(name,"fine") == 0) {
    using namespace Gran_Fine;
    return gran_fine_load(active_gran);
  }

  if(strcmp(name,"bins") == 0) {
    using namespace Gran_Bins;
    return gran_bins_load(active_gran);
  }

  if(strcmp(name,"bins_tmp") == 0) {
    using namespace Gran_Bins_Nminus;
    return gran_bins_nminus_load(active_gran);
  }

  if(strcmp(name,"offload") == 0) {
    using namespace Gran_Offload;
    return gran_offload_load(active_gran);
  }

  if(strcmp(name,"offload_bins") == 0) {
    using namespace Gran_Offload_Bins;
    return gran_offload_bins_load(active_gran);
  }

  if(strcmp(name,"offload_bins_nminus") == 0) {
    using namespace Gran_Offload_Bins_Nminus;
    return gran_offload_bins_nminus_load(active_gran);
  }

  
  return -2;
}

int gran_init() {
  return active_gran->init();
}

int gran_thread_req() {
  return active_gran->thread_req();
}

int gran_stencil_init(gran_request** gran_req, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int stencil) {
  return active_gran->stencil_init(gran_req, recv_bufs, send_bufs, recv_ids, send_ids, num_entries, entry_size, stencil);
}

int gran_stencil_start(gran_request* gran_req) {
  return active_gran->stencil_start(gran_req);
}

int gran_stencil_ready(gran_request* gran_req, int send_dir, int index) {
  return active_gran->stencil_ready(gran_req, send_dir, index);
}

int gran_stencil_end(gran_request* gran_req) {
  return active_gran->stencil_end(gran_req);
}

int gran_stencil_finish(gran_request* gran_req) {
  return active_gran->stencil_finish(gran_req);
}
