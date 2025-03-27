/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_GRANULARITY_OFFLOAD_BINS
#define MINIMOD_GRANULARITY_OFFLOAD_BINS

#include <stdio.h>
#include <atomic>

#include "gran.h"
#include "util_ring.h"

namespace Gran_Offload_Bins
{
  typedef struct
  {
    comm_request** comm_req;
    int stencil;

    std::atomic<int>** bin_ready;
    int* bin_thresh;
    std::atomic<int> ready_count;
    int part_count;
  } gran_offload_bins_req;

  typedef struct
  {
    gran_offload_bins_req* req;
    int send_dir;
    int index;
  } gran_offload_bins_info;

  int gran_offload_bins_load(gran_funcs*);
} 
#endif
