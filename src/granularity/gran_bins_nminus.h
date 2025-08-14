/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_GRANULARITY_BINS_TMP
#define MINIMOD_GRANULARITY_BINS_TMP

#include <stdio.h>
#include <atomic>
#include "../gran.h"

#define GRAN_BINS_NUMBER 2

namespace Gran_Bins_Nminus
{
  typedef struct
  {
    comm_request** comm_req;
    int stencil;

    std::atomic<int>** bin_ready;
    int* bin_thresh;

    std::atomic<int>* bin_count;
  } gran_bins_nminus_req;

  int gran_bins_nminus_load(gran_funcs*);
} 
#endif
