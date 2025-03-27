/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_GRANULARITY_FUNN
#define MINIMOD_GRANULARITY_FUNN

#include <stdio.h>
#include <atomic>

#include "../gran.h"
#include "../utility/util_queue.h"

namespace Gran_Funn
{
  typedef struct
  {
    int send_dir;
    int index;
  } gran_funn_info;

  typedef struct
  {
    comm_request** comm_req;
    int stencil;
    Util_Queue::util_queue_data* queue;
  } gran_funn_req;

  int gran_funn_load(gran_funcs*);
} 
#endif
