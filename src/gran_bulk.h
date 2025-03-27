/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_GRANULARITY_BULK
#define MINIMOD_GRANULARITY_BULK

#include <stdio.h>

#include "gran.h"

namespace Gran_Bulk
{
  typedef struct
  {
    comm_request** comm_req;
    int stencil;
    int* to_send; 
    int* to_send_size;
  } gran_bulk_req;

  int gran_bulk_load(gran_funcs*);
} 
#endif
