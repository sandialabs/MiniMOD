/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_GRANULARITY_FINE
#define MINIMOD_GRANULARITY_FINE

#include <stdio.h>

#include "../gran.h"

namespace Gran_Fine
{
  typedef struct
  {
    comm_request** comm_req;
    int stencil;
  } gran_fine_req;

  int gran_fine_load(gran_funcs*);
} 
#endif
