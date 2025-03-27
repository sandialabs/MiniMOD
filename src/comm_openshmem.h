/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_COMM_OPENSHMEM
#define MINIMOD_COMM_OPENSHMEM

#ifdef HAVE_OSHMEM

#include <stdio.h>
#include "comm.h"

namespace Comm_Openshmem
{
  typedef struct
  {
    void* recv_buf;
    void* send_buf;
    int recv_dir;
    int send_dir;
    int recv_id;
    int send_id;
    int num_entries;
    int entry_size;
    int ready;
  } comm_openshmem_req;

  int comm_openshmem_load(comm_funcs*);
} 
#endif

#endif
