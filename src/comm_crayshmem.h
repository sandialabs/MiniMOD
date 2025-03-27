/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_COMM_CRAYSHMEM
#define MINIMOD_COMM_CRAYSHMEM

#ifdef HAVE_CRAY_SHMEM

#include <stdio.h>
#include "comm.h"

namespace Comm_Crayshmem
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
  } comm_crayshmem_req;

  int comm_crayshmem_load(comm_funcs*);
} 
#endif

#endif
