/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifdef HAVE_MPI

#ifndef MINIMOD_COMM_MPI_PERS
#define MINIMOD_COMM_MPI_PERS

#include <stdio.h>
#include <mpi.h>
#include "comm.h"

namespace Comm_MPI_Pers
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
    MPI_Request* recv_reqs;
    MPI_Request* send_reqs;
  } comm_mpi_pers_req;

  int comm_mpi_pers_load(comm_funcs*);
} 
#endif

#endif
