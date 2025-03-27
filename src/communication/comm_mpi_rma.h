/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_COMM_MPI_RMA
#define MINIMOD_COMM_MPI_RMA

#ifdef HAVE_MPI

#include <stdio.h>
#include <mpi.h>
#include "../comm.h"

namespace Comm_MPI_RMA
{
  typedef struct
  {
    void* recv_buf;
    void* send_buf;
    MPI_Win win; 
    int recv_dir;
    int send_dir;
    int recv_id;
    int send_id;
    int my_id;
    int num_entries;
    int entry_size;
    int ready;
  } comm_mpi_rma_req;

  int comm_mpi_rma_load(comm_funcs*);
} 
#endif
#endif
