/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifdef HAVE_MPI

#ifndef MINIMOD_COMM_MPI_PCL
#define MINIMOD_COMM_MPI_PCL

#include <stdio.h>
#include <mpi.h>
#include "comm.h"
#include "mpipcl.h"

namespace Comm_MPI_PCL
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
    MPIX_Request recv_req;
    MPIX_Request send_req;
  } comm_mpi_pcl_req;

  int comm_mpi_pcl_load(comm_funcs*);
} 
#endif

#endif
