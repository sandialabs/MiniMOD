/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "comm.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


comm_funcs* active_comm;

int comm_load(char* name)
{
  active_comm = (comm_funcs*) malloc(sizeof(comm_funcs));
  
  
  
#ifdef HAVE_MPI
  if(strcmp(name,"mpi_pt2pt") == 0)
  {
    using namespace Comm_MPI_PT2PT;
    return comm_mpi_pt2pt_load(active_comm);
  }
#endif
#ifdef HAVE_MPI
  if(strcmp(name,"mpi_pers") == 0)
  {
    using namespace Comm_MPI_Pers;
    return comm_mpi_pers_load(active_comm);
  }
#endif
#ifdef HAVE_MPI
  if(strcmp(name,"mpi_rma") == 0)
  {
    using namespace Comm_MPI_RMA;
    return comm_mpi_rma_load(active_comm);
  }
#endif
#ifdef HAVE_MPI_PART
  if(strcmp(name,"mpi_part") == 0)
  {
    using namespace Comm_MPI_Part;
    return comm_mpi_part_load(active_comm);
  }
#endif  
#ifdef HAVE_MPI_PCL
  if(strcmp(name,"mpi_pcl") == 0)
  {
    using namespace Comm_MPI_PCL;
    return comm_mpi_pcl_load(active_comm);
  }
#endif
#ifdef HAVE_OSHMEM
  if(strcmp(name,"openshmem") == 0)
  {
    using namespace Comm_Openshmem;
    return comm_openshmem_load(active_comm);
  }
#endif
#ifdef HAVE_CRAY_SHMEM
  if(strcmp(name,"crayshmem") == 0)
  {
    using namespace Comm_Crayshmem;
    return comm_crayshmem_load(active_comm);
  }
#endif

  
  return -2;
}

int comm_init(int threads)
{
  return active_comm->init(threads);
}

int comm_get_id(int* id)
{
  return active_comm->get_id(id); 
}

int comm_get_job_size (int* size)
{
  return active_comm->job_size(size);
}

int comm_channel_init(comm_request** comm_req, void** recv_buf, void* send_buf, int recv_id, int send_id, int num_entries, int entry_size)
{
  return active_comm->channel_init(comm_req, recv_buf, send_buf, recv_id, send_id, num_entries, entry_size);
}

int comm_channel_start(comm_request* comm_req)
{
  return active_comm->channel_start(comm_req);
}

int comm_channel_send(comm_request* comm_req, int index)
{
  return active_comm->channel_send(comm_req, index);
}

int comm_channel_send_range(comm_request* comm_req, int index_start, int index_end)
{
  return active_comm->channel_send_range(comm_req, index_start, index_end);
}

int comm_channel_end(comm_request* comm_req)
{
  return active_comm->channel_end(comm_req);
}

int comm_channel_finalize (comm_request* comm_req)
{
  return active_comm->channel_finalize(comm_req);
}

int comm_barrier()
{
  return active_comm->barrier();
}

int comm_broadcast(void* buffer, int count, comm_datatype datatype, int root)
{
  return active_comm->broadcast(buffer, count, datatype, root);
}

int comm_reduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op, int root)
{
  return active_comm->reduce(send_buf, recv_buf, count, datatype, op, root);
}

int comm_allreduce(void* send_buf, void* recv_buf, int count, comm_datatype datatype, comm_reduce_op op)
{
  return active_comm->allreduce(send_buf, recv_buf, count, datatype, op);
}

int comm_allgather(void* send_buf, int send_count, comm_datatype send_datatype, void* recv_buf, int recv_count, comm_datatype recv_datatype)
{
  return active_comm->allgather(send_buf, send_count, send_datatype, recv_buf, recv_count, recv_datatype);
}

int comm_clean()
{
  int err = active_comm->clean();
  free(active_comm);
  return err;
}
