/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "threads.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


threads_funcs* active_threads;

int threads_load(char* name)
{
  active_threads = (threads_funcs*) malloc(sizeof(threads_funcs));
  if(strcmp(name,"pthreads") == 0)
  {
    using namespace Threads_Pthreads;
    return threads_pthreads_load(active_threads);
  }
#ifdef HAVE_OMP
  if(strcmp(name,"openmp") == 0)
  {
    using namespace Threads_OpenMP;
    return threads_openmp_load(active_threads);
  }
#endif

  printf(name);

  printf("COMM - Not found\n");
  return -2;
}


int threads_init (threads_request** request, int num_threads)
{
  return active_threads->init(request, num_threads);
}

int threads_run (threads_request* req, int (*func) (int thread_num, int num_threads, void* arg_struct), void* args)
{
  return active_threads->run(req, func, args);
}

int threads_synch (threads_request* req)
{
  return active_threads->synch(req);
}

int threads_begin_critical (threads_request* req)
{
  return active_threads->begin_critical(req);
}

int threads_end_critical (threads_request* req)
{
  return active_threads->end_critical(req);
}

int threads_clean()
{
  int err = active_threads->clean();
  free(active_threads);
  return err;
}
