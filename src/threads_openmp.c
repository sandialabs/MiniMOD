/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "threads_openmp.h"
#include <stdlib.h>
#include <omp.h>


#define DEBUG 0
#define VERBOSE 0

#ifdef HAVE_OMP

namespace Threads_OpenMP 
{
  int threads_openmp_clean()
  {
    return -1; 
  }

  int threads_openmp_init(threads_request** req, int num_threads)
  {
    threads_openmp_request* _req = (threads_openmp_request*)malloc(sizeof(threads_openmp_request));
    _req->num_threads     = num_threads;
    omp_init_lock (&_req->lck);
    *req = (threads_request*) _req;

    return 0;
  }

  int threads_openmp_run(threads_request* request, int (*func)(int thread_num, int num_threads, void* arg_struct), void* args) 
  {
    int ret = 0;
    threads_openmp_request* req = (threads_openmp_request*) request;
    if(req->num_threads) omp_set_num_threads(req->num_threads);
    else omp_set_num_threads(omp_get_max_threads());  
    #pragma omp parallel
    {
       int _ret = (*func)(omp_get_thread_num(), omp_get_num_threads(), args);
       if(_ret) ret = _ret;
    }
    return ret;
  }

  int threads_openmp_synch(threads_request* req)
  {
    #pragma omp barrier
    return 0; 
  }

  int threads_openmp_begin_critical(threads_request* request)
  {
    threads_openmp_request* req = (threads_openmp_request*) request;
    omp_set_lock(&req->lck);
    return 0;
  }

  int threads_openmp_end_critical(threads_request* request)
  {
    threads_openmp_request* req = (threads_openmp_request*) request;
    omp_unset_lock(&req->lck);
    return 0;
  }


  int threads_openmp_load(threads_funcs* funcs)
  {
    funcs->clean = &threads_openmp_clean;
    funcs->init  = &threads_openmp_init;
    funcs->run   = &threads_openmp_run;
    funcs->synch = &threads_openmp_synch;
    funcs->begin_critical = &threads_openmp_begin_critical;
    funcs->end_critical = &threads_openmp_end_critical;
    return 0;
  }
} 

#endif
