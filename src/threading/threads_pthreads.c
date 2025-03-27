/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "threads_pthreads.h"
#include <stdlib.h>
#include <omp.h>


#define DEBUG 0
#define VERBOSE 0

namespace Threads_Pthreads 
{
  int threads_pthreads_clean()
  {
    return -1; 
  }

  int threads_pthreads_init(threads_request** req, int num_threads)
  { 
    threads_pthreads_request* _req = (threads_pthreads_request*)malloc(sizeof(threads_pthreads_request));
    _req->num_threads     = num_threads;
    _req->threads         = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
    _req->thread_arg_bins = (threads_pthreads_arg_wrapper*) malloc(sizeof(threads_pthreads_arg_wrapper)*num_threads);
    pthread_barrier_init (&_req->barrier, NULL, num_threads);
    pthread_mutex_init (&_req->mutex, NULL);
    *req = (threads_request*) _req; 
    return 0;
  }

  void * threads_pthreads_func_wrapper(void* args)
  {
    threads_pthreads_arg_wrapper* _args = (threads_pthreads_arg_wrapper*)args;
    (*(_args->func))(_args->thread_num, _args->num_threads, _args->args);
    return 0;
  }

  int threads_pthreads_run(threads_request* request, int (*func)(int thread_num, int num_threads, void* arg_struct), void* args) 
  {
    int ret = 0;
    threads_pthreads_request* req = (threads_pthreads_request*) request;
    for(int i = 0; i < req->num_threads; i++)
    {
      req->thread_arg_bins[i].thread_num = i;
      req->thread_arg_bins[i].num_threads = req->num_threads;
      req->thread_arg_bins[i].func = func;
      req->thread_arg_bins[i].args = args;
      if(i == req->num_threads - 1)
      {
        threads_pthreads_func_wrapper((void*)&req->thread_arg_bins[i]);
      }
      else
      {
        int _ret = pthread_create(&req->threads[i], NULL, threads_pthreads_func_wrapper, (void*) &req->thread_arg_bins[i]); 
        if(_ret) { 
          printf("pthread error\n"); exit(-1);
        }
      }
    } 

    for (int i = 0; i < req->num_threads - 1; i++)
    {
      pthread_join(req->threads[i], NULL);
    }
    return ret;
  }

  int threads_pthreads_synch(threads_request* request)
  {
    threads_pthreads_request* req = (threads_pthreads_request*) request;
    return pthread_barrier_wait(&req->barrier);
  }

  int threads_pthreads_begin_critical(threads_request* req)
  {
    threads_pthreads_request* _req = (threads_pthreads_request*) req;
    pthread_mutex_lock(&_req->mutex);
    return 0;
  }

  int threads_pthreads_end_critical(threads_request* req)
  {
    threads_pthreads_request* _req = (threads_pthreads_request*) req;
    pthread_mutex_unlock(&_req->mutex);
    return 0;
  }


  int threads_pthreads_load(threads_funcs* funcs)
  {
    funcs->clean = &threads_pthreads_clean;
    funcs->init  = &threads_pthreads_init;
    funcs->run   = &threads_pthreads_run;
    funcs->synch = &threads_pthreads_synch;
    funcs->begin_critical = &threads_pthreads_begin_critical;
    funcs->end_critical = &threads_pthreads_end_critical;
    return 0;
  }
} 
