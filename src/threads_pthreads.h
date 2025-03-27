/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_THREADS_PTHREADS
#define MINIMOD_THREADS_PTHREADS

#include <stdio.h>
#include <pthread.h>

#include "threads.h"


namespace Threads_Pthreads
{

  typedef struct {
    int thread_num;
    int num_threads;
    void* args;
    int (*func)(int thread_num, int num_threads, void* arg_struct);
  } threads_pthreads_arg_wrapper;

  typedef struct
  {
    int num_threads;
    pthread_t* threads;
    pthread_barrier_t   barrier; 
    pthread_mutex_t mutex;
    threads_pthreads_arg_wrapper* thread_arg_bins;
    
  } threads_pthreads_request;

  int threads_pthreads_load(threads_funcs*);
} 
#endif
