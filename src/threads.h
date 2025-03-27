/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_THREADS
#define MINIMOD_THREADS

#include <stddef.h>

typedef void threads_request;

typedef struct
{
  int (*clean)       ();
  int (*init)        (threads_request** request, int num_threads);
  int (*run)         (threads_request*, int (*func)(int thread_num, int num_threads, void* arg_struct), void* args);
  int (*synch)       (threads_request*);
  int (*begin_critical) (threads_request*);
  int (*end_critical)   (threads_request*);
} threads_funcs;

int threads_load       (char* name);
int threads_init       (threads_request** request, int num_threads);
int threads_run        (threads_request*, int (*func)(int thread_num, int num_threads, void* arg_struct), void* args);
int threads_synch      (threads_request*);
int threads_begin_critical (threads_request*);
int threads_end_critical   (threads_request*);

#include"threading/threads_openmp.h" 
#include"threading/threads_pthreads.h"



#endif
