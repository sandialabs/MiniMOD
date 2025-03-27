/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_THREADS_OPENMP
#define MINIMOD_THREADS_OPENMP

#include <stdio.h>
#include <omp.h>
#include "threads.h"

#ifdef HAVE_OMP

namespace Threads_OpenMP
{
  typedef struct
  {
    
    int num_threads;
    omp_lock_t lck;
  } threads_openmp_request;

  int threads_openmp_load(threads_funcs*);
} 
#endif

#endif
