/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_KERNEL
#define MINIMOD_KERNEL

#include "comm.h"
#include "gran.h"
#include "threads.h"

typedef struct
{
  int (*init) (void**, int, int);
  int (*parse_args) (void*, int, char**);
  int (*run)  (void*);
  int (*clean)(void*);
  int id;
  int job_size;
  void* data;
} kernelfuncs;


int kernel_load(char* name, int id, int job_size);

int kernel_parse_args(int, char**);

int kernel_run();

int kernel_clean();


#include"kernel/kernel_sample.h" 
#include"kernel/kernel_sample_thread.h"
#include"kernel/kernel_sample_threads.h"
#include"kernel/kernel_cmb.h"
#include"kernel/kernel_collective_bench.h"
#include"kernel/kernel_channel_bench.h"
#include"kernel/kernel_gran_bench.h"
#include"kernel/kernel_partsend.h"

#endif
