/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "kernel_sample.h"
#include <string.h>
#include <sys/stat.h>
#include <ctime> 


#define DEBUG 0 
#define WRITE 1
#define TIMING 1
#define WINDOW 1
#define STENCIL 0

/*
 * Please note the naming convetion of <module_type>_<module_name>_<function_name>: when building a new module 
 * please follow this naming convetion to avoid conflicts. 
 */

typedef struct
{
} kernel_channel_bench_data;

int kernel_channel_bench_init(void** kernel_data, int id, int job_size)
{
  return 0; 
}

int kernel_channel_bench_run(void* _kernel_data)
{
  int iters  = 10;
  int verify = 1;

  int max_exp = 16; 
  int min_exp = 0; 

  int max_size = 2;
  int min_size = 2;

  for(int i = 0; i < max_exp; i++) max_size *=2;

  if(min_exp == 0) min_size = 1;
  else for(int i = 0; i < min_exp; i++) min_size *=2;

  long* recv_buf = 0; 
  long* send_buf = (long*) malloc(sizeof(long)*max_size);
  
  int me = -1;
  comm_get_id(&me);

  int world = -1;
  comm_get_job_size(&world);
  
  if(me == 0)
  {
    printf("Starting MiniMod Collective Benchmark Tests\n");
    printf("Running for %d iterations ", iters);
    if(verify)
    {
      printf("with correctness verification "); 
    }
    printf("from 2^%d elements to 2^%d elements (max_size=%d) ", min_exp, max_exp, max_size);
    printf("over %d processes\n", world);

  }


  if(me == 0) printf("Starting Fine Grained Test\n");

  for(int size = min_size; size <= max_size; size *= 2)
  {
    std::clock_t t_start;
    if(me == 0) t_start = std::clock();

    comm_request* comm_req = 0;
    comm_channel_init(&comm_req, (void**) &recv_buf, (void*)send_buf, (me + 1) % world, (me + world -1) % world, size, sizeof(long));
    for(int iter = 0; iter < iters; iter++)
    {
      if(verify)
      {
        for(int element = 0; element < size; element++)
        {
          ((long*)send_buf)[element] = iter + element + me;
        }
      }
      comm_channel_start(comm_req);
      for(int element = 0; element < size; element++)
      {
        comm_channel_send(comm_req, element);
      }  
      comm_channel_end(comm_req);
      if(verify)
      {
        for(int element = 0; element < size; element++)
        {
          if(((long*)recv_buf)[element] != iter + element + (me + 1) % world)
          {
            printf("Comm Error, output doesn't match expected\n");
          }
        }
      }
    }
    comm_channel_finalize(comm_req); 

    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("FINE_SEND size: %d time: %.5e Per call: %.5e\n", size, t_total, t_total/iters);
    }
  }

  if(me == 0) printf("Starting Bulk Test\n");

  for(int size = min_size; size <= max_size; size *= 2)
  {
    std::clock_t t_start;
    if(me == 0) t_start = std::clock();

    comm_request* comm_req = 0;
    comm_channel_init(&comm_req, (void**) &recv_buf, (void*)send_buf, (me + 1) % world, (me + world -1) % world, 1, size * sizeof(long));
    for(int iter = 0; iter < iters; iter++)
    {
      if(verify)
      {
        for(int element = 0; element < size; element++)
        {
          ((long*)send_buf)[element] = iter + element + me;
        }
      }
      comm_channel_start(comm_req);
      
      
        comm_channel_send(comm_req, 0);;
      
      comm_channel_end(comm_req);
      if(verify)
      {
        for(int element = 0; element < size; element++)
        {
          if(((long*)recv_buf)[element] != iter + element + (me + 1) % world)
          {
            printf("Comm Error, output doesn't match expected\n");
          }
        }
      }
    }
    comm_channel_finalize(comm_req); 

    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("BULK_SEND - size: %d time: %.5e Per call: %.5e\n", size, t_total, t_total/iters);
    }
  }



  return 0;
}

int kernel_channel_bench_clean(void* _kernel_data)
{
  return 0; 
}

int kernel_channel_bench_load(kernelfuncs* funcs)
{
  funcs->init  = &kernel_channel_bench_init;
  funcs->run   = &kernel_channel_bench_run;
  funcs->clean = &kernel_channel_bench_clean;
  return 0;
}

