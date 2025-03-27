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
} kernel_collective_bench_data;

int kernel_collective_bench_init(void** kernel_data, int id, int job_size)
{
  return 0; 
}

int kernel_collective_bench_run(void* _kernel_data)
{
  int iters  = 10000;
  int verify = 0;

  int max_exp = 20; 
  int min_exp = 0; 

  int max_size = 2;
  int min_size = 2;

  for(int i = 0; i < max_exp; i++) max_size *=2;

  if(min_exp == 0) min_size = 1;
  else for(int i = 0; i < min_exp; i++) min_size *=2;

  long* recv_buf = (long*) malloc(sizeof(long)*max_size);
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

    printf("\n\nBegining barrier test\n");
  }


  std::clock_t t_start;
  if(me == 0) t_start = std::clock();
  for(int i = 0; i < iters; i++)
  {
    comm_barrier();
  }
  double t_total;
  if(me == 0)
  {
    t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
    printf("Barrier - Total time: %.5e Per call: %.5e\n", t_total, t_total/iters);
  }



  
  for(int i = 0; i < max_size; i++)
  {
    ((int*)send_buf)[i] = me;
    ((int*)recv_buf)[i] = 0;
  }
  
  int size = min_size; 
  for(int j = 0; j < max_exp; j++)
  {
    if(me == 0) t_start = std::clock();
    for(int i = 0; i < iters; i++)
    {
      int turn=i%world;
      if(turn == me)
      {
        comm_broadcast(send_buf, size, COMM_DATATYPE_INT, turn);
      }
      else
      {
        comm_broadcast(recv_buf, size, COMM_DATATYPE_INT, turn);
      }
      if(verify && turn != me)
      {
        for(int k = 0; k < max_size; k++)
        {
          if(k < size)
          {
            if(((int*)recv_buf)[k]!=turn)
            {
              printf("Error! Expected value %d recived value %d. Element %d collective size %d\n", turn, ((int*)recv_buf)[k], k, size);
              exit(-1);
            } 
          }
          else
          {
            if(((int*)recv_buf)[k]!=0)       
            {
              printf("Error! Expected value %d (overwrite check sent values are %d) recived value %d. Element %d collective size %d\n", 0, turn, ((int*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
        }
      }
    }
    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("Broadcast - Type: int Size: %d bytes Total time: %.5e Per call: %.5e\n", size*4, t_total, t_total/iters);
    }
    size *= 2;
  }



  
  for(int i = 0; i < max_size; i++)
  {
    ((float*)send_buf)[i] = me;
    ((float*)recv_buf)[i] = 0;
  }
  
  size = min_size; 
  for(int j = 0; j < max_exp; j++)
  {
    if(me == 0) t_start = std::clock();
    for(int i = 0; i < iters; i++)
    {
      int turn=i%world;
      if(turn == me)
      {
        comm_broadcast(send_buf, size, COMM_DATATYPE_FLOAT, turn);
      }
      else
      {
        comm_broadcast(recv_buf, size, COMM_DATATYPE_FLOAT, turn);
      }
      if(verify && turn != me)
      {
        for(int k = 0; k < max_size; k++)
        {
          if(k < size)
          {
            if(((float*)recv_buf)[k]!=turn)
            {
              printf("Error! Expected value %d recived value %f. Element %d collective size %d\n", turn, ((float*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
          else
          {
            if(((float*)recv_buf)[k]!=0)
            {
              printf("Error! Expected value %d (overwrite check sent values are %d) recived value %f. Element %d collective size %d\n", 0, turn, ((float*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
        }
      }
    }
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("Broadcast - Type: float Size: %d bytes Total time: %.5e Per call: %.5e\n", size*4, t_total, t_total/iters);
    }
    size *= 2;
  }



  
  for(int i = 0; i < max_size; i++)
  {
    ((double*)send_buf)[i] = me;
    ((double*)recv_buf)[i] = 0;
  }
  
  size = min_size; 
  for(int j = 0; j < max_exp; j++)
  {
    if(me == 0) t_start = std::clock();
    for(int i = 0; i < iters; i++)
    {
      int turn=i%world;
      if(turn == me)
      {
        comm_broadcast(send_buf, size, COMM_DATATYPE_DOUBLE, turn);
      }
      else
      {
        comm_broadcast(recv_buf, size, COMM_DATATYPE_DOUBLE, turn);
      }
      if(verify && turn != me)
      {
        for(int k = 0; k < max_size; k++)
        {
          if(k < size)
          {
            if(((double*)recv_buf)[k]!=turn)
            {
              printf("Error! Expected value %d recived value %f. Element %d collective size %d\n", turn, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
          else
          {
            if(((double*)recv_buf)[k]!=0)
            {
              printf("Error! Expected value %d (overwrite check sent values are %d) recived value %f. Element %d collective size %d\n", 0, turn, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
        }
      }
    }
    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("Broadcast - Type: double Size: %d bytes Total time: %.5e Per call: %.5e\n", size*8, t_total, t_total/iters);
    }
    size *= 2;
  }



  
  for(int i = 0; i < max_size; i++)
  {
    ((double*)send_buf)[i] = ((me+i)%world);
    ((double*)recv_buf)[i] = -1;
  }
  
  size = min_size; 
  for(int j = 0; j < max_exp; j++)
  {
    if(me == 0) t_start = std::clock();
    for(int i = 0; i < iters; i++)
    {
      comm_allreduce(send_buf, recv_buf, size, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_MAX);
      if(verify)
      {
        for(int k = 0; k < max_size; k++)
        {
          if(k < size)
          {
            if(((double*)recv_buf)[k]!=world-1)
            {
              printf("Error! Expected value %d recived value %f. Element %d collective size %d\n", world-1, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
          else
          {
            if(((double*)recv_buf)[k]!=-1)
            {
              printf("Error! Expected value %d (overwrite check sent values are %d) recived value %f. Element %d collective size %d\n", 0, world-1, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
        }
      }
    }
    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("AllReduce - Op: MAX Type: double Size: %d bytes Total time: %.5e Per call: %.5e\n", size*8, t_total, t_total/iters);
    }
    size *= 2;
  }

  
  for(int i = 0; i < max_size; i++)
  {
    ((double*)send_buf)[i] = -1 * ((me+i)%world);
    ((double*)recv_buf)[i] = 1;
  }
  
  size = min_size; 
  for(int j = 0; j < max_exp; j++)
  {
    if(me == 0) t_start = std::clock();
    for(int i = 0; i < iters; i++)
    {
      comm_allreduce(send_buf, recv_buf, size, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_MIN);
      if(verify)
      {
        for(int k = 0; k < max_size; k++)
        {
          if(k < size)
          {
            if(((double*)recv_buf)[k]!=0-(world-1))
            {
              printf("Error! Expected value %d recived value %f. Element %d collective size %d\n", world-1, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
          else
          {
            if(((double*)recv_buf)[k]!=1)
            {
              printf("Error! Expected value %d (overwrite check sent values are %d) recived value %f. Element %d collective size %d\n", 0, world-1, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
        }
      }
    }
    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("AllReduce - Op: MIN Type: double Size: %d bytes Total time: %.5e Per call: %.5e\n", size*8, t_total, t_total/iters);
    }
    size *= 2;
  }

  
  for(int i = 0; i < max_size; i++)
  {
    ((double*)send_buf)[i] = i;
    ((double*)recv_buf)[i] = 0;
  }
  
  size = min_size; 
  for(int j = 0; j < max_exp; j++)
  {
    if(me == 0) t_start = std::clock();
    for(int i = 0; i < iters; i++)
    {
      comm_allreduce(send_buf, recv_buf, size, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_SUM);
      if(verify)
      {
        for(int k = 0; k < max_size; k++)
        {
          if(k < size)
          {
            if(((double*)recv_buf)[k]!=world*k)
            {
              printf("Error! Expected value %d recived value %f local value %f. Element %d collective size %d\n", world*k, ((double*)recv_buf)[k], ((double*)send_buf)[k], k, size);
              exit(-1);
            }
          }
          else
          {
            if(((double*)recv_buf)[k]!=0)
            {
              printf("Error! Expected value %d (overwrite check sent values are %d) recived value %f. Element %d collective size %d\n", 0, world-1, ((double*)recv_buf)[k], k, size);
              exit(-1);
            }
          }
        }
      }
    }
    double t_total;
    if(me == 0)
    {
      t_total = (std::clock() - t_start)/(double) CLOCKS_PER_SEC;
      printf("AllReduce - Op: SUM Type: double Size: %d bytes Total time: %.5e Per call: %.5e\n", size*8, t_total, t_total/iters);
    }
    size *= 2;
  }

  


  return 0;
}

int kernel_collective_bench_clean(void* _kernel_data)
{
  return 0; 
}

int kernel_collective_bench_load(kernelfuncs* funcs)
{
  funcs->init  = &kernel_collective_bench_init;
  funcs->run   = &kernel_collective_bench_run;
  funcs->clean = &kernel_collective_bench_clean;
  return 0;
}

