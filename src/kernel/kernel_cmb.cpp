/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "kernel_cmb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <stdint.h> 
#include <time.h>  
#include <inttypes.h> 

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#define BILLION 1000000000L



typedef struct
{
  
  int id;
  int job_size;

  
  int npeers;
  int niters;
  int nmsgs;
  int nbytes;
  int msg_size;
  int cache_size;
  int ppn;
  int machine_output;
  int threads;
  int validate;
  int laggard;
  char filename[64];
  uint64_t comp_time;
  uint64_t* comp_thread;
  double noise;

  
  char* cache_buf;
  void* send_buf;
  void** recv_bufs;

  
  int* peer_send_id;
  int* peer_recv_id;

  
  gran_request* gran_req;
} kernel_cmb_data;

typedef struct
{
  
  int npeers;
  int niters;
  int nmsgs;
  int nbytes;
  int msg_size;
  int cache_size;
  int ppn;
  int machine_output;
  int threads;
  int validate;
  int laggard;
  char filename[64];
  uint64_t comp_time;
  double noise;
} kernel_cmb_input_arg_data;


kernel_cmb_input_arg_data* TEMP_INPUT_ARG_DATA;

void usage(void)
{
    fprintf(stderr, "Usage: msgrate [OPTION]...\n\n");
    fprintf(stderr, "  -p  <num>     Number of peers used in communication\n");
    fprintf(stderr, "  -i  <num>     Number of iterations per test\n");
    fprintf(stderr, "  -m  <num>     Number of messages per peer per iteration\n");
    fprintf(stderr, "  -s  <size>    Number of bytes per message\n");
    fprintf(stderr, "  -b  <size>    Size of buffer sent to each peer in number of bytes\n");
    fprintf(stderr, "  -c  <size>    Cache size in bytes\n");
    fprintf(stderr, "  -n  <num>     Number of procs per node\n");
    fprintf(stderr, "  -t  <num>     Number of threads per process\n");
    fprintf(stderr, "  -ct <ns>      Nanoseconds of simulated comp time per thread\n");
    fprintf(stderr, "  -cn <num>     Fractional amount slower the laggered thread is from the rest\n");
    fprintf(stderr, "  -f  <path>    Read comptime from an inputfile\n");
    fprintf(stderr, "  -l            Use a single laggard thread\n");
    fprintf(stderr, "  -o            Format output to be machine readable\n");
    fprintf(stderr, "  -v            Perform communication correctness validation\n");
    fprintf(stderr, "  -h            Display this help message and exit\n");
}

void cache_invalidate(kernel_cmb_data* kernel_data)
{
  kernel_data->cache_buf[0] = 1;
  for (int i = 1 ; i < kernel_data->cache_size ; ++i) {
    kernel_data->cache_buf[i] = kernel_data->cache_buf[i - 1];
  }
}

void display_result(const char *test, double result, kernel_cmb_data* kernel_data)
{
  if (0 == kernel_data->id) {
    if (kernel_data->machine_output) {
      printf("%.16f ", result);
    } else {
      printf("%10s: %.16f\n", test, result);
    }
  }
}

int kernel_cmb_parse_args(void* _data, int argc, char** argv)
{
  kernel_cmb_input_arg_data* arg = (kernel_cmb_input_arg_data*) malloc(sizeof(kernel_cmb_input_arg_data));

  
  arg->npeers = 6;
  arg->niters = 512;
  arg->nmsgs = 2;
  arg->nbytes = -1;
  arg->msg_size = -1;
  arg->cache_size = (8 * 1024 * 1024 / sizeof(int));
  arg->ppn = 1;
  arg->machine_output = 0;
  arg->laggard = 0;
  arg->threads = 2;
  arg->comp_time = 1000;
  arg->noise = 0.04;
  arg->validate = 0;
  for (int i = 0; i < 64; ++i)
  {
    arg->filename[i] = 0;
  }

  
  
  char **foobar = (char **)malloc(sizeof(char*) * (argc+1));
  foobar[0] = (char*)"foo"; 
  for (int i = 1;i <= argc; ++i) {
    foobar[i] = argv[i-1];
  }
  
  #define OPTIONAL 2
  static struct option long_options[] = {
          {"p", OPTIONAL, NULL, 0 },
          {"i", OPTIONAL, NULL, 1 },
          {"m", OPTIONAL, NULL, 2 },
          {"s", OPTIONAL, NULL, 3 },
          {"c", OPTIONAL, NULL, 4 },
          {"n", OPTIONAL, NULL, 5 },
          {"o", OPTIONAL, NULL, 6 },
          {"t", OPTIONAL, NULL, 7 },
          {"ct", OPTIONAL, NULL, 8 },
          {"cn", OPTIONAL, NULL, 9 },
          {"v", OPTIONAL, NULL, 10 },
          {"b", OPTIONAL, NULL, 11 },
          {"f", OPTIONAL, NULL, 12 },
          {"l", OPTIONAL, NULL, 13 },
          {"h", OPTIONAL, NULL, 14 },
          {0,0,0,0}
  };

  while (1) {
    int option_index = 0;
    int c = getopt_long_only(argc+1, foobar, "", long_options, &option_index );
    if ( c == -1) break;
    switch (c) {
      case 0:
        arg->npeers = atoi(optarg);
        break;
      case 1:
        arg->niters = atoi(optarg);
        break;
      case 2:
        arg->nmsgs = atoi(optarg);
        break;
      case 3:
        arg->msg_size = atoi(optarg);
        break;
      case 4:
        arg->cache_size = atoi(optarg) / sizeof(int);
        break;
      case 5:
        arg->ppn = atoi(optarg);
        break;
      case 6:
        arg->machine_output = 1;
        break;
      case 7:
        arg->threads = atoi(optarg);
        break;
      case 8:
        arg->comp_time = atoll(optarg);
        break;
      case 9:
        arg->noise = atof(optarg);
        break;
      case 10:
        arg->validate = 1;
        break;
      case 11:
        arg->nbytes = atoi(optarg);
        break;
      case 12:
        strcpy(arg->filename, optarg); 
        break;
      case 13:
        arg->laggard = 1;
        break;
      case 14:
      default:
        usage();
        return(-1);
    }  
  }

  TEMP_INPUT_ARG_DATA = arg;
  free(foobar);

  return(0);
}

int kernel_cmb_init(void** kernel_data, int id, int job_size)
{
  kernel_cmb_data* new_data = (kernel_cmb_data*) malloc(sizeof(kernel_cmb_data));

  new_data->id = id;
  new_data->job_size = job_size;

  
  new_data->npeers = TEMP_INPUT_ARG_DATA->npeers;
  new_data->niters = TEMP_INPUT_ARG_DATA->niters;
  new_data->nbytes = TEMP_INPUT_ARG_DATA->nbytes;
  new_data->msg_size = TEMP_INPUT_ARG_DATA->msg_size;
  new_data->nmsgs = TEMP_INPUT_ARG_DATA->nmsgs;
  new_data->cache_size = TEMP_INPUT_ARG_DATA->cache_size;
  new_data->ppn = TEMP_INPUT_ARG_DATA->ppn;
  new_data->machine_output = TEMP_INPUT_ARG_DATA->machine_output;
  new_data->threads = TEMP_INPUT_ARG_DATA->threads;
  new_data->validate = TEMP_INPUT_ARG_DATA->validate;
  new_data->laggard = TEMP_INPUT_ARG_DATA->laggard;
  new_data->comp_time = TEMP_INPUT_ARG_DATA->comp_time;
  new_data->noise = TEMP_INPUT_ARG_DATA->noise;
  strcpy(new_data->filename, TEMP_INPUT_ARG_DATA->filename);
  free(TEMP_INPUT_ARG_DATA);

  
  if (-1 != new_data->nbytes && -1 != new_data->msg_size) {
    fprintf(stderr, "Error: Must set either buffer size or message size, not both.\n");
    return(-1);
  } else if (-1 == new_data->nbytes && -1 == new_data->msg_size) {
    new_data->msg_size = 1024;
  }
  if (-1 != new_data->nbytes) {
    new_data->msg_size = new_data->nbytes / new_data->nmsgs;
  } else {
    new_data->nbytes = new_data->msg_size * new_data->nmsgs;
  }

  /* sanity check */
  if (job_size <= new_data->npeers) {
    fprintf(stderr, "Error: job size (%d) <= number of peers (%d)\n", new_data->job_size, new_data->npeers);
    return(-1);
  } else if (new_data->ppn < 1) {
    fprintf(stderr, "Error: process per node (-n #) must be positive\n");
    return(-1);
  } else if ((double)job_size / (double)new_data->ppn <= (double)new_data->npeers) {
    fprintf(stderr, "Error: node count <= number of peers %i / %i <= %i\n",job_size,new_data->ppn,new_data->npeers);
    return(-1);
  } else if (job_size % 2 == 1) {
    fprintf(stderr, "Error: node count of %d isn't even.\n", job_size);
    return(-1);
  } else if (new_data->nbytes % new_data->nmsgs != 0) {
    fprintf(stderr, "Error: number of messages must be an even divisor of the buffer size.\n");
    return(-1);
  } else if (new_data->msg_size % 8 != 0) {
    fprintf(stderr, "Error: message size must be an evenly divisible into doubles.\n");
    return(-1);
  } else if (new_data->nmsgs != new_data->threads) {
    fprintf(stderr, "Error: number of messages must be equal to number of threads.\n");
    return(-1);
  }

  /* allocate buffers */
  new_data->peer_send_id = (int*)malloc(sizeof(int) * new_data->npeers);
  if (NULL == new_data->peer_send_id) exit(-1);
  new_data->peer_recv_id = (int*)malloc(sizeof(int) * new_data->npeers);
  if (NULL == new_data->peer_recv_id) exit(-1);
  new_data->cache_buf = (char*)malloc(sizeof(char) * new_data->cache_size);
  if (NULL == new_data->cache_buf) exit(-1);
  new_data->send_buf = (void*)malloc(new_data->nbytes);
  if (NULL == new_data->send_buf) exit(-1);
  new_data->recv_bufs = (void**)malloc(sizeof(void*) * new_data->npeers);
  if (NULL == new_data->recv_bufs) exit(-1);

  new_data->comp_thread = (uint64_t*)calloc(sizeof(uint64_t) , new_data->threads);
  if (NULL == new_data->comp_thread) exit(-1);

  if (new_data->validate)
  {
    if (!id) printf("- - -  DO NOT USE THIS DATA  - - - \n- - -   validation enabled   - - -\n");
    int* send_buf_int = (int*)new_data->send_buf;
    int num_ints = new_data->nbytes / sizeof(int);
    for (int i = 0; i < num_ints; ++i)
    {
      send_buf_int[i] = i;
    }
  }

  /* calculate peers */
  for (int i = 0 ; i < new_data->npeers ; ++i) 
  {
    new_data->peer_send_id[i] = (id + i + 1) % job_size;
    new_data->peer_recv_id[i] = (id - i - 1 + job_size) % job_size;
  }

  
  void** peer_send_buf = (void**)malloc(sizeof(void*) * new_data->npeers);
  void*** peer_recv_buf = (void***)malloc(sizeof(void**) * new_data->npeers);
  int* entry_num = (int*)malloc(sizeof(int) * new_data->npeers);
  int* entry_size = (int*)malloc(sizeof(int) * new_data->npeers);
  for(int i = 0; i < new_data->npeers; i++)
  {
    peer_send_buf[i] = new_data->send_buf;
    peer_recv_buf[i] = (void**)&new_data->recv_bufs[i];
    entry_num[i] = new_data->nmsgs;
    entry_size[i] = new_data->msg_size;
  }

  int err = gran_stencil_init(&new_data->gran_req, peer_recv_buf, peer_send_buf, new_data->peer_recv_id, new_data->peer_send_id, entry_num, entry_size, new_data->npeers);
  if (err) 
  {
    printf("error is gran_stencil_init!\n");
    return err;
  }

  free(peer_send_buf);
  free(peer_recv_buf);
  free(entry_num);
  free(entry_size);

  *kernel_data = new_data;
  return 0;
}

uint64_t timespec_to_nano(struct timespec* ts)
{
  uint64_t tmp = 0;
  tmp += ts->tv_nsec;
  tmp += ts->tv_sec * BILLION;
  
  return tmp;
}

typedef struct {
  kernel_cmb_data* kernel_data;
} kernel_cmb_args;

int kernel_cmb_func(int thread_num, int num_threads, void* args)
{
  kernel_cmb_args* kernel_args = (kernel_cmb_args*) args;
  kernel_cmb_data* kernel_data = kernel_args->kernel_data;

  struct timespec req,rem;










  
  req.tv_nsec = kernel_data->comp_thread[thread_num] % BILLION;
  req.tv_sec = kernel_data->comp_thread[thread_num] / BILLION;
  
  clock_nanosleep(CLOCK_REALTIME,0,&req, &rem);

  
  






















 
  













  for (int i = 0; i < kernel_data->npeers; i++)
  {
    int err = gran_stencil_ready(kernel_data->gran_req, i, thread_num);
    if (err) 
    {
      printf("(%d)error in gran_stencil_ready(gran_req<%p>, dir<%d>, index<%d>\n",kernel_data->id,kernel_data->gran_req,i,thread_num);
      fflush(stdout);
      return err;
    }
  }
  return 0;
}

int kernel_cmb_run(void* _kernel_data)
{
  kernel_cmb_data* kernel_data = (kernel_cmb_data*)_kernel_data;
  kernel_cmb_args* args = (kernel_cmb_args*)malloc(sizeof(kernel_cmb_args));
  args->kernel_data = kernel_data;

  double tmp, iter, total_largest = 0;
  uint64_t iter_t, total_t = 0;
  struct timespec start_time;
  struct timespec end_time;

  threads_request* t_req;
  threads_init(&t_req, kernel_data->threads);

  int err = 0;
  
  int file_input = kernel_data->filename[0] ? 1 : 0;
  FILE* fp;
  if (file_input)
  {
    fp = fopen(kernel_data->filename , "r");
    if(fp == NULL) {
      printf("Error opening provided comptime file.");
      return(-1);
    }
  }

  boost::mt19937 rng; 
  boost::normal_distribution<> nd(kernel_data->comp_time, kernel_data->noise);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> vg(rng, nd);

  
  for (int i = 0 ; i <kernel_data->niters ; ++i) {
    
    if (!(i % 20) && !kernel_data->id) printf("(%d)iter:%d\n",kernel_data->id,i);

    
    cache_invalidate(kernel_data);

    
    if (kernel_data->validate)
    {
      int** recv_bufs_int = (int**)kernel_data->recv_bufs;
      int num_ints = kernel_data->nbytes / sizeof(int);
      for (int p = 0; p < kernel_data->npeers; ++p)
      {
        for (int j = 0; j < num_ints; ++j)
        {
          recv_bufs_int[p][j] = 0;
        }
      }
    }

    
    if(file_input)
    {
      int t = 0;
      char tmp_buffer[256];
      
      fgets(tmp_buffer, INT_MAX, fp);
      
      char *token = strtok(tmp_buffer, ",");
      while (token) { 
        uint64_t n = (uint64_t)atoll(token);
        kernel_data->comp_thread[t++] = n;

        token = strtok(NULL, ",");
      }
    }
    else if (kernel_data->laggard)
    {
      for (int t = 0; t < kernel_data->threads; ++t)
      {
        
        kernel_data->comp_thread[t] = (t == 0) ? (uint64_t)(kernel_data->comp_time * (kernel_data->noise + 1)) : kernel_data->comp_time;

      }
    }
    else
    {  
      for (int t = 0; t < kernel_data->threads; ++t)
      {
        
        double sample = vg();
        kernel_data->comp_thread[t] = (sample > 0) ? (uint64_t)sample : 0;

      }
    }

    
    comm_barrier();
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    gran_stencil_start(kernel_data->gran_req);
    err = threads_run(t_req, &kernel_cmb_func, args);
    if (err) printf("(%d)iter:%d (threads_run error!)\n",kernel_data->id,i);
    gran_stencil_end(kernel_data->gran_req);
    clock_gettime(CLOCK_MONOTONIC, &end_time);

    
    iter_t = timespec_to_nano(&end_time) - timespec_to_nano(&start_time);
    
    total_t += iter_t;
    
    iter = ((double)iter_t / BILLION); 
    comm_allreduce(&iter, &tmp, 1, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_MAX);
    
    total_largest += tmp;

    
    if (kernel_data->validate)
    {
      int** recv_bufs_int = (int**)kernel_data->recv_bufs;
      int num_ints = kernel_data->nbytes / sizeof(int);
      for(int p = 0; p < kernel_data->npeers; ++p)
      {
        for(int j = 0; j < num_ints; ++j)
        {
          if (recv_bufs_int[p][j] != j)
          {
            printf("(%d)VALIDATION FAILED: recv_buf from peer:%d incorrect durring iter:%d! (%d != %d)\n",kernel_data->id,p,i,recv_bufs_int[p][j],j);
            exit(1);
          }
        }
      }
    }
  }

  
  double total = ((double)total_t / BILLION); 
  comm_allreduce(&total, &tmp, 1, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_MAX);
  double max_local_time = tmp; 
  double max_local_time_iter = max_local_time / kernel_data->niters;
  double max_global_time = total_largest; 
  double max_global_time_iter = max_global_time / kernel_data->niters;
  comm_allreduce(&total, &tmp, 1, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_SUM);
  double avg_time = tmp / kernel_data->job_size; 
  double avg_time_iter = avg_time / kernel_data->niters; 
  double min_local_msgrate = (kernel_data->niters * kernel_data->npeers * kernel_data->nmsgs * 2) / max_local_time; 
  double min_global_msgrate = (kernel_data->niters * kernel_data->npeers * kernel_data->nmsgs * 2) / max_global_time; 
  double avg_msgrate = (kernel_data->niters * kernel_data->npeers * kernel_data->nmsgs * 2) / avg_time; 
  double min_local_bandwidth = min_local_msgrate * kernel_data->nbytes; 
  double min_global_bandwidth = min_global_msgrate * kernel_data->nbytes; 
  double avg_bandwidth = avg_msgrate * kernel_data->nbytes; 
 
  if (!kernel_data->id)printf("#STARTRESULTS\n"); 
  display_result("#max_local_time        ", max_local_time, kernel_data);
  display_result("#max_global_time       ", max_global_time, kernel_data);
  display_result("#avg_time              ", avg_time, kernel_data);
  display_result("#max_local_time_iter   ", max_local_time_iter, kernel_data);
  display_result("#max_global_time_iter  ", max_global_time_iter, kernel_data);
  display_result("#avg_time_iter         ", avg_time_iter, kernel_data);
  display_result("#min_local_msgrate     ", min_local_msgrate, kernel_data);
  display_result("#min_global_msgrate    ", min_global_msgrate, kernel_data);
  display_result("#avg_msgrate           ", avg_msgrate, kernel_data);
  display_result("#min_global_bandwidth  ", min_local_bandwidth, kernel_data);
  display_result("#min_local_bandwidth   ", min_global_bandwidth, kernel_data);
  display_result("#avg_bandwidth         ", avg_bandwidth, kernel_data);
  if (!kernel_data->id)printf("#ENDRESULTS\n"); 

  fflush(stdout);
  comm_barrier();
 
  free(args);

  return 0;
}

int kernel_cmb_clean(void* _kernel_data)
{
  kernel_cmb_data* kernel_data = (kernel_cmb_data*)_kernel_data;

  
  gran_stencil_finish(kernel_data->gran_req);

  
  free(kernel_data->peer_send_id);
  free(kernel_data->peer_recv_id);
  free(kernel_data->cache_buf);
  free(kernel_data->send_buf);
  free(kernel_data->recv_bufs);
  free(kernel_data->comp_thread);
  free(kernel_data);
  return 0;
}

int kernel_cmb_load(kernelfuncs* funcs)
{
  funcs->init  = &kernel_cmb_init;
  funcs->parse_args = &kernel_cmb_parse_args;
  funcs->run   = &kernel_cmb_run;
  funcs->clean = &kernel_cmb_clean;
  return 0;
}
