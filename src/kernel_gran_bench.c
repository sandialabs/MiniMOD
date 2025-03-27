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

#ifdef TIME_CALLS
#include <time.h>

#define NSEC_PER_SEC (1000000000)

/* calculate difference between two timestamps (stop - start) in nsecs */
static inline uint64_t timespec_diff(struct timespec *start, struct timespec *stop) {
    struct timespec result;
    if ((stop->tv_nsec - start->tv_nsec) < 0) {
        result.tv_sec = stop->tv_sec - start->tv_sec - 1;
        result.tv_nsec = (stop->tv_nsec - start->tv_nsec) + NSEC_PER_SEC; 
    } else {
        result.tv_sec = stop->tv_sec - start->tv_sec;
        result.tv_nsec = stop->tv_nsec - start->tv_nsec;
    }
    return ((uint64_t)(result.tv_sec) * NSEC_PER_SEC) + (uint64_t)result.tv_nsec;
}
#endif

typedef struct
{
} kernel_gran_bench_data;

int kernel_gran_bench_init(void** kernel_data, int id, int job_size)
{
  return 0; 
}

int kernel_gran_bench_run(void* _kernel_data)
{
  printf("Executing: kernel_gran_bench\n");
  fflush(stdout);
  
  
  int iters  = 5000;
  int verify = 0;
  int err;

  
  
  
  
  
  
  

  
  int min_exp = 11;
  int max_exp = 21; 
  int long_size = sizeof(long);
  int min_buf_size = long_size*pow(2,min_exp); 
  int max_buf_size = long_size*pow(2,max_exp);
  int cur_buf_size;
   
  
  
  int min_chunks = 1;
  int max_chunks = 512;
  int num_chunks;
  int chunk_size;
  #ifdef TIME_CALLS
    iters++; 
  #endif

  long* recv_buf = 0; 
  long* send_buf = (long*) malloc(max_buf_size);
  long* dummy_buf = 0;
  
  long** dummy_buf_ptr = &dummy_buf;
  long** recv_buf_ptr = &recv_buf;
  long** send_buf_ptr = &send_buf;

  int me = -1;
  comm_get_id(&me);

  int world = -1;
  comm_get_job_size(&world);

  int send_id = (me + world -1) % world;
  int recv_id = (me + 1) % world;

  int dead_id = -1;
  
  #ifdef TIME_CALLS
    timespec t_starts_start, t_starts_end;
    timespec t_startr_start, t_startr_end;
    timespec t_ends_start, t_ends_end;
    timespec t_endr_start, t_endr_end;
    timespec t_readys_start, t_readys_end; 

    uint64_t *call_starts_times = (uint64_t *)malloc(iters * sizeof(uint64_t));
    if (call_starts_times == NULL) {
      printf("%s:%d: Unable to allocate call_starts_times array\n",__FILE__,__LINE__);
      exit(-1);
    } 
    uint64_t *call_startr_times = (uint64_t *)malloc(iters * sizeof(uint64_t));
    if (call_startr_times == NULL) {
      printf("%s:%d: Unable to allocate call_startr_times array\n",__FILE__,__LINE__);
      exit(-1);
    } 
    uint64_t *call_ends_times = (uint64_t *)malloc(iters * sizeof(uint64_t));
    if (call_ends_times == NULL) {
      printf("%s:%d: Unable to allocate call_ends_times array\n",__FILE__,__LINE__);
      exit(-1);
    } 
    uint64_t *call_endr_times = (uint64_t *)malloc(iters * sizeof(uint64_t));
    if (call_endr_times == NULL) {
      printf("%s:%d: Unable to allocate call_endr_times array\n",__FILE__,__LINE__);
      exit(-1);
    } 
    uint64_t **call_readys_times = (uint64_t **)malloc(iters * sizeof(uint64_t *));
    if (call_readys_times == NULL) {
      printf("%s:%d: Unable to allocate call_readys_times array\n",__FILE__,__LINE__);
      exit(-1);
    } 
  #endif

  #ifdef DO_WORK
    #define I (512)
    #define J (512)
    #define K (512)
    #define NUM_MM (100) 

    double A[I][J];
    double B[J][K];
    double C[I][K];

    printf("Initializing work arrays ...\n");
    fflush(stdout);
  
    for (int i = 0; i < I; ++i) {
      for (int j = 0; j < J; ++j) {
        A[i][j] = (float)rand()/(float)(RAND_MAX/10);
      }
    }

    for (int j = 0; j < J; ++j) {
      for (int k = 0; k < K; ++k) {
        B[j][k] = (float)rand()/(float)(RAND_MAX/10);
      }
    }

    for (int i = 0; i < I; ++i) {
      for (int k = 0; k < K; ++k) {
        C[i][k] = 0.0;
      }
    }

    printf(".... Done.\n");
    fflush(stdout);
  #endif

  if(me == 0)
  {
    printf("MiniMod Granularity Benchmark\n");
    printf("Iterations        : %d\n", iters);
    printf("Verification      : ");
    if(verify)
    { 
      printf("True\n");
    } else { printf("False\n");}
    printf("Min buf size      : %d bytes\n", min_buf_size);
    printf("Max buf size      : %d bytes\n", max_buf_size);
    printf("Min chunks        : %d\n", min_chunks);
    printf("Max chunks        : %d\n", max_chunks);
    printf("Processes         : %d\n", world);

    #ifdef TIME_CALLS
    printf("Timing communication calls\n");
    #endif

    #ifdef DO_WORK
    printf("Doing work between each call to ready.\n");
    printf("Matrices: i = %d, j = %d, k = %d\n", I, J, K);
    printf("MMs per ready call: %d\n", NUM_MM);
    #endif
  }

  std::clock_t t_start, t_end;
  double t_total;

  printf("\n");
  fflush(stdout); 


  for(cur_buf_size = min_buf_size; cur_buf_size <= max_buf_size; cur_buf_size*=2)
  {
    #ifdef TIME_CALLS
      
      int chunk_multiplier = 1;
      for (num_chunks = min_chunks; num_chunks <= max_chunks && num_chunks <= (cur_buf_size/long_size); num_chunks = chunk_multiplier * chunk_multiplier * 2, chunk_multiplier *=2)
    #else
      for (num_chunks = min_chunks; num_chunks <= max_chunks && num_chunks <= (cur_buf_size/long_size); num_chunks *= 2)
    #endif
    {
      gran_request* gran_req = 0;
      gran_request* gran_req2 = 0;

      chunk_size = cur_buf_size/num_chunks;

      if(me == 0)
      {
        
        
        err = gran_stencil_init(&gran_req, (void***) &dummy_buf_ptr, (void**)send_buf_ptr, &dead_id, &send_id, &num_chunks, &chunk_size, 1);
        if (err != 0) {
          printf("%s:%d: ERROR: gran_stencil_init returned error %d\n", __FILE__,__LINE__,err);
          exit(-1);
        }
        err = gran_stencil_init(&gran_req2, (void***) &recv_buf_ptr, (void**)send_buf_ptr, &recv_id, &dead_id, &num_chunks, &chunk_size, 1);
        if (err != 0) {
          printf("%s:%d: ERROR: gran_stencil_init returned error %d\n", __FILE__,__LINE__,err);
          exit(-1);
        }
      }
      else
      {
        
        
        err = gran_stencil_init(&gran_req, (void***) &recv_buf_ptr, (void**)send_buf_ptr, &recv_id, &dead_id, &num_chunks, &chunk_size, 1);
        if (err != 0) {
          printf("%s:%d: ERROR: gran_stencil_init returned error %d\n", __FILE__,__LINE__,err);
          exit(-1);
        }
        err = gran_stencil_init(&gran_req2, (void***) &dummy_buf_ptr, (void**)send_buf_ptr, &dead_id, &send_id, &num_chunks, &chunk_size, 1);
        if (err != 0) {
          printf("%s:%d: ERROR: gran_stencil_init returned error %d\n", __FILE__,__LINE__,err);
          exit(-1);
        }
      }

      t_start = std::clock();

      #ifdef TIME_CALLS
        for (int i = 0; i < iters; ++i) {
          call_starts_times[i] = 0;
          call_startr_times[i] = 0;
          call_ends_times[i] = 0;
          call_endr_times[i] = 0;
        }
        for (int i = 0; i < iters; ++i) {
          call_readys_times[i] = (uint64_t *)malloc(num_chunks * sizeof(uint64_t));
          if (call_readys_times[i] == NULL) {
            printf("%s:%d: Unable to allocate call_readys_times[i] array for i = %d\n",__FILE__,__LINE__,i);
            exit(-1);
          }
        } 
      #endif
          
      for(int iter = 0; iter < iters; iter++)
      {
        if(verify)
        {
          for(int element = 0; element < num_chunks; element++)
          {
            ((long*)send_buf)[element] = iter + element + me;
          }
        }
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_starts_start);
        #endif
        gran_stencil_start(gran_req);
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_starts_end);
          call_starts_times[iter] = timespec_diff(&t_starts_start, &t_starts_end);
        #endif

        if (me == 0) {
          for(int k = 0; k < num_chunks; ++k) {
            #ifdef DO_WORK
              for (int f = 0; f < NUM_MM; ++f) {
                for (int i = 0; i < I; ++i) {
                  for (int k = 0; k < K; ++k) {
                    for (int j = 0; j < J; ++j) {
                      C[i][k] += A[i][j] * B[j][k];
                    }
                  }
                }
              }
            #endif
            #ifdef TIME_CALLS
              clock_gettime(CLOCK_MONOTONIC, &t_readys_start);
            #endif
            gran_stencil_ready(gran_req, 0, k);
            #ifdef TIME_CALLS
              clock_gettime(CLOCK_MONOTONIC, &t_readys_end);
              call_readys_times[iter][k] = timespec_diff(&t_readys_start, &t_readys_end);
            #endif
          }
        }  
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_ends_start);
        #endif
        gran_stencil_end(gran_req);
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_ends_end);
          call_ends_times[iter] = timespec_diff(&t_ends_start, &t_ends_end);
        #endif

        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_startr_start);
        #endif
        gran_stencil_start(gran_req2);
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_startr_end);
          call_startr_times[iter] = timespec_diff(&t_startr_start, &t_startr_end);
        #endif
        if (me == 1) {
          for(int k = 0; k < num_chunks; ++k) {
            gran_stencil_ready(gran_req2, 0, k);
          }
        }
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_endr_start);
        #endif
        gran_stencil_end(gran_req2);
        #ifdef TIME_CALLS
          clock_gettime(CLOCK_MONOTONIC, &t_endr_end);
          call_endr_times[iter] = timespec_diff(&t_endr_start, &t_endr_end);
        #endif

        if(verify)
        {
          for(int element = 0; element < num_chunks; element++)
          {
            if(((long*)recv_buf)[element] != iter + element + (me + 1) % world)
            {
              printf("Comm Error, output doesn't match expected\n");
            }
          }
        }
      } 
      t_end = std::clock();

      #ifdef TIME_CALLS
        if (0 == me) {
          for (int i = 0; i < iters; ++i) {
            printf("buf_B: %d chunks: %d chunk_B: %d iter: %d starts_t_ns: %lu ends_t_ns: %lu startr_t_ns: %lu endr_t_ns: %lu readys_t_ns:",cur_buf_size, num_chunks, chunk_size, i, call_starts_times[i], call_ends_times[i], call_startr_times[i], call_endr_times[i]);
            for (int k=0; k < num_chunks; ++k) {
              printf(" %lu", call_readys_times[i][k]);
            }
            printf("\n");
          }
          fflush(stdout);
        }
        for (int i = 0; i < iters; ++i ) {
          free(call_readys_times[i]);
          call_starts_times[i] = 0;
          call_startr_times[i] = 0;
          call_ends_times[i] = 0;
          call_endr_times[i] = 0;
        }
      #endif

      gran_stencil_finish(gran_req); 
      gran_stencil_finish(gran_req2);

      if(me == 0)
      {
        t_total = (t_end-t_start)/(double) CLOCKS_PER_SEC;
        #ifdef TIME_CALLS
          {}
        #else
        
          printf("iterations: %d, buf_size: %d B num_chunks: %d bytes_per_chunk: %d total_time_msecs: %.6f msecs_per_iteration: %.6f bytes_per_second: %.5f MiB_per_second: %.5f\n", iters, cur_buf_size, num_chunks, chunk_size, t_total*1000, (t_total*1000)/iters, (cur_buf_size*(double)iters)/t_total, ((cur_buf_size*(double)iters)/t_total)/1048576);
        #endif      
        fflush(stdout);
      }
    } 
  } 

  #ifdef TIME_CALLS
    free(call_starts_times);
    free(call_startr_times);
    free(call_ends_times);
    free(call_endr_times);
    free(call_readys_times); 
  #endif

  return 0;
}

int kernel_gran_bench_clean(void* _kernel_data)
{
  return 0; 
}

int kernel_gran_bench_load(kernelfuncs* funcs)
{
  funcs->init  = &kernel_gran_bench_init;
  funcs->run   = &kernel_gran_bench_run;
  funcs->clean = &kernel_gran_bench_clean;
  return 0;
}

