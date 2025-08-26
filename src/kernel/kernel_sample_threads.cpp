/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifdef HAVE_OMP

#include <omp.h>
#include "kernel_sample_threads.h"
#include <string.h>
#include <sys/stat.h>
#include <ctime> 
#include <ratio> 
#include <chrono>
#include <iostream>
#include <atomic>


#define DEBUG 0 
#define WRITE 0
#define TIMING 1
#define WINDOW 1
#define STENCIL 0
/*
 * In this sample kernel, we run a simple 2D diffusion loop, using a 5-point stencil.
 *
 * This code is based on https:
 * 
 * Please note the naming convetion of <module_type>_<module_name>_<function_name>: when building a new module 
 * please follow this naming convetion to avoid conflicts. 
 */






typedef struct
{
  
  int id;
  int job_size;
  int proc_x;
  int proc_y;
  double** main_array; 
  double** next_array; 
  int halo_x;
  int halo_y;
  int iters;

  
  int total_halo_x;
  int parts_x;
  double dx;
  double dx2;
  int total_halo_y;
  int parts_y;
  std::atomic<int>* part_north_bins;
  std::atomic<int>* part_south_bins;
  int part_y_thresh;
  double dy2;
  double dy;
  double dt;
  double diff;

  
  int north_id;
  double* north_recv;
  double* north_send;
  int east_id;
  double* east_recv;
  double* east_send;
  int south_id;
  double* south_recv;
  double* south_send;
  int west_id;
  double* west_recv;
  double* west_send;

  double* prev_north_recv;
  double* prev_east_recv;
  double* prev_south_recv;
  double* prev_west_recv;


  
  gran_request* gran_req;
  gran_request* prev_gran_req;
} kernel_sample_threads_data;



int kernel_sample_threads_halo_x;
int kernel_sample_threads_halo_y;
int kernel_sample_threads_threads;
int kernel_sample_threads_part_size;

int kernel_sample_threads_parse_args(void* _data, int argc, char** argv)
{

  int c;
  int x = 0;
  int y = 0;
  int t = 0;
  int p = 0;
  
  
  char **foobar = (char **)malloc(sizeof(char*) * (argc+1));
  foobar[0] = (char*)"foo"; 
  for (int i = 1;i <= argc; ++i) {
    foobar[i] = argv[i-1];
  }
  
  while ((c = getopt(argc+1, foobar, "x:y:t:p:")) != -1) {
      
      
    switch (c) {
      case 'x':
        x = atoi(optarg);
        break;
      case 'y':
        y = atoi(optarg);
        break;
      case 't':
        t = atoi(optarg);
        break;
      case 'p':
        p = atoi(optarg);
        break;
      case '?':
        printf("Unknown command line argument (%c)\n", c);
        exit(-1);
      default:
        printf("Kernel requires three arguments: -x <halo_x> -y <halo_x> -t <threads>\n");
        exit(-1);
    }
  }

  if (0 == x) {
    printf("ERROR: no halo_x value provided, or value is 0 (use -x option)\n");
    exit(-1);
  }
  if (0 == y) {
    printf("ERROR: no halo_y value provided, or value is 0 (use -y option)\n");
    exit(-1);
  }
  if (0 == t) {
    printf("ERROR: no halo_t value provided, or value is 0 (use -t option)\n");
    exit(-1);
  }
  if (0 == p) {
    printf("ERROR: no halo_p value provided, or value is 0 (use -p option)\n");
    exit(-1);
  }

  kernel_sample_threads_halo_x = x;
  kernel_sample_threads_halo_y = y;
  kernel_sample_threads_threads = t;
  kernel_sample_threads_part_size = p;

  
  
  









  
  free(foobar);

  return(0);
}

/*
 * Useful debug function to printout local region 
 */
void kernel_sample_threads_print_matrix(kernel_sample_threads_data* kernel_data)
{
  int i;
  for(i = 0; i < kernel_data->job_size; i++)
  {
    if(kernel_data->id == i)
    {
      int j, k;
      for(j = 0; j < kernel_data->halo_x; j++) 
      {
        for(k = 0; k < kernel_data->halo_y; k++) 
        {
          if(kernel_data->main_array[j][k] == 300)
          {
            printf("   ");
          }
          else
          {
            printf("%.0f ",kernel_data->main_array[j][k]/10);
          }
        }
        printf("\n");
      }
    }
    comm_barrier();
  }
}

void kernel_sample_threads_write(kernel_sample_threads_data* kernel_data, int iter)                                                             
{
  
  char iter_s[20];
  sprintf(iter_s, "%d", iter);
  char rank_s[20];
  sprintf(rank_s, "%d", kernel_data->id);

  mkdir("out",S_IRWXU);
  char fn[44]; 
  strcat(strcat(strcat(strcpy(fn,"out/"), rank_s), "_"), iter_s);
  FILE* fp = fopen(fn,"w+");
  if (fp == NULL)
  {
    printf("BLOODY MURDER : %s\n", fn);
  }
  else 
  {
    int i;
    for(i = 0; i < kernel_data->job_size; i++)
    {
      if(kernel_data->id == i)
      {
        int j, k;
        for(j = 0; j < kernel_data->halo_x; j++)
        {
          for(k = 0; k < kernel_data->halo_y; k++)
          {
            if(k == (kernel_data->halo_y - 1))
            {
              fprintf(fp, "%.0f", kernel_data->main_array[j][k]);
            }
            else
            {
              fprintf(fp, "%.0f,", kernel_data->main_array[j][k]);
            }
          }
          fprintf(fp,"\n");
        }
      }
      comm_barrier();
    }
    fclose(fp);
  }
}


/*
 * This function takes an input of num_procs and returns the rectangle with the least perimeter
 */
int kernel_sample_threads_process_map(int num_procs, int* x, int* y)
{
  int i;
  int best = 1;
  
  for(i = 1; i < num_procs/i; i++)
  {
    if(num_procs % i == 0)
    {
      best = i;
    }
  }
  if(i == num_procs/i && num_procs % i == 0) best = i; 
  *x = num_procs/best;
  *y = best;
  return 0;
}

int kernel_sample_threads_init(void** kernel_data, int id, int job_size)
{
  int i, j;
  kernel_sample_threads_data* new_data = (kernel_sample_threads_data*) malloc(sizeof(kernel_sample_threads_data));

  new_data->id = id;
  new_data->job_size = job_size;

  new_data->iters = 100;

  
  kernel_sample_threads_process_map(job_size, &new_data->proc_x, &new_data->proc_y); 
  if(id == 0)
  {
    printf("For job of size %d we will use a process map of  %d x %d\n", job_size, new_data->proc_x, new_data->proc_y);
    
    printf("halo_x      : %d\n",kernel_sample_threads_halo_x);
    printf("halo_y      : %d\n", kernel_sample_threads_halo_y);
    printf("threads     : %d\n", kernel_sample_threads_threads);
  }

  
  if(id < new_data->proc_x)
    new_data->north_id = -1;
  else
    new_data->north_id = id - new_data->proc_x;

  if(id % new_data->proc_x == new_data->proc_x -1)
    new_data->east_id = -1;
  else
    new_data->east_id = id + 1;

  if(id % new_data->proc_x == 0)
    new_data->west_id = -1;
  else
    new_data->west_id = id - 1;

  if(id/new_data->proc_x == new_data->proc_y - 1)
    new_data->south_id = -1;
  else
    new_data->south_id = id + new_data->proc_x;

  
  new_data->halo_x = kernel_sample_threads_halo_x;
  if(new_data->halo_x % kernel_sample_threads_part_size != 0)
  {
    printf("Error: halo_x doesn't divide evenly into part size");
  }
  new_data->parts_x = kernel_sample_threads_halo_x/kernel_sample_threads_part_size;
  
  new_data->halo_y = kernel_sample_threads_halo_y;
  if(new_data->halo_y % kernel_sample_threads_part_size != 0)
  { 
    printf("Error: halo_y doesn't divide evenly into part size");
  }
  new_data->parts_y = kernel_sample_threads_halo_y/kernel_sample_threads_part_size;
  new_data->part_y_thresh = kernel_sample_threads_part_size;
  new_data->part_south_bins = (std::atomic<int>*)calloc(new_data->parts_y, sizeof(std::atomic<int>));
  new_data->part_north_bins = (std::atomic<int>*)calloc(new_data->parts_y, sizeof(std::atomic<int>));

  


  
  
  double width  = 10.0;
  double height = 10.0;
  
  new_data->dx = width/new_data->halo_x;
  new_data->dx2 = new_data->dx * new_data->dx;
  new_data->dy = height/new_data->halo_y;
  new_data->dy2 = new_data->dy * new_data->dy;
  
  new_data->diff = 4.0;
  
  new_data->dt = new_data->dx2 * new_data->dy2 / (2 * new_data->diff * (new_data->dx2 + new_data->dy2)); 

  

  
  new_data->main_array = (double**)calloc(new_data->halo_x,sizeof(double*));
  new_data->next_array = (double**)calloc(new_data->halo_x,sizeof(double*));
  for(i = 0; i < new_data->halo_x; i++)
  {
    new_data->main_array[i] = (double*)calloc(new_data->halo_y,sizeof(double));
    new_data->next_array[i] = (double*)calloc(new_data->halo_y,sizeof(double));
  }

  
  new_data->north_send = (double*)calloc(new_data->halo_y,sizeof(double));
  new_data->east_send = (double*)calloc(new_data->halo_x,sizeof(double));
  new_data->south_send = (double*)calloc(new_data->halo_y,sizeof(double));
  new_data->west_send = (double*)calloc(new_data->halo_x,sizeof(double));

  
  #if WINDOW == 1
    new_data->gran_req = NULL;
  #endif
 
  
  int cool = 300;
  int hot  = 500;

  int r = 2;
  int r2 = r*r;
  int cx = 2; 
  int cy = 2;



  for(i = 0; i < new_data->halo_x; i++)
  {
    for(j = 0; j < new_data->halo_y; j++)
    {
      int dist2 = (i*new_data->dx-cx) * (i*new_data->dx-cx) + (j*new_data->dy-cy) * (j*new_data->dy-cy);
      if(dist2 < r2)
        new_data->main_array[i][j] = hot;
      else
        new_data->main_array[i][j] = cool;
    }
  }
 
  
  if(DEBUG) kernel_sample_threads_print_matrix(new_data);
  if(WRITE) kernel_sample_threads_write(new_data, -1);

  
  
  
 
  
  void*** recv_bufs = (void***)malloc(sizeof(void**) * 4);
  recv_bufs[0] = (void**)&new_data->north_recv; 
  recv_bufs[1] = (void**)&new_data->east_recv; 
  recv_bufs[2] = (void**)&new_data->south_recv; 
  recv_bufs[3] = (void**)&new_data->west_recv; 
  void*** prev_recv_bufs = (void***)malloc(sizeof(void**) * 4);
  prev_recv_bufs[0] = (void**)&new_data->prev_north_recv;
  prev_recv_bufs[1] = (void**)&new_data->prev_east_recv;
  prev_recv_bufs[2] = (void**)&new_data->prev_south_recv;
  prev_recv_bufs[3] = (void**)&new_data->prev_west_recv;
  void** send_bufs = (void**)malloc(sizeof(void*) * 4);
  send_bufs[0] = (void*)new_data->south_send; 
  send_bufs[1] = (void*)new_data->west_send; 
  send_bufs[2] = (void*)new_data->north_send; 
  send_bufs[3] = (void*)new_data->east_send; 
  int* recv_ids = (int*)malloc(sizeof(int) * 4);
  recv_ids[0] = new_data->north_id;
  recv_ids[1] = new_data->east_id;
  recv_ids[2] = new_data->south_id;
  recv_ids[3] = new_data->west_id;
  int* send_ids = (int*)malloc(sizeof(int) * 4);
  send_ids[0] = new_data->south_id;
  send_ids[1] = new_data->west_id;
  send_ids[2] = new_data->north_id;
  send_ids[3] = new_data->east_id;
  int* num_entries = (int*)malloc(sizeof(int) * 4);
  num_entries[0] = new_data->parts_y;
  num_entries[1] = new_data->parts_x;
  num_entries[2] = new_data->parts_y;
  num_entries[3] = new_data->parts_x;
  int* entry_sizes = (int*)malloc(sizeof(int) * 4);
  entry_sizes[0] = kernel_sample_threads_part_size * sizeof(double);
  entry_sizes[1] = kernel_sample_threads_part_size * sizeof(double);
  entry_sizes[2] = kernel_sample_threads_part_size * sizeof(double);
  entry_sizes[3] = kernel_sample_threads_part_size * sizeof(double);


 
 
  
  int err = gran_stencil_init(&new_data->gran_req, recv_bufs, send_bufs, recv_ids, send_ids, num_entries, entry_sizes, 4);
  if (err) return err;
  err = gran_stencil_init(&new_data->prev_gran_req, prev_recv_bufs, send_bufs, recv_ids, send_ids, num_entries, entry_sizes, 4);

  free(recv_bufs);
  free(prev_recv_bufs);
  free(send_bufs);
  free(recv_ids);
  free(send_ids);
  free(num_entries);
  free(entry_sizes);

  if (err) return err;
 

  for(j = 0; j < new_data->halo_y; j++)
  {
    new_data->prev_north_recv[j] = new_data->main_array[new_data->halo_x - 1][j];
    new_data->prev_south_recv[j] = new_data->main_array[0][j];
  }

  for(i = 0; i < new_data->halo_x; i++)
  {
    new_data->prev_east_recv[i] = new_data->main_array[i][0];
    new_data->prev_west_recv[i] = new_data->main_array[i][new_data->halo_y - 1];
  }

  *kernel_data = new_data;
  return 0; 
}


typedef struct {
  kernel_sample_threads_data* kernel_data;
} kernel_sample_threads_args;



int kernel_sample_threads_func(int thread_num, int num_threads, void* args)
{
    kernel_sample_threads_args* kernel_args = (kernel_sample_threads_args*) args;
    kernel_sample_threads_data* kernel_data = kernel_args->kernel_data;
    
    
    

    for(int i = thread_num;  i < kernel_data->halo_x; i += num_threads) 
    {
      
      

      int north_ready, east_ready, south_ready, west_ready;
      for (int j = 0; j < kernel_data->halo_y; j++)
      {
        north_ready = east_ready = south_ready = west_ready = 0;
        double uxx = 0;
        double uyy = 0;

        
        
        if(i == 0)
        {
          if(kernel_data->north_id == -1)
          {
            
          }
          else
          {
            uxx += kernel_data->prev_north_recv[j] - kernel_data->main_array[i][j];
            north_ready = 1;
          }
        }
        else
        {
          uxx += kernel_data->main_array[i-1][j] - kernel_data->main_array[i][j];
        }
        
        if(j == kernel_data->halo_y - 1)
        {
          if(kernel_data->east_id == -1)
          {
            
          }
          else
          {
            uyy += kernel_data->prev_east_recv[i] - kernel_data->main_array[i][j];
            east_ready = 1;
          }
        }
        else
        {
          uyy += kernel_data->main_array[i][j+1] - kernel_data->main_array[i][j];
        }
        
        if(i == kernel_data->halo_x - 1)
        {
          if(kernel_data->south_id == -1)
          {
            
          }
          else
          {
            uxx += kernel_data->prev_south_recv[j] - kernel_data->main_array[i][j];
            south_ready = 1;
          }
        }
        else
        {
          uxx += kernel_data->main_array[i+1][j] - kernel_data->main_array[i][j];
        }
        
        if(j == 0)
        {
          if(kernel_data->west_id == -1)
          {
            
          }
          else
          {
            uyy += kernel_data->prev_west_recv[i] - kernel_data->main_array[i][j];
            west_ready = 1;
          }
        }
        else
        {
          uyy += kernel_data->main_array[i][j-1] - kernel_data->main_array[i][j];
        }
        uxx = uxx / kernel_data->dx2;
        uyy = uyy / kernel_data->dy2;
        kernel_data->next_array[i][j] = kernel_data->main_array[i][j] + kernel_data->dt * kernel_data->diff * (uxx + uyy);

        
        if(north_ready)
        {
          kernel_data->north_send[j] = kernel_data->next_array[i][j];
          #if WINDOW == 1
          
          int my_bin = j/kernel_sample_threads_part_size;
          int my_num = ++(kernel_data->part_north_bins[my_bin]);
          if(my_num == kernel_data->part_y_thresh)
          {

            gran_stencil_ready(kernel_data->gran_req, 2, my_bin);
          }
          #endif
        }
        if(east_ready)
        {
          kernel_data->east_send[i] = kernel_data->next_array[i][j];
          #if WINDOW == 1
            if(0 == (i+1) % (kernel_sample_threads_part_size))
            {
 
              gran_stencil_ready(kernel_data->gran_req, 3, i/kernel_sample_threads_part_size);
            }
          #endif
        }
        if(south_ready)
        {
          kernel_data->south_send[j] = kernel_data->next_array[i][j];
          #if WINDOW == 1
          int my_bin = j/kernel_sample_threads_part_size;
          int my_num = ++(kernel_data->part_south_bins[my_bin]);
  
          if(my_num == kernel_data->part_y_thresh)
          {

            gran_stencil_ready(kernel_data->gran_req, 0, my_bin);
          }
          #endif
        }
        if(west_ready)
        {
          kernel_data->west_send[i] = kernel_data->next_array[i][j];
          #if WINDOW == 1
            if(0 == (i+1) % (kernel_sample_threads_part_size))
            {

              gran_stencil_ready(kernel_data->gran_req, 1, i/kernel_sample_threads_part_size);
            }
          #endif
        }
      }
    }

  return 0;
}

int kernel_sample_threads_run(void* _kernel_data)
{
  kernel_sample_threads_data* kernel_data = (kernel_sample_threads_data*)_kernel_data;
  kernel_sample_threads_args* args = (kernel_sample_threads_args*)malloc(sizeof(kernel_sample_threads_args));
  comm_barrier();

  if(kernel_data->id == 0) printf("Starting computation\n");

  #if TIMING == 1
    std::chrono::steady_clock::time_point t_start, t_end, t_comm_start, t_comm_end;
    
    
    
    
    std::chrono::duration<double> t_total; 
    std::chrono::duration<double> t_comm; 
    t_total = std::chrono::steady_clock::duration::zero();
    t_comm = std::chrono::steady_clock::duration::zero();
    t_start = std::chrono::steady_clock::now();
  #endif

  threads_request* t_req;

  
  threads_init(&t_req, kernel_sample_threads_threads); 

  
  

  for (int t = 0; t < kernel_data->iters; t++)
  {
    

    
    
    {
      gran_stencil_start(kernel_data->gran_req);
    }

    args->kernel_data = kernel_data;
 
    for(int i = 0; i < kernel_data->parts_y; i++)
    {
      kernel_data->part_north_bins[i] = 0;
      kernel_data->part_south_bins[i] = 0;
    }

    
    threads_run(t_req, &kernel_sample_threads_func, args);

    
    double** swap = kernel_data->next_array;
    kernel_data->next_array = kernel_data->main_array;
    kernel_data->main_array = swap;

    
    if(DEBUG) kernel_sample_threads_print_matrix(kernel_data);
    if(WRITE) kernel_sample_threads_write(kernel_data, t);

    #if TIMING == 1
      
      t_comm_start = std::chrono::steady_clock::now();
    #endif

    
    #if WINDOW == 1
      gran_stencil_end(kernel_data->gran_req);
    #endif

    
    double* buff_swp;
    buff_swp = kernel_data->prev_north_recv;
    kernel_data->prev_north_recv = kernel_data->north_recv;
    kernel_data->north_recv = buff_swp;

    buff_swp = kernel_data->prev_south_recv;
    kernel_data->prev_south_recv = kernel_data->south_recv;
    kernel_data->south_recv = buff_swp;

    buff_swp = kernel_data->prev_east_recv;
    kernel_data->prev_east_recv = kernel_data->east_recv;
    kernel_data->east_recv = buff_swp;

    buff_swp = kernel_data->prev_west_recv;
    kernel_data->prev_west_recv = kernel_data->west_recv;
    kernel_data->west_recv = buff_swp;

    gran_request* gran_swap = kernel_data->prev_gran_req;
    kernel_data->prev_gran_req = kernel_data->gran_req;
    kernel_data->gran_req = gran_swap;

    comm_barrier();

    #if TIMING == 1
      
      t_comm_end = std::chrono::steady_clock::now();
      t_comm += (std::chrono::duration_cast<std::chrono::duration<double>>(t_comm_end - t_comm_start));
    #endif
  }
  if(kernel_data->id == 0) printf("finished computation\n");

  
  #if TIMING == 1
    if(kernel_data->id == 0)
    {
      
      t_end = std::chrono::steady_clock::now();
      t_total = (std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start));
    }

    double in = t_comm.count();
    double out = 0;

    comm_allreduce(&in, &out, 1, COMM_DATATYPE_DOUBLE, COMM_REDUCE_OP_MIN);

    if(kernel_data->id == 0)
    {
      
      printf("HEAT_DIFFUSION per-node-halo_x: %d per-node-halo_x: %d total_time: %.5e secs comm_wait_time: %.5e secs\n", kernel_data->halo_x, kernel_data->halo_y, t_total.count(), out);
      fflush(stdout);
      
    }
  #endif
  return 0;
}

int kernel_sample_threads_clean(void* _kernel_data)
{
  
 
  kernel_sample_threads_data* kernel_data = (kernel_sample_threads_data*)_kernel_data;
  int i;
  for(i = 0; i < kernel_data->halo_x; i++)
  {
    
    free(kernel_data->main_array[i]);
    free(kernel_data->next_array[i]);
  }
  free(kernel_data->main_array);
  free(kernel_data->next_array);
  free(kernel_data); 
  return 0; 
}

int kernel_sample_threads_load(kernelfuncs* funcs)
{
  funcs->init  = &kernel_sample_threads_init;
  funcs->parse_args = &kernel_sample_threads_parse_args;
  funcs->run   = &kernel_sample_threads_run;
  funcs->clean = &kernel_sample_threads_clean;
  return 0;
}

#endif
