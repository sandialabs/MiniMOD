/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifdef HAVE_OMP

#include <omp.h>
#include "kernel_config.h"
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


#define FIVE_POINT 4
#define SEVEN_POINT 6
#define NINE_POINT 8
#define TWENTYSEVEN_POINT 26

#define TWO-DIMENSION 2
#define THREE-DIMENSION 3

#define HAS-DIAGONALS 1
#define NO-DIAGONALS 0
/*
 * In this configurable kernel, we run an abstracted computation loop, using a variety of stencils.
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
  int proc_z;
  double*** main_array; 
  double*** next_array; 
  int halo_x; 
  int halo_y; 
  int halo_z; 
  int threads; 
  int part_size; 
  int iters;
  int stencil;
  int three_dimensional;
  int has_diagonals;
  int dim_scaling;
  int stencil_length;

  
  int parts_x;
  int parts_y;
  int parts_z;

  
  int* neighbor_send_id;
  int* neighbor_recv_id;
  int* neighbor_entries;
  double** neighbor_recv;
  double** neighbor_prev_recv;
  double** neighbor_send;
  int**** dir_of;
  int**** elem_ready_of;

  
  gran_request* gran_req;
  gran_request* prev_gran_req;
} kernel_config_data;



typedef struct
{
  int halo_x;
  int halo_y;
  int halo_z;
  int threads;
  int part_size;
} kernel_config_input_arg_data;


kernel_config_input_arg_data* TEMP_INPUT_ARG_DATA;

int kernel_config_parse_args(void* _data, int argc, char** argv)
{
  kernel_config_input_arg_data* new_input_arg_data = (kernel_config_input_arg_data*) malloc(sizeof(kernel_config_input_arg_data));

  int c;
  int x = 0;
  int y = 0;
  int z = 0;
  int t = 0;
  int p = 0;
  
  
  char **foobar = (char **)malloc(sizeof(char*) * (argc+1));
  foobar[0] = (char*)"foo"; 
  for (int i = 1;i <= argc; ++i) {
    foobar[i] = argv[i-1];
  }
  
  while ((c = getopt(argc+1, foobar, "x:y:z:t:p:")) != -1) {
      
      
    switch (c) {
      case 'x':
        x = atoi(optarg);
        break;
      case 'y':
        y = atoi(optarg);
        break;
      case 'z':
        z = atoi(optarg);
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
  if (0 == z) {
    printf("ERROR: no halo_z value provided, or value is 0 (use -z option)\n");
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

  new_input_arg_data->halo_x = x;
  new_input_arg_data->halo_y = y;
  new_input_arg_data->halo_z = z;
  new_input_arg_data->threads = t;
  new_input_arg_data->part_size = p;

  TEMP_INPUT_ARG_DATA = new_input_arg_data;
  free(foobar);

  return(0);
}

/*
 * Useful debug function to printout local region 
 */
void kernel_config_print_matrix(kernel_config_data* kernel_data)
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
          if(kernel_data->main_array[j][k][1] == 300)
          {
            printf("   ");
          }
          else
          {
            printf("%.0f ",kernel_data->main_array[j][k][1]/10);
          }
        }
        printf("\n");
      }
    }
    comm_barrier();
  }
}

void kernel_config_write(kernel_config_data* kernel_data, int iter)                                                             
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
int kernel_config_process_map(int num_procs, int* x, int* y)
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

int kernel_config_binning_map(void* _kernel_data)
{
  kernel_config_data* kernel_data = (kernel_config_data*)_kernel_data;

  
  kernel_data->dir_of = (int****)malloc(sizeof(int***)*3);
  kernel_data->elem_ready_of = (int****)malloc(sizeof(int***)*3);
  for (int x = -1; x <= 1; x++)
  {
    kernel_data->dir_of[x+1] = (int***)malloc(sizeof(int**)*3);
    kernel_data->elem_ready_of[x+1] = (int***)malloc(sizeof(int**)*3);
    for (int y = -1; y <= 1; y++)
    {
      kernel_data->dir_of[x+1][y+1] = (int**)malloc(sizeof(int*)*3);
      kernel_data->elem_ready_of[x+1][y+1] = (int**)malloc(sizeof(int*)*3);
      for (int z = -1; z <= 1; z++)
      {
        kernel_data->dir_of[x+1][y+1][z+1] = (int*)malloc(sizeof(int)*kernel_data->stencil_length);
        kernel_data->elem_ready_of[x+1][y+1][z+1] = (int*)calloc(sizeof(int),kernel_data->stencil_length);
      }
    }
  }

  
  int count = 0;
  for (int k = -1; k <= 1; k++)
  {
    for (int j = -1; j <= 1; j++)
    {
      for (int i = -1; i <= 1; i++)
      {
        
        int dim = abs(i) + abs(j) + abs(k);
        if (dim == 0 || (abs(k) && !kernel_data->three_dimensional) || (dim > 1 && !kernel_data->has_diagonals))
        {
          for(int d = 1; d <= kernel_data->stencil_length; d++)
          {
            kernel_data->dir_of[i+1][j+1][k+1][d-1] = -1;
          }
        }
	else
        {
          for(int d = 1; d <= kernel_data->stencil_length; d++)
          {
            kernel_data->dir_of[i+1][j+1][k+1][d-1] = count++;
          }
        }
      }
    }
  }
  return 0;
}

int kernel_config_get_neighbor(void* _kernel_data, int id, int dist, int x, int y, int z)
{
  kernel_config_data* kernel_data = (kernel_config_data*)_kernel_data;
  int ret = id;

  if(z > 0)
  {
    if((ret + ((kernel_data->proc_x * kernel_data->proc_y) * dist)) >= kernel_data->job_size)
    {
      return -1;
    }
    else
    {
      ret = ret + ((kernel_data->proc_x * kernel_data->proc_y) * dist);
    }
  }
  else if(z < 0)
  {
    if((ret - ((kernel_data->proc_x * kernel_data->proc_y) * dist)) < 0)
    {
      return -1;
    }
    else
    {
      ret = ret - ((kernel_data->proc_x * kernel_data->proc_y) * dist);
    }
  }

  if(y > 0)
  {
    if((((kernel_data->proc_x * kernel_data->proc_y) - (ret % (kernel_data->proc_x * kernel_data->proc_y))) <= (kernel_data->proc_x * dist)) || ((ret + (kernel_data->proc_x * dist)) >= kernel_data->job_size))
    {
      return -1;
    }
    else
    {
      ret = ret + (kernel_data->proc_x * dist);
    }
  }
  else if(y < 0)
  {
    if((ret % (kernel_data->proc_x * kernel_data->proc_y)) < (kernel_data->proc_x * dist))
    {
      return -1;
    }
    else
    {
      ret = ret - (kernel_data->proc_x * dist);
    }
  }

  if(x > 0)
  {
    if(((((ret % (kernel_data->proc_x * kernel_data->proc_y)) % kernel_data->proc_x) + dist) >= kernel_data->proc_x) || ((ret + dist) >= kernel_data->job_size))
    {
      return -1;
    }
    else
    {
      ret = ret + dist;
    }
  }
  else if(x < 0)
  {
    if(((ret % (kernel_data->proc_x * kernel_data->proc_y)) % kernel_data->proc_x) < dist)
    {
      return -1;
    }
    else
    {
      ret = ret - dist;
    }
  }

  return ret;
}

int kernel_config_init(void** kernel_data, int id, int job_size)
{
  kernel_config_data* new_data = (kernel_config_data*) malloc(sizeof(kernel_config_data));

  new_data->id = id;
  new_data->job_size = job_size;

  new_data->iters = 100;
  new_data->stencil = FIVE_POINT;
  new_data->stencil_length = 1;
  new_data->has_diagonals = 0;
  new_data->dim_scaling = 1;

  
  new_data->halo_x = TEMP_INPUT_ARG_DATA->halo_x;
  new_data->halo_y = TEMP_INPUT_ARG_DATA->halo_y;
  new_data->halo_z = TEMP_INPUT_ARG_DATA->halo_z;
  new_data->threads = TEMP_INPUT_ARG_DATA->threads;
  new_data->part_size = TEMP_INPUT_ARG_DATA->part_size;
  free(TEMP_INPUT_ARG_DATA);

  
  
  kernel_config_process_map(job_size, &new_data->proc_x, &new_data->proc_y); 
  if(id == 0)
  {
    printf("For job of size %d we will use a process map of  %d x %d\n", job_size, new_data->proc_x, new_data->proc_y);
    
    printf("halo_x      : %d\n",new_data->halo_x);
    printf("halo_y      : %d\n", new_data->halo_y);
    printf("threads     : %d\n", new_data->threads);
  }
  new_data->proc_z = 1;
  new_data->three_dimensional = new_data->proc_z > 1;

  
  new_data->neighbor_send_id = (int*)malloc(sizeof(int)*new_data->stencil);
  new_data->neighbor_recv_id = (int*)malloc(sizeof(int)*new_data->stencil);
  new_data->neighbor_entries = (int*)malloc(sizeof(size_t)*new_data->stencil);

  int count = 0;


    for (int z = -1; z <= 1; z++)
    {
      for (int y = -1; y <= 1; y++)
      {
        for (int x = -1; x <= 1; x++)
        {
          
          int dim = abs(x) + abs(y) + abs(z);
	  if (dim == 0 || (abs(z) && !new_data->three_dimensional) || (dim > 1 && !new_data->has_diagonals)) continue;
          
          int tmp_size =  new_data->halo_x * new_data->halo_y * new_data->halo_z;
          if(new_data->dim_scaling)
          {
	    if(abs(x)) tmp_size = tmp_size / new_data->halo_x;
	    if(abs(y)) tmp_size = tmp_size / new_data->halo_y;
	    if(abs(z)) tmp_size = tmp_size / new_data->halo_z;
          }
	  
	  for(int dist = 1; dist <= new_data->stencil_length; dist++)
          {
	    new_data->neighbor_send_id[count] = kernel_config_get_neighbor(new_data, id, dist, x, y, z);
	    new_data->neighbor_recv_id[count] = kernel_config_get_neighbor(new_data, id, dist, -x, -y, -z);
            printf("(%d)Howdy! send:(%d)\n",id,kernel_config_get_neighbor(new_data, id, 1, x, y, z));
            new_data->neighbor_entries[count] = tmp_size;
	    count++;
          }
        }
      }
    }


  
  kernel_config_binning_map(kernel_data);

  
  
  if(new_data->halo_x % new_data->part_size != 0)
  {
    printf("Error: halo_x doesn't divide evenly into part size");
  }
  new_data->parts_x = new_data->halo_x/new_data->part_size;
  
  if(new_data->halo_y % new_data->part_size != 0)
  { 
    printf("Error: halo_y doesn't divide evenly into part size");
  }
  new_data->parts_y = new_data->halo_y/new_data->part_size;

  if(new_data->three_dimensional && (new_data->halo_z % new_data->part_size != 0))
  { 
    printf("Error: halo_z doesn't divide evenly into part size");
  }
  new_data->parts_z = new_data->halo_z/new_data->part_size;

  

  
  new_data->main_array = (double***)calloc(new_data->halo_x,sizeof(double**));
  new_data->next_array = (double***)calloc(new_data->halo_x,sizeof(double**));
  for(int i = 0; i < new_data->halo_x; i++)
  {
    new_data->main_array[i] = (double**)calloc(new_data->halo_y,sizeof(double*));
    new_data->next_array[i] = (double**)calloc(new_data->halo_y,sizeof(double*));
    for(int j = 0; j < new_data->halo_y; j++)
    {
      new_data->main_array[i][j] = (double*)calloc(new_data->halo_z,sizeof(double));
      new_data->next_array[i][j] = (double*)calloc(new_data->halo_z,sizeof(double));
    }
  }

  
  new_data->neighbor_send = (double**)calloc(new_data->stencil,sizeof(double*));
  new_data->neighbor_recv = (double**)calloc(new_data->stencil,sizeof(double*));
  new_data->neighbor_prev_recv = (double**)calloc(new_data->stencil,sizeof(double*));
  for (int i = 0; i < new_data->stencil; i++)
  {
    new_data->neighbor_send[i] = (double*)calloc(new_data->neighbor_entries[i],sizeof(double));
  }

  
  #if WINDOW == 1
    new_data->gran_req = NULL;
  #endif

  for(int i = 0; i < new_data->halo_x; i++)
  {
    for(int j = 0; j < new_data->halo_y; j++)
    {
      for(int k = 0; k < new_data->halo_z; k++)
      {
        new_data->main_array[i][j][k] = rand();
      }
    }
  }
 
  
  if(DEBUG) kernel_config_print_matrix(new_data);
  if(WRITE) kernel_config_write(new_data, -1);

  
  
  
  int* entry_num = (int*)malloc(sizeof(int) * new_data->stencil);
  int* entry_size = (int*)malloc(sizeof(int) * new_data->stencil);
  for (int i = 0; i < new_data->stencil; i++)
  {
    entry_num[i] = (new_data->neighbor_entries[i] == 1 ? 1 : new_data->neighbor_entries[i] / new_data->part_size);
    entry_size[i] = (new_data->neighbor_entries[i] == 1 ? sizeof(double) : new_data->part_size * sizeof(double));
  }

  int err;
  err = gran_stencil_init(&new_data->gran_req, (void***)&new_data->neighbor_recv, (void**)new_data->neighbor_send, new_data->neighbor_recv_id, new_data->neighbor_send_id, entry_num, entry_size, new_data->stencil);
  if (err) return err;
  err = gran_stencil_init(&new_data->prev_gran_req, (void***)&new_data->neighbor_recv, (void**)new_data->neighbor_send, new_data->neighbor_recv_id, new_data->neighbor_send_id, entry_num, entry_size, new_data->stencil);
  if (err) return err;

  free(entry_num);
  free(entry_size);

  *kernel_data = new_data;
  return 0; 
}


typedef struct {
  kernel_config_data* kernel_data;
} kernel_config_args;



int kernel_config_func(int thread_num, int num_threads, void* args)
{
  kernel_config_args* kernel_args = (kernel_config_args*) args;
  kernel_config_data* kernel_data = kernel_args->kernel_data;

  int dist = 1;
  int x, y, z;

  struct timespec req,rem;
  req.tv_nsec = 100;
  if (thread_num == 0) req.tv_nsec = 200;

  for(int i = thread_num;  i < kernel_data->halo_x; i += num_threads) 
  {
    for (int j = 0; j < kernel_data->halo_y; j++)
    {
      for (int k = 0; k < kernel_data->halo_y; k++)
      {
        
        clock_nanosleep(CLOCK_REALTIME,0,&req, &rem);

	
	
        
        x = y = z = 0;
        if(i == 0) x = -1;
        if(i == kernel_data->halo_x - 1) x = 1;
        if(j == 0) y = -1;
        if(j == kernel_data->halo_y - 1) y = 1;
        if(k == 0) z = -1;
        if(k == kernel_data->halo_z - 1) z = 1;
       
	
        if(!x && !y & !z) continue;

        if(x)
        {
          int d = kernel_data->dir_of[x+1][1][1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[x+1][1][1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }

        if(y)
        {
          int d = kernel_data->dir_of[1][y+1][1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[1][y+1][1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }

        if(z)
        {
          int d = kernel_data->dir_of[1][1][z+1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[1][1][z+1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }

        if(x && y)
        {
          int d = kernel_data->dir_of[x+1][y+1][1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[x+1][y+1][1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }

        if(x && z)
        {
          int d = kernel_data->dir_of[x+1][1][z+1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[x+1][1][z+1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }

        if(y && z)
        {
          int d = kernel_data->dir_of[1][y+1][z+1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[1][y+1][z+1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }

        if(x && y && z)
        {
          int d = kernel_data->dir_of[x+1][y+1][z+1][dist - 1];
          if(d != -1)
          {
            int r = kernel_data->elem_ready_of[x+1][y+1][z+1][dist - 1]++;
            if (kernel_data->neighbor_entries[d] == 1 || ((kernel_data->neighbor_entries[d] / kernel_data->part_size) % r) == 0)
            {
              gran_stencil_ready(kernel_data->gran_req, d, ((kernel_data->neighbor_entries[d] / r) - 1));
            }
          }
        }
      }
    }
  }

  return 0;
}

int kernel_config_run(void* _kernel_data)
{
  kernel_config_data* kernel_data = (kernel_config_data*)_kernel_data;
  kernel_config_args* args = (kernel_config_args*)malloc(sizeof(kernel_config_args));
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

  
  threads_init(&t_req, kernel_data->threads); 

  
  

  for (int t = 0; t < kernel_data->iters; t++)
  {
    

    
    
    {
      gran_stencil_start(kernel_data->gran_req);
    }

    args->kernel_data = kernel_data;

    
    for (int x = -1; x <= 1; x++)
    {
      for (int y = -1; y <= 1; y++)
      {
        for (int z = -1; z <= 1; z++)
        {
          for (int d = 1; d <= kernel_data->stencil_length; d++)
          {
            kernel_data->elem_ready_of[x+1][y+1][z+1][d-1] = 0;
          }
        }
      }
    }

    
    threads_run(t_req, &kernel_config_func, args);

    
    double*** swap = kernel_data->next_array;
    kernel_data->next_array = kernel_data->main_array;
    kernel_data->main_array = swap;

    
    if(DEBUG) kernel_config_print_matrix(kernel_data);
    if(WRITE) kernel_config_write(kernel_data, t);

    #if TIMING == 1
      
      t_comm_start = std::chrono::steady_clock::now();
    #endif

    
    #if WINDOW == 1
      gran_stencil_end(kernel_data->gran_req);
    #endif

    
    





















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

int kernel_config_clean(void* _kernel_data)
{
  
 
  kernel_config_data* kernel_data = (kernel_config_data*)_kernel_data;
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

int kernel_config_load(kernelfuncs* funcs)
{
  funcs->init  = &kernel_config_init;
  funcs->parse_args = &kernel_config_parse_args;
  funcs->run   = &kernel_config_run;
  funcs->clean = &kernel_config_clean;
  return 0;
}

#endif
