/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "kernel.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

kernelfuncs * active_kernel;


int kernel_parse_args_default(void* data,int argc, char** argv)
{
  printf("Not parsing args today!");
  
  return 0;
}

int kernel_load(char* name, int id, int job_size)
{
  active_kernel = (kernelfuncs*) malloc(sizeof(kernelfuncs));
  active_kernel->parse_args = &kernel_parse_args_default; 
  active_kernel->id = id;
  active_kernel->job_size = job_size;
  if(id == 0) 
  {
    printf("KERNEL - Searching for ");
    printf("%s\n",name);
  }
  if(strcmp(name,"sample") == 0) 
  {
    if(id == 0) printf("KERNEL - Loading 'sample'\n");
    return kernel_sample_load(active_kernel);
  }

#ifdef HAVE_OMP
  if(strcmp(name,"sample_thread") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'sample'\n");
    return kernel_sample_thread_load(active_kernel);
  }
#endif

  if(strcmp(name,"sample_threads") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'sample'\n");
    return kernel_sample_threads_load(active_kernel);
  }

  if(strcmp(name,"cmb") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'cmb'\n");
    return kernel_cmb_load(active_kernel);
  }

  if(strcmp(name,"collective_bench") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'collective_bench'\n");
    return kernel_collective_bench_load(active_kernel);
  }

  if(strcmp(name,"channel_bench") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'channel_bench'\n");
    return kernel_channel_bench_load(active_kernel);
  }

  if(strcmp(name,"gran_bench") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'gran_bench'\n");
    return kernel_gran_bench_load(active_kernel);
  }

  if(strcmp(name,"partsend") == 0)
  {
    if(id == 0) printf("KERNEL - Loading 'partsend'\n");
    return kernel_partsend_load(active_kernel);
  }

#ifdef HAVE_MINIMD
  if(strcmp(name,"minimd") == 0) 
  {
    using namespace Kernel_MiniMD;
    if(id == 0) printf("KERNEL - Loading 'minimd'\n");
    return kernel_minimd_load(active_kernel);
  }
#endif

#ifdef HAVE_MINIFE
  if(strcmp(name,"minife") == 0) 
  {

    if(id == 0) printf("KERNEL - Loading 'minife'\n");
    return kernel_minife_load(active_kernel);
  }
#endif

  if(id == 0) printf("KERNEL - Not found\n");
  return -2;
}

int kernel_parse_args(int argc, char** argv)
{
  return active_kernel->parse_args(&active_kernel->data, argc, argv);
}

int kernel_run()
{
  int err = 0;

  err = active_kernel->init(&active_kernel->data, active_kernel->id, active_kernel->job_size);
  if(err) return err;
  err = active_kernel->run(active_kernel->data);
  if(err) return err;
  err = active_kernel->clean(active_kernel->data);


  return err;
}

int kernel_clean()
{
  free(active_kernel);
  active_kernel = 0;
  return 0;
}
