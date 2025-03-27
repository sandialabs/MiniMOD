/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include<stdio.h>
#include "kernel.h"
#include "comm.h"
#include "gran.h"
#include "threads.h"
#include <stdlib.h>

void help()
{
  printf("Run as ./minimod.x <kernel_module> <comm_module> <gran_module> <thread_module> <kernel_args>\n");
  printf("kernel_args are kernel module spesific and may be reqired or optional.\n"); 
  exit(0);
}


void print_err(const char* message, int err)
{
  printf("%s\n",message);
  if(err == -1)
  {
    printf("Funtionality Under Development\n");
  } 
  else if (err == -2)
  {
    printf("Module not found\n");
  }
  else if (err == -9)
  {
    printf("Loaded comm module doesn't support threads\n");
  }
  exit(err);
}



int main(int argc, char** argv)
{
  if(argc < 4)
  {
    help();
  }
  
  int err = 0;      
  int my_id = -1;
  int job_size = -1;



 


  err = threads_load(argv[4]);
  if(err) print_err("Error in thread_load\n", err);

  err = gran_load(argv[3]);
  if(err) print_err("Error in gran_load\n", err);

  err = comm_load(argv[2]);
  if(err) print_err("Error in comm_load\n", err);

 
  err = gran_init();
  if(err) print_err("Error in gran_init\n", err);

  err = comm_init(gran_thread_req());
  if(err) print_err("Error in comm_init\n", err);

  err = comm_get_id(&my_id);  
  if(err) print_err("Error in comm_get_id\n", err);

  err = comm_get_job_size(&job_size);
  if(err) print_err("Error in comm_get_job_size\n", err);


  if(my_id == 0) 
  {
    printf("-- Starting MiniMOD --\n");
    printf("Run Command: ");
    for(int i = 0; i < argc; ++i) printf("%s ",argv[i]);
    printf("\n");
  }

  err = kernel_load(argv[1], my_id, job_size);
  if(err) print_err("Error in kernel_load\n", err);

  err = kernel_parse_args(argc - 5, argv + 5); 
  if(err) print_err("Error in kernel_parse_args\n", err);


  if(my_id == 0) printf("KERNEL Start\n");
  err = kernel_run();
  if(err) print_err("Error found in running kernel\n",err);





  err = kernel_clean();
  if(err) print_err("Error found in cleaning kernel\n", err);
  if(my_id == 0) printf("KERNEL Finish\n");



  err = comm_clean();
  if(err) print_err("Error found in cleaning comm\n", err);

  return 0;
}
