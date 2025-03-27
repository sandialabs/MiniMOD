/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_GRANULARITY
#define MINIMOD_GRANULARITY

#include <stddef.h>
#include "comm.h"

typedef void gran_request;

typedef struct
{
  int (*init)               ();
  int (*clean)              ();
  int (*thread_req)         ();
  int (*stencil_init)       (gran_request** request, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int num_directions);
  int (*stencil_start)      (gran_request*);
  int (*stencil_ready)      (gran_request*, int, int);
  int (*stencil_end)        (gran_request*);
  int (*stencil_finish)     (gran_request*);
  int (*stencil_start_u)      (gran_request*);
  int (*stencil_ready_u)      (gran_request*, int direction, int entry, int size);
  int (*stencil_match_u)      (gran_request*, int direction, int n, int size);
  int (*stencil_end_u)        (gran_request*);

} gran_funcs;

int gran_init               ();
int gran_load               (char* name);
int gran_thread_req         (); 
int gran_stencil_init       (gran_request** request, void*** recv_bufs, void** send_bufs, int* recv_ids, int* send_ids, int* num_entries, int* entry_size, int num_directions);
int gran_stencil_start      (gran_request*);
int gran_stencil_ready      (gran_request*, int direction, int entry);
int gran_stencil_end        (gran_request*);
int gran_stencil_finish     (gran_request*); 

#include"granularity/gran_bulk.h" 
#include"granularity/gran_fine.h"
#include"granularity/gran_bins.h"
#include"granularity/gran_bins_tmp.h"
#include"granularity/gran_funn.h"
#include"granularity/gran_offload.h"
#include"granularity/gran_offload_bins.h"
#include"granularity/gran_offload_bins_nminus.h"

#endif
