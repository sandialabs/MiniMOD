/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_UTILITY_RING
#define MINIMOD_UTILITY_RING

#include <stdio.h>
#include <atomic>

#define RING_START_LENGTH 512 

namespace Util_Ring
{
  typedef struct _util_ring_event
  {
    volatile void* info;
    volatile int state;
    _util_ring_event* next;
  } util_ring_event;

  typedef struct
  {
    util_ring_event* ring; 
    int length;
    int head;
    std::atomic<int> tail;
  } util_ring_data;

  int util_ring_init(util_ring_data**  data);
  int util_ring_add_event(util_ring_data* data, void* info);
  int util_ring_get_event(util_ring_data* data, void** info);
} 
#endif
