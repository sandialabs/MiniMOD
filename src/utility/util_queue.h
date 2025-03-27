/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef MINIMOD_UTILITY_QUEUE
#define MINIMOD_UTILITY_QUEUE

#include <stdio.h>
#include <atomic>

namespace Util_Queue
{
  typedef struct _util_queue_event
  {
    void* info;
    _util_queue_event* next;
  } util_queue_event;

  typedef struct
  {
    util_queue_event* head;
    std::atomic<util_queue_event*> tail;
  } util_queue_data;

  int util_queue_init(util_queue_data**  data);
  int util_queue_create_event(util_queue_event** event, void* info);
  int util_queue_add_event(util_queue_data* data, util_queue_event* event);
  int util_queue_get_event(util_queue_data* data, util_queue_event** event);
} 
#endif
