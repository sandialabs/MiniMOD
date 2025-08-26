/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "util_queue.h"
#include <stdlib.h>
#include <mutex>

namespace Util_Queue
{
  
  std::mutex head_mutex;

  int util_queue_init(util_queue_data**  data)
  {
    util_queue_data* queue = (util_queue_data*)malloc(sizeof(util_queue_data));
    if(queue == NULL) return -4; 

    queue->head = NULL;
    queue->tail.store(NULL);

    *data = queue;

    return 0;
  }
  
  int util_queue_create_event(util_queue_event** event, void* info)
  {
    *event = (util_queue_event*)malloc(sizeof(util_queue_event));
    if (*event == NULL) return -1;
    (*event)->info = info;

    return 0;
  }

  int util_queue_add_event(util_queue_data* data, util_queue_event* event)
  {
    const std::lock_guard<std::mutex> lock(head_mutex);
    
    util_queue_event* prev_tail = data->tail.exchange(event);
    
    
    if (prev_tail == NULL)
    {

      data->head = event;
    } 
    
    else
    {
      prev_tail->next = event;
    }

    return 0;
  }

  int util_queue_get_event(util_queue_data* data, util_queue_event** event)
  {

    
    if(data->head == NULL)
    {
      *event = NULL;
    } 
    
    else 
    {
      const std::lock_guard<std::mutex> lock(head_mutex);
      *event = data->head;
      util_queue_event* tmp = data->head;
      
      
      
      if(data->tail.compare_exchange_strong(tmp, NULL))
      {
        data->head = NULL;
      }
      
      else
      {
        data->head = data->head->next;
      }
    }

    return 0;
  }

} 
