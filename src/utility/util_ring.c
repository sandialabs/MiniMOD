/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "util_ring.h"
#include <stdlib.h>

namespace Util_Ring
{
  
  int util_ring_init(util_ring_data**  data)
  {
    util_ring_data* _data = (util_ring_data*)malloc(sizeof(util_ring_data));
    if(NULL == _data) return -1;

    util_ring_event* _ring = (util_ring_event*)calloc(RING_START_LENGTH, sizeof(util_ring_event));
    if(NULL == _ring) return -1;

    _data->length = RING_START_LENGTH;
    _data->head = 0;
    _data->tail.store(0);

    _data->ring = _ring;
    *data = _data;

    return 0;
  } 

  int util_ring_add_event(util_ring_data* data, void* info)
  {
    
    int tmp = data->tail.fetch_add(1) % data->length;
    
    if (1 == data->ring[tmp].state) 
    {
      
      printf("overwrote my ring buffer!\n");
      return -1;
    }
    data->ring[tmp].info = info;
    data->ring[tmp].state = 1;

    return 0;
  }

  int util_ring_get_event(util_ring_data* data, void** info)
  {
    
    if (0 == data->ring[data->head].state)
    {
      *info = NULL;
    } 
    
    else
    {
      *info = (void*)data->ring[data->head].info;
      data->ring[data->head].state = 0;
      data->head = (data->head + 1) % data->length;
    }

    

    return 0;
  }

} 
