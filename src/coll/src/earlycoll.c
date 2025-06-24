/* MiniMOD  1.0 - A modular communication benchmark 
 * Copyright (2025) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */


#include"earlycoll.h"

void MPIX_Start(MPIX_Request *request) {
    request->operation_func(request, request->sendbuf, request->rank, request->chunk_size, request->world_size);
}
