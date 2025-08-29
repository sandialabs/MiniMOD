#include"earlycoll.h"

void MPIX_Start(MPIX_Request *request) {
    request->operation_func(request, request->sendbuf, request->rank, request->chunk_size, request->world_size);
}
