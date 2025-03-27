CXX=mpicxx
#CXX=oshcxx
#CXX=CC  # cray-shmem compiler
CXXFLAGS := -g -Wall #-O3

#Comment/uncomment these lines to change which modules compile TODO: configure
CXXFLAGS += -DHAVE_MPI -pthread -lpthread

#mpi_part = true
#open_shmem = true
#cray_shmem = true

ifdef mpi_part
  CXXFLAGS += -DHAVE_MPI_PART
endif
ifdef open_shmem
  CXXFLAGS += -DHAVE_OSHMEM -loshmem
endif
ifdef cray_shmem
  CXXFLAGS += -DHAVE_CRAY_SHMEM 
endif

CXXFLAGS += -DHAVE_OMP -fopenmp #Assuming GNU compilers

SRC := src
OBJ := obj

kernel_objects = $(OBJ)/kernel_sample.o $(OBJ)/kernel_sample_thread.o $(OBJ)/kernel_collective_bench.o $(OBJ)/kernel_channel_bench.o $(OBJ)/kernel_gran_bench.o $(OBJ)/kernel_sample_threads.o $(OBJ)/kernel_cmb.o $(OBJ)/kernel_partsend.o
kernel_headers = kernel_sample.h kernel_sample_thread.h kernel_collective_bench.h kernel_channel_bench.h kernel_gran_bench.h kernel_sample_threads.h kernel_cmb.h kernel_partsend.h

granularity_objects = $(OBJ)/gran_bulk.o $(OBJ)/gran_fine.o $(OBJ)/gran_bins.o $(OBJ)/gran_bins_tmp.o $(OBJ)/gran_funn.o $(OBJ)/gran_offload.o $(OBJ)/gran_offload_bins.o $(OBJ)/gran_offload_bins_nminus.o
granularity_headers = gran_bulk.h gran_fine.h gran_bins.h gran_bins_tmp.h gran_funn.h gran_offload.h gran_offload_bins.h gran_offload_bins_nminus.h

communication_objects = $(OBJ)/comm_mpi_pt2pt.o $(OBJ)/comm_mpi_rma.o  $(OBJ)/comm_mpi_pers.o
communication_headers = comm_mpi_pt2pt.h comm_mpi_rma.h comm_mpi_pers.h

threadding_objects = $(OBJ)/threads_openmp.o $(OBJ)/threads_pthreads.o
threadding_headers = threads_openmp.h threads_pthreads.h

utility_objects = $(OBJ)/util_queue.o $(OBJ)/util_ring.o
utility_headers = util_queue.h util_ring.h

ifdef mpi_part
  communication_objects += $(OBJ)/comm_mpi_part.o
  communication_headers += comm_mpi_part.h
endif
ifdef open_shmem
  communication_objects += $(OBJ)/comm_openshmem.o
  communication_headers += comm_openshmem.h
endif
ifdef cray_shmem
  communication_objects += $(OBJ)/comm_crayshmem.o
  communication_headers += comm_crayshmem.h
endif

objects = $(OBJ)/driver.o $(OBJ)/kernel.o $(OBJ)/gran.o $(OBJ)/comm.o $(OBJ)/threads.o $(kernel_objects) $(granularity_objects) $(communication_objects) $(threadding_objects) $(utility_objects)
headers = kernel.h gran.h comm.h threads.h $(kernel_headers) $(granularity_headers) $(communication_headers) $(threadding_headers) $(utility_headers)

minimod.x: $(objects)
	$(CXX) $(CXXFLAGS) -o minimod.x $^

$(OBJ)/%.o: %.c $(headers)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm $(objects) minimod.x

