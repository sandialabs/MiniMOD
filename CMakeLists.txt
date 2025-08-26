# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(minimod LANGUAGES CXX)

# Default options
option(ENABLE_OPEN_SHMEM "Enable OpenSHMEM" OFF)
option(ENABLE_CRAY_SHMEM "Enable Cray SHMEM" OFF)

# Define a custom build type variable (default to Release)
if(NOT BUILD)
    set(BUILD Release)  # Default to optimized HPC build
endif()

# Print the selected build type
message(STATUS "Build type: ${BUILD}")

# Define compiler flags for different build types
if(BUILD STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "-g -Wall -pthread -lpthread -DHAVE_OMP -fopenmp")  # Debug flags
elseif(BUILD STREQUAL Release)
    set(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native -funroll-loops -fopenmp -DNDEBUG")  # Optimized HPC flags
else()
    message(FATAL_ERROR "Unknown build type: ${BUILD}. Valid options are Debug or Release.")
endif()


# Check MPI version
include(CheckCXXSourceCompiles)

set(MPI_VERSION_TEST_CODE "
#include <mpi.h>
#if MPI_VERSION >= 4
int main() { return 0; }
#else
#error \"MPI version is less than 4.0\"
#endif
")

check_cxx_source_compiles("${MPI_VERSION_TEST_CODE}" MPI_SUPPORTS_PARTITIONED_COMM)

if(MPI_SUPPORTS_PARTITIONED_COMM)
    message(STATUS "MPI version supports partitioned communication (MPI 4.0 or greater).")
    add_definitions(-DHAVE_MPI_PARTITIONED_COMM)
    set(PARTITIONED_COMM_SOURCES
        ${SRC_DIR}/communication/comm_mpi_part.c
    )
else()
    message(STATUS "MPI version does not support partitioned communication (requires MPI 4.0 or greater).")
    set(PARTITIONED_COMM_SOURCES "")
endif()

# Set default compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pthread -lpthread -DHAVE_OMP -fopenmp")

# Check for MPI support
find_package(MPI REQUIRED)
if(MPI_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_MPI")
endif()

# Add conditional flags based on options
if(ENABLE_MPI_PART)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_MPI_PART")
endif()

if(ENABLE_OPEN_SHMEM)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_OSHMEM -loshmem")
endif()

if(ENABLE_CRAY_SHMEM)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CRAY_SHMEM")
endif()

# Define source and header directories
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(OBJ_DIR "${CMAKE_BINARY_DIR}/obj")

# Create necessary build directories
file(MAKE_DIRECTORY "${OBJ_DIR}/kernel")
file(MAKE_DIRECTORY "${OBJ_DIR}/granularity")
file(MAKE_DIRECTORY "${OBJ_DIR}/communication")
file(MAKE_DIRECTORY "${OBJ_DIR}/threading")
file(MAKE_DIRECTORY "${OBJ_DIR}/utility")

# Define source files
set(KERNEL_SOURCES
    ${SRC_DIR}/kernel/kernel_sample.cpp
    ${SRC_DIR}/kernel/kernel_sample_thread.cpp
    ${SRC_DIR}/kernel/kernel_collective_bench.cpp
    ${SRC_DIR}/kernel/kernel_channel_bench.cpp
    ${SRC_DIR}/kernel/kernel_gran_bench.cpp
    ${SRC_DIR}/kernel/kernel_sample_threads.cpp
    ${SRC_DIR}/kernel/kernel_cmb.cpp
    ${SRC_DIR}/kernel/kernel_partsend.cpp
)

set(GRANULARITY_SOURCES
    ${SRC_DIR}/granularity/gran_bulk.cpp
    ${SRC_DIR}/granularity/gran_fine.cpp
    ${SRC_DIR}/granularity/gran_bins.cpp
    ${SRC_DIR}/granularity/gran_bins_nminus.cpp
    ${SRC_DIR}/granularity/gran_offload.cpp
    ${SRC_DIR}/granularity/gran_offload_bins.cpp
    ${SRC_DIR}/granularity/gran_offload_bins_nminus.cpp
)

set(COMMUNICATION_SOURCES
    ${SRC_DIR}/communication/comm_mpi_pt2pt.cpp
    ${SRC_DIR}/communication/comm_mpi_rma.cpp
    ${SRC_DIR}/communication/comm_mpi_pers.cpp
    ${PARTITIONED_COMM_SOURCES}
)

set(THREADING_SOURCES
    ${SRC_DIR}/threading/threads_openmp.cpp
    ${SRC_DIR}/threading/threads_pthreads.cpp
)

set(UTILITY_SOURCES
    ${SRC_DIR}/utility/util_queue.cpp
    ${SRC_DIR}/utility/util_ring.cpp
)

# Add conditional sources based on options
if(ENABLE_OPEN_SHMEM)
    list(APPEND COMMUNICATION_SOURCES ${SRC_DIR}/communication/comm_openshmem.cpp)
endif()

if(ENABLE_CRAY_SHMEM)
    list(APPEND COMMUNICATION_SOURCES ${SRC_DIR}/communication/comm_crayshmem.cpp)
endif()

# Combine all sources
set(SOURCES
    ${SRC_DIR}/driver.cpp
    ${SRC_DIR}/kernel.cpp
    ${SRC_DIR}/gran.cpp
    ${SRC_DIR}/comm.cpp
    ${SRC_DIR}/threads.cpp
    ${KERNEL_SOURCES}
    ${GRANULARITY_SOURCES}
    ${COMMUNICATION_SOURCES}
    ${THREADING_SOURCES}
    ${UTILITY_SOURCES}
)

# Add executable target
add_executable(minimod ${SOURCES})

# Link libraries
target_link_libraries(minimod MPI::MPI_CXX)

# Clean target
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E remove -f ${OBJ_DIR}/*.o minimod
    COMMENT "Cleaning build files"
)
