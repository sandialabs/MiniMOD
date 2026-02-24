/* Remote Access Memory Channels (RAMC) 0.5
   Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
   Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
   certain rights in this software. */
/*! \file 

 *  \brief Exercises the lazy synchronization between initiator and target 
 *
 *  The initiator cannot write until the target's status is ready-to-write. This 
 *  test makes sure that mechanism works by 
 *  having the target sleep by +1 seconds per iteration to add delay before the 
 *  target enters a ready-to-write status. Outputs the delay time measured by 
 *  initiator (i.e., how long it had to wait to do the put).
 *
 *  Command line options:
 * 
 *  -i : Number of iterations \n 
 *  -b : Message size \n 
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include "../src/ramc_api.h"

#define MAX_ATTEMPTS 5000

#define INITIAL_STATUS_VALUE 2

void opt_error(char *msg) {
  fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Usage: put -i <iterations> -b <buffer_size_bytes>\n");
  exit(1);
}

int main(int argc, char **argv) {

  int opt;
  int iterations = 0;
  size_t size = 0;
  int num_procs = 0;
  int myrank;
  int numranks;
  int ret;
  int attempts;

  // struct to hold information for my own target window
  struct ramc_target_win_info_s my_target_info;
  // struct to hold information for the target buffer of the other proc 
  struct ramc_init_win_info_s peer_info;
  uint64_t tag = 777;

  while ((opt = getopt(argc, argv, "i:b:")) != -1) {
    switch (opt) {
      case 'i':
        iterations = atoi(optarg);
        break;
      case 'b':
        size = atoi(optarg);
        break;
      case '?':
        opt_error("Unknown argument");
        break;
      default:
        opt_error("Shouldn't have gotten to this point");
        break;
    }
  }

  if (iterations > 10) {
    fprintf(stderr, "Number of iterations should be 10 or less\n");
    exit(1);
  }

  ramc_init();

  myrank = ramc_get_rank();
  numranks = ramc_get_numranks();

  if (0 == myrank && 2 != numranks) {
    printf("Test requires 2 and only 2 processes (current = %d)\n", numranks);
    exit(1);
  }

  if (0 == myrank) {
    if (0 == iterations) opt_error("Missing number of iterations");
  }

  int peer_rank;
  if (0 == myrank) {
    peer_rank = 1;
  } else {
    peer_rank = 0;
  }

  // allocate buffer
  // will be used as source on initiator and target on target
  uint8_t *buffer = (uint8_t *)malloc(size * sizeof(uint8_t));

  // initialize buffer
  for (int j = 0; j < size; ++j) {
    if (0 == myrank) {
      buffer[j] = 1;
    } else {
      buffer[j] = 0;
    }
  }

  // create target window on target
  if (1 == myrank) {
    ret = ramc_tgt_create_window(buffer, size * sizeof(uint8_t), tag, INITIAL_STATUS_VALUE, &my_target_info);
    if (ret != RAMC_SUCCESS) {
      fprintf(stderr, "Rank %d: Error creating target window\n", myrank);
      exit(1);
    }
    // post target window info to BB
    ret = ramc_tgt_post_window(&my_target_info);
    // set BB status to active
    ret = ramc_tgt_activate_bb();
    // block on 1 BB reads
    ret = ramc_tgt_await_bb_reads(1);
    // set BB to inactive
    ret = ramc_tgt_deactivate_bb(); 
  } else {
    // check if target's target buffer is active (could block)
    attempts = 0;
    do {
      ret = ramc_init_check_bb_status(peer_rank, tag);
      attempts++;
      } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
    if (attempts >= MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for peer\n", myrank);
      exit(1);
    }
    // get target buffer info
    ret = ramc_init_get_bb_posting(peer_rank, INITIAL_STATUS_VALUE, &peer_info);
  }

  // initiator should have target's buffer info now

  if (1 == myrank) {
    printf("# Put test using %zu uint8_t\n", size);
    printf("Initial buffer contents on target: ");
    for (int j = 0; j < size; ++j) {
      printf("%hhu ", buffer[j]);
    }
    printf("\n");
    fflush(stdout);
  }

  ret = ramc_barrier_binary();

  struct timeval tv;
  uint64_t start, end;
  for (int i = 0; i < iterations; ++i) {
    if (0 == myrank) {
      // update my value to indicate I want to write
      ret = ramc_init_increment_status(&peer_info, 1);

      gettimeofday(&tv, NULL);
      start = (uint64_t)((tv.tv_sec * 1000000) + (tv.tv_usec));

      // when target buffer is write enabled, do the put
      do {
        ret = ramc_init_check_win_status(&peer_info);
      } while (ret != RAMC_SUCCESS);

      gettimeofday(&tv, NULL);
      end = (uint64_t)((tv.tv_sec * 1000000) + (tv.tv_usec));
      printf("Rank 0: Waited %lu secs until doing put\n", (end-start)/1000000);

      ret = ramc_put(buffer, size, &peer_info, 0);
      // update my value to indicate I won't write 
      ret = ramc_init_increment_status(&peer_info, 1);
      for (int j = 0; j < size; ++j) {
        buffer[j]++;
      }
      // Don't need barrier because using status/epoch check
      //ret = ramc_barrier_binary();
    } else {
      // delay
      sleep(i+1);
      // write enable my target buffer
      ret = ramc_tgt_increment_win_status(&my_target_info, 1);
      // blocking wait on getting a put
      ret = ramc_tgt_await_win_ops(&my_target_info, 1); // checks MR remote_write counter
      // make my buffer read only
      ret = ramc_tgt_increment_win_status(&my_target_info, 1);
      printf("Rank 1: My target buffer was updated to: ");
      for (int j = 0; j < size; ++j) {
        printf("%hhu ", buffer[j]);
      }
      printf("\n");
      fflush(stdout);
      // using the status/epoch check removes the need for a barrier
      //ret = ramc_barrier_binary();
    }
  }
  
  ret = ramc_barrier_binary();
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error encountered in final barrier\n", myrank);
    exit(1);
  }

  if (1 == myrank) {
    ret = ramc_tgt_destroy_window(&my_target_info);
  }
  ret = ramc_finalize();
  
  return 0;
}
