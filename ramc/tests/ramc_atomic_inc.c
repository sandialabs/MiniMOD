/* Remote Access Memory Channels (RAMC) 0.5
   Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
   Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
   certain rights in this software. */

/*! \file 
 *  \brief Exercises atomic increment.
 *
 *  This code was written to simply exercise the experimental atomic increment 
 *  over a RAMC channel. 
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "../src/ramc_api.h"

#define MAX_ATTEMPTS 1000
#define INITIAL_STATUS_VALUE 2

#define TGT_BUF_SIZE 4

void opt_error(char *msg) {
  fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Usage: put -i <iterations>\n");
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

  // target buffer
  uint64_t target_buf[TGT_BUF_SIZE];
  for (int i = 0; i < TGT_BUF_SIZE; ++i) target_buf[i] = 0;

  // struct to hold information for my own target window
  struct ramc_target_win_info_s my_target_info;
  // struct to hold information for the target buffer of the other proc 
  struct ramc_init_win_info_s peer_info;
  uint64_t tag = 777;

  while ((opt = getopt(argc, argv, "i:")) != -1) {
    switch (opt) {
      case 'i':
        iterations = atoi(optarg);
        break;
      case '?':
        opt_error("Unknown argument");
        break;
      default:
        opt_error("Shouldn't have gotten to this point");
        break;
    }
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

  if (1 == myrank) {
    ret = ramc_tgt_create_window(&target_buf, TGT_BUF_SIZE*sizeof(uint64_t), tag, INITIAL_STATUS_VALUE, &my_target_info);
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
      } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
    if (attempts > MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for peer\n", myrank);
      exit(1);
    }
    // get target buffer info
    ret = ramc_init_get_bb_posting(peer_rank, INITIAL_STATUS_VALUE, &peer_info);  
  }

  // initiator should have target's buffer info now

  if (1 == myrank) {
    printf("# Atomic increment test using target buffer of %d uint64_t\n", TGT_BUF_SIZE);
    printf("Initial buffer contents on target: ");
    for (int j = 0; j < TGT_BUF_SIZE; ++j) {
      printf("%lu ", target_buf[j]);
    }
    printf("\n");
    fflush(stdout);
  }

  ret = ramc_barrier_binary();

  for (int i = 0; i < iterations; ++i) {
    if (0 == myrank) {
      // update my epoch to indicate I want to write
      ret = ramc_init_increment_status(&peer_info, 1);

      // when target buffer is write enabled, do the atomic increment 
      do {
        ret = ramc_init_check_win_status(&peer_info);
      } while (ret != RAMC_SUCCESS);

      for (int j = 0; j < TGT_BUF_SIZE; ++j) {
        ret = ramc_atomic_inc(&peer_info, j);
      }
      // update my epoch to indicate I won't write 
      ret = ramc_init_increment_status(&peer_info, 1);
    } else {
      // write enable my target buffer
      ret = ramc_tgt_increment_win_status(&my_target_info, 1);
      // blocking wait on getting TGT_BUF_SIZE atomic increments
      ret = ramc_tgt_await_win_ops(&my_target_info, TGT_BUF_SIZE); // checks MR remote_write counter
      // make my buffer read only
      ret = ramc_tgt_increment_win_status(&my_target_info, 1);
      printf("Rank 1: My target buffer was updated to: ");
      for (int j = 0; j < TGT_BUF_SIZE; ++j) {
        printf("%lu ", target_buf[j]);
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
