/* Remote Access Memory Channels (RAMC) 0.5
   Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
   Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
   certain rights in this software. */

/*! \file 
 *  \brief Ping-pong latency micro-benchmark.
 *
 *  Modeled after the MPI P2P latency test in the OSU message becnhmarks (OMB).
 *  https://mvapich.cse.ohio-state.edu/benchmarks/
 *
 *  Command line options:
 * 
 *  -i : Number of iterations \n 
 *  -m : maximum message size bytes exponent (i.e., 2^m) \n 
 *
 *  This version uses a traditional explicit notification method to let 
 *  targets know data has been written (so they can return the message). After 
 *  putting the data, an initiator issues a second put writing the next iteration 
 *  number; the target waits for this value to increment before doing its put.
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
  fprintf(stderr, "Usage: ramc_pp_latency -i <iterations> -m <max_message_size_exponent>\n");
  exit(1);
}

int main(int argc, char **argv) {

  int opt;
  int iterations = 0;
  int max_exp = 0;
  int num_procs = 0;
  int myrank;
  int numranks;
  int ret;
  int attempts;
  uint64_t cur_msg_size = 0;
  uint64_t max_msg_size = 0;

  int warmups = 50;

  // each process has the buffer beign ping ponged and a second 
  // buffer used for atomic increments to notify the target of 
  // completion
  // The first is the target buffer (tgt_buf)
  // The second is the target notification (tgt_not)
  // struct to hold information for my own target window
  struct ramc_target_win_info_s tgt_buf_info;
  struct ramc_target_win_info_s tgt_not_info;
  
  // struct to hold information for the target buffer of the other proc 
  struct ramc_init_win_info_s peer_buf_info;
  struct ramc_init_win_info_s peer_not_info;
  
  // tag used during channel creation
  uint64_t tag = 777;

  while ((opt = getopt(argc, argv, "i:m:")) != -1) {
    switch (opt) {
      case 'i':
        iterations = atoi(optarg);
        break;
      case 'm':
        max_exp = atoi(optarg);
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

  max_msg_size = (uint64_t)pow(2, max_exp);

  // allocate target buffer and notificaton buffer 
  uint8_t *tgt_buf = (uint8_t *)malloc(max_msg_size * sizeof(uint8_t));
  uint64_t tgt_not = 0; 

  // create windows for the buffer and notification
  ret = ramc_tgt_create_window(tgt_buf, max_msg_size * sizeof(uint8_t), tag, INITIAL_STATUS_VALUE, &tgt_buf_info);
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error creating window for tgt_buf\n", myrank);
    exit(1);
  }
  // Note we use a different tag to avoid a peer getting the wrong addressing info
  ret = ramc_tgt_create_window(&tgt_not, sizeof(uint64_t), tag+1, INITIAL_STATUS_VALUE, &tgt_not_info);
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error creating window for tgt_not\n", myrank);
    exit(1);
  }

  // post tgt_buf info to BB
  ret = ramc_tgt_post_window(&tgt_buf_info);

  // set BB status to active
  ret = ramc_tgt_activate_bb();

  // check if peer's target buffer is active (could block)
  attempts = 0;
  do {
    ret = ramc_init_check_bb_status(peer_rank, tag);
    } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
  if (attempts > MAX_ATTEMPTS) {
    fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for peer\n", myrank);
    exit(1);
  }
  // get peer target buffer info
  ret = ramc_init_get_bb_posting(peer_rank, INITIAL_STATUS_VALUE, &peer_buf_info);  
  
  // block on 1 BB reads
  ret = ramc_tgt_await_bb_reads(1);

  // set BB to inactive
  ret = ramc_tgt_deactivate_bb(); 
  
  // post tgt_not info to BB
  ret = ramc_tgt_post_window(&tgt_not_info);

  // set BB status to active
  ret = ramc_tgt_activate_bb();

  // check if peer's target buffer is active (could block)
  attempts = 0;
  do {
    ret = ramc_init_check_bb_status(peer_rank, tag+1);
    } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
  if (attempts > MAX_ATTEMPTS) {
    fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for peer\n", myrank);
    exit(1);
  }
  // get peer target buffer info
  ret = ramc_init_get_bb_posting(peer_rank, INITIAL_STATUS_VALUE, &peer_not_info);  
  
  // block on 1 BB reads
  ret = ramc_tgt_await_bb_reads(1);

  // set BB to inactive
  ret = ramc_tgt_deactivate_bb(); 

  //
  // everyone should have the addressing info for their peer's target buffer and notification buffer
  //

  if (0 == myrank) {
    printf("# Pingpong latency test, using explicit notification (no MR counters)\n");
    printf("# Iterations per message size : %d\n", iterations);
    printf("# Warmups per message size    : %d\n", warmups);
    printf("msg_size_bytes,total_pp_latency_usecs,one_way_pp_latency\n");
    fflush(stdout);
  }

  ret = ramc_barrier_binary();

  struct timeval tv;
  uint64_t start;
  double result = 0.0;
  uint64_t expected_val = 0;
  for (cur_msg_size = 1; cur_msg_size <= max_msg_size; cur_msg_size=cur_msg_size*2) {
    // Do warmups
    for (int i = 0; i < warmups; ++i) {
      expected_val++;
      if (0 == myrank) {
        ret = ramc_put(tgt_buf, cur_msg_size, &peer_buf_info, 0);
        //ret = ramc_atomic_inc(&peer_not_info, 0); // notification with atomic
        ret = ramc_put(&expected_val, sizeof(uint64_t), &peer_not_info, 0); // notification without atomic
        do {;} while (tgt_not < expected_val);
      } else {
        do {;} while (tgt_not < expected_val);
        ret = ramc_put(tgt_buf, cur_msg_size, &peer_buf_info, 0);
        //ret = ramc_atomic_inc(&peer_not_info, 0); // notification with atomic
        ret = ramc_put(&expected_val, sizeof(uint64_t), &peer_not_info, 0); // notification without atomic
      }
    }
    // Get times
    uint64_t accumulate = 0;
    for (int i = 0; i < iterations; ++i) {
      expected_val++;
      if (0 == myrank) {
        gettimeofday(&tv, NULL);
        start = (uint64_t)((tv.tv_sec * 1000000) + (tv.tv_usec));
        ret = ramc_put(tgt_buf, cur_msg_size, &peer_buf_info, 0);
        ret = ramc_atomic_inc(&peer_not_info, 0);
        do {;} while (tgt_not < expected_val);
        gettimeofday(&tv, NULL);
        accumulate += (uint64_t)((tv.tv_sec * 1000000) + (tv.tv_usec)) - start;
      } else {
        do {;} while (tgt_not < expected_val);
        ret = ramc_put(tgt_buf, cur_msg_size, &peer_buf_info, 0);
        ret = ramc_atomic_inc(&peer_not_info, 0);
      }
    }
    if (0 == myrank) {
      result = accumulate/(double)iterations;
      printf("%lu,%f,%f\n", cur_msg_size, result, result/2.0);
    }
  }

  ret = ramc_barrier_binary();
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error encountered in final barrier\n", myrank);
    exit(1);
  }

  ret = ramc_tgt_destroy_window(&tgt_buf_info);
  ret = ramc_tgt_destroy_window(&tgt_not_info);
  ret = ramc_finalize();
  
  return 0;
}
