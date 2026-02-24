/* Remote Access Memory Channels (RAMC) 0.5
   Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
   Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
   certain rights in this software. */

/*! \file 
 *  \brief Bandwidth micro-benchmark.
 *
 *  Modeled after the MPI P2P BW test in the OSU message becnhmarks (OMB).
 *
 *  Command line options:
 * 
 *  -i : Number of iterations \n 
 *  -m : maximum message size bytes exponent (i.e., 2^m) \n 
 * 
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include "../src/ramc_api.h"

#define MAX_ATTEMPTS 1000
#define MSGS_PER_ITERATION 2000
#define MEGA 1000000.0

#define INITIAL_STATUS_VALUE 2

void opt_error(char *msg) {
  fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Usage: ramc_pp_bw -i <iterations> -m <max_message_size_exponent>\n");
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
  uint64_t ack = 1;
  double bytes = 0;
  double mbytes = 0;
  double seconds = 0.0;
  double mbytes_per_second = 0;

  // struct to hold information for my own target window
  struct ramc_target_win_info_s my_target_info;
  // struct to hold information for the target buffer of the other proc 
  struct ramc_init_win_info_s peer_info;
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

  // allocate target window
  // initiator uses the same window for ACKs (yes it is overkill but whatever)
  uint8_t *buffer = (uint8_t *)malloc(max_msg_size * sizeof(uint8_t));

  // create target window
  ret = ramc_tgt_create_window(buffer, max_msg_size * sizeof(uint8_t), tag, INITIAL_STATUS_VALUE, &my_target_info);
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error creating target window\n", myrank);
    exit(1);
  }

  // post target window info to BB
  ret = ramc_tgt_post_window(&my_target_info);

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
  ret = ramc_init_get_bb_posting(peer_rank, INITIAL_STATUS_VALUE, &peer_info);  
  
  // block on 1 BB reads
  ret = ramc_tgt_await_bb_reads(1);

  // set BB to inactive
  ret = ramc_tgt_deactivate_bb(); 

  // everyone should have their peer's target buffer info now

  if (0 == myrank) {
    printf("# BW test\n");
    printf("# Trials per message size: %d\n", iterations);
    printf("# NO WARMUP\n");
    printf("msg_size_bytes,bw_mb_s\n");
    fflush(stdout);
  }

  ret = ramc_barrier_binary();

  struct timeval tv_start, tv_end;
  uint64_t start;
  double result = 0.0;
  for (cur_msg_size = 1; cur_msg_size <= max_msg_size; cur_msg_size=cur_msg_size*2) {
    uint64_t accumulate = 0;
    for (int i = 0; i < iterations; ++i) {
      ret = ramc_barrier_binary();

      // we use a non-blocking put (i.e., it does not wait on the local counter to 
      // test for local put completion); the ACK from the target ensures the local puts 
      // have completed.
      if (0 == myrank) {
        gettimeofday(&tv_start, NULL);
        for (int j = 0; j < MSGS_PER_ITERATION; ++j) {
//          ret = ramc_put(buffer, cur_msg_size, &peer_info, 0);
          ret = ramc_put_nb(buffer, cur_msg_size, &peer_info, 0);
        }
        // wait for explicit ACK from target
        ret = ramc_tgt_await_win_ops(&my_target_info, 1);
        gettimeofday(&tv_end, NULL);
        start = (uint64_t)((tv_start.tv_sec * 1000000) + (tv_start.tv_usec));
        accumulate += (uint64_t)((tv_end.tv_sec * 1000000) + (tv_end.tv_usec)) - start;
      } else {
        ret = ramc_tgt_await_win_ops(&my_target_info, MSGS_PER_ITERATION);
        ret = ramc_put(&ack, sizeof(uint64_t), &peer_info, 0);
        ++ack;
      }
    }
    if (0 == myrank) {
      bytes = MSGS_PER_ITERATION * cur_msg_size * iterations; 
      mbytes = (float)bytes/MEGA;
      seconds = (float)accumulate/MEGA;
      mbytes_per_second = mbytes/seconds;
      printf("%lu,%f,%f\n", cur_msg_size, seconds, mbytes_per_second);
      fflush(stdout);
    } 
  }

  ret = ramc_barrier_binary();
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error encountered in final barrier\n", myrank);
    exit(1);
  }

  ret = ramc_finalize();

  fprintf(stderr, "Rank %d: DONE!\n", myrank);
  fflush(stderr);
  
  return 0;
}
