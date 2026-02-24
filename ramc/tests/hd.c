/* Remote Access Memory Channels (RAMC) 0.5
   Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
   Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
   certain rights in this software. */

/*! \file 
 *  \brief A simple heat diffusion RAMC example.
 * 
 *  Code does a simple heat diffusion with a 5 point stencil. This is strictly 
 *  weak scaling: Each process computes a single temperature based on its 
 *  neighbors, so adding more processes increases the problem size. (Think 
 *  game of life except with heat diffusion equations.) 
 *
 *  Code based on: https://medium.com/data-science/simulation-101-conductive-heat-transfer-a4f09b3e16b4
 *  
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../src/ramc_api.h"

#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3

#define MAX_ATTEMPTS 1000

#define INITIAL_STATUS_VALUE 2

#define MIN_TEMP 20.0
#define MAX_TEMP 1000.0
#define DX 0.1
#define DY 0.1
#define ALPHA 5.0 // thermal diffusion

void opt_error(char *app, char *msg) {
  fprintf(stderr, "%s\n", msg);
  fprintf(stderr, "Usage: %s -i <iterations>\n", app);
  exit(1);
}

int check_square(int n) {
    if (n < 0) {
        return 0; 
    }
    if (n == 0) {
        return 1;
    }
    double sqrt_n = sqrt((double)n);
    if (floor(sqrt_n) == ceil(sqrt_n)) {
        return 1;
    } else {
        return 0;
    }
}

int main(int argc, char **argv) {

  int opt;
  int iterations = 0;
  int num_procs = 0;
  int x, y;
  int myrank;
  int numranks;
  float mytemp = MIN_TEMP;
  float mytemp_old = MIN_TEMP;
  int rank_n, rank_s, rank_e, rank_w;
  int ret;
  int expected_ops = 4; // non-boundary procs expect 4 bb reads and 4 writes
  int attempts;
  // struct to hold information for my own target window
  struct ramc_target_win_info_s my_target_info;
  // structs to hold addressing information for the target buffers of my four neighbors
  struct ramc_init_win_info_s info_n, info_e, info_s, info_w;
  uint64_t tag = 777;
  // bools controlling whether to write in a direction
  int write_n = 0, write_e = 0, write_s = 0, write_w = 0;
  int idx_n, idx_e, idx_s, idx_w; // indices into neighbor temps

  // calculate delta in time
  float dt = (DX*DX)/(2*ALPHA);

  while ((opt = getopt(argc, argv, "i:")) != -1) {
    switch (opt) {
      case 'i':
        iterations = atoi(optarg);
        break;
      case '?':
        opt_error(argv[0], "Unknown argument");
        break;
      default:
        opt_error(argv[0], "Shouldn't have gotten to this point");
        break;
    }
  }

  ramc_init();

  myrank = ramc_get_rank();
  numranks = ramc_get_numranks();
  num_procs = numranks;

  char outfn[256];
  snprintf(outfn, 256, "out.%d.txt", myrank);
  FILE *fp = fopen(outfn, "w");

  if (0 == myrank) {
    if (0 == iterations) opt_error(argv[0], "Missing number of iterations");

    if (!check_square(num_procs)) {
      fprintf(stderr, "Number of processes must be a perfect square (current: %d)\n", num_procs);
      exit(1);
    }
  }

  x = (int)sqrt(num_procs);
  y = x;
  if (x < 7) {
    fprintf(stderr, "Please make dimensions > 6x6\n");
  }

  rank_n = (myrank >= x) ? (myrank-x) : -1;
  rank_e = (myrank % x == (x-1)) ? -1 : (myrank + 1);
  rank_s = (myrank < (num_procs - x)) ? myrank + x : -1;
  rank_w = (myrank % x == 0) ? -1 : (myrank - 1);

  // initialize hot spot
  // one big hot spot filling up everything except the outer two rows/cols
  //if ((myrank >= 2*x) && (myrank <= num_procs-(2*x)) && (myrank % x > 1) && (myrank %x < x-2)) {
  //  mytemp = MAX_TEMP;
  //  mytemp_old = MAX_TEMP;
  //} 

  // spot in the center but for use when there are more procs
//  if ((myrank >= 5*x) && (myrank <= num_procs-(5*x)) && (myrank % x > 4) && (myrank %x < x-4)) {
//    mytemp = MAX_TEMP;
//    mytemp_old = MAX_TEMP;
//  } 

  // ad hoc for 140x140 (19600 procs)
  if ((myrank >= 34*x) && (myrank <= num_procs-(34*x)) && (myrank % x > 33) && (myrank %x < x-33)) {
    mytemp = MAX_TEMP;
    mytemp_old = MAX_TEMP;
  } 

  if (0 == myrank) {
    fprintf(fp, "# Grid size  : %d by %d (%d procs)\n", x, y, x*y);
    fprintf(fp, "# Iterations : %d\n", iterations);
  } 

  // adjust indices into neighbor temps array
  idx_n = NORTH * sizeof(float);
  idx_e = EAST * sizeof(float);
  idx_s = SOUTH * sizeof(float);
  idx_w = WEST * sizeof(float);

  // allocate target window
  float *neighbor_temps = (float *)malloc(4 * sizeof(float));

  // set MIN_TEMP values for nodes on boundaries because these will not be written to
  if (rank_n < 0) neighbor_temps[NORTH] = MIN_TEMP;
  if (rank_e < 0) neighbor_temps[EAST] = MIN_TEMP;
  if (rank_s < 0) neighbor_temps[SOUTH] = MIN_TEMP;
  if (rank_w < 0) neighbor_temps[WEST] = MIN_TEMP;

  fprintf(fp, "# Rank %d: Starting: mytemp: %f  Neighbor Temps: N: %f E:%f S: %f W: %f\n", myrank, mytemp, neighbor_temps[NORTH], neighbor_temps[EAST], neighbor_temps[SOUTH], neighbor_temps[WEST]);
  
  // adjust expected number of bb reads and temp writes to take into account boundaries
  if (myrank >= 0 && myrank < x) expected_ops --;
  if (myrank >= num_procs-x && myrank < num_procs) expected_ops --;
  if (myrank % x == 0) expected_ops--;
  if (myrank % x == (x-1)) expected_ops--;
  assert(expected_ops > 0);

  fprintf(fp, "# Rank %d: expected ops: %d\n", myrank, expected_ops);

  // create target window
  ret = ramc_tgt_create_window(neighbor_temps, 4 * sizeof(float), tag, INITIAL_STATUS_VALUE, &my_target_info);
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error creating target window\n", myrank);
    exit(1);
  }

  // post target window info to BB
  ret = ramc_tgt_post_window(&my_target_info);

  // set BB status to active
  ret = ramc_tgt_activate_bb();

  if (rank_n >= 0) {
    // check if N is active (could block)
    attempts = 0;
    do {
      ret = ramc_init_check_bb_status(rank_n, tag);
      } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
    if (attempts > MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for north rank (%d)\n", myrank, rank_n);
      exit(1);
    }
    // get N info
    ret = ramc_init_get_bb_posting(rank_n, INITIAL_STATUS_VALUE, &info_n);  
  }
  
  if (rank_e >= 0) {
    // check if E is active
    attempts = 0;
    do {
      ret = ramc_init_check_bb_status(rank_e, tag);
      } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
    if (attempts > MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for east rank (%d)\n", myrank, rank_e);
      exit(1);
    }
    // get E info
    ret = ramc_init_get_bb_posting(rank_e, INITIAL_STATUS_VALUE, &info_e);
  }

  if (rank_s >= 0) {
    // check if S is active
    attempts = 0;
    do {
      ret = ramc_init_check_bb_status(rank_s, tag);
      } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
    if (attempts > MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for south rank (%d)\n", myrank, rank_s);
      exit(1);
    }
    // get S info
    ret = ramc_init_get_bb_posting(rank_s, INITIAL_STATUS_VALUE, &info_s);
  }
  
  if (rank_w >= 0) {
    // check if W is active
    attempts = 0;
    do {
      ret = ramc_init_check_bb_status(rank_w, tag);
      } while (ret != RAMC_SUCCESS && attempts < MAX_ATTEMPTS);
    if (attempts > MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded max attempts getting BB status for west rank (%d)\n", myrank, rank_w);
      exit(1);
    }
    // get W info
    ret = ramc_init_get_bb_posting(rank_w, INITIAL_STATUS_VALUE, &info_w);
  }

  // block on expected number of BB reads
  ret = ramc_tgt_await_bb_reads(expected_ops);

  // set BB to inactive
  ret = ramc_tgt_deactivate_bb(); 

  // everyone should have neighbor target buffer info now

  fprintf(fp, "Epoch: %lu Rank: %d Temp: %f\n", ramc_tgt_get_win_status_raw(&my_target_info), myrank, mytemp);

  ret = ramc_barrier_binary();
  
  for (int i = 0; i < iterations; ++i) {
    // copy current temp
    mytemp_old = mytemp;

    // make my target write enabled
    ret = ramc_tgt_increment_win_status(&my_target_info, 1);

    // make my neighbor info have the write value
    if (rank_n >= 0) {
      ret = ramc_init_increment_status(&info_n, 1);
      write_n = 1;
    }
    if (rank_e >= 0) {
      ret = ramc_init_increment_status(&info_e, 1);
      write_e = 1;
    }
    if (rank_s >= 0) {
      ret = ramc_init_increment_status(&info_s, 1);
      write_s = 1;
    }
    if (rank_w >= 0) {
      ret = ramc_init_increment_status(&info_w, 1);
      write_w = 1;
    }

    // TODO: replace with non-blocking puts
    attempts = 0;
    do {
      if (write_n) {
        if (ramc_init_check_win_status(&info_n) == RAMC_SUCCESS) {
          // do the put
          // Note this is the proc to the south of the target so the value goes in the SOUTH location at the target
          ret = ramc_put(&mytemp, sizeof(float), &info_n, idx_s);
          assert(ret == RAMC_SUCCESS);
          // mark that the write is done
          write_n = 0;
        }
      }
      if (write_e) {
        if (ramc_init_check_win_status(&info_e) == RAMC_SUCCESS) {
          // do the put 
          ret = ramc_put(&mytemp, sizeof(float), &info_e, idx_w);
          assert(ret == RAMC_SUCCESS);
          // mark that the write is done
          write_e = 0;
        }
      }
      if (write_s) {
        if (ramc_init_check_win_status(&info_s) == RAMC_SUCCESS) {
          // do the put 
          ret = ramc_put(&mytemp, sizeof(float), &info_s, idx_n);
          assert(ret == RAMC_SUCCESS);
          // mark that the write is done
          write_s = 0;
        }
      }
      if (write_w) {
        if (ramc_init_check_win_status(&info_w) == RAMC_SUCCESS) {
          // do the put 
          ret = ramc_put(&mytemp, sizeof(float), &info_w, idx_e);
          assert(ret == RAMC_SUCCESS);
          // mark that the write is done
          write_w = 0;
        }
      }
      attempts++;
    } while ((write_n || write_e || write_s || write_w) && attempts < MAX_ATTEMPTS);

    if (attempts >= MAX_ATTEMPTS) {
      fprintf(stderr, "Rank %d: Exceeded maximum attempts at updating neighbors (N: %d E: %d S: %d W: %d)\n", myrank, write_n, write_e, write_s, write_w);
      fflush(stderr);
      exit(1);
    }

    // block on receiving expected number of writes
    ret = ramc_tgt_await_win_ops(&my_target_info, expected_ops);
    assert(ret == RAMC_SUCCESS);

    // switch to next target epoch to prevent writes
    ret = ramc_tgt_increment_win_status(&my_target_info, 1);
    
    // make my neighbor info have the right epoch (no longer writable)
    ret = ramc_init_increment_status(&info_n, 1);
    ret = ramc_init_increment_status(&info_e, 1);
    ret = ramc_init_increment_status(&info_s, 1);
    ret = ramc_init_increment_status(&info_w, 1);

    //fprintf(fp, "Epoch: %lu temps for calculation: self: %f N: %f E: %f S: %f W: %f\n", my_target_status, mytemp, neighbor_temps[NORTH], neighbor_temps[EAST], neighbor_temps[SOUTH], neighbor_temps[WEST]);

    // update my temp 
    //mytemp = mytemp + TIME_STEP*THERMAL_DIFFUSION*((neighbor_temps[EAST] - (2*mytemp) + neighbor_temps[WEST]) + (neighbor_temps[NORTH] - (2*mytemp) + neighbor_temps[SOUTH]));
    mytemp += (dt * ((neighbor_temps[EAST] - 2*mytemp_old  + neighbor_temps[WEST])/(DX*DX) + (neighbor_temps[NORTH]- 2*mytemp_old + neighbor_temps[SOUTH])/(DY*DY)) + mytemp_old) - mytemp; 

    fprintf(fp, "Epoch: %lu Rank: %d Temp: %f\n", ramc_tgt_get_win_status_raw(&my_target_info), myrank, mytemp);

  }

  fprintf(fp, "# COMPLETE\n");
  fclose(fp);
  ret = ramc_barrier_binary();
  if (ret != RAMC_SUCCESS) {
    fprintf(stderr, "Rank %d: Error encountered in final barrier\n", myrank);
    exit(1);
  }

  ret = ramc_tgt_destroy_window(&my_target_info);
  ret = ramc_finalize();
  
  return 0;
}
