/* Remote Access Memory Channels (RAMC) 0.5
 * Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef __RAMC_CORE_H__
#define __RAMC_CORE_H__

#include "pmi_stuff.h"

#include <stdint.h>
#include <stdatomic.h>
#include <rdma/fi_rma.h>
#include <stdio.h>
#include <stdlib.h>

// provider and nic to use
// TODO: The following are used to select the appropriate NIC exclusively on MI300A
// TODO: Generalize (non-trivial!)
#define LIBFABRIC_PROVIDER "cxi"
#define LIBFABRIC_DOMAIN_0 "cxi0"
#define LIBFABRIC_DOMAIN_1 "cxi1"
#define LIBFABRIC_DOMAIN_2 "cxi2"
#define LIBFABRIC_DOMAIN_3 "cxi3"
#define NUM_NUMA 4 // MI300A
#define CORES_PER_NUMA 24 // MI300A

// Initial call to fi_getname returns the address and the length
// This value is used to allocate an array to store the adddress and check if 
// there is sufficient space.
#define ADDR_LENGTH (256)

// timeout (ms) for waiting on CQ and counter completions
#define COMPLETION_TIMEOUT (10000)

// The Bulletin Board
// The following two structs are posted to the BB
// The first is REMOTE_READ only, and is used to indicate 
// whether the BB is active and the tag 
// for that information. The second contains the addressing 
// information for the target window.
// TODO: we have the status as an atomic but does it need to be? 
struct ramcc_bb_status_s {
  uint64_t tag;
  atomic_uchar status;
};

// Stores addressing information for a target's (i) status value and (ii) buffer
// struct fi_rma_iov is from libfabric, and stores an address and memory key
struct ramcc_target_addr_info_s{
  struct fi_rma_iov status; 
  struct fi_rma_iov target;
};

// General RAMC return codes
#define RAMC_SUCCESS 0
#define RAMC_FAILURE 1
#define RAMC_ERROR 2
#define RAMC_TIMEOUT 3

// BB return codes
#define RAMC_BB_STATUS_INACTIVE 0
#define RAMC_BB_STATUS_ACTIVE 1

// Status value return codes
#define RAMC_TARGET_READ 0
#define RAMC_TARGET_WRITE 1
#define RAMC_TARGET_BEHIND 2
#define RAMC_TARGET_AHEAD 3

// reserved target status values
#define RAMC_TARGET_STATUS_DESTROYED 0
#define RAMC_TARGET_STATUS_ABORT 1
#define RAMC_TARGET_STATUS_INITIALIZE 2

// transport information
struct ramcc_tp_s {
    struct fi_info *provider;       // fi_info struct with selected provider info
    struct fid_fabric *fabric;
    struct fid_domain *domain;
    struct fid_av *av;              // address vector
    struct fid_ep *ep;              // endpoint
    struct fid_eq *eq;              // event queue
    struct fid_cq *rxcq;            // recv completion queue
    struct fid_cq *txcq;            // transmit completion queue
    void *cqe;                      // a completion queue entry to be used for waiting
    size_t inject_size;             // maximum size for fi_inject (from provider)

    // counters for endpoints
    struct fid_cntr *write_cntr_ep; // FI_WRITE
    uint64_t write_cntr_ep_val;
    uint64_t write_cntr_ep_expected;
    struct fid_cntr *read_cntr_ep; // FI_READ
    uint64_t read_cntr_ep_val;
    uint64_t read_cntr_ep_expected;

    // counters for memory regions
    // TODO: do we use this?
    //struct fid_cntr *remote_write_cntr_mr;
    //uint64_t remote_write_cntr_mr_val;

    // PMI information 
    int pmi_myrank;
    int pmi_numranks;
    int pmi_max_val_len;
    int pmi_max_key_len;

    fi_addr_t *addresses; // addresses of all PEs, indexed by PMI rank
    size_t addr_len;      // length of address used by selected provider
    char my_addr[ADDR_LENGTH]; // address of this PE

    // for FI_MSG based barrier
    uint64_t barrier_send[1];
    uint64_t barrier_recv[1];
    struct fid_mr *barrier_send_mr;
    struct fid_mr *barrier_recv_mr;

    // The bulletin board
    struct fi_rma_iov *bb_statuses; // will store addresses of BB status
    struct fi_rma_iov *bb_postings; // will store addresses of BB postings
    struct ramcc_bb_status_s *bb_status;
    struct fid_mr *bb_status_mr;
    struct ramcc_target_addr_info_s *bb_posting; 
    struct fid_mr *bb_posting_mr;
    struct fid_cntr *bb_remote_read_cntr; // for counting number of reads of BB addr info
    uint64_t bb_remote_read_cntr_val; // stores last read counter value
};


// When a target creates a window, the information for that window is stored in 
// ramc_target_win_info_s struct; this struct is then used with the BB to post 
// addressing information; it also contains the counter attached to the target's window
//
// When a channel is created, the initiator stores target addressing information; 
// this struct is the used for subsequent communication operations
//
// The target information held by the target (creator) of a window
// Note that the current window status is part of the target buffer structure, so is accessed through the addr_info field
struct ramc_target_win_info_s {
  int rank;                     // rank of target owning this buffer
  uint64_t tag;                 // tag of target buffer
  void *buffer;                 // buffer
  struct fid_mr * mr_buf;           // pointer to OFI-internal memory region struct
  struct fid_mr * mr_status;           // pointer to OFI-internal memory region struct
  uint64_t status_val;          // status value of target buffer
  struct ramcc_target_addr_info_s addr_info; // libfabric addressing info for status_val and buffer
  struct fid_cntr *window_cntr; // for counting FI_REMOTE_* ops on this buffer
  uint64_t window_cntr_val;     // last read counter value
};

// The target information held by the initiator using a target's window
// tag is not included because matching already occurred and there no good way to get it
struct ramc_init_win_info_s {
  int rank; // rank of target
  struct ramcc_target_addr_info_s addr_info; // addressing info for target status value and buffer
  uint64_t status_val;    // expected status value of target 
};

//////////////////////////////////////////////////////////////
// Error reporting stuff
//////////////////////////////////////////////////////////////

static void handle_fi_error(int errcode, char *file, int line, char *fn)
{
    fprintf(stderr, "%s:%d: %s: %s\n", file, line, fn, fi_strerror(-errcode));
    abort();
}

#define RAMC_ERR(fmt, ...) do { \
    fprintf(stderr, "***RAMC ERROR*** %s:%d:%s(): Rank %d: " fmt, \
      __FILE__, __LINE__, __func__, ramcc_tp_get_myrank(), ##__VA_ARGS__); \
} while (0)

#ifdef __RAMC_DEBUG__
#define RAMC_DEBUG(fmt, ...) do { \
    fprintf(stderr, "[RAMC DEBUG] %s:%d:%s(): Rank %d: " fmt, \
      __FILE__, __LINE__, __func__, ramcc_tp_get_myrank(), ##__VA_ARGS__); \
} while (0)
#else
#define RAMC_DEBUG(fmt, args...)
#endif

#define FI_CHECK(fn) { \
    int errcode;\
    errcode = (fn);\
    if (errcode < 0) handle_fi_error(errcode, __FILE__, __LINE__, #fn); }

// Transport initialization, operation, and teardown
int ramcc_tp_init(void);
int ramcc_tp_finalize(void);
int ramcc_tp_get_myrank(void);
int ramcc_tp_get_numranks(void);
void ramcc_tp_init_ep_counters();
void ramcc_tp_close_ep_counters();
int ramcc_tp_init_barrier();
int ramcc_tp_close_barrier();
int ramcc_tp_setup_av();
int ramcc_tp_await_completion(struct fid_cq *, void *, uint64_t);
int ramcc_tp_await_completion_multi(struct fid_cq *cq, struct fi_cq_entry *cqe, const int num_entries);
int ramcc_tp_await_completion_cntr(struct fid_cntr *, uint64_t, uint64_t *, uint64_t);
int ramcc_tp_check_for_error(struct fid_cq *);

// Transport utility
int ramcc_tp_barrier_linear();
int ramcc_tp_barrier_binary();
void ramcc_tp_print_ep_cntrs();
int ramcc_tp_get_hwt();

int ramcc_target_create_window(void *target_buf, size_t target_len, uint64_t tag, uint64_t init_status_val, struct ramc_target_win_info_s *target_info);
int ramcc_target_destroy_window(struct ramc_target_win_info_s *target_info);
int ramcc_target_post_window(struct ramc_target_win_info_s *);
int ramcc_target_activate_bb(void);
int ramcc_target_deactivate_bb(void);
int ramcc_target_await_bb_reads(uint64_t expected_num_reads); // BLOCKING
int ramcc_target_test_bb_reads(uint64_t expected_num_reads); // NON-BLOCKING
int ramcc_target_await_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops); // BLOCKING
int ramcc_target_test_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops); // NON-BLOCKING

int ramcc_init_get_bb_status(int target_rank, uint8_t *status, uint64_t *tag); // BLOCKING (on read cntr)
// TODO: Assumes getting a bb posting will be for a newly initialized target window!!! 
int ramcc_init_get_bb_posting(int target_rank, uint64_t init_status_val, struct ramc_init_win_info_s *); // BLOCKING (on endpoint read counter)
int ramcc_init_get_win_status(struct ramc_init_win_info_s *target_info, uint64_t *tgt_value); // BLOCKING
int ramcc_init_put(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset); // BLOCKING
// TODO: Currently await_all_* blocks on ALL ops of that type completing, regardless of what the target is
// I.e., it does not track operations targeting a specific buffer
int ramcc_init_put_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset); // NON-BLOCKING
int ramcc_init_get(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset);
int ramcc_init_get_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset);
int ramcc_init_await_all_puts(); // BLOCKING
int ramcc_init_await_all_gets(); // BLOCKING

// ---- Experimental ----
int ramcc_init_atomic_inc(struct ramc_init_win_info_s *target_info, size_t target_offset);


#endif // __RAMC_CORE_H__
