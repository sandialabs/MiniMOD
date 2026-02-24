/* Remote Access Memory Channels (RAMC) 0.5
 * Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "ramc_core.h"

#include <stdbool.h>
#include <rdma/fabric.h>
#include <rdma/fi_cm.h>
#include <rdma/fi_domain.h>
#include <rdma/fi_endpoint.h>
#include <rdma/fi_atomic.h>
#include <rdma/fi_cxi_ext.h> // For CXI
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <sched.h>
#include <string.h>
#include <sys/time.h>

///////////////////////////////////////////////////////////////////////////
//
// Transport (OFI) info and functions
// 
///////////////////////////////////////////////////////////////////////////

struct ramcc_tp_s tp;

int ramcc_tp_get_myrank()
{
  return tp.pmi_myrank;
}

int ramcc_tp_get_numranks()
{
  return tp.pmi_numranks;
}

// Initialize the transport
int ramcc_tp_init()
{
  struct fi_info* hints;
  struct fi_info* all_providers;
  int err, pmi_err;
  tp.addr_len = ADDR_LENGTH;

  // Note: moved this initialization from below so as to have the rank 
  // available for purposes of selecting a NIC.
  pmi_err = pmi_kvs_setup();
  if (pmi_err != 0) {
    printf("Error setting up PMI KVS\n");
    exit(1);
  }

  tp.pmi_myrank = pmi_get_myrank();
  tp.pmi_numranks = pmi_get_numranks();
	tp.pmi_max_key_len = pmi_get_max_key_len();
	tp.pmi_max_val_len = pmi_get_max_val_len();
  
  /////////////////////////////////////////////////////////////////////
  // Select an appropriate NIC
  // TODO: The current method assumes we are running on El Cap or similar (MI300A).
  // TODO: Generalize (this will be non-trivial!)
  // TODO: See header file for system-specific defines
  /////////////////////////////////////////////////////////////////////
  // On other systems, we were able to include sched.h and then use:
  // int core = sched_getcpu();
  // However, this function is still not found on Eldorado
  // So here we use an alternative method to figure out which numa domain the 
  // process is on
  int hwt = ramcc_tp_get_hwt();
  int numa = (hwt / CORES_PER_NUMA) % NUM_NUMA;

  RAMC_DEBUG("Using hardware thread: %d In NUMA domain: %d\n", hwt, numa);

  hints = fi_allocinfo();

  // These hints are based on working through the CXI tests (specifically the code 
  // path invoked for the simple_write rma test) included with libfabric v 2.1
  // cxit_allocinfo_common: calls fi_allocinfo and sets the first bits of the hints:
  hints->fabric_attr->prov_name = strdup(LIBFABRIC_PROVIDER); // limits returned providers to cxi only
  hints->domain_attr->mr_mode = FI_MR_ENDPOINT | FI_MR_ALLOCATED | FI_MR_PROV_KEY; 
  // cxit_setup_enabled_ep: adds the following to the hints
  hints->domain_attr->data_progress = FI_PROGRESS_MANUAL;
  hints->domain_attr->control_progress = FI_PROGRESS_MANUAL; 
  // TODO: figure out what this should be set to, or what the default is
  hints->tx_attr->size = 8192; // will impact a fan-in / fan-out broacast, for example; not sure what the max is default; 131072 (the default CRAYPICH CQ size) is too large
  // on fi_display_info, using hints of provider name and mr_mode results in a tx_attr->msg_order being 
  // set to YES for many ordering options. 
  // Individual ordering bits can be set by e.g. hints->tx_attr->msg_order |= FI_ORDER_ATOMIC_RAR
  // For good measure, lets explicitly make them 0
  hints->tx_attr->msg_order = 0;

  // Select the local NIC, assuming numbering corresponds
  // TODO: This is specific to MI300A
  switch (numa) {
    case 0:
      hints->domain_attr->name = strdup(LIBFABRIC_DOMAIN_0);
      break;
    case 1:
      hints->domain_attr->name = strdup(LIBFABRIC_DOMAIN_1);
      break;
    case 2:
      hints->domain_attr->name = strdup(LIBFABRIC_DOMAIN_2);
      break;
    case 3:
      hints->domain_attr->name = strdup(LIBFABRIC_DOMAIN_3);
      break;
    defualt:
      RAMC_ERR("Unkown NUMA domain");
      return RAMC_FAILURE;
  }

  // Now the CXI test code does the getinfo
  FI_CHECK(fi_getinfo(FI_VERSION(1,15),0,0,0,hints,&all_providers));

  // just use the first one for the time being
  tp.provider = all_providers;

  // at this point fi_getinfo will return cxi providers with all of the capabilities (including RMA) except FI_SOURCE and FI_SOURCE_ERR
  // cxit_create_fabric_info: performs getinfo using current hints (so just the stuff above); then adds to the provider:
  tp.provider->ep_attr->tx_ctx_cnt = tp.provider->domain_attr->tx_ctx_cnt;
  tp.provider->ep_attr->rx_ctx_cnt = tp.provider->domain_attr->rx_ctx_cnt;
  // For CXI, the provider sets these to 'no' so if we want them we have to set them here
  tp.provider->caps |= FI_SOURCE | FI_SOURCE_ERR;
  tp.provider->rx_attr->caps |= FI_SOURCE | FI_SOURCE_ERR;

  fi_freeinfo(hints);

  // grab some info from the provider to reduce pointer chasing later
  tp.inject_size = tp.provider->tx_attr->inject_size;
    
  FI_CHECK(fi_fabric(tp.provider->fabric_attr, &tp.fabric, NULL));
  FI_CHECK(fi_domain(tp.fabric, tp.provider, &tp.domain, NULL));
  FI_CHECK(fi_endpoint(tp.domain, tp.provider, &tp.ep, NULL));

  // copied from cxit_create_eq
  struct fi_eq_attr eq_attr = {
    .size = 32,
    .flags = FI_WRITE,
    .wait_obj = FI_WAIT_NONE
  };
	FI_CHECK(fi_eq_open(tp.fabric, &eq_attr, &(tp.eq), NULL));

  FI_CHECK(fi_ep_bind(tp.ep, &(tp.eq->fid), 0));

  // CQ_FORMAT_CONTEXT is the simplest CQ entry format
  // CQs are bound with SELECTIVE_COMPLETION so entries are not generated unless explicitly requested 
  // through fi_*msg operations or there is an error
  struct fi_cq_attr tx_cq_attr = {
    .format = FI_CQ_FORMAT_CONTEXT,
  };
  struct fi_cq_attr rx_cq_attr = { 
      .format = FI_CQ_FORMAT_CONTEXT
  };
  FI_CHECK(fi_cq_open(tp.domain, &tx_cq_attr, &(tp.txcq), NULL));
  FI_CHECK(fi_cq_open(tp.domain, &rx_cq_attr, &(tp.rxcq), NULL));
  FI_CHECK(fi_ep_bind(tp.ep, &(tp.txcq->fid), FI_TRANSMIT | FI_SELECTIVE_COMPLETION));
  FI_CHECK(fi_ep_bind(tp.ep, &(tp.rxcq->fid), FI_RECV | FI_SELECTIVE_COMPLETION));

  struct fi_av_attr av_attr = {};
  // TODO: Should this be AV_MAP or AV_TABLE? 
  //av_attr.type = FI_AV_TABLE;
  av_attr.type = FI_AV_MAP;
  FI_CHECK(fi_av_open(tp.domain, &av_attr, &(tp.av), NULL));
  FI_CHECK(fi_ep_bind(tp.ep, &(tp.av->fid), 0));

  FI_CHECK(fi_enable(tp.ep));

//  pmi_err = pmi_kvs_setup();
//  if (pmi_err != 0) {
//    printf("Error setting up PMI KVS\n");
//    exit(1);
//  }
//
//  tp.pmi_myrank = pmi_get_myrank();
//  tp.pmi_numranks = pmi_get_numranks();
//	tp.pmi_max_key_len = pmi_get_max_key_len();
//	tp.pmi_max_val_len = pmi_get_max_val_len();

  err = ramcc_tp_init_barrier();

  // TODO: Do we have any cleanup?
  tp.cqe = (void *)malloc(sizeof(struct fi_cq_entry));

  FI_CHECK(fi_getname(&(tp.ep->fid), tp.my_addr, &(tp.addr_len)));

  char *pmi_key = (char *)malloc(tp.pmi_max_key_len);
  snprintf(pmi_key, tp.pmi_max_key_len, "%d", tp.pmi_myrank); 

  // insert my address into the PMI KV store
  // TODO: define PMI error macro in pmi_stuff.h
  err = pmi_kvs_insert_value(pmi_key, tp.my_addr, tp.addr_len);
  if (err != 0) {
    printf("%s:%s:%d: Rank %d: Could not insert my address in PMI KVS\n", __FILE__, __FUNCTION__, __LINE__, tp.pmi_myrank);
    exit(1);
  }
  err = pmi_commit_kvs();
  if (err != 0) {
    printf("%s:%s:%d: Rank %d: Could not commit KVS\n", __FILE__, __FUNCTION__, __LINE__, tp.pmi_myrank); 
    exit(1);
  }

  tp.addresses = (fi_addr_t *)malloc(tp.pmi_numranks * sizeof(fi_addr_t));
  if (NULL == tp.addresses) {
    RAMC_ERR("Could not allocate array to store addresses to insert into AV\n");
  }

  err = ramcc_tp_setup_av();
  if (err != 0) {
      RAMC_ERR("Error setting up AV\n");
      exit(1);
  }

  ramcc_tp_init_ep_counters();

  ///////////////////////////////////////////////////////////////////
  //
  // Bulletin Board Setup
  //
  // TODO: Should the BB share the buffer size (in bytes), too?
  // 
  ///////////////////////////////////////////////////////////////////

  // Allocate space for the two structs that comprise the BB
  // TODO: free this
  tp.bb_status = (struct ramcc_bb_status_s *)malloc(sizeof(struct ramcc_bb_status_s));
  tp.bb_posting = (struct ramcc_target_addr_info_s *)malloc(sizeof(struct ramcc_target_addr_info_s));

  // Initialize the BB status
  atomic_store(&(tp.bb_status->status), RAMC_BB_STATUS_INACTIVE);

  // Arrays, indexed by rank, with address info for each rank's BB
  // TODO: free this
  tp.bb_statuses = (struct fi_rma_iov *)malloc(tp.pmi_numranks * sizeof(struct fi_rma_iov));
  tp.bb_postings = (struct fi_rma_iov *)malloc(tp.pmi_numranks * sizeof(struct fi_rma_iov));
  
  // Register the BB memory regions
  // NOTE: If one tries to register (without bind and enable) more than one MR in a row, 
  // they will end up with the same memory key.
  //struct fid_mr *bb_status_mr = NULL;
  FI_CHECK(fi_mr_reg(tp.domain, tp.bb_status, sizeof(struct ramcc_bb_status_s), FI_REMOTE_READ, 0, 0, 0, &tp.bb_status_mr, NULL));
  FI_CHECK(fi_mr_bind(tp.bb_status_mr, &(tp.ep->fid), 0));
  FI_CHECK(fi_mr_enable(tp.bb_status_mr));

  //struct fid_mr *bb_posting_mr = NULL;
  FI_CHECK(fi_mr_reg(tp.domain, tp.bb_posting, sizeof(struct ramcc_target_addr_info_s), FI_REMOTE_READ, 0, 0, 0, &tp.bb_posting_mr, NULL));
  struct fi_cntr_attr cntr_attr = {};
  cntr_attr.events = FI_CNTR_EVENTS_COMP;
  cntr_attr.wait_obj = FI_WAIT_NONE;
  FI_CHECK(fi_mr_bind(tp.bb_posting_mr, &(tp.ep->fid), 0));
  FI_CHECK(fi_cntr_open(tp.domain, &cntr_attr, &(tp.bb_remote_read_cntr), NULL));
  // This should also count REMOTE_READS on CXI
  FI_CHECK(fi_mr_bind(tp.bb_posting_mr, &(tp.bb_remote_read_cntr->fid), FI_REMOTE_WRITE));
  FI_CHECK(fi_mr_enable(tp.bb_posting_mr));
  tp.bb_remote_read_cntr_val = fi_cntr_read(tp.bb_remote_read_cntr);

  // Exchange addressing info for BBs using PMI
  struct fi_rma_iov rma_iov_bb_status;
  struct fi_rma_iov rma_iov_bb_posting;
  rma_iov_bb_status.addr = (uint64_t)tp.bb_status;
  rma_iov_bb_status.key = fi_mr_key(tp.bb_status_mr);
  rma_iov_bb_posting.addr = (uint64_t)tp.bb_posting;
  rma_iov_bb_posting.key = fi_mr_key(tp.bb_posting_mr);

  // create key for bb status
  snprintf(pmi_key, tp.pmi_max_key_len, "bbs-%d", tp.pmi_myrank);
  err = pmi_kvs_insert_value(pmi_key, &rma_iov_bb_status, sizeof(struct fi_rma_iov));
  if (err != 0) {
    printf("%s:%s:%d: Rank %d: Could not insert my BB status info in PMI KVS\n", __FILE__, __FUNCTION__, __LINE__, tp.pmi_myrank);
    exit(1);
  }
  snprintf(pmi_key, tp.pmi_max_key_len, "bbp-%d", tp.pmi_myrank);
  err = pmi_kvs_insert_value(pmi_key, &rma_iov_bb_posting, sizeof(struct fi_rma_iov));
  if (err != 0) {
    printf("%s:%s:%d: Rank %d: Could not insert my BB status info in PMI KVS\n", __FILE__, __FUNCTION__, __LINE__, tp.pmi_myrank);
    exit(1);
  }
  err = pmi_commit_kvs();
  if (err != 0) {
    printf("%s:%s:%d: Rank %d: Could not commit KVS\n", __FILE__, __FUNCTION__, __LINE__, tp.pmi_myrank); 
    exit(1);
  }

  for (int i = 0; i < tp.pmi_numranks; ++i) {
    snprintf(pmi_key, tp.pmi_max_key_len, "bbs-%d", i);
    err = pmi_kvs_retrieve_value(pmi_key, &(tp.bb_statuses[i]), sizeof(struct fi_rma_iov));
    if (err != 0) {
      RAMC_ERR("Could not retrieve BB status address for key %s\n", pmi_key);
      return RAMC_FAILURE;
    }
    snprintf(pmi_key, tp.pmi_max_key_len, "bbp-%d", i);
    err = pmi_kvs_retrieve_value(pmi_key, &(tp.bb_postings[i]), sizeof(struct fi_rma_iov));
    if (err != 0) {
      RAMC_ERR("Could not retrieve BB posting address for key %s\n", pmi_key);
      return RAMC_FAILURE;
    }
  }

  ///////////////////////////////////////////////////////////////////
  //
  // DONE Bulletin Board Setup
  // 
  ///////////////////////////////////////////////////////////////////

  free(pmi_key);
  return RAMC_SUCCESS;
}

// initialize the barrier buffers
int ramcc_tp_init_barrier()
{
  tp.barrier_send[0] = 0;
  tp.barrier_recv[0] = 0;

  tp.barrier_send_mr = NULL;
  tp.barrier_recv_mr = NULL;
  
  return RAMC_SUCCESS;
}

int ramcc_tp_close_barrier()
{
  // Since we are not using VERBS nothing needs to be done
  return RAMC_SUCCESS;
}

int ramcc_tp_setup_av() 
{
    int err;
    char *pmi_key = (char *)malloc(tp.pmi_max_key_len);
    if (NULL == pmi_key) {
        RAMC_ERR("Could not allocate array to store keys\n");
        return RAMC_FAILURE;
    }
    char *addr_to_insert = (char *)malloc(tp.pmi_numranks * tp.addr_len);
    if (NULL == addr_to_insert) {
        RAMC_ERR("Could not allocate array to store addresses to insert into AV\n");
        return RAMC_FAILURE;
    }

    for (int i = 0; i < tp.pmi_numranks; ++i) {
        snprintf(pmi_key, tp.pmi_max_key_len, "%d", i);
        err = pmi_kvs_retrieve_value(pmi_key, &(addr_to_insert[i * tp.addr_len]), tp.addr_len);
        if (err != 0) {
            RAMC_ERR("Could not retrieve address for key %s when initializing AV\n", pmi_key);
            return RAMC_FAILURE;
       }
    }
  
    err = fi_av_insert(tp.av, addr_to_insert, tp.pmi_numranks, tp.addresses, 0, NULL);
    if (err < 0) {
        RAMC_ERR("Failed to insert array of addresses into AV\n");
        return RAMC_FAILURE;
    }
    if (err != tp.pmi_numranks) {
        RAMC_ERR("Failed to insert the correct number of addresses into AV (actual: %d inserted: %d)\n", tp.pmi_numranks, err);
        return RAMC_FAILURE;
    }

    free(pmi_key);
    free(addr_to_insert);

    return RAMC_SUCCESS;
}

int ramcc_tp_finalize() {
  free(tp.addresses);
  free(tp.cqe);
  free(tp.bb_statuses);
  free(tp.bb_postings);
  free(tp.bb_status);
  free(tp.bb_posting);
  FI_CHECK(fi_close(&(tp.bb_status_mr->fid)));
  FI_CHECK(fi_close(&(tp.bb_posting_mr->fid)));
  FI_CHECK(fi_close(&(tp.fabric->fid)));
  FI_CHECK(fi_close(&(tp.ep->fid)));
  FI_CHECK(fi_close(&(tp.bb_remote_read_cntr->fid)));
  ramcc_tp_close_ep_counters();
  FI_CHECK(fi_close(&(tp.av->fid)));
  FI_CHECK(fi_close(&(tp.rxcq->fid)));
  FI_CHECK(fi_close(&(tp.txcq->fid)));
  FI_CHECK(fi_close(&(tp.domain->fid)));

  return RAMC_SUCCESS;

//For reference: the order of closing in the minimal cxi example:
//    FI_CHECK(fi_close(&(recv_mr->fid)));
//    FI_CHECK(fi_close(&(tp.ep->fid)));
//    FI_CHECK(fi_close(&(tp.remote_write_cntr_mr->fid)));
//    close_ep_counters(&tp);
//    FI_CHECK(fi_close(&(tp.av->fid)));
//    FI_CHECK(fi_close(&(tp.rxcq->fid)));
//    FI_CHECK(fi_close(&(tp.txcq->fid)));
//    FI_CHECK(fi_close(&(tp.domain->fid)));
//    FI_CHECK(fi_close(&(tp.fabric->fid)));

}


// Wait for specified number of entries to appear on the specified CQ
// TODO: currently does not have a timeout
int ramcc_tp_await_completion_multi(struct fid_cq *cq, struct fi_cq_entry *cqe, const int num_entries) {
  int ret;
  int count = 0;

  while (count < num_entries) {
    do {
      ret = fi_cq_read(cq, cqe, 1);
    } while (ret == -FI_EAGAIN);

    if (ret != 1) {
      RAMC_ERR("fi_cq_read: %d", ret);
      struct fi_cq_err_entry ebuf = {0};
      int ret = fi_cq_readerr(cq, (void *)&ebuf, 0);
      if (ret > 0) {
        const char *errmsg = fi_cq_strerror(cq, ebuf.prov_errno, ebuf.err_data, NULL, 0);
        RAMC_ERR("Error from fi_cq_strerror: %s\n", errmsg);
        return RAMC_ERROR;
      }
    }
    count++;
  }
  return RAMC_SUCCESS;
}

// Wait for single entry to appear on specified CQ
int ramcc_tp_await_completion(struct fid_cq *cq, void *cqe, uint64_t timeout_ms)
{
	  int ret;
    struct timeval tv;
    uint64_t start;
    gettimeofday(&tv, NULL);
    start = (uint64_t)((tv.tv_sec * 1000) + (tv.tv_usec / 1000));

    if (0 == timeout_ms) {
      do {
        ret = fi_cq_read(cq, cqe, 1);
      //} while (ret != 1);
      } while (ret == -FI_EAGAIN);
    } else {
	    do {
		    ret = fi_cq_read(cq, cqe, 1);
        gettimeofday(&tv, NULL);
        if ((((uint64_t)(tv.tv_sec * 1000) + (tv.tv_usec/1000)) - start) > timeout_ms) {
          RAMC_ERR("Timed out (%lu ms timeout)\n", timeout_ms);
          return RAMC_TIMEOUT;
        }
	    } while (ret == -FI_EAGAIN);
    }

    if (ret != 1) {
      struct fi_cq_err_entry err_entry = {};
      ret = fi_cq_readerr(cq, &err_entry, 0);
      RAMC_ERR("fi_cq_read error: %s Provider error: %s\n", fi_strerror(err_entry.err), fi_cq_strerror(cq, err_entry.prov_errno, err_entry.err_data, NULL, 0));
      return RAMC_ERROR;
    }

  return RAMC_SUCCESS;
}

// When using SELECTIVE COMPLETION only errors will show up on a queue unless we 
// explicitly say otherwise. 
// Checks for unexpected entry on the specified CQ
int ramcc_tp_check_for_error(struct fid_cq *cq) {
  int ret;
  ret = fi_cq_read(cq, tp.cqe, 1);
  if (1 == ret) {
    printf("Rank %d: Unexpected CQ event!\n", tp.pmi_myrank);
  }
  return 0;
}

// Wait on or check counter value for expected value
// Function obeys the following rules:
//    If timeout_ms == 0: wait indefinitely; return SUCCESS if it gets through
//    If timeout_ms == 1: check counter value once; return SUCCESS if counter value >= expected, FAILURE otherwise
//    If timeout_ms > 1: check value until timeout reached; return SUCESS if counter value >= expected, TIMEOUT otherwise
int ramcc_tp_await_completion_cntr(struct fid_cntr *cntr, uint64_t expected_value, uint64_t *cntr_val, uint64_t timeout_ms)
{
    int ret;
    struct timeval tv;
    uint64_t start;

    // TODO: find a better place (?) for checking for errors
    ret = ramcc_tp_check_for_error(tp.txcq);
    ret = ramcc_tp_check_for_error(tp.rxcq);

    if (0 == timeout_ms) {
      do {
        *cntr_val = fi_cntr_read(cntr);
	    } while (*cntr_val < expected_value);
    } else if (0 == timeout_ms) {
      *cntr_val = fi_cntr_read(cntr);
      if (*cntr_val >= expected_value) {
        return RAMC_SUCCESS;
      } else {
        return RAMC_FAILURE;
      }
	  } else {
      gettimeofday(&tv, NULL);
      start = (uint64_t)((tv.tv_sec * 1000) + (tv.tv_usec / 1000));
      do {
		    *cntr_val = fi_cntr_read(cntr);
        gettimeofday(&tv, NULL);
        if ((((uint64_t)(tv.tv_sec * 1000) + (tv.tv_usec/1000)) - start) > timeout_ms) {
          RAMC_DEBUG("Timed out (%lu ms timeout, last read counter value: %lu)\n", timeout_ms, *cntr_val);
          return RAMC_TIMEOUT;
        }
	    } while (*cntr_val < expected_value);
    }
    return RAMC_SUCCESS;
}


// do a linear barrier (fan in, fan out)
// This version uses iovec structs so that fi_sndmsg/fi_recvmsg can be used with FI_COMPLETION
// TODO: check to see if this is generating events on the Tx CQ, and if so, whether they are being 
// cleaned up
int ramcc_tp_barrier_linear()
{
  int ret;
  struct iovec iov_send = {tp.barrier_send, sizeof(uint64_t)};
  struct fi_msg msg_send = {};
  msg_send.msg_iov = &iov_send;
  msg_send.iov_count = 1;

  struct iovec iov_recv = {tp.barrier_recv, sizeof(uint64_t)};
  struct fi_msg msg_recv = {};
  msg_recv.msg_iov = &iov_recv;
  msg_recv.iov_count = 1;

  if (0 == tp.pmi_myrank) {
    for (int p = 1; p < tp.pmi_numranks; ++p) {
      msg_recv.addr = tp.addresses[p];
      FI_CHECK(fi_recvmsg(tp.ep, &msg_recv, FI_COMPLETION));
    }
    ret = ramcc_tp_await_completion_multi(tp.rxcq, tp.cqe, tp.pmi_numranks-1);
    for (int p = 1; p < tp.pmi_numranks; ++p) {
      msg_send.addr = tp.addresses[p];
      // should this generate CQ entries and block on them?
      FI_CHECK(fi_sendmsg(tp.ep, &msg_send, 0));
    }
  } else {
    msg_send.addr = tp.addresses[0];
    FI_CHECK(fi_sendmsg(tp.ep, &msg_send, 0));
    msg_recv.addr = tp.addresses[0];
    FI_CHECK(fi_recvmsg(tp.ep, &msg_recv, FI_COMPLETION));
    ret = ramcc_tp_await_completion_multi(tp.rxcq, tp.cqe, 1);
  }
  return RAMC_SUCCESS;
}

// do a binary tree barrier (both up and down)
// breadth first numbering:
// children: Left: 2r+1, Right: 2(r+1)
// parent: floor((r-1)/2)
int ramcc_tp_barrier_binary()
{
  int ret;
  struct iovec iov_send = {tp.barrier_send, sizeof(uint64_t)};
  struct fi_msg msg_send = {};
  msg_send.msg_iov = &iov_send;
  msg_send.iov_count = 1;

  struct iovec iov_recv = {tp.barrier_recv, sizeof(uint64_t)};
  struct fi_msg msg_recv = {};
  msg_recv.msg_iov = &iov_recv;
  msg_recv.iov_count = 1;

  int parent = floor((tp.pmi_myrank-1)/2);
  int left_child = (2 * tp.pmi_myrank) + 1;
  int right_child = 2 * (tp.pmi_myrank + 1);

  int num_waits = 0;

  if (right_child < tp.pmi_numranks) {
    num_waits++;
    msg_recv.addr = tp.addresses[right_child];
    FI_CHECK(fi_recvmsg(tp.ep, &msg_recv, FI_COMPLETION));
  }
  if (left_child < tp.pmi_numranks) {
    num_waits++;
    msg_recv.addr = tp.addresses[left_child];
    FI_CHECK(fi_recvmsg(tp.ep, &msg_recv, FI_COMPLETION));
  }
  if (num_waits > 0) {
    ret = ramcc_tp_await_completion_multi(tp.rxcq, tp.cqe, num_waits);
  }
  if (tp.pmi_myrank > 0) {
    msg_send.addr = tp.addresses[parent];
    FI_CHECK(fi_sendmsg(tp.ep, &msg_send, 0));
  }

  if (tp.pmi_myrank > 0) {
    msg_recv.addr = tp.addresses[parent];
    FI_CHECK(fi_recvmsg(tp.ep, &msg_recv, FI_COMPLETION));
    ret = ramcc_tp_await_completion_multi(tp.rxcq, tp.cqe, 1);
  }
  if (left_child < tp.pmi_numranks) {
    msg_send.addr = tp.addresses[left_child];
    FI_CHECK(fi_sendmsg(tp.ep, &msg_send, 0));
  }
  if (right_child < tp.pmi_numranks) {
    msg_send.addr = tp.addresses[right_child];
    FI_CHECK(fi_sendmsg(tp.ep, &msg_send, 0));
  }
  return RAMC_SUCCESS;
}

void ramcc_tp_init_ep_counters()
{
  struct fi_cntr_attr cntr_attr = {};
  cntr_attr.events = FI_CNTR_EVENTS_COMP; // counter counts completion events
  cntr_attr.wait_obj = FI_WAIT_NONE;
  cntr_attr.flags = 0;

  tp.write_cntr_ep = NULL;
  tp.read_cntr_ep = NULL;

  tp.write_cntr_ep_val = 0;
  tp.write_cntr_ep_expected = 0;
  tp.read_cntr_ep_val = 0;
  tp.read_cntr_ep_expected = 0;

  FI_CHECK(fi_cntr_open(tp.domain, &cntr_attr, &(tp.write_cntr_ep), NULL));
  FI_CHECK(fi_cntr_open(tp.domain, &cntr_attr, &(tp.read_cntr_ep), NULL));

  FI_CHECK(fi_ep_bind(tp.ep, &(tp.write_cntr_ep->fid), FI_WRITE));
  FI_CHECK(fi_ep_bind(tp.ep, &(tp.read_cntr_ep->fid), FI_READ));
}

void ramcc_tp_close_ep_counters()
{
  FI_CHECK(fi_close(&(tp.write_cntr_ep->fid)));
  FI_CHECK(fi_close(&(tp.read_cntr_ep->fid)));
}

void ramcc_tp_print_ep_cntrs()
{
    RAMC_DEBUG("write_cntr_ep: %lu read_cntr_ep: %lu\n", fi_cntr_read(tp.write_cntr_ep), fi_cntr_read(tp.read_cntr_ep));
}

// returns the hardware thread the process is running on at this moment
// (it could migrate, but it should not migrate off the current NUMA domain (?)
int ramcc_tp_get_hwt() {
  char command[256];
  char output[256];
  char line[256];
  FILE *fp;
  pid_t pid;

  for (int i = 0; i < 256; ++i) command[i] = 32;
  pid = getpid();
  snprintf(command, 256, "ps -o psr %d", pid);
  fp = popen(command, "r");
  if (NULL == fp) {
    printf("Could not execute command\n");
    exit(1);
  }
  for (int i = 0; i < 2; ++i) {
    fgets(line, sizeof(line), fp);
  }
  return atoi(line);
}

//////////////////////////////////////////////////////////////////////
// End of Transport Routines
// Begin Core RAMC Routines
//
// TARGET
//
//////////////////////////////////////////////////////////////////////

// 
// Args:
// domain from transport struct
// buffer address
// size in bytes
// access options
// offset (must be 0; reserved for future use)
// MR key requested (since we use PROV_KEY it is ignored)
// Additional flags (none in this case)
// struct fid_mr **
// context (here we do not need it)

// Called by target to create target window
// TODO: Add access flags (REMOTE_READ, REMOTE_WRITE, WRITE, etc.)
// target_buf : the buffer to use as the target window
// target_len : length of buffer in bytes
// tag  : the tag
// target_info : struct to store the resulting addressing information
//
// returns: 0 if success, otherwise will fail FI_CHECK in some fashion and everything will crash
// TODO: fix this error behavior
int ramcc_target_create_window(void *target_buf, size_t target_len, uint64_t tag, uint64_t init_status_val, struct ramc_target_win_info_s *target_info) {

  // TARGET_STATUS_INITIALIZE is an even value; values less than it are reserved for e.g. abort
  //target_info->status_val = RAMC_TARGET_STATUS_INITIALIZE;
  target_info->status_val = init_status_val;
  target_info->buffer = target_buf; // keep a record of the location of the buffer

  target_info->rank = ramcc_tp_get_myrank();
  target_info->tag = tag;

  // register target window
  FI_CHECK(fi_mr_reg(tp.domain, target_buf, target_len, FI_REMOTE_WRITE | FI_REMOTE_READ, 0, 0, 0, &(target_info->mr_buf), NULL));
  FI_CHECK(fi_mr_bind(target_info->mr_buf, &(tp.ep->fid), 0));
  // open counter and bind to this memory region
  struct fi_cntr_attr cntr_attr = {};
  cntr_attr.events = FI_CNTR_EVENTS_COMP;
  cntr_attr.wait_obj = FI_WAIT_NONE;
  FI_CHECK(fi_cntr_open(tp.domain, &cntr_attr, &(target_info->window_cntr), NULL));
  FI_CHECK(fi_mr_bind(target_info->mr_buf, &(target_info->window_cntr->fid), FI_REMOTE_WRITE));
  FI_CHECK(fi_mr_enable(target_info->mr_buf));
  target_info->window_cntr_val = 0; // counter starts at zero

  // fill in the address information for the buffer
  target_info->addr_info.target.key = fi_mr_key(target_info->mr_buf);
  target_info->addr_info.target.addr = (uint64_t)target_buf;
  target_info->addr_info.target.len = target_len;

  // register status buffer 
  FI_CHECK(fi_mr_reg(tp.domain, (void *)&(target_info->status_val), sizeof(uint64_t), FI_REMOTE_READ, 0, 0, 0, &(target_info->mr_status), NULL));
  FI_CHECK(fi_mr_bind(target_info->mr_status, &(tp.ep->fid), 0));
  FI_CHECK(fi_mr_enable(target_info->mr_status));
  
  target_info->addr_info.status.key = fi_mr_key(target_info->mr_status);
  target_info->addr_info.status.addr = (uint64_t)&(target_info->status_val);
  target_info->addr_info.status.len = sizeof(uint64_t);
  
  RAMC_DEBUG("Created target window (target_addr: %lu target_size: %zu target_mr_key: %lu status_addr: %lu status_size: %zu status_mr_key: %lu tag: %lu)\n", target_info->addr_info.target.addr, target_info->addr_info.target.len, target_info->addr_info.target.key, target_info->addr_info.status.addr, target_info->addr_info.status.len, target_info->addr_info.status.key, target_info->tag);

  return RAMC_SUCCESS;
}

int ramcc_target_destroy_window(struct ramc_target_win_info_s *target_info) {
  // the buffer is bound to the counter so you have to close the MR before the counter
  target_info->status_val = RAMC_TARGET_STATUS_DESTROYED;
  FI_CHECK(fi_close(&(target_info->mr_buf->fid)));
  FI_CHECK(fi_close(&(target_info->mr_status->fid)));
  FI_CHECK(fi_close(&(target_info->window_cntr->fid)));

  return RAMC_SUCCESS;
}

// Posts target window information to bulletin board
// Does not activate the posting
int ramcc_target_post_window(struct ramc_target_win_info_s *target_info) {

  tp.bb_status->tag = target_info->tag;
  tp.bb_status->status = RAMC_BB_STATUS_INACTIVE;

  tp.bb_posting->target.addr = target_info->addr_info.target.addr;
  tp.bb_posting->target.len = target_info->addr_info.target.len;
  tp.bb_posting->target.key = target_info->addr_info.target.key;
  tp.bb_posting->status.addr = target_info->addr_info.status.addr;
  tp.bb_posting->status.len = target_info->addr_info.status.len;
  tp.bb_posting->status.key = target_info->addr_info.status.key;

  RAMC_DEBUG("Posted target addressing info to BB (tag: %lu, win_addr: %lu, win_size %zu, win_mr_key %lu, status_addr %lu, status_size %zu, status_mr_key %lu)\n", tp.bb_status->tag, \
                              tp.bb_posting->target.addr, \
                              tp.bb_posting->target.len, \
                              tp.bb_posting->target.key, \
                              tp.bb_posting->status.addr, \
                              tp.bb_posting->status.len, \
                              tp.bb_posting->status.key);

  return RAMC_SUCCESS;
}

//// Activates bulletin board posting
int ramcc_target_activate_bb(void) {

  atomic_store(&(tp.bb_status->status), RAMC_BB_STATUS_ACTIVE);

  RAMC_DEBUG("BB activated (tp.bb_status->status = %d)\n", tp.bb_status->status);

  return RAMC_SUCCESS;
}

//// Deactivates bulletin board posting
int ramcc_target_deactivate_bb(void) {

  atomic_store(&(tp.bb_status->status), RAMC_BB_STATUS_INACTIVE);

  RAMC_DEBUG("BB deactivated (tp.bb_status->status = %d)\n", tp.bb_status->status);

  return RAMC_SUCCESS;
}

// Wait for specified number of reads of BB posting to occur
// BLOCKING
// Note: the counter will count all operatiosn (writes, too), so the name of this function 
// is indicative of the fact only REMOTE_READS should be done with the BB buffer
// TODO: Confirm the BB target window info buffer is set to only accept REMOTE_READs
// Returns RAMC_SUCCESS or never returns
int ramcc_target_await_bb_reads(uint64_t expected_num_reads) {

  int ret;
  // Note: we don't want to use the current counter value here because a read may have 
  // arrived between setting the BB status to active and the wait started; so we use 
  // the cached last known value
  uint64_t target_value = tp.bb_remote_read_cntr_val + expected_num_reads;

  RAMC_DEBUG("Wating for %lu reads to occur on BB (current cntr_val: %lu expected cntr_val: %lu)\n", expected_num_reads, tp.bb_remote_read_cntr_val, target_value);

  // TODO: using timeout_ms = 0 will result in the following blocking
  ret = ramcc_tp_await_completion_cntr(tp.bb_remote_read_cntr, target_value, &(tp.bb_remote_read_cntr_val), 0);

  RAMC_DEBUG("BB remote read counter reached expected value (current = %lu, expected = %lu)\n", tp.bb_remote_read_cntr_val, target_value);

  // if it reaches this far the wait could only have succeeded
  return RAMC_SUCCESS;
}

// Test if the number of BB reads is where it is expected
// Non-blocking
// Note: the counter will count all operatiosn (writes, too), so the name of this function 
// is indicative of the fact only REMOTE_READS should be done with the BB buffer
// TODO: Confirm the BB target window info buffer is set to only accept REMOTE_READs
// Returns RAMC_SUCCESS or RAMC_FAILURE
int ramcc_target_test_bb_reads(uint64_t expected_num_reads) {

  int ret;
  uint64_t tmp_cntr_val;
  // Note: we don't want to use the current counter value here because a read may have 
  // arrived between setting the BB status to active and the wait started; so we use 
  // the cached last known value
  uint64_t target_value = tp.bb_remote_read_cntr_val + expected_num_reads;

  // using timeout_ms = 1 means check the counter once
  ret = ramcc_tp_await_completion_cntr(tp.bb_remote_read_cntr, target_value, &tmp_cntr_val, 1);

  if (ret == RAMC_SUCCESS) {
    tp.bb_remote_read_cntr_val = tmp_cntr_val;
    return RAMC_SUCCESS;
  } else {
    return RAMC_FAILURE;
  }
}

// Update the status for a specified target window
// This version takes the simple approach of just incrementing the status value by 1
// The idea is that READ ONLY is defined as even, and WRITE ONLY is defined as odd
//
// TODO: what about closing down a connection? We could reserve 0 and 1 to be signals to initiators that something has 
// changed, so start status with 2 (= read only). Setting a target win status to 0 would tell the initiator to discard the channel.
// This would cost an additional step when checking status from an initiator
//int ramcc_target_update_win_status(struct ramc_target_win_info_s *target_info) {
  //++(*(uint64_t *)target_info->addr_info.status.addr);
//  ++(target_info->status_val);
//  return RAMC_SUCCESS;
//}

//int ramcc_target_increment_win_status(struct ramc_target_win_info_s *target_info, uint64_t value) {
//  target_info->status_val = target_info->status_val + value;
//  return RAMC_SUCCESS;
//}

//int ramcc_target_set_win_status(struct ramc_target_win_info_s *target_info, uint64_t value) {
//  target_info->status_val = value;
//  return RAMC_SUCCESS;
//}

// Wait for specified number of operations (remote reads or remote writes) to target window 
// Blocking
int ramcc_target_await_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops) {

  int ret;
  uint64_t target_value = target_info->window_cntr_val + expected_num_ops;

  RAMC_DEBUG("Wating for %lu ops to occur on target buffer (current cntr_val: %lu expected cntr_val: %lu)\n", expected_num_ops, target_info->window_cntr_val, target_value);

  ret = ramcc_tp_await_completion_cntr(target_info->window_cntr, target_value, &(target_info->window_cntr_val), 0);

  RAMC_DEBUG("Target buffer counter remote operation counter reached expected value (current = %lu, expected = %lu)\n", target_info->window_cntr_val, target_value);

  // if it gets this far it must have suceeded
  return RAMC_SUCCESS;
}

// Test for specified number of operations (remote reads or remote writes) to target window 
// Checks the counter ONCE
int ramcc_target_test_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops) {

  int ret;
  uint64_t tmp_cntr_val;
  uint64_t target_value = target_info->window_cntr_val + expected_num_ops;

  RAMC_DEBUG("Wating for %lu ops to occur on target buffer (current cntr_val: %lu expected cntr_val: %lu)\n", expected_num_ops, target_info->window_cntr_val, target_value);
  
  ret = ramcc_tp_await_completion_cntr(target_info->window_cntr, target_value, &(target_info->window_cntr_val), 1);

  if (ret == RAMC_SUCCESS) {
    target_info->window_cntr_val = tmp_cntr_val;
    return RAMC_SUCCESS;
  } else {
    return RAMC_FAILURE;
  }
}

// returns BB status and tag
// BLOCKS on read endpoint counter
int ramcc_init_get_bb_status(int target_rank, uint8_t *status, uint64_t *tag) {
  
  int ret;
  struct ramcc_bb_status_s tmp_status;

  FI_CHECK(fi_read(tp.ep, &tmp_status, sizeof(struct ramcc_bb_status_s), NULL, (tp.addresses)[target_rank], 0, (tp.bb_statuses)[target_rank].key, NULL));
  tp.read_cntr_ep_expected++;
  // TODO: perhaps we need to do an error check here too (on CQ)?
  ret = ramcc_tp_await_completion_cntr(tp.read_cntr_ep, tp.read_cntr_ep_expected, &(tp.read_cntr_ep_val), 0);

  *status = tmp_status.status;
  *tag = tmp_status.tag;
  return ret;
}

// Read current BB post from target rank 
// Assumes status was previously checked
// Populates provided ramc_init_win_info_s struct with target window addressing information
// BLOCKs on single read of read endpoint counter
int ramcc_init_get_bb_posting(int target_rank, uint64_t init_status_val, struct ramc_init_win_info_s *target_info) {

  int ret;
  target_info->rank = target_rank;
  // TODO: assumes that the status of the target window does not need to be read, as it always starts in 
  // read-only mode. Might instead want to read the current status from the actual target window?
  target_info->status_val = init_status_val;

  //RAMR_DEBUG("Reading BB posting info: target rank %d key: %lu\n", target_rank, (tp.bb_postings)[target_rank].key);

  FI_CHECK(fi_read(tp.ep, &(target_info->addr_info), sizeof(struct ramcc_target_addr_info_s), NULL, (tp.addresses)[target_rank], 0, (tp.bb_postings)[target_rank].key, NULL));
  tp.read_cntr_ep_expected++;
  ret = ramcc_tp_await_completion_cntr(tp.read_cntr_ep, tp.read_cntr_ep_expected, &(tp.read_cntr_ep_val), 0); // BLOCK on read counter

  RAMC_DEBUG("Retrieved target info from BB: target_addr: %lu target_size: %zu target_mr_key: %lu status_addr: %lu status_size: %zu status_mr_key: %lu\n", target_info->addr_info.target.addr, target_info->addr_info.target.len, target_info->addr_info.target.key, target_info->addr_info.status.addr, target_info->addr_info.status.len, target_info->addr_info.status.key);

  return RAMC_SUCCESS;
}

// Returns the status value of the target window
// Blocking on successful read of target status
int ramcc_init_get_win_status(struct ramc_init_win_info_s *target_info, uint64_t *tgt_value) {

  int ret;

  FI_CHECK(fi_read(tp.ep, tgt_value, sizeof(uint64_t), NULL, (tp.addresses)[target_info->rank], 0, target_info->addr_info.status.key, NULL));
  tp.read_cntr_ep_expected++;
  ret = ramcc_tp_await_completion_cntr(tp.read_cntr_ep, tp.read_cntr_ep_expected, &(tp.read_cntr_ep_val), 0); // BLOCKS on read counter

  //RAMR_DEBUG("Returning target status (key: %lu) of window (key: %lu) on rank %d: %lu\n", target_info->addr_info.status.key, target_info->addr_info.target.key, target_info->rank, *status);

  return RAMC_SUCCESS;
}

int ramcc_init_put(void *source_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  int ret;
  if (length > tp.inject_size) {
    FI_CHECK(fi_write(tp.ep, source_buf, length, NULL, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key, NULL));
    tp.write_cntr_ep_expected++;
    RAMC_DEBUG("Used fi_write. Waiting for PUT completion using write counter (current: %lu expected: %lu)\n", tp.write_cntr_ep_val, tp.write_cntr_ep_expected);
    ret = ramcc_tp_await_completion_cntr(tp.write_cntr_ep, tp.write_cntr_ep_expected, &(tp.write_cntr_ep_val), 0);
  } else {
    RAMC_DEBUG("Used fi_inject_write. Do not need to wait on write endpoint counter\n");
    FI_CHECK(fi_inject_write(tp.ep, source_buf, length, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key));
    tp.write_cntr_ep_expected++;
  }
  
  //RAMC_DEBUG("Wrote %zu bytes to target rank : %d key: %lu offset: %zu\n", length, target_info->rank, target_info->addr_info.target.key, target_offset);

  return RAMC_SUCCESS;
}

// TODO: test the non-blocking versions
// NOTE: With inject_write the point is that returning from the call is sufficient, so no additional wait is requried. Using the NB version 
// introduces a superfluous wait unless one takes care they know the message size and can skip it
int ramcc_init_put_nb(void *source_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  int ret;
  if (length > tp.inject_size) {
    FI_CHECK(fi_write(tp.ep, source_buf, length, NULL, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key, NULL));
  } else {
    FI_CHECK(fi_inject_write(tp.ep, source_buf, length, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key));
  }
  tp.write_cntr_ep_expected++;

  return RAMC_SUCCESS;
}

int ramcc_init_await_all_puts() {
  int ret;
  RAMC_DEBUG("Waiting for PUT completion using write counter (current: %lu expected: %lu)\n", tp.write_cntr_ep_val, tp.write_cntr_ep_expected);
  ret = ramcc_tp_await_completion_cntr(tp.write_cntr_ep, tp.write_cntr_ep_expected, &(tp.write_cntr_ep_val), 0);
  
  //RAMC_DEBUG("Wrote %zu bytes to target rank : %d key: %lu offset: %zu\n", length, target_info->rank, target_info->addr_info.target.key, target_offset);

  return RAMC_SUCCESS;
}

int ramcc_init_get(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  int ret;
  FI_CHECK(fi_read(tp.ep, init_buf, length, NULL, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key, NULL));
  tp.read_cntr_ep_expected++;
  RAMC_DEBUG("Waiting for GET completion using write counter (current: %lu expected: %lu)\n", tp.read_cntr_ep_val, tp.read_cntr_ep_expected);
  ret = ramcc_tp_await_completion_cntr(tp.read_cntr_ep, tp.read_cntr_ep_expected, &(tp.read_cntr_ep_val), 0);
  
  //RAMC_DEBUG("Wrote %zu bytes to target rank : %d key: %lu offset: %zu\n", length, target_info->rank, target_info->addr_info.target.key, target_offset);

  return RAMC_SUCCESS;
}

int ramcc_init_get_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  int ret;
  FI_CHECK(fi_read(tp.ep, init_buf, length, NULL, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key, NULL));
  tp.read_cntr_ep_expected++;
  return RAMC_SUCCESS;
}

int ramcc_init_await_all_gets() {
  int ret;
  RAMC_DEBUG("Waiting for GET completion using write counter (current: %lu expected: %lu)\n", tp.read_cntr_ep_val, tp.read_cntr_ep_expected);
  ret = ramcc_tp_await_completion_cntr(tp.read_cntr_ep, tp.read_cntr_ep_expected, &(tp.read_cntr_ep_val), 0);
  return RAMC_SUCCESS;
}

//------- Atomics ------

// simplified to just do a fi_atomic with FI_SUM as the op and uint64_t as the type
// we don't need a source buf because it will always contain 1
// we don't need length (which in this case is a count) because that will always be 1
// TODO: is not waiting on completion counter OK? When does fi_atomic return?
int ramcc_init_atomic_inc(struct ramc_init_win_info_s *target_info, size_t target_offset) {
  int ret;
  uint64_t one = 1;
  target_offset = target_offset * sizeof(uint64_t); // convert count offset to byte offset
  FI_CHECK(fi_atomic(tp.ep, &one, 1, NULL, (tp.addresses)[target_info->rank], target_offset, target_info->addr_info.target.key, FI_UINT64, FI_SUM, NULL));
  tp.write_cntr_ep_expected++; 
  //ret = ramcc_tp_await_completion_cntr(tp.write_cntr_ep, tp.write_cntr_ep_expected, &(tp.write_cntr_ep_val), 0);
  
  //RAMC_DEBUG("Wrote %zu bytes to target rank : %d key: %lu offset: %zu\n", length, target_info->rank, target_info->addr_info.target.key, target_offset);

  return RAMC_SUCCESS;
}
