/* Remote Access Memory Channels (RAMC) 0.5
 * Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include "ramc_api.h"

int ramc_init() {
  return ramcc_tp_init();
}

int ramc_finalize() {
  return ramcc_tp_finalize();
}

//////////////////////////////////////////////////////////////
// Target
//////////////////////////////////////////////////////////////

int ramc_tgt_create_window(void *target_buf, size_t target_len, uint64_t tag, uint64_t init_status_val, struct ramc_target_win_info_s *target_info) {
  return ramcc_target_create_window(target_buf, target_len, tag, init_status_val, target_info);
}

int ramc_tgt_destroy_window(struct ramc_target_win_info_s *target_info) {
  return ramcc_target_destroy_window(target_info);
}

int ramc_tgt_post_window(struct ramc_target_win_info_s *target_info) {
  return ramcc_target_post_window(target_info);
}

int ramc_tgt_activate_bb() {
  return ramcc_target_activate_bb();
}

int ramc_tgt_deactivate_bb() {
  return ramcc_target_deactivate_bb();
}

int ramc_tgt_await_bb_reads(uint64_t expected_num_reads) {
  return ramcc_target_await_bb_reads(expected_num_reads);
}

int ramc_tgt_test_bb_reads(uint64_t expected_num_reads) {
  return ramcc_target_test_bb_reads(expected_num_reads);
}

int ramc_tgt_await_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops) {
  return ramcc_target_await_win_ops(target_info, expected_num_ops);
}

int ramc_tgt_test_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops) {
  return ramcc_target_test_win_ops(target_info, expected_num_ops);
}

//int ramc_tgt_update_win_status(struct ramc_target_win_info_s *target_info) {
//  return ramcc_target_update_win_status(target_info);
//}

int ramc_tgt_increment_win_status(struct ramc_target_win_info_s *target_info, uint64_t value) {
  target_info->status_val = target_info->status_val + value;
  return RAMC_SUCCESS;
}

int ramc_tgt_set_win_status(struct ramc_target_win_info_s *target_info, uint64_t value) {
  target_info->status_val = value;
  return RAMC_SUCCESS;
}

//int ramc_tgt_get_win_status(struct ramc_target_win_info_s *target_info) {
//  if (target_info->status_val & 0x1) {
//    return RAMC_TARGET_WRITE;
//  } else {
//    return RAMC_TARGET_READ;
//  }
//}

uint64_t ramc_tgt_get_win_status_raw(struct ramc_target_win_info_s *target_info) {
  return target_info->status_val;
}

//////////////////////////////////////////////////////////////
// Initiator
//////////////////////////////////////////////////////////////


// Check BB status
// Returns RAMC_SUCCESS if BB is active and tag matches
// Returns RAMC_FAILURE otherwise
int ramc_init_check_bb_status(int target_rank, uint64_t tag) {
  int ret;
  uint8_t bb_status;
  uint64_t bb_tag;
  ret = ramcc_init_get_bb_status(target_rank, &bb_status, &bb_tag);

  if (RAMC_BB_STATUS_ACTIVE == bb_status && bb_tag == tag) {
    return RAMC_SUCCESS;
  } else {
    return RAMC_FAILURE;
  }
}

// get status and tag from target BB
int ramc_init_get_bb_status(int target_rank, uint8_t *bb_status, uint64_t *bb_tag) {
  return ramcc_init_get_bb_status(target_rank, bb_status, bb_tag);
}

int ramc_init_get_bb_posting(int target_rank, uint64_t init_status_val, struct ramc_init_win_info_s *target_info) {
  return ramcc_init_get_bb_posting(target_rank, init_status_val, target_info);
}

// returns the target's status value
// and the state (calculated from that value)
int ramc_init_get_win_status(struct ramc_init_win_info_s *target_info, uint64_t *tgt_value, uint64_t *tgt_state) {
  int ret;
  ret = ramcc_init_get_win_status(target_info, tgt_value);
  *tgt_state = *tgt_value & 0x1;
  return ret;
}

int ramc_init_increment_status(struct ramc_init_win_info_s *target_info, uint64_t value) {
  target_info->status_val = target_info->status_val + value;
  return RAMC_SUCCESS;
}

int ramc_init_set_status(struct ramc_init_win_info_s *target_info, uint64_t value) {
  target_info->status_val = value;
  return RAMC_SUCCESS;
}

// Return RAMC_SUCCESS if initiator's status value == target's status value
// Return RAMC_TARGET_BEHIND if target status is behind
// Return RAMC_TARGET_AHEAD if targhet status is ahead
int ramc_init_check_win_status(struct ramc_init_win_info_s *target_info) {
  int ret;
  uint64_t tgt_value;

  ret = ramcc_init_get_win_status(target_info, &tgt_value);

  if (tgt_value == target_info->status_val) {
    return RAMC_SUCCESS;
  } else if (tgt_value < target_info->status_val) {
    return RAMC_TARGET_BEHIND;
  } else {
    return RAMC_TARGET_AHEAD;
  }
}

int ramc_put(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  return ramcc_init_put(init_buf, length, target_info, target_offset);
}

int ramc_put_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  return ramcc_init_put_nb(init_buf, length, target_info, target_offset);
}
  
int ramc_await_all_puts() {
  return ramcc_init_await_all_puts();
}

int ramc_get(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  return ramcc_init_get(init_buf, length, target_info, target_offset);
}

int ramc_get_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset) {
  return ramcc_init_get(init_buf, length, target_info, target_offset);
}

int ramc_await_all_gets() {
  return ramcc_init_await_all_gets();
}

//////////////////////////////////////////////////////////////
// Experimental 
//////////////////////////////////////////////////////////////

int ramc_atomic_inc(struct ramc_init_win_info_s *target_info, size_t target_offset) {
  return ramcc_init_atomic_inc(target_info, target_offset);
}


//////////////////////////////////////////////////////////////
// Utility
//////////////////////////////////////////////////////////////

int ramc_get_rank() {
  return ramcc_tp_get_myrank();
}

int ramc_get_numranks() {
  return ramcc_tp_get_numranks();
}

int ramc_barrier_linear() {
  return ramcc_tp_barrier_linear();
}

int ramc_barrier_binary() {
  return ramcc_tp_barrier_binary();
}
