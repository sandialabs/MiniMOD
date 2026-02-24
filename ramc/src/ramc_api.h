/* Remote Access Memory Channels (RAMC) 0.5
 * Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

/**
* @file
* @author Whit Schonbein <wwschon@sandia.gov>
* @version 0.5
*
* @section LICENSE
*
* License to be determined.
*
* @section DESCRIPTION
*
* The RAMC user-level API.
*/

#ifndef __RAMC_API_H__
#define __RAMC_API_H__

#include "ramc_core.h"

// /////////////////////////////////////////////////////////
// Setup and Teardown
// /////////////////////////////////////////////////////////
/**
* @brief Initializes RAMC
* @return RAMC_SUCCESS if initialization succeeds; RAMC_FAILURE otherwise.
*
*
* All RAMC processes call this function before any others. Uses PMI to 
* exchange libfabric endpoint addresses, sets up a bulletin board on each 
* process and exchanges addressing information for each processes' bulletin 
* board, among other things. 
*
*/
int ramc_init(void);

/**
* @brief Finalizes RAMC
* @return RAMC_SUCCESS if finalization succeeds; RAMC_FAILURE otherwise.
*
*
* The last RAMC function called by all RAMC processes. Cleans up libfabric 
* resources.
*
*/
int ramc_finalize(void);

// /////////////////////////////////////////////////////////
// Target API
// /////////////////////////////////////////////////////////

/**
* @brief Creates window that will be the target of remote operations.
*
* @param[in] target_buf The buffer that is the ultimate target of remote operations
* @param[in] target_len Length (in bytes) of the target buffer
* @param[in] tag Tag for matching during sharing addressing information via the bulletin board
* @param[in] init_status_val Value to initialize target window status value
* @param[out] target_info Structure that will hold target window information
* @return RAMC_SUCCESS if finalization succeeds; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_create_window(void *target_buf, size_t target_len, uint64_t tag, uint64_t init_status_val, struct ramc_target_win_info_s *target_info);

/**
* @brief Frees resources associated with specified target window
*
* @param[in] target_info Structure specifying target window 
* @return RAMC_SUCCESS if finalization succeeds; RAMC_FAILURE otherwise.
* 
* \remark This function closes (fi_close) (i) the data memory region, (ii) the data memory 
* region counter, and (iii) the status memory region. It does NOT free the data buffer. 
*
* \todo Currently an attempt to remotely access a destroyed target will result in 
* OFI throwing an error. Leaving the status live (i.e., not closing it) and setting it to 
* RAMC_TARGET_STATUS_DESTROYED is not an easy option because then the user-provided target 
* window structure has to be kept around until the end of the program. We should check to see 
* if we can handle the OFI error in a way that would return RAMC_FAILURE (e.g.).
*
*/
int ramc_tgt_destroy_window(struct ramc_target_win_info_s *target_info);

/**
* @brief Posts target window addressing information to bulletin board
*
* @param[in] target_info Structure specifying target window 
* @return RAMC_SUCCESS if finalization succeeds; RAMC_FAILURE otherwise.
*
* \remark This function does not set the bulletin board status to ACTIVE
*
*/
int ramc_tgt_post_window(struct ramc_target_win_info_s *target_info);

/**
* @brief Sets bulletin board status to ACTIVE
*
* @return RAMC_SUCCESS if finalization succeeds; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_activate_bb(void);

/**
* @brief Sets bulletin board status to INACTIVE
*
* @return RAMC_SUCCESS if finalization succeeds; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_deactivate_bb(void);

/**
* @brief Waits on the bulletin board remote read counter to reach the specified value
*
* @param[in] expected_num_reads Expected number of reads
*
* \remark This call is blocking, and will not return until the specified number of 
* reads has been reached.
*
* @return RAMC_SUCCESS when the counter has reached the specified value.
*
*/
int ramc_tgt_await_bb_reads(uint64_t expected_num_reads);

/**
* @brief Tests if the bulletin board remote read counter has reached the specified value
*
* @param[in] expected_num_reads Expected number of reads
*
* \remark This call is non-blocking
*
* @return RAMC_SUCCESS if the bulletin board counter has reached the specified value; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_test_bb_reads(uint64_t expected_num_reads);

/**
* @brief Waits on the target buffer memory region counter to reach the specified value 
*
* @param[in] target_info Pointer to target window
* @param[in] expected_num_reads Expected number of operations
*
* \remark This call is blocking, and will not return until the expected number 
* of operations has been reached.  
*
* \todo Target windows are currently set up to count REMOTE_READs and REMOTE_WRITEs; we can make 
* this more flexible by allowing the user to specify the types of permitted operations during 
* window creation. 
* This will also require test and await functions to also take an argument 
* indicating what to wait for.
*
* @return RAMC_SUCCESS when the bulletin board counter has reached the specified value
*
*/
int ramc_tgt_await_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops);

/**
* @brief Tests for the target buffer memory region counter to reach the specified value 
*
* @param[in] target_info Pointer to target window
* @param[in] expected_num_reads Expected number of operations
*
* \remark This call is NON-BLOCKING
*
* @return RAMC_SUCCESS if the bulletin board counter has reached the specified value; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_test_win_ops(struct ramc_target_win_info_s *target_info, uint64_t expected_num_ops);

/**
* @brief Increments the target window status by the specified amount
*
* @param[in] target_info Pointer to target window
* @param[in] value Increment by this amount
*
* @return RAMC_SUCCESS if the status is incremented; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_increment_win_status(struct ramc_target_win_info_s *target_info, uint64_t value);

/**
* @brief Sets the target window status to the specified value 
*
* @param[in] target_info Pointer to target window
* @param[in] value Set status to this value
*
* @return RAMC_SUCCESS if the status is set; RAMC_FAILURE otherwise.
*
*/
int ramc_tgt_set_win_status(struct ramc_target_win_info_s *target_info, uint64_t value);

/**
* @brief Retrieves the target window status value
*
* @param[in] target_info Pointer to target window
*
* @return Target window status value
*
* \remark This function is provided so a user can construct their own logic around 
* the target window status value. The value can also be accessed directly from the 
* target_info structure (target_info->status_val).
*
*/
uint64_t ramc_tgt_get_win_status_raw(struct ramc_target_win_info_s *target_info);

// /////////////////////////////////////////////////////////
//
// Initiator API
//
// /////////////////////////////////////////////////////////

/**
* @brief Check whether the target's bulletin board status is active and the 
* specified tag matches.
*
* @param[in] target_rank Rank whose bulletin board is to be checked
* @param[in] tag The tag that must match if the bulletin board is active
*
* @return RAMC_SUCCESS if that target's bulletin board is active and the tag 
* matches; RAMC_FAILURE otherwise
*
* \remark This call is BLOCKING insofar as the call blocks until the read 
* completes (i.e., until the target's status and tag has been received as indicated 
* by the local endpoint FI_READ counter). A non-blocking version would not wait 
* for the endpoint counter to increment.
*
* \todo Do we want a non-blocking version?
*
*/
int ramc_init_check_bb_status(int target_rank, uint64_t tag);

/**
* @brief Retrieve the target's bulletin board status and tag.
*
* @param[in] target_rank Rank whose bulletin board is to be checked
* @param[out] status The target's bulletin board status
* @param[out] tag The target's bulletin board tag
*
* @return RAMC_SUCCESS if that target's bulletin board is active and the tag 
* matches; RAMC_FAILURE otherwise
*
* This function is provided in the event the user wants to inspect the 
* status and tag of a target. 
*
* \remark This call is BLOCKING insofar as the call blocks until the read 
* completes (i.e., until the target's status and tag has been received as indicated 
* by the local endpoint FI_READ counter). A non-blocking version would not wait 
* for the endpoint counter to increment.
*
* \todo Do we want a non-blocking version?
*
*/
int ramc_init_get_bb_status(int target_rank, uint8_t *status, uint64_t *tag); 

/**
* @brief Retrieve the target's bulletin board posting (i.e., target window and 
* status addressing information).
*
* @param[in] target_rank Rank whose bulletin board is to be read 
* @param[in] init_status_value The local initiator status value for the channel should be set to this value
* @param[out] target_info Target addressing information
*
* @return RAMC_SUCCESS if target bulletin board information is retrieved successfully
*
* \remark This call is BLOCKING insofar as it waits for the endpoint read counter to 
* increment before returning; a non-blocking version would not wait for the data 
* to be delivered to the initiator.
*
* \todo We reserve status values of 0 and 1 for other purposes, so the minimum 
* value should be 2. Add a check.
* \todo This call is based on the assumption the initiator is joining the channel 
* before the target has engaged in any other channels. How could we let an initiator 
* create a channel with a target that has already been engaged with other initiators? E.g., 
* add a flag that tells the initiator to GET the status value from the target and use 
* that as the initial value.
*
*/
int ramc_init_get_bb_posting(int target_rank, uint64_t init_status_val, struct ramc_init_win_info_s *target_info);

/**
* @brief Increment the initiator's channel status by the specified value
*
* @param[in] target_info Pointer to structure encapsulating information about the target of the channel
* @param[in] value Amount by which to increment the initiator channel status value
*
* @return RAMC_SUCCESS if initiator status value is updated succesfully
*
*/
int ramc_init_increment_status(struct ramc_init_win_info_s *target_info, uint64_t value);

/**
* @brief Set the initiator's channel status to the specified value
*
* @param[in] target_info Pointer to structure encapsulating information about the target of the channel
* @param[in] value Value to set the initiator status value 
*
* @return RAMC_SUCCESS if initiator status value is updated succesfully
*
*/
int ramc_init_set_status(struct ramc_init_win_info_s *target_info, uint64_t value);

/**
* @brief Compares initiator's status value with that of the target.
*
* @return RAMC_SUCCESS if initiator's status value == target's status value; 
* RAMC_TARGET_BEHIND if target status is behind (less than initiator's status);
* RAMC_TARGET_AHEAD if target status is ahead (greater than initiator's status).
*
* \remark This call is BLOCKING in the sense it does not return until the 
* target's status value has been retrieved.
*/
int ramc_init_check_win_status(struct ramc_init_win_info_s *target_info);

/**
* @brief Retrieves target's status value.
*
* \remark This call is BLOCKING in the sense it does not return until the 
* target's status value has been retrieved.
*
* \todo This function assumes odd/even state. It needs to be updated to just 
* return the raw target state
*/
int ramc_init_get_win_status(struct ramc_init_win_info_s *target_info, uint64_t *tgt_value, uint64_t *tgt_state);

// /////////////////////////////////////////////////////////
// Communication operations (all called by initiator)
// /////////////////////////////////////////////////////////

/**
* @brief Put data to target (BLOCKING)
*
* @param[in] init_buf Pointer to data to be written
* @param[in] length Number of bytes to be written
* @param[in] target_info Structure with target addressing information (obtained from bulletin board)
* @param[in] target_offset Offset into target buffer in bytes (beginning at zero)
*
* \remark This call is BLOCKING in the sense it does not return until the 
* initiator's write endpoint counter has been updated (i.e., until an ACK is received from the 
* target indicating the data was received by the target NIC).
*
* \remark Put uses fi_inject_write if the length is less than or equal to the 
* maximum injection size indicated by the CXI provider (192 bytes the last time 
* we checked); otherwise, fi_write is used.
*
*/
int ramc_put(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset);

/**
* @brief Put data to target (NON-BLOCKING)
*
* @param[in] init_buf Pointer to data to be written
* @param[in] length Number of bytes to be written
* @param[in] target_info Structure with target addressing information (obtained from bulletin board)
* @param[in] target_offset Offset into target buffer in bytes (beginning at zero)
*
* \remark This call is NON-BLOCKING in the sense it returns without waiting 
* for the initiator's write endpoint counter to be updated. The application 
* should wait (ramc_await_all_puts) at some point in the future.
*
*/
int ramc_put_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset);

/**
* @brief Wait on all outstanding put operations (BLOCKING)
*
* \remark Each call to a RAMC put operation increments an internal value tracking the 
* number of puts that have been issued by the application. This value is 
* shared across all channels. This function waits (BLOCKING) for the remote 
* write endpoint counter to reach this internal value.
*
* \todo See if there is a way to add context to endpoint counters so that an 
* initiator can distinguish between puts issued on different channels. 
*/
int ramc_await_all_puts(); 

/**
* @brief Get data from target (BLOCKING)
*
* @param[in] init_buf Pointer to where data from target should be placed 
* @param[in] length Number of bytes to be retrieved 
* @param[in] target_info Structure with target addressing information (obtained from bulletin board)
* @param[in] target_offset Offset into target buffer in bytes (beginning at zero)
*
* \remark This call is BLOCKING in the sense it does not return until the 
* initiator's read endpoint counter has incremented (indicating the data is 
* visible in the initiator's buffer).
*
*/
int ramc_get(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset);

/**
* @brief Get data from target (NON-BLOCKING)
*
* @param[in] init_buf Pointer to where data from target should be placed 
* @param[in] length Number of bytes to be retrieved 
* @param[in] target_info Structure with target addressing information (obtained from bulletin board)
* @param[in] target_offset Offset into target buffer in bytes (beginning at zero)
*
* \remark This call is NON-BLOCKING in the sense it may return before the 
* initiator's read counter has updated. That is, returning from this call is not 
* a guarantee the data to be retrieved from the target is visible to the caller.
*
*/
int ramc_get_nb(void *init_buf, size_t length, struct ramc_init_win_info_s *target_info, size_t target_offset);

/**
* @brief Wait on all outstanding get operations (BLOCKING)
*
* \remark Each call to a RAMC get operation increments an internal value tracking the 
* number of gets issued by the application. This value is 
* shared across all channels. This function waits (BLOCKING) for the remote 
* write endpoint counter to reach this internal value.
*
* \todo See if there is a way to add context to endpoint counters so that an 
* initiator can distinguish between puts issued on different channels. 
*
*/
int ramc_await_all_gets();

// /////////////////////////////////////////////////////////
// Utilities
// /////////////////////////////////////////////////////////

/**
* @brief Get the rank of the calling RAMC process.
*
* \remark Ranks are assigned by PMI 
*
*/
int ramc_get_rank();

/**
* @brief Get the total number of processes in the RAMC job.
*
* \remark Nothing remarkable
*/
int ramc_get_numranks();

/**
* @brief Perform a linear barrier.
*
* \remark This is a simple N-to-1 linear barrier, included as a utility. It uses 
* the libfabric send and receive message calls, with buffers included in the 
* RAMC transport structure.
*
*/
int ramc_barrier_linear();

/**
* @brief Perform a linear barrier.
*
* \remark This is a simple binary-tree barrier, included as a utility. It uses 
* the libfabric send and receive message calls, with buffers included in the 
* RAMC transport structure.
*
*/
int ramc_barrier_binary();

// /////////////////////////////////////////////////////////
// Experimental
// /////////////////////////////////////////////////////////
int ramc_atomic_inc(struct ramc_init_win_info_s *target_info, size_t target_offset);

#endif // __RAMC_API_H__
