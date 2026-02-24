/* Remote Access Memory Channels (RAMC) 0.5
 * Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pmi.h>

#define PMI_ERR(...) do { \
    printf ("@ %s (%d): ", __FILE__, __LINE__); \
    printf (__VA_ARGS__); \
} while (0)

// TODO: Can the pmi_kvs_key and pmi_kvs_value variables be made non-global?
static char *pmi_kvs_name; 
static char *pmi_kvs_key; 
static char *pmi_kvs_value; 
static int pmi_max_key_len;
static int pmi_max_val_len;
static int pmi_numranks;
static int pmi_myrank;

// PMI keys and values are strings
// These helper functions encode a binary address to string and vice versa
int encode_to_string(const void *in, int in_len, char *out, int out_len)
{
  // confirm sufficient space for encoding
  if ((in_len * 2) + 1 > out_len) return 1;

  for (size_t i = 0; i < in_len; ++i) {
    snprintf(out + (2 * i), 3, "%02x", ((unsigned char *)in)[i]);
  }

  // don't need to add string termination because that is done by snprintf

  return 0;
}

int decode_to_binary(const char *in, void *out, size_t out_len) 
{
  if (out_len != strlen(in) / 2) return 1;

  for (size_t i = 0; i < out_len; ++i) {
    sscanf(in + (2 * i), "%02hhx", &((unsigned char*)out)[i]);
  }
  return 0;
}

int pmi_get_max_key_len()
{
  return pmi_max_key_len;
}

int pmi_get_max_val_len()
{
  return pmi_max_val_len;
}

int pmi_get_myrank()
{
  return pmi_myrank;
}

int pmi_get_numranks() 
{
  return pmi_numranks;
}

char *pmi_get_kvs_name()
{
  return pmi_kvs_name;
}

int pmi_do_barrier()
{
  int err = PMI_Barrier();
  return err;
}

// sets up the PMI KVS, gets number of ranks, my rank, etc.
int pmi_kvs_setup() {

  int pmi_initialized;
  int pmi_err;
  int pmi_max_name_len; // max length of the name of the PMI KVS store

  pmi_err = PMI_Initialized(&pmi_initialized);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Already initialized\n");
    return 1;
  }

  pmi_err = PMI_Init(&pmi_initialized);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not initialize PMI\n");
    return 1;
  }
  
  pmi_err = PMI_Get_size(&pmi_numranks);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not get number of ranks\n");
    return 1;
  }

  pmi_err = PMI_Get_rank(&pmi_myrank);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not get my rank\n");
    return 1;
  }

  pmi_err = PMI_KVS_Get_name_length_max(&pmi_max_name_len);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not get maximum KVS name length\n");
    return 1;
  }

  // allocate and store name of KVS
  // note this pointer is global
  pmi_kvs_name = (char *)malloc(pmi_max_name_len);
  if (NULL == pmi_kvs_name) {
    PMI_ERR("Could not malloc space for the PMI KVS name\n");
    return 1;
  }

  pmi_err = PMI_KVS_Get_my_name(pmi_kvs_name, pmi_max_name_len);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not get KVS name\n");
    return 1;
  }

  pmi_err = PMI_KVS_Get_key_length_max(&pmi_max_key_len);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not get maximum KVS key length\n");
    return 1;
  }

  pmi_err = PMI_KVS_Get_value_length_max(&pmi_max_val_len);
  if (pmi_err != PMI_SUCCESS) {
    PMI_ERR("Could not get maximum KVS value length\n");
    return 1;
  }

  pmi_kvs_key = (char *)malloc(pmi_max_key_len);
  if (NULL == pmi_kvs_key) {
    PMI_ERR("Could not malloc space for PMI KVS key\n");
    return 1;
  }
  
  pmi_kvs_value = (char *)malloc(pmi_max_val_len);
  memset(pmi_kvs_value, 0, pmi_max_val_len * sizeof(char));
  if (NULL == pmi_kvs_value) {
    PMI_ERR("Could not malloc space for PMI KVS value\n");
    return 1;
  }

#ifdef LIBFABRIC_DEBUG
  if (0 == my_info->pmi_myrank) {
    printf("Rank 0: KVS name                  : %s\n", pmi_kvs_name);
    printf("Rank 0: Maximum KVS name length   : %d\n", pmi_max_name_len);
    printf("Rank 0: Maximum KVS key length    : %d\n", pmi_max_key_len);
    printf("Rank 0: Maximum KVS value length  : %d\n", pmi_max_val_len);
  } 

#endif
  return 0;
}



// encode value into a string and insert it into the PMI KVS
// char *key: the key to use (string)
// void *value: the value, in binary, to insert into the KVS
// size_t value_len: the length (bytes) of the value to be inserted
// pmi_max_val_len: maximum length of value in chars
int pmi_kvs_insert_value(char *key, void *value, size_t value_len) {

  int err;

  // encode the value into a string (pmi_kvs_value) with a maximum length of pmi_max_val_len chars
  err = encode_to_string(value, value_len, pmi_kvs_value, pmi_max_val_len);
  if (err != 0) {
    PMI_ERR("Error encoding raw address into hex string\n");
    return 1;
  }

  // insert in PMI KVS
  err = PMI_KVS_Put(pmi_kvs_name, key, pmi_kvs_value);
  if (err != PMI_SUCCESS) {
    PMI_ERR("Error inserting address %s into PMI KVS\n", pmi_kvs_value);
    return 1;
  }

  return 0;
}

// char *key: the key to use to lookup value
// void *value: the value to be returned (gets inserted at this address)
// size_t value_len: the length of the value in bytes
int pmi_kvs_retrieve_value(char *key, void *value, size_t value_len) {
  int err;

  // get the string value
  err = PMI_KVS_Get(pmi_kvs_name, key, pmi_kvs_value, pmi_max_val_len);
  if (err != PMI_SUCCESS) {
    PMI_ERR("Could not retrieve key = %s from PMI KVS\n", key);
    return 1;
  }

  // decode into binary to return
  err = decode_to_binary(pmi_kvs_value, value, value_len);
  if (err != 0) {
    PMI_ERR("Could not decode value = %s to binary\n", pmi_kvs_value);
    return 1;
  }

  return 0;
}

// commit local KVs to other ranks
int pmi_commit_kvs() {
  int err;

  err = PMI_KVS_Commit(pmi_kvs_name);
  err = PMI_Barrier();
  if (err != PMI_SUCCESS) {
    PMI_ERR("PMI Barrier failed\n");
    return 1;
  }

  return 0;
}
