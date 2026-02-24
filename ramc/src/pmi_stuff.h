/* Remote Access Memory Channels (RAMC) 0.5
 * Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
 * certain rights in this software. */

#ifndef PMI_STUFF_H
#define PMI_STUFF_H

#include <stddef.h>

int encode_to_string(const void *, int, char *, int);
int decode_to_binary(const char *, void *, size_t); 
int pmi_get_max_key_len();
int pmi_get_max_val_len();
int pmi_get_myrank();
int pmi_get_numranks();
char *pmi_get_kvs_name();
int pmi_kvs_setup();
int pmi_kvs_insert_value(char *, void *, size_t);
int pmi_kvs_retrieve_value(char *, void *, size_t);
int pmi_commit_kvs();
int pmi_do_barrier();

#endif
