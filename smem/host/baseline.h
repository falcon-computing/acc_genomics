#ifndef BASELINE_H
#define BASELINE_H
#include <stdio.h>
#include "host/host_types.h"
#include "host/ocl.h"

#define bwt_set_intv1(c, ik) ((ik).x[0] = L2_baseline[(int)(c)]+1, (ik).x[2] = L2_baseline[(int)(c)+1]-L2_baseline[(int)(c)], (ik).x[1] = L2_baseline[3-(c)]+1, (ik).info = 0)

void mem_collect_intv_new(uint32_t *bwt, int len, const uint8_t *seq, smem_aux_t *a);
int smem_baseline(const uint32_t* bwt, const uint64_t* bwt_para, uint8_t* seq, uint8_t* seq_len, int batch_size, bwtintv_t* mem_output, int* mem_num, double mem_request_size[BANK_NUM]);

#endif
