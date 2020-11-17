#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ap_int.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <dirent.h>
#include "common/common.h"

typedef ap_uint<4> uint4_t; 
typedef ap_uint<512> uint512_t; 

typedef struct {
	bwtint_t x0, x1, x2, info;
} bwtintv_t_fpga;

#define MAX_INTV 0
#define MAX_INTV_ALLOC 200
#define BATCH_SIZE 20
#define SEQ_LENGTH 150

#ifdef HLS_
void mem_collect_intv_fpga(
uint512_t *bwt_0, uint512_t *bwt_1, uint512_t *bwt_2, uint512_t *bwt_3, uint8_t *seq, bwtintv_t_fpga *mem_output, int* mem_num);
#endif

 
