#ifndef OCL_H
#define OCL_H

#include "host/host_types.h"
#include "CL/cl.h"
#include <CL/cl_ext.h>
#include "host/xcl.h"

#define BWT_SIZE 775451201

struct smemOclCtx{
    cl_platform_id platform_id;         // platform id
    cl_device_id device_id;             // compute device id 
    cl_context context;                 // compute context
    cl_program program;
    cl_command_queue command;          // compute command queue
    cl_kernel kernel_smem[BANK_NUM];
    cl_kernel kernel_dram[BANK_NUM];
    cl_mem bwt_buffer[BANK_NUM];
    cl_mem bwt_para_buffer[BANK_NUM];
    cl_mem seq_buffer[BANK_NUM];
    cl_mem seq_len_buffer[BANK_NUM];
    cl_mem mem_buffer[BANK_NUM];
    cl_mem mem_num_buffer[BANK_NUM];
    bwtintv_t* host_mem_buffer[BANK_NUM];
};


void ocl_init(char* btsm, const uint32_t* bwt, const uint64_t* bwt_para, uint64_t bwt_size, bwtintv_t* mem, int batch_size);

int smem_ocl(char* btsm, const uint32_t* bwt, const uint64_t* bwt_para, uint8_t* seq, uint8_t* seq_len, \
        int batch_size, bwtintv_t* mem_output, int* mem_num, double kernel_time[BANK_NUM]);
#endif
