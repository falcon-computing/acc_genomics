/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: OpenCL function headers for smith-waterman FPGA kernel call 
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/
#ifndef SMITHWATERMANHOST_H
#define SMITHWATERMANHOST_H
#include <CL/opencl.h>
#include <CL/cl_ext.h>
int _init_opencl(const char * bitstream);
int _init_kernel_buffer();
void _smithWatermanRun(char *inputs,int refLength,int batchSize,int overhang_strategy,int w_match,int w_mismatch,int w_open,int w_extend,short *outputs);
int _release_smithWaterman();
#endif
