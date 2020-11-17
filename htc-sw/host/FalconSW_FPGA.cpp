/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: Wrapper Functions for Xilinx Smith-Waterman FPGA-Host interface, 
 *              Mostly data format transformation between host and FPGA
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/
#include <time.h>
#include "common.h"

#include "smithWatermanHost.h"

#define MAX_FPGA_SEQ_LENGTH 512

bool FalconSWFPGA_init(char* bitstream){
    static bool init = false;
    if(init)
        return true;
    if(!init){
        if(!_init_opencl(bitstream))
            return false;
        _init_kernel_buffer();
        init = true;
    }
    return true;
}
double FalconSWFPGA_run(char* ref, int refLength, char alts[][MAX_SEQ_LENGTH], int* altLengths, int batchSize, int overhang_strategy, int w_match, int w_mismatch, int w_open, int w_extend, struct Cigar* cigarResults, int* alignmentOffsets, bool isFPGA){
    struct timespec time1, time2, time_diff;
    double pure_kernel_time = 0;
    int max_length = -1;
    int min_length = 2 * MAX_SEQ_LENGTH;
    if(refLength > max_length)
        max_length = refLength;
    if(refLength < min_length)
        min_length = refLength;
    for(int i = 0; i < batchSize; i++){
        if(altLengths[i] > max_length)
            max_length = altLengths[i];
        if(altLengths[i] < min_length)
            min_length = altLengths[i];
    }
    if((!isFPGA) || batchSize <= 0 || max_length >= MAX_FPGA_SEQ_LENGTH - 1 || min_length <= 0){
        //if FPGA does not exist or sequence length is over 512 then use AVX impl
        clock_gettime(CLOCK_REALTIME, &time1);
        SWPairwiseAlignmentMultiBatch(ref, refLength, alts, batchSize, altLengths, cigarResults, alignmentOffsets, overhang_strategy, 0);
        clock_gettime(CLOCK_REALTIME, &time2);
        time_diff = diff_time(time1, time2);
        pure_kernel_time += (double)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);
        return pure_kernel_time;
    }
    
    char inputs[MAX_FPGA_SEQ_LENGTH * (MAX_BATCH_SIZE + 1) + MAX_BATCH_SIZE];
    for(int i = 0; i < batchSize; i++){
        inputs[i * 2] = altLengths[i] & 0xff;
        inputs[i * 2 + 1] = ((signed short)(altLengths[i] & 0xff00) >> 8);
    }
    for(int i = 0; i < MAX_FPGA_SEQ_LENGTH; i++){
        inputs[i + 2 * batchSize] = ref[i];
    }
    for(int i = 0; i < batchSize; i++){
        for(int j = 0; j < MAX_FPGA_SEQ_LENGTH; j++){
            inputs[2 * batchSize + MAX_FPGA_SEQ_LENGTH + i * MAX_FPGA_SEQ_LENGTH + j] = alts[i][j];
        }
    }
    
    short outputs[(2 * MAX_FPGA_SEQ_LENGTH + 2) * MAX_BATCH_SIZE];

    clock_gettime(CLOCK_REALTIME, &time1);
    _smithWatermanRun(&inputs[0], refLength, batchSize, overhang_strategy, w_match, w_mismatch, w_open, w_extend, &outputs[0]);
    clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    pure_kernel_time += (double)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);
     
    int offset_ptr = batchSize + 2;

    for(int i = 0; i < batchSize; i++){
        //the first 2 shorts value are cigarInfoLength
        cigarResults[i].CigarElementNum = outputs[i + 2];
    }
    for(int i = 0; i < batchSize; i++){
        offset_ptr += 2 * (cigarResults[i].CigarElementNum) + 1;
        alignmentOffsets[i] = outputs[offset_ptr - 1];
        for(int j = 0; j < cigarResults[i].CigarElementNum; j++){
            cigarResults[i].cigarElements[j].length = outputs[offset_ptr - 3 - 2 * j];
            cigarResults[i].cigarElements[j].state = outputs[offset_ptr - 2 - 2 * j];
        }
    }
    return pure_kernel_time;
}

void FalconSWFPGA_release(){
    _release_smithWaterman();
}
