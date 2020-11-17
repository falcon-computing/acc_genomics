#ifndef PAIRHMM_H
#define PAIRHMM_H
#define __constant
#define __kernel
#define __global
#include <string.h> 
#include <stdio.h>
#include <stdlib.h>
#include "hls_stream.h"
#include "ap_int.h"
#include "ap_utils.h"
#include "hls_half.h"

#include <math.h>
#include "common/ph2pr.h"
#include "common/ph2pr_sub1.h"
#include "common/ph2pr_div3.h"


#include "common/common.h"

#define PH2PR_SUB1_10 float(0.9000000059604644775390625)
#define PH2PR_10 float(0.0999999940395355224609375)

const int read_fifo_depth = 2 * (1 + MAX_READ_LEN * READ_BLOCK_SIZE);
const int read_block_size = READ_BLOCK_SIZE;
const int hap_width = 4 * HAP_BLOCK_SIZE;

//this is the pipeline stages of the computeEngine pipeline, modify this based on that shown in HLS report

#define READ_INFO_OFFSET 0
#define READ_DATA_OFFSET ((MAX_RSDATA_NUM / READ_BLOCK_SIZE) * 64  / 512)
#define HAP_LEN_OFFSET READ_DATA_OFFSET + (MAX_RSDATA_NUM * MAX_READ_LEN * 64 / 512)
#define HAP_DATA_OFFSET HAP_LEN_OFFSET + (MAX_HAPDATA_NUM * 64 / 512)
#define NUMREADPU_OFFSET HAP_DATA_OFFSET + (16 * MAX_HAPDATA_NUM * MAX_HAP_LEN / HAP_BLOCK_SIZE / 512)
#define ITERNUM_OFFSET NUMREADPU_OFFSET + (16 * 64 / 512)
#define END_OFFSET ITERNUM_OFFSET + (64 * 64 / 512)

/*
typedef struct{
    float M;
    float X;
    float Y;
}cellData;*/

typedef float custom_float;

typedef struct{
    custom_float M;
    custom_float X;
    custom_float Y;
}cellData;

typedef union{
    int hex;
    float val;
}floatUnion;

typedef struct{
    float val;
    int index;
}accData;

typedef struct{
    accData pack[8];
}accDataBundle;

typedef struct{
    int hapLen;
    float oneDivHapLen;
}HapLenPack;

typedef ap_uint<64> hapLenType;
typedef ap_uint<hap_width> hapType;
typedef ap_uint<64> hapTypeBundle;
typedef ap_uint<512> busType;
typedef ap_uint<64> readType;
typedef ap_uint<64> longlong;

typedef struct{
    ap_uint<512> reads[MAX_READ_LEN / 8];
}readBundle;

typedef struct{
    ap_uint<8> readLen0;
    ap_uint<8> readLen1;
    ap_uint<21> resultOffset;
    ap_uint<21> iterNum;
    bool oneOrTwo;
}readInfo;
//The first data of readData is readInfo, bit(7, 0) is readLen0, bit(15, 8) is readLen1, bit (36, 16) is resultOffset, bit(57, 37) is iterNum, bit 58 is oneOrTwo

template <class T>  T reg(T x) {
#pragma HLS pipeline
#pragma HLS inline self off
#pragma HLS interface ap_ctrl_none register port=return
    return x;
}
template<unsigned int PU_NUM>
void pmm_core(ap_uint<512> *input_data, int numRead, int numHap, float* output_data);
#endif
