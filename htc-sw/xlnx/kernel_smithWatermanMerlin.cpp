/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: Xilinx Smith-Waterman Vivado HLS implementation 
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/

#define __constant
#define __kernel
#define __global
#include <string.h> 
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "ap_int.h"

#define INIT_BASIC(offset, modVal) \
    a_base_local[offset] = ref[32 * i_j_gap_div32 - offset + modVal]; \
    sw_firstCol_cur[offset] = sw_firstCol[32 * i_j_gap_div32 - offset + modVal]; \
    sw_firstCol_nxt[offset] = sw_firstCol[32 * i_j_gap_div32 - offset + 1 + modVal]; 

#define INIT_CUSTOM(offset, modVal) \
    a_base_local[offset] = ref[32 * i_j_gap_div32 - offset + modVal]; \
    sw_firstCol_cur[offset] = sw_firstCol[32 * i_j_gap_div32 - offset + modVal]; \
    sw_firstCol_nxt[offset] = sw_firstCol[32 * i_j_gap_div32 - offset + 1 + modVal]; \


#define INIT_2(modVal) \
    INIT_CUSTOM(0, modVal) \
    INIT_CUSTOM(1, modVal) \
    INIT_CUSTOM(2, modVal) \
    INIT_CUSTOM(3, modVal) \
    INIT_CUSTOM(4, modVal) \
    INIT_CUSTOM(5, modVal) \
    INIT_CUSTOM(6, modVal) \
    INIT_CUSTOM(7, modVal) \
    INIT_CUSTOM(8, modVal) \
    INIT_CUSTOM(9, modVal) \
    INIT_CUSTOM(10, modVal) \
    INIT_CUSTOM(11, modVal) \
    INIT_CUSTOM(12, modVal) \
    INIT_CUSTOM(13, modVal) \
    INIT_CUSTOM(14, modVal) \
    INIT_CUSTOM(15, modVal) \
    INIT_CUSTOM(16, modVal) \
    INIT_CUSTOM(17, modVal) \
    INIT_CUSTOM(18, modVal) \
    INIT_CUSTOM(19, modVal) \
    INIT_CUSTOM(20, modVal) \
    INIT_CUSTOM(21, modVal) \
    INIT_CUSTOM(22, modVal) \
    INIT_CUSTOM(23, modVal) \
    INIT_CUSTOM(24, modVal) \
    INIT_CUSTOM(25, modVal) \
    INIT_CUSTOM(26, modVal) \
    INIT_CUSTOM(27, modVal) \
    INIT_CUSTOM(28, modVal) \
    INIT_CUSTOM(29, modVal) \
    INIT_CUSTOM(30, modVal) \
    INIT_CUSTOM(31, modVal) 

#define INIT_0(modVal) \
    INIT_BASIC(0, modVal) \
    INIT_BASIC(1, modVal) \
    INIT_BASIC(2, modVal) \
    INIT_BASIC(3, modVal) \
    INIT_BASIC(4, modVal) \
    INIT_BASIC(5, modVal) \
    INIT_BASIC(6, modVal) \
    INIT_BASIC(7, modVal) \
    INIT_BASIC(8, modVal) \
    INIT_BASIC(9, modVal) \
    INIT_BASIC(10, modVal) \
    INIT_BASIC(11, modVal) \
    INIT_BASIC(12, modVal) \
    INIT_BASIC(13, modVal) \
    INIT_BASIC(14, modVal) \
    INIT_BASIC(15, modVal) \
    INIT_BASIC(16, modVal) \
    INIT_BASIC(17, modVal) \
    INIT_BASIC(18, modVal) \
    INIT_BASIC(19, modVal) \
    INIT_BASIC(20, modVal) \
    INIT_BASIC(21, modVal) \
    INIT_BASIC(22, modVal) \
    INIT_BASIC(23, modVal) \
    INIT_BASIC(24, modVal) \
    INIT_BASIC(25, modVal) \
    INIT_BASIC(26, modVal) \
    INIT_BASIC(27, modVal) \
    INIT_BASIC(28, modVal) \
    INIT_BASIC(29, modVal) \
    INIT_BASIC(30, modVal) \
    INIT_BASIC(31, modVal) 

#define INIT_1(baseVal, modVal) \
    lastDiag_pre[baseVal] = lastDiag[64 * jdiv64 - 1 + modVal]; \
    lastDiag_cur[baseVal] = lastDiag[64 * jdiv64 + modVal]; \
    lastLastDiag_pre[baseVal] = lastLastDiag[64 * jdiv64 - 1 + modVal]; \
    lastLastDiag_cur[baseVal] = lastLastDiag[64 * jdiv64 + modVal]; \
    curDiag_pre[baseVal] = curDiag[64 * jdiv64 - 1 + modVal]; \
    curDiag_cur[baseVal] = curDiag[64 * jdiv64 + modVal];


#ifdef CUSTOM
static int swCellWisePipeline(int lastLastDiag_pre,int lastLastDiag_cur,int lastDiag_pre,int lastDiag_cur,int curDiag_pre,int curDiag_cur,int sw_firstCol_cur,int sw_firstCol_nxt,int sw_firstRow_cur,int sw_firstRow_nxt,ap_int<11> i_j_gap, ap_int<11> j,ap_uint<2> imod3,char a_base,char b_base,short w_match,short w_mismatch,short w_open,short w_extend, int MATRIX_MIN_CUTOFF, int best_gap_v_in,short gap_size_v_in,int best_gap_h_in, short gap_size_h_in, char imod2, int* best_gap_v_out, short* gap_size_v_out, int* best_gap_h_out, short* gap_size_h_out, short* btrack_out)
#else
static int swCellWisePipeline(int lastLastDiag_pre,int lastLastDiag_cur,int lastDiag_pre,int lastDiag_cur,int curDiag_pre,int curDiag_cur,int sw_firstCol_cur,int sw_firstCol_nxt,int sw_firstRow_cur,int sw_firstRow_nxt,int i_j_gap,int j,int imod3,char a_base,char b_base,short w_match,short w_mismatch,short w_open,short w_extend, int MATRIX_MIN_CUTOFF, int best_gap_v_in,short gap_size_v_in,int best_gap_h_in, short gap_size_h_in, char imod2, int* best_gap_v_out, short* gap_size_v_out, int* best_gap_h_out, short* gap_size_h_out, short* btrack_out)
#endif
{  
#pragma HLS inline
    int step_diag;
    int step_down;
    int step_right;
    short kd;
    short ki;
    int sw_lastRow_lastCol;
    int sw_lastRow_curCol;
    int sw_curRow_lastCol;
    int prev_gap_v;
    int prev_gap_h;
    bool diagHighestOrEqual;
    signed int lowInitValue = - 1073741824;
    if (j == 0) {
        sw_lastRow_lastCol = sw_firstCol_cur;
    }
    else if (i_j_gap == 0) {
        sw_lastRow_lastCol = sw_firstRow_cur;
    }   
    else {
#ifdef CUSTOM
        if (imod3 == (ap_uint<2>)(0)) {
#else
        if (imod3 == 0) {
#endif
            sw_lastRow_lastCol = lastLastDiag_pre;
        }
#ifdef CUSTOM
        else if (imod3 == (ap_uint<2>)(1)) {
#else
        else if (imod3 == 1) {
#endif
            sw_lastRow_lastCol = lastDiag_pre;
        }
        else {
            sw_lastRow_lastCol = curDiag_pre;
        }
    }
    if (j == 0) {
        sw_curRow_lastCol = sw_firstCol_nxt;
    }
    else {
#ifdef CUSTOM
        if (imod3 == (ap_uint<2>)(0)) {
#else
        if (imod3 == 0) {
#endif
            sw_curRow_lastCol = lastDiag_pre;
        }
#ifdef CUSTOM
        else if (imod3 == (ap_uint<2>)(1)) {
#else
        else if (imod3 == 1) {
#endif
            sw_curRow_lastCol = curDiag_pre;
        }
        else {
            sw_curRow_lastCol = lastLastDiag_pre;
        }
    }
    if (i_j_gap == 0) {
        sw_lastRow_curCol = sw_firstRow_nxt;
    }
    else {
#ifdef CUSTOM
        if (imod3 == (ap_uint<2>)(0)) {
#else
        if (imod3 == 0) {
#endif
            sw_lastRow_curCol = lastDiag_cur;
        }
#ifdef CUSTOM
        else if (imod3 == (ap_uint<2>)(1)) {
#else
        else if (imod3 == 1) {
#endif
            sw_lastRow_curCol = curDiag_cur;
        }
        else {
            sw_lastRow_curCol = lastLastDiag_cur;
        }
    }
    int rose_temp;
    if (a_base == b_base) {
        rose_temp = w_match;
    }
    else {
        rose_temp = w_mismatch;
    }
    step_diag = sw_lastRow_lastCol + rose_temp;
    int cur_best_gap_v;
    short cur_gap_size_v;
    prev_gap_v = sw_lastRow_curCol + w_open;
    cur_best_gap_v = best_gap_v_in + w_extend;
    int cur_best_gap_v_1;
    if (prev_gap_v > cur_best_gap_v) {
        cur_best_gap_v_1 = prev_gap_v;
        cur_gap_size_v = 1;
    }
    else {
        cur_best_gap_v_1 = cur_best_gap_v;
        cur_gap_size_v = (short)gap_size_v_in + 1;
    }
    {
#pragma HLS latency min=1
        *best_gap_v_out = cur_best_gap_v_1;
        *gap_size_v_out = cur_gap_size_v;
    }
    {
#pragma HLS latency min=1
        step_down = cur_best_gap_v_1;
        kd = cur_gap_size_v;
    }
    int cur_best_gap_h;
    short cur_gap_size_h;
    prev_gap_h = sw_curRow_lastCol + w_open;
    if((bool)j)
        cur_best_gap_h = best_gap_h_in + w_extend;
    else
        cur_best_gap_h = - 1073741824 + w_extend;

    int cur_best_gap_h_1;
    if (prev_gap_h > cur_best_gap_h) {
        cur_best_gap_h_1 = prev_gap_h;
        cur_gap_size_h = 1;
    }
    else {
        cur_best_gap_h_1 = cur_best_gap_h;
        if ((bool)j) {
            cur_gap_size_h = (short)gap_size_h_in + 1;
        }
        else {
            cur_gap_size_h = 1;
        }
    }
    {
#pragma HLS latency min=1
        *best_gap_h_out = cur_best_gap_h_1;
        *gap_size_h_out = cur_gap_size_h;
    }
    {
#pragma HLS latency min=1
        step_right = cur_best_gap_h_1;
        ki = cur_gap_size_h;
    }
    diagHighestOrEqual = ((step_diag >= step_down && step_diag >= step_right));
    int curDiagVal;
    short curBtrack;
    if (diagHighestOrEqual) {
        if (MATRIX_MIN_CUTOFF > step_diag) {
            curDiagVal = MATRIX_MIN_CUTOFF;
        }
        else {
            curDiagVal = step_diag;
        }
        curBtrack = 0;
    }
    else {
        if (step_right >= step_down) {
            if (MATRIX_MIN_CUTOFF > step_right) {
                curDiagVal = MATRIX_MIN_CUTOFF;
            }
            else {
                curDiagVal = step_right;
            }
            curBtrack = -ki;
        }
        else {
            if (MATRIX_MIN_CUTOFF > step_down) {
                curDiagVal = MATRIX_MIN_CUTOFF;
            }
            else {
                curDiagVal = step_down;
            }
            curBtrack = kd;
        }
    }
    {
#pragma HLS latency min=1
    *btrack_out = (short)curBtrack;
    }
    return curDiagVal;
}

static void merlin_memcpy_00(char* dst,int dst_idx_0,char *src,int src_idx_0,unsigned int len,unsigned int max_len)
{
#pragma HLS inline off
#pragma HLS function_instantiate variable=dst_idx_0,src_idx_0
    long long i;
    long long total_offset1 = dst_idx_0;
    long long total_offset2 = src_idx_0;
    for (i = 0; i < len; ++i) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT max=512
        dst[total_offset1 + i] = src[total_offset2 + i];
    }
}

static int smithWatermanKernel(char *ref_global,char *alt_global,short refLength,short altLength,ap_uint<2> overhang_strategy,int cutoff,short w_match,short w_mismatch,short w_open,short w_extend,short *lengths,short *states,short *cigarElementNum,short *alignment_offset_ptr)
{
#pragma HLS inline off
#pragma HLS function_instantiate variable=altLength

    signed int lowInitValue = - 1073741824;
    int MATRIX_MIN_CUTOFF;
    if ((bool )cutoff) {
        MATRIX_MIN_CUTOFF = 0;
    }
    else {
        MATRIX_MIN_CUTOFF = ((int )(- 1e8));
    }
    short ncol = altLength + 1;
    short nrow = refLength + 1;
    nrow > 0 && nrow < 512?((void )0) : __assert_fail("nrow > 0 && nrow < 512","com_falconcomputing_genomics_haplotypecaller_FalconSWPairwiseAlignmentHLS.cpp",((unsigned int )357),__PRETTY_FUNCTION__);
    ncol > 0 && ncol < 512?((void )0) : __assert_fail("ncol > 0 && ncol < 512","com_falconcomputing_genomics_haplotypecaller_FalconSWPairwiseAlignmentHLS.cpp",((unsigned int )358),__PRETTY_FUNCTION__);
    char ref[512];
#pragma HLS array_partition variable=ref cyclic factor=32 dim=1
    char alt[512];
#pragma HLS array_partition variable=alt cyclic factor=32 dim=1
    short btrack[512 * 512];
#pragma HLS array_partition variable=btrack cyclic factor=32 dim=1
    
#ifdef CUSTOM
    ap_uint<11> i = 0;
    ap_uint<11> j = 0;
#else
    short i = 0;
    short j = 0;
#endif

    short k = 0;
    int lastLastDiag[512];
#pragma HLS array_partition variable=lastLastDiag cyclic factor=64 dim=1
    int lastDiag[512];
#pragma HLS array_partition variable=lastDiag cyclic factor=64 dim=1
    int curDiag[512];
#pragma HLS array_partition variable=curDiag cyclic factor=64 dim=1
    int best_gap_v[512];
#pragma HLS array_partition variable=best_gap_v cyclic factor=32 dim=1
    int best_gap_h_last[512];
#pragma HLS array_partition variable=best_gap_h_last cyclic factor=32 dim=1
    int best_gap_h_curr[512];
#pragma HLS array_partition variable=best_gap_h_curr cyclic factor=32 dim=1
    short gap_size_v[512];
#pragma HLS array_partition variable=gap_size_v cyclic factor=32 dim=1
    short gap_size_h_last[512];
#pragma HLS array_partition variable=gap_size_h_last cyclic factor=32 dim=1
    short gap_size_h_curr[512];
#pragma HLS array_partition variable=gap_size_h_curr cyclic factor=32 dim=1
    int sw_firstRow[512];
#pragma HLS array_partition variable=sw_firstRow cyclic factor=32 dim=1
    int sw_firstCol[512];
#pragma HLS array_partition variable=sw_firstCol cyclic factor=32 dim=1
    int sw_finalRow[512];
    int sw_finalCol[512];
    short sw_finalCol_index = 0;

    merlin_memcpy_00(ref,0,ref_global,0,((unsigned long )512) * 1UL,512UL);
    merlin_memcpy_00(alt,0,alt_global,0,((unsigned long )512) * 1UL,512UL);
    int i_sub_3;
   
    for (i = 0; i < 16; ++i) {
#pragma HLS pipeline
        for (i_sub_3 = 0; i_sub_3 < 32; ++i_sub_3) {
#pragma HLS unroll
            best_gap_v[i * 32 + i_sub_3] = - 1073741824;
            gap_size_v[i * 32 + i_sub_3] = 0;
            best_gap_h_last[i * 32 + i_sub_3] = - 1073741824;
            gap_size_h_last[i * 32 + i_sub_3] = 0;
            best_gap_h_curr[i * 32 + i_sub_3] = - 1073741824;
            gap_size_h_curr[i * 32 + i_sub_3] = 0;
            sw_firstRow[i * 32 + i_sub_3] = 0;
            sw_firstCol[i * 32 + i_sub_3] = 0;
        }
    }
    
    for(i = 0; i < ncol; ++i){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=512
        sw_finalRow[i] = 0;
    }

    for(i = 0; i < nrow; ++i){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=512
        sw_finalCol[i] = 0;
    }
  
    if (overhang_strategy == 1 || overhang_strategy == 2) {
        sw_firstRow[1] = w_open;
        int currentValue = w_open;
        for (i = 2; i < 512; i++) {
#pragma HLS pipeline
            if (i < ncol) {
                currentValue += w_extend;
                sw_firstRow[i] = currentValue;
            }
        }
        sw_finalCol[0] = currentValue;
        sw_firstCol[1] = w_open;
        currentValue = w_open;
    
        for (i = 2; i < 512; i++) {
#pragma HLS pipeline
            if (i < nrow) {
                currentValue += w_extend;
                sw_firstCol[i] = currentValue;
            }
        }
        sw_finalRow[0] = currentValue;
    }
    
    i = 0;
    j = 0;
#ifdef CUSTOM
    ap_uint<2> imod3 = 0;
    ap_uint<5> j_PE = 0;
#else
    int imod3 = 0;
    char j_PE = 0;
#endif
    for (k = 0; k < (2 * 512 - 1) * 512 / 32; k++) {
    
#pragma HLS pipeline II=1
#pragma HLS dependence variable=lastLastDiag inter false
#pragma HLS dependence variable=lastDiag inter false
#pragma HLS dependence variable=curDiag inter false
#pragma HLS dependence variable=best_gap_v inter false
#pragma HLS dependence variable=gap_size_v inter false
#pragma HLS dependence variable=best_gap_h_last inter false
#pragma HLS dependence variable=best_gap_h_curr inter false
#pragma HLS dependence variable=gap_size_h_last inter false
#pragma HLS dependence variable=gap_size_h_curr inter false
#pragma HLS dependence variable=best_gap_h_last intra false
#pragma HLS dependence variable=best_gap_h_curr intra false
#pragma HLS dependence variable=gap_size_h_last intra false
#pragma HLS dependence variable=gap_size_h_curr intra false
#pragma HLS dependence variable=sw_finalCol inter false
#pragma HLS dependence variable=sw_finalRow inter false
#pragma HLS dependence variable=btrack inter false

#ifdef CUSTOM
        j = (((ap_uint<11>)j_PE) << 5);
#else
        j = ((short)j_PE << 5);
#endif
       
        char imod2 = (char )(i % 2);
//working on j + 1, i - j + 1 element in matrix
//b_base = alt[j];
//a_base = ref[i_j_gap];
        int lastDiag_pre[32];
#pragma HLS array_partition variable=lastDiag_pre complete dim=1
        int lastDiag_cur[32];
#pragma HLS array_partition variable=lastDiag_cur complete dim=1
        int lastLastDiag_pre[32];
#pragma HLS array_partition variable=lastLastDiag_pre complete dim=1
        int lastLastDiag_cur[32];
#pragma HLS array_partition variable=lastLastDiag_cur complete dim=1
        int curDiag_pre[32];
#pragma HLS array_partition variable=curDiag_pre complete dim=1
        int curDiag_cur[32];
#pragma HLS array_partition variable=curDiag_cur complete dim=1
        int sw_firstCol_cur[32];
#pragma HLS array_partition variable=sw_firstCol_cur complete dim=1
        int sw_firstCol_nxt[32];
#pragma HLS array_partition variable=sw_firstCol_nxt complete dim=1
        int sw_firstRow_cur[32];
#pragma HLS array_partition variable=sw_firstRow_cur complete dim=1
        int sw_firstRow_nxt[32];
#pragma HLS array_partition variable=sw_firstRow_nxt complete dim=1
        int curDiagVal[32];
#pragma HLS array_partition variable=curDiagVal complete dim=1
#ifdef CUSTOM        
        ap_int<11> i_j_gap_local[32];
#else
        int i_j_gap_local[32];
#endif
#pragma HLS array_partition variable=i_j_gap_local complete dim=1
#ifdef CUSTOM        
        ap_int<11> j_local[32];
#else
        int j_local[32];
#endif
#pragma HLS array_partition variable=j_local complete dim=1
        char a_base_local[32];
#pragma HLS array_partition variable=a_base_local complete dim=1
        char b_base_local[32];
#pragma HLS array_partition variable=b_base_local complete dim=1
        short w_match_local[32];
#pragma HLS array_partition variable=w_match_local complete dim=1
        short w_mismatch_local[32];
#pragma HLS array_partition variable=w_mismatch_local complete dim=1
        short w_open_local[32];
#pragma HLS array_partition variable=w_open_local complete dim=1
        short w_extend_local[32];
#pragma HLS array_partition variable=w_extend_local complete dim=1

        //char i_j_gap_div32 = ((i - j) >> 5) & 0xf;
        char i_j_gap_div32 = ((i - j) >> 5);
#ifdef CUSTOM
        ap_uint<5> i_j_gap_mod32 = (ap_uint<5>)((i - j) & 0x1f);
#else
        int i_j_gap_mod32 = ((i - j) & 0x1f);
#endif

        for (int index = 0; index < 32; index++) {
#pragma HLS unroll
#pragma HLS latency min=1
            i_j_gap_local[index] = i - j - index;
            j_local[index] = j + index;
            b_base_local[index] = alt[j + index];
            sw_firstRow_cur[index] = sw_firstRow[j + index];
            sw_firstRow_nxt[index] = sw_firstRow[j + index + 1];
            w_match_local[index] = w_match;
            w_mismatch_local[index] = w_mismatch;
            w_open_local[index] = w_open;
            w_extend_local[index] = w_extend;
        }

       #ifdef CUSTOM
        if(i_j_gap_mod32 == (ap_uint<5>)(0)){
            INIT_2(0)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(1)){
            INIT_2(1)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(2)){
            INIT_2(2)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(3)){
            INIT_2(3)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(4)){
            INIT_2(4)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(5)){
            INIT_2(5)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(6)){
            INIT_2(6)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(7)){
            INIT_2(7)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(8)){
            INIT_2(8)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(9)){
            INIT_2(9)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(10)){
            INIT_2(10)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(11)){
            INIT_2(11)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(12)){
            INIT_2(12)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(13)){
            INIT_2(13)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(14)){
            INIT_2(14)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(15)){
            INIT_2(15)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(16)){
            INIT_2(16)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(17)){
            INIT_2(17)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(18)){
            INIT_2(18)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(19)){
            INIT_2(19)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(20)){
            INIT_2(20)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(21)){
            INIT_2(21)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(22)){
            INIT_2(22)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(23)){
            INIT_2(23)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(24)){
            INIT_2(24)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(25)){
            INIT_2(25)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(26)){
            INIT_2(26)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(27)){
            INIT_2(27)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(28)){
            INIT_2(28)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(29)){
            INIT_2(29)
        }
        else if(i_j_gap_mod32 == (ap_uint<5>)(30)){
            INIT_2(30)
        }
        else{
            INIT_2(31)
        }
#else
        if(i_j_gap_mod32 == (0)){
            INIT_0(0)
        }
        else if(i_j_gap_mod32 == (1)){
            INIT_0(1)
        }
        else if(i_j_gap_mod32 == (2)){
            INIT_0(2)
        }
        else if(i_j_gap_mod32 == (3)){
            INIT_0(3)
        }
        else if(i_j_gap_mod32 == (4)){
            INIT_0(4)
        }
        else if(i_j_gap_mod32 == (5)){
            INIT_0(5)
        }
        else if(i_j_gap_mod32 == (6)){
            INIT_0(6)
        }
        else if(i_j_gap_mod32 == (7)){
            INIT_0(7)
        }
        else if(i_j_gap_mod32 == (8)){
            INIT_0(8)
        }
        else if(i_j_gap_mod32 == (9)){
            INIT_0(9)
        }
        else if(i_j_gap_mod32 == (10)){
            INIT_0(10)
        }
        else if(i_j_gap_mod32 == (11)){
            INIT_0(11)
        }
        else if(i_j_gap_mod32 == (12)){
            INIT_0(12)
        }
        else if(i_j_gap_mod32 == (13)){
            INIT_0(13)
        }
        else if(i_j_gap_mod32 == (14)){
            INIT_0(14)
        }
        else if(i_j_gap_mod32 == (15)){
            INIT_0(15)
        }
        else if(i_j_gap_mod32 == (16)){
            INIT_0(16)
        }
        else if(i_j_gap_mod32 == (17)){
            INIT_0(17)
        }
        else if(i_j_gap_mod32 == (18)){
            INIT_0(18)
        }
        else if(i_j_gap_mod32 == (19)){
            INIT_0(19)
        }
        else if(i_j_gap_mod32 == (20)){
            INIT_0(20)
        }
        else if(i_j_gap_mod32 == (21)){
            INIT_0(21)
        }
        else if(i_j_gap_mod32 == (22)){
            INIT_0(22)
        }
        else if(i_j_gap_mod32 == (23)){
            INIT_0(23)
        }
        else if(i_j_gap_mod32 == (24)){
            INIT_0(24)
                }
        else if(i_j_gap_mod32 == (25)){
            INIT_0(25)
        }
        else if(i_j_gap_mod32 == (26)){
            INIT_0(26)
        }
        else if(i_j_gap_mod32 == (27)){
            INIT_0(27)
        }
        else if(i_j_gap_mod32 == (28)){
            INIT_0(28)
        }
        else if(i_j_gap_mod32 == (29)){
            INIT_0(29)
        }
        else if(i_j_gap_mod32 == (30)){
            INIT_0(30)
        }
        else{
            INIT_0(31)
        }
#endif

        bool j_PEmod2;
        {
#pragma HLS latency min=1
            j_PEmod2 = (bool)(j_PE & 0x1);
        }
        
        //char jdiv64 = (j >> 6) & 0x7;
        char jdiv64 = (j >> 6);
        if (!j_PEmod2) {
            /*if (j != 0) {
                lastDiag_pre[0] = lastDiag[64 * jdiv64 - 1];
            }
            lastDiag_cur[0] = lastDiag[64 * jdiv64];
            if (j != 0) {
                lastLastDiag_pre[0] = lastLastDiag[64 * jdiv64 - 1];
            }
            lastLastDiag_cur[0] = lastLastDiag[64 * jdiv64];
            if (j != 0) {
                curDiag_pre[0] = curDiag[64 * jdiv64 - 1];
            }
            curDiag_cur[0] = curDiag[64 * jdiv64];*/
            
            lastDiag_pre[0] = lastDiag[64 * jdiv64 - 1];
            lastDiag_cur[0] = lastDiag[64 * jdiv64];
            lastLastDiag_pre[0] = lastLastDiag[64 * jdiv64 - 1];
            lastLastDiag_cur[0] = lastLastDiag[64 * jdiv64];
            curDiag_pre[0] = curDiag[64 * jdiv64 - 1];
            curDiag_cur[0] = curDiag[64 * jdiv64];


        
            INIT_1(1, 1)
            INIT_1(2, 2)
            INIT_1(3, 3)
            INIT_1(4, 4)
            INIT_1(5, 5)
            INIT_1(6, 6)
            INIT_1(7, 7)
            INIT_1(8, 8)
            INIT_1(9, 9)
            INIT_1(10, 10)
            INIT_1(11, 11)
            INIT_1(12, 12)
            INIT_1(13, 13)
            INIT_1(14, 14)
            INIT_1(15, 15)
            INIT_1(16, 16)
            INIT_1(17, 17)
            INIT_1(18, 18)
            INIT_1(19, 19)
            INIT_1(20, 20)
            INIT_1(21, 21)
            INIT_1(22, 22)
            INIT_1(23, 23)
            INIT_1(24, 24)
            INIT_1(25, 25)
            INIT_1(26, 26)
            INIT_1(27, 27)
            INIT_1(28, 28)
            INIT_1(29, 29)
            INIT_1(30, 30)
            INIT_1(31, 31)

            for (int index = 0; index < 32; index++) {
#pragma HLS unroll
                if (i < refLength + altLength - 1 && i >= j_local[index] && i <= j_local[index] + refLength - 1 && j_local[index] < altLength) {
                    int best_gap_v_in = best_gap_v[j_local[index]];
                    short gap_size_v_in = gap_size_v[j_local[index]];
                    int best_gap_h_in;
                    short gap_size_h_in;
                    int best_gap_v_out;
                    short gap_size_v_out;
                    int best_gap_h_out;
                    short gap_size_h_out;
                    short btrack_out;
                    if((bool)imod2){
                        best_gap_h_in = best_gap_h_curr[j_local[index] - 1];
                        gap_size_h_in = gap_size_h_curr[j_local[index] - 1];
                    }
                    else{
                        best_gap_h_in = best_gap_h_last[j_local[index] - 1];
                        gap_size_h_in = gap_size_h_last[j_local[index] - 1];
                    }

                    curDiagVal[index] = swCellWisePipeline(lastLastDiag_pre[index],lastLastDiag_cur[index],lastDiag_pre[index],lastDiag_cur[index],curDiag_pre[index],curDiag_cur[index],sw_firstCol_cur[index],sw_firstCol_nxt[index],sw_firstRow_cur[index],sw_firstRow_nxt[index], i_j_gap_local[index], j_local[index], imod3, a_base_local[index],b_base_local[index],w_match_local[index],w_mismatch_local[index],w_open_local[index],w_extend_local[index], MATRIX_MIN_CUTOFF, best_gap_v_in, gap_size_v_in, best_gap_h_in, gap_size_h_in, imod2, &best_gap_v_out, &gap_size_v_out, &best_gap_h_out, &gap_size_h_out, &btrack_out);

                    btrack[i_j_gap_local[index] * 512 + j_local[index]] = btrack_out;
                    if(imod2){
                        best_gap_h_last[j_local[index]] = best_gap_h_out;
                        gap_size_h_last[j_local[index]] = gap_size_h_out;
                    }
                    else{
                        best_gap_h_curr[j_local[index]] = best_gap_h_out;
                        gap_size_h_curr[j_local[index]] = gap_size_h_out;
                    }
                    best_gap_v[j_local[index]] = best_gap_v_out;
                    gap_size_v[j_local[index]] = gap_size_v_out;
                    
                }
            }
#ifdef CUSTOM
            if (imod3 == (ap_uint<2>)(0)) {
#else
            if (imod3 == 0) {
#endif
                curDiag[64 * jdiv64] = curDiagVal[0];
                curDiag[64 * jdiv64 + 1] = curDiagVal[1];
                curDiag[64 * jdiv64 + 2] = curDiagVal[2];
                curDiag[64 * jdiv64 + 3] = curDiagVal[3];
                curDiag[64 * jdiv64 + 4] = curDiagVal[4];
                curDiag[64 * jdiv64 + 5] = curDiagVal[5];
                curDiag[64 * jdiv64 + 6] = curDiagVal[6];
                curDiag[64 * jdiv64 + 7] = curDiagVal[7];
                curDiag[64 * jdiv64 + 8] = curDiagVal[8];
                curDiag[64 * jdiv64 + 9] = curDiagVal[9];
                curDiag[64 * jdiv64 + 10] = curDiagVal[10];
                curDiag[64 * jdiv64 + 11] = curDiagVal[11];
                curDiag[64 * jdiv64 + 12] = curDiagVal[12];
                curDiag[64 * jdiv64 + 13] = curDiagVal[13];
                curDiag[64 * jdiv64 + 14] = curDiagVal[14];
                curDiag[64 * jdiv64 + 15] = curDiagVal[15];
                curDiag[64 * jdiv64 + 16] = curDiagVal[16];
                curDiag[64 * jdiv64 + 17] = curDiagVal[17];
                curDiag[64 * jdiv64 + 18] = curDiagVal[18];
                curDiag[64 * jdiv64 + 19] = curDiagVal[19];
                curDiag[64 * jdiv64 + 20] = curDiagVal[20];
                curDiag[64 * jdiv64 + 21] = curDiagVal[21];
                curDiag[64 * jdiv64 + 22] = curDiagVal[22];
                curDiag[64 * jdiv64 + 23] = curDiagVal[23];
                curDiag[64 * jdiv64 + 24] = curDiagVal[24];
                curDiag[64 * jdiv64 + 25] = curDiagVal[25];
                curDiag[64 * jdiv64 + 26] = curDiagVal[26];
                curDiag[64 * jdiv64 + 27] = curDiagVal[27];
                curDiag[64 * jdiv64 + 28] = curDiagVal[28];
                curDiag[64 * jdiv64 + 29] = curDiagVal[29];
                curDiag[64 * jdiv64 + 30] = curDiagVal[30];
                curDiag[64 * jdiv64 + 31] = curDiagVal[31];
            }
#ifdef CUSTOM
            else if (imod3 == (ap_uint<2>)(1)) {
#else
            else if (imod3 == 1) {
#endif
                lastLastDiag[64 * jdiv64] = curDiagVal[0];
                lastLastDiag[64 * jdiv64 + 1] = curDiagVal[1];
                lastLastDiag[64 * jdiv64 + 2] = curDiagVal[2];
                lastLastDiag[64 * jdiv64 + 3] = curDiagVal[3];
                lastLastDiag[64 * jdiv64 + 4] = curDiagVal[4];
                lastLastDiag[64 * jdiv64 + 5] = curDiagVal[5];
                lastLastDiag[64 * jdiv64 + 6] = curDiagVal[6];
                lastLastDiag[64 * jdiv64 + 7] = curDiagVal[7];
                lastLastDiag[64 * jdiv64 + 8] = curDiagVal[8];
                lastLastDiag[64 * jdiv64 + 9] = curDiagVal[9];
                lastLastDiag[64 * jdiv64 + 10] = curDiagVal[10];
                lastLastDiag[64 * jdiv64 + 11] = curDiagVal[11];
                lastLastDiag[64 * jdiv64 + 12] = curDiagVal[12];
                lastLastDiag[64 * jdiv64 + 13] = curDiagVal[13];
                lastLastDiag[64 * jdiv64 + 14] = curDiagVal[14];
                lastLastDiag[64 * jdiv64 + 15] = curDiagVal[15];
                lastLastDiag[64 * jdiv64 + 16] = curDiagVal[16];
                lastLastDiag[64 * jdiv64 + 17] = curDiagVal[17];
                lastLastDiag[64 * jdiv64 + 18] = curDiagVal[18];
                lastLastDiag[64 * jdiv64 + 19] = curDiagVal[19];
                lastLastDiag[64 * jdiv64 + 20] = curDiagVal[20];
                lastLastDiag[64 * jdiv64 + 21] = curDiagVal[21];
                lastLastDiag[64 * jdiv64 + 22] = curDiagVal[22];
                lastLastDiag[64 * jdiv64 + 23] = curDiagVal[23];
                lastLastDiag[64 * jdiv64 + 24] = curDiagVal[24];
                lastLastDiag[64 * jdiv64 + 25] = curDiagVal[25];
                lastLastDiag[64 * jdiv64 + 26] = curDiagVal[26];
                lastLastDiag[64 * jdiv64 + 27] = curDiagVal[27];
                lastLastDiag[64 * jdiv64 + 28] = curDiagVal[28];
                lastLastDiag[64 * jdiv64 + 29] = curDiagVal[29];
                lastLastDiag[64 * jdiv64 + 30] = curDiagVal[30];
                lastLastDiag[64 * jdiv64 + 31] = curDiagVal[31];
            }
            else {
                lastDiag[64 * jdiv64] = curDiagVal[0];
                lastDiag[64 * jdiv64 + 1] = curDiagVal[1];
                lastDiag[64 * jdiv64 + 2] = curDiagVal[2];
                lastDiag[64 * jdiv64 + 3] = curDiagVal[3];
                lastDiag[64 * jdiv64 + 4] = curDiagVal[4];
                lastDiag[64 * jdiv64 + 5] = curDiagVal[5];
                lastDiag[64 * jdiv64 + 6] = curDiagVal[6];
                lastDiag[64 * jdiv64 + 7] = curDiagVal[7];
                lastDiag[64 * jdiv64 + 8] = curDiagVal[8];
                lastDiag[64 * jdiv64 + 9] = curDiagVal[9];
                lastDiag[64 * jdiv64 + 10] = curDiagVal[10];
                lastDiag[64 * jdiv64 + 11] = curDiagVal[11];
                lastDiag[64 * jdiv64 + 12] = curDiagVal[12];
                lastDiag[64 * jdiv64 + 13] = curDiagVal[13];
                lastDiag[64 * jdiv64 + 14] = curDiagVal[14];
                lastDiag[64 * jdiv64 + 15] = curDiagVal[15];
                lastDiag[64 * jdiv64 + 16] = curDiagVal[16];
                lastDiag[64 * jdiv64 + 17] = curDiagVal[17];
                lastDiag[64 * jdiv64 + 18] = curDiagVal[18];
                lastDiag[64 * jdiv64 + 19] = curDiagVal[19];
                lastDiag[64 * jdiv64 + 20] = curDiagVal[20];
                lastDiag[64 * jdiv64 + 21] = curDiagVal[21];
                lastDiag[64 * jdiv64 + 22] = curDiagVal[22];
                lastDiag[64 * jdiv64 + 23] = curDiagVal[23];
                lastDiag[64 * jdiv64 + 24] = curDiagVal[24];
                lastDiag[64 * jdiv64 + 25] = curDiagVal[25];
                lastDiag[64 * jdiv64 + 26] = curDiagVal[26];
                lastDiag[64 * jdiv64 + 27] = curDiagVal[27];
                lastDiag[64 * jdiv64 + 28] = curDiagVal[28];
                lastDiag[64 * jdiv64 + 29] = curDiagVal[29];
                lastDiag[64 * jdiv64 + 30] = curDiagVal[30];
                lastDiag[64 * jdiv64 + 31] = curDiagVal[31];
            }
        }
        else {
            /*if (j != 0) {
                lastDiag_pre[0] = lastDiag[64 * jdiv64 + 31];
            }
            lastDiag_cur[0] = lastDiag[64 * jdiv64 + 32];
            if (j != 0) {
                lastLastDiag_pre[0] = lastLastDiag[64 * jdiv64 + 31];
            }
            lastLastDiag_cur[0] = lastLastDiag[64 * jdiv64 + 32];
            if (j != 0) {
                curDiag_pre[0] = curDiag[64 * jdiv64 + 31];
            }
            curDiag_cur[0] = curDiag[64 * jdiv64 + 32];*/

            lastDiag_pre[0] = lastDiag[64 * jdiv64 + 31];
            lastDiag_cur[0] = lastDiag[64 * jdiv64 + 32];
            lastLastDiag_pre[0] = lastLastDiag[64 * jdiv64 + 31];
            lastLastDiag_cur[0] = lastLastDiag[64 * jdiv64 + 32];
            curDiag_pre[0] = curDiag[64 * jdiv64 + 31];
            curDiag_cur[0] = curDiag[64 * jdiv64 + 32];
        

            INIT_1(1, 33)
            INIT_1(2, 34)
            INIT_1(3, 35)
            INIT_1(4, 36)
            INIT_1(5, 37)
            INIT_1(6, 38)
            INIT_1(7, 39)
            INIT_1(8, 40)
            INIT_1(9, 41)
            INIT_1(10, 42)
            INIT_1(11, 43)
            INIT_1(12, 44)
            INIT_1(13, 45)
            INIT_1(14, 46)
            INIT_1(15, 47)
            INIT_1(16, 48)
            INIT_1(17, 49)
            INIT_1(18, 50)
            INIT_1(19, 51)
            INIT_1(20, 52)
            INIT_1(21, 53)
            INIT_1(22, 54)
            INIT_1(23, 55)
            INIT_1(24, 56)
            INIT_1(25, 57)
            INIT_1(26, 58)
            INIT_1(27, 59)
            INIT_1(28, 60)
            INIT_1(29, 61)
            INIT_1(30, 62)
            INIT_1(31, 63)

            for (int index = 0; index < 32; index++) {
#pragma HLS unroll
                if (i < refLength + altLength - 1 && i >= j_local[index] && i <= j_local[index] + refLength - 1 && j_local[index] < altLength) {
                    int best_gap_v_in = best_gap_v[j_local[index]];
                    short gap_size_v_in = gap_size_v[j_local[index]];
                    int best_gap_h_in;
                    short gap_size_h_in;
                    int best_gap_v_out;
                    short gap_size_v_out;
                    int best_gap_h_out;
                    short gap_size_h_out;
                    short btrack_out;
                    if((bool)imod2){
                        best_gap_h_in = best_gap_h_curr[j_local[index] - 1];
                        gap_size_h_in = gap_size_h_curr[j_local[index] - 1];
                    }
                    else{
                        best_gap_h_in = best_gap_h_last[j_local[index] - 1];
                        gap_size_h_in = gap_size_h_last[j_local[index] - 1];
                    }

                    curDiagVal[index] = swCellWisePipeline(lastLastDiag_pre[index],lastLastDiag_cur[index],lastDiag_pre[index],lastDiag_cur[index],curDiag_pre[index],curDiag_cur[index],sw_firstCol_cur[index],sw_firstCol_nxt[index],sw_firstRow_cur[index],sw_firstRow_nxt[index], i_j_gap_local[index], j_local[index], imod3, a_base_local[index],b_base_local[index],w_match_local[index],w_mismatch_local[index],w_open_local[index],w_extend_local[index], MATRIX_MIN_CUTOFF, best_gap_v_in,gap_size_v_in,best_gap_h_in, gap_size_h_in, imod2, &best_gap_v_out, &gap_size_v_out, &best_gap_h_out, &gap_size_h_out, &btrack_out);

                    btrack[i_j_gap_local[index] * 512 + j_local[index]] = btrack_out;
                    if(imod2){
                        best_gap_h_last[j_local[index]] = best_gap_h_out;
                        gap_size_h_last[j_local[index]] = gap_size_h_out;
                    }
                    else{
                        best_gap_h_curr[j_local[index]] = best_gap_h_out;
                        gap_size_h_curr[j_local[index]] = gap_size_h_out;
                    }
                    best_gap_v[j_local[index]] = best_gap_v_out;
                    gap_size_v[j_local[index]] = gap_size_v_out;
                    
                }
            }
#ifdef CUSTOM
            if (imod3 == (ap_uint<2>)(0)) {
#else
            if (imod3 == 0) {
#endif
                curDiag[64 * jdiv64 + 32] = curDiagVal[0];
                curDiag[64 * jdiv64 + 33] = curDiagVal[1];
                curDiag[64 * jdiv64 + 34] = curDiagVal[2];
                curDiag[64 * jdiv64 + 35] = curDiagVal[3];
                curDiag[64 * jdiv64 + 36] = curDiagVal[4];
                curDiag[64 * jdiv64 + 37] = curDiagVal[5];
                curDiag[64 * jdiv64 + 38] = curDiagVal[6];
                curDiag[64 * jdiv64 + 39] = curDiagVal[7];
                curDiag[64 * jdiv64 + 40] = curDiagVal[8];
                curDiag[64 * jdiv64 + 41] = curDiagVal[9];
                curDiag[64 * jdiv64 + 42] = curDiagVal[10];
                curDiag[64 * jdiv64 + 43] = curDiagVal[11];
                curDiag[64 * jdiv64 + 44] = curDiagVal[12];
                curDiag[64 * jdiv64 + 45] = curDiagVal[13];
                curDiag[64 * jdiv64 + 46] = curDiagVal[14];
                curDiag[64 * jdiv64 + 47] = curDiagVal[15];
                curDiag[64 * jdiv64 + 48] = curDiagVal[16];
                curDiag[64 * jdiv64 + 49] = curDiagVal[17];
                curDiag[64 * jdiv64 + 50] = curDiagVal[18];
                curDiag[64 * jdiv64 + 51] = curDiagVal[19];
                curDiag[64 * jdiv64 + 52] = curDiagVal[20];
                curDiag[64 * jdiv64 + 53] = curDiagVal[21];
                curDiag[64 * jdiv64 + 54] = curDiagVal[22];
                curDiag[64 * jdiv64 + 55] = curDiagVal[23];
                curDiag[64 * jdiv64 + 56] = curDiagVal[24];
                curDiag[64 * jdiv64 + 57] = curDiagVal[25];
                curDiag[64 * jdiv64 + 58] = curDiagVal[26];
                curDiag[64 * jdiv64 + 59] = curDiagVal[27];
                curDiag[64 * jdiv64 + 60] = curDiagVal[28];
                curDiag[64 * jdiv64 + 61] = curDiagVal[29];
                curDiag[64 * jdiv64 + 62] = curDiagVal[30];
                curDiag[64 * jdiv64 + 63] = curDiagVal[31];
             }
#ifdef CUSTOM
            else if (imod3 == (ap_uint<2>)(1)) {
#else
            else if (imod3 == 1) {
#endif
                lastLastDiag[64 * jdiv64 + 32] = curDiagVal[0];
                lastLastDiag[64 * jdiv64 + 33] = curDiagVal[1];
                lastLastDiag[64 * jdiv64 + 34] = curDiagVal[2];
                lastLastDiag[64 * jdiv64 + 35] = curDiagVal[3];
                lastLastDiag[64 * jdiv64 + 36] = curDiagVal[4];
                lastLastDiag[64 * jdiv64 + 37] = curDiagVal[5];
                lastLastDiag[64 * jdiv64 + 38] = curDiagVal[6];
                lastLastDiag[64 * jdiv64 + 39] = curDiagVal[7];
                lastLastDiag[64 * jdiv64 + 40] = curDiagVal[8];
                lastLastDiag[64 * jdiv64 + 41] = curDiagVal[9];
                lastLastDiag[64 * jdiv64 + 42] = curDiagVal[10];
                lastLastDiag[64 * jdiv64 + 43] = curDiagVal[11];
                lastLastDiag[64 * jdiv64 + 44] = curDiagVal[12];
                lastLastDiag[64 * jdiv64 + 45] = curDiagVal[13];
                lastLastDiag[64 * jdiv64 + 46] = curDiagVal[14];
                lastLastDiag[64 * jdiv64 + 47] = curDiagVal[15];
                lastLastDiag[64 * jdiv64 + 48] = curDiagVal[16];
                lastLastDiag[64 * jdiv64 + 49] = curDiagVal[17];
                lastLastDiag[64 * jdiv64 + 50] = curDiagVal[18];
                lastLastDiag[64 * jdiv64 + 51] = curDiagVal[19];
                lastLastDiag[64 * jdiv64 + 52] = curDiagVal[20];
                lastLastDiag[64 * jdiv64 + 53] = curDiagVal[21];
                lastLastDiag[64 * jdiv64 + 54] = curDiagVal[22];
                lastLastDiag[64 * jdiv64 + 55] = curDiagVal[23];
                lastLastDiag[64 * jdiv64 + 56] = curDiagVal[24];
                lastLastDiag[64 * jdiv64 + 57] = curDiagVal[25];
                lastLastDiag[64 * jdiv64 + 58] = curDiagVal[26];
                lastLastDiag[64 * jdiv64 + 59] = curDiagVal[27];
                lastLastDiag[64 * jdiv64 + 60] = curDiagVal[28];
                lastLastDiag[64 * jdiv64 + 61] = curDiagVal[29];
                lastLastDiag[64 * jdiv64 + 62] = curDiagVal[30];
                lastLastDiag[64 * jdiv64 + 63] = curDiagVal[31];
            }
            else {
                lastDiag[64 * jdiv64 + 32] = curDiagVal[0];
                lastDiag[64 * jdiv64 + 33] = curDiagVal[1];
                lastDiag[64 * jdiv64 + 34] = curDiagVal[2];
                lastDiag[64 * jdiv64 + 35] = curDiagVal[3];
                lastDiag[64 * jdiv64 + 36] = curDiagVal[4];
                lastDiag[64 * jdiv64 + 37] = curDiagVal[5];
                lastDiag[64 * jdiv64 + 38] = curDiagVal[6];
                lastDiag[64 * jdiv64 + 39] = curDiagVal[7];
                lastDiag[64 * jdiv64 + 40] = curDiagVal[8];
                lastDiag[64 * jdiv64 + 41] = curDiagVal[9];
                lastDiag[64 * jdiv64 + 42] = curDiagVal[10];
                lastDiag[64 * jdiv64 + 43] = curDiagVal[11];
                lastDiag[64 * jdiv64 + 44] = curDiagVal[12];
                lastDiag[64 * jdiv64 + 45] = curDiagVal[13];
                lastDiag[64 * jdiv64 + 46] = curDiagVal[14];
                lastDiag[64 * jdiv64 + 47] = curDiagVal[15];
                lastDiag[64 * jdiv64 + 48] = curDiagVal[16];
                lastDiag[64 * jdiv64 + 49] = curDiagVal[17];
                lastDiag[64 * jdiv64 + 50] = curDiagVal[18];
                lastDiag[64 * jdiv64 + 51] = curDiagVal[19];
                lastDiag[64 * jdiv64 + 52] = curDiagVal[20];
                lastDiag[64 * jdiv64 + 53] = curDiagVal[21];
                lastDiag[64 * jdiv64 + 54] = curDiagVal[22];
                lastDiag[64 * jdiv64 + 55] = curDiagVal[23];
                lastDiag[64 * jdiv64 + 56] = curDiagVal[24];
                lastDiag[64 * jdiv64 + 57] = curDiagVal[25];
                lastDiag[64 * jdiv64 + 58] = curDiagVal[26];
                lastDiag[64 * jdiv64 + 59] = curDiagVal[27];
                lastDiag[64 * jdiv64 + 60] = curDiagVal[28];
                lastDiag[64 * jdiv64 + 61] = curDiagVal[29];
                lastDiag[64 * jdiv64 + 62] = curDiagVal[30];
                lastDiag[64 * jdiv64 + 63] = curDiagVal[31];
            }
        }
       

        if(i - j - refLength + 1 >= 0 && i - j - refLength + 1 < 32 && i - refLength + 1 >= 0 && i - refLength + 1 < altLength){
            sw_finalRow[i + 2 - refLength] = curDiagVal[i - j - refLength + 1];
        }
    
        if(altLength - 1 - j >= 0 && altLength - 1 - j < 32 && i - altLength + 1 >= 0 && i - altLength + 1 < refLength){
            sw_finalCol[i + 2 - altLength] = curDiagVal[altLength - 1 - j];
        }

        
        if (j_PE == 512 / 32 - 1) {
            i = i + 1;
        }
        if (j_PE == 512 / 32 - 1) {
#ifdef CUSTOM
            ap_uint<2> rose_temp;
            if (imod3 == (ap_uint<2>)(2)) {
                rose_temp = 0;
            }
#else
            int rose_temp;
            if (imod3 == 2) {
                rose_temp = 0;
            }
#endif
            else {
                rose_temp = imod3 + 1;
            }
            imod3 = rose_temp;
        }
        if (j_PE == 512 / 32 - 1) {
            j_PE = 0;
        }
        else {
            j_PE = j_PE + 1;
        }
        if(k > (refLength + altLength - 1) * 512 / 32){
            break;
        }
    }
    signed short p1 = 0;
    signed short p2 = 0;
    int maxscore = (int )(- 2147483648);
    int segment_length = 0;
    int curScore = 0;
    if (overhang_strategy == 1) {
        p1 = refLength;
        p2 = altLength;
    }
    else {
        p2 = altLength;
        for (i = 1; i < nrow; ++i) {
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=510
            curScore = sw_finalCol[i];
            if (curScore >= maxscore) {
                p1 = i;
                maxscore = curScore;
            }
        }
//now look for a larger score on the bottom-most row
        if (overhang_strategy != 2) {
            for (j = 1; j < ncol; ++j){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=510
                curScore = sw_finalRow[j];
                short abs1 = refLength - j;
                if (abs1 < 0) {
                    abs1 = -abs1;
                }
                short abs2 = p1 - p2;
                if (abs2 < 0) {
                    abs2 = -abs2;
                }
                if (curScore > maxscore || (curScore == maxscore && abs1 < abs2)) {
                    p1 = refLength;
                    p2 = j;
                    maxscore = curScore;
                    segment_length = altLength - j;
                }
            }
        }
    }
    p1--;
    p2--;
    short cigarSize = (short )0;
    if (segment_length > 0 && overhang_strategy == 0) {
        lengths[0] = ((short )segment_length);
        states[0] = ((short )4);
        cigarSize++;
        segment_length = 0;
    }
    short state = (short )0;
    short btr = 0;
    char new_state = 0;
    short step_length = 0;
    do {
        btr = (btrack[p1 * 512 + p2]);
        step_length = 1;
        if (btr > 0) {
            new_state = 2;
            step_length = btr;
            p1 -= btr;
        }
        else {
            if (btr < 0) {
                new_state = 1;
                step_length = -btr;
                p2 += btr;
            }
            else {
                new_state = 0;
                p1--;
                p2--;
            }
        }
    
        if (new_state == ((char)state)) {
            segment_length += step_length;
        }
        else {
            lengths[cigarSize] = ((short )segment_length);
            states[cigarSize] = state;
            cigarSize++;
            segment_length = step_length;
            state = ((short )new_state);
        }
    }while (p1 >= 0 && p2 >= 0);
    p1++;
    p2++;
    if (overhang_strategy == 0) {
        lengths[cigarSize] = ((short )segment_length);
        states[cigarSize] = state;
        cigarSize++;
        if (p2 > 0) {
            lengths[cigarSize] = ((short )p2);
            states[cigarSize] = ((short )4);
            cigarSize++;
        }
        alignment_offset_ptr[0] = ((short )p1);
    }
    else {
        if (overhang_strategy == 3) {
            lengths[cigarSize] = ((short )(segment_length + p2));
            states[cigarSize] = state;
            cigarSize++;
            alignment_offset_ptr[0] = ((short )(p1 - p2));
        }
        else {
            lengths[cigarSize] = ((short )segment_length);
            states[cigarSize] = state;
            cigarSize++;
            if (p1 > 0) {
                lengths[cigarSize] = ((short )p1);
                states[cigarSize] = ((short )2);
                cigarSize++;
            }
            else {
                if (p2 > 0) {
                    lengths[cigarSize] = ((short )p2);
                    states[cigarSize] = ((short )1);
                    cigarSize++;
                }
            }
            alignment_offset_ptr[0] = ((short )0);
        }
    }
    cigarElementNum[0] = cigarSize;
    return 0;
}

static int smithWatermanMerlin_dummy_pos;

static void merlin_memcpy_2(char dst[512],int dst_idx_0,char *src,int src_idx_0,unsigned long len,unsigned long max_len)
{
#pragma HLS inline off
#pragma HLS function_instantiate variable=dst_idx_0,src_idx_0
    long long i;
    long long total_offset1 = dst_idx_0;
    long long total_offset2 = src_idx_0;
    for (i = 0; i < len; ++i) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT max=512
        dst[total_offset1 + i] = src[total_offset2 + i];
    }
}

static void merlin_memcpy_3(char dst[512],int dst_idx_1,char *src,int src_idx_0,unsigned long len,unsigned long max_len)
{
#pragma HLS inline off
#pragma HLS function_instantiate variable=dst_idx_1,src_idx_0
    long long i;
    long long total_offset1 = dst_idx_1;
    long long total_offset2 = src_idx_0;
    for (i = 0; i < len; ++i) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT max=512
        dst[total_offset1 + i] = src[total_offset2 + i];
    }
}

extern "C" { 
__kernel void smithWatermanMerlin(char *inputs,int refLength,int batchSize,int overhang_strategy,int w_match,int w_mismatch,int w_open,int w_extend,short *outputs)
{
#pragma HLS INTERFACE m_axi port=inputs offset=slave bundle=inputs depth=134152
#pragma HLS INTERFACE m_axi port=outputs offset=slave bundle=outputs depth=266764
#pragma HLS INTERFACE s_axilite port=batchSize bundle=control
#pragma HLS INTERFACE s_axilite port=inputs bundle=control
#pragma HLS INTERFACE s_axilite port=outputs bundle=control
#pragma HLS INTERFACE s_axilite port=overhang_strategy bundle=control
#pragma HLS INTERFACE s_axilite port=refLength bundle=control
#pragma HLS INTERFACE s_axilite port=w_extend bundle=control
#pragma HLS INTERFACE s_axilite port=w_match bundle=control
#pragma HLS INTERFACE s_axilite port=w_mismatch bundle=control
#pragma HLS INTERFACE s_axilite port=w_open bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
 
/*
 * in inputs, the first 2*batchSize bytes are altLengths, each 2 bytes is a short
 * The rest data is ref and alts
 * in outputs, the first 2 short values are cigarInfoLength, then batchSize of short value is cigarSize, the rest value is cigarInfo.
 * */
  
  
    short batchCount = 0;
    char ref[512];
    char altLengths_tmp[260 * 2];
    short altLengths[260];
    for(short i = 0; i < 2 * batchSize; i++){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=520
        altLengths_tmp[i] = inputs[i];
    }

    for (short i = 0; i < batchSize; i++){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT max=260
#pragma HLS latency min=3
        short tmp1 = (short)(altLengths_tmp[2 * i + 1] << 8);
        short tmp2 = (signed short)(altLengths_tmp[2 * i] & 0xff);
        altLengths[i] = tmp1 | tmp2;
    }

    merlin_memcpy_2(ref,0,inputs,2 * batchSize,((unsigned long )512) * 1UL,512UL);
    int cigarElementPtr = 0;
#define KERNEL_NUM 3
    short w_match_local[KERNEL_NUM];
#pragma HLS array_partition variable=w_match_local complete dim=1
    short w_mismatch_local[KERNEL_NUM];
#pragma HLS array_partition variable=w_mismatch_local complete dim=1
    short w_open_local[KERNEL_NUM];
#pragma HLS array_partition variable=w_open_local complete dim=1
    short w_extend_local[KERNEL_NUM];
#pragma HLS array_partition variable=w_extend_local complete dim=1
    short refLength_local[KERNEL_NUM];
#pragma HLS array_partition variable=refLength_local complete dim=1
    ap_uint<2> overhang_strategy_local[KERNEL_NUM];
#pragma HLS array_partition variable=overhang_strategy_local complete dim=1

    for(char i = 0; i < KERNEL_NUM; i++){
#pragma HLS latency min=1
        refLength_local[i] = refLength;
        w_match_local[i] = (short)w_match;
        w_mismatch_local[i] = (short)w_mismatch;
        w_open_local[i] = (short)w_open;
        w_extend_local[i] = (short)w_extend;
        overhang_strategy_local[i] = (ap_uint<2>)overhang_strategy;
    }

    batchSize >= 0 && batchSize < 260?((void )0) : __assert_fail("batchSize >= 0 && batchSize < 260","smithWatermanMerlin.cpp",((unsigned int )47),__PRETTY_FUNCTION__);
    for (batchCount = 0; batchCount < batchSize; batchCount += KERNEL_NUM) {
#pragma HLS LOOP_TRIPCOUNT max=87
        short altLengthsLocal[KERNEL_NUM];
#pragma HLS array_partition variable=altLengthsLocal complete dim=1
        char alt[KERNEL_NUM][512];
#pragma HLS array_partition variable=alt complete dim=1
        short lengths[KERNEL_NUM][512];
#pragma HLS array_partition variable=lengths complete dim=1
        short states[KERNEL_NUM][512];
#pragma HLS array_partition variable=states complete dim=1
        short cigarElementNum[KERNEL_NUM];
#pragma HLS array_partition variable=cigarElementNum complete dim=1
        short alignment_offset[KERNEL_NUM];
#pragma HLS array_partition variable=alignment_offset complete dim=1
        for (char i = 0; i < KERNEL_NUM; ++i) {
            if (batchCount + i < batchSize) {
                merlin_memcpy_3(alt[i],0,inputs,2 * batchSize + (batchCount + i + 1) * 512,((unsigned long )512) * 1UL,512UL);
                altLengthsLocal[i] = ((short )altLengths[batchCount + i]);
            }
        }
        for (char i = 0; i < KERNEL_NUM; ++i){
#pragma HLS unroll
            if (batchCount + i < batchSize) {
                smithWatermanKernel(ref,alt[i],refLength_local[i],altLengthsLocal[i],overhang_strategy_local[i],0,w_match_local[i],w_mismatch_local[i],w_open_local[i],w_extend_local[i],lengths[i],states[i],&cigarElementNum[i],&alignment_offset[i]);
            }
        }
        for (char i = 0; i < KERNEL_NUM; ++i) {
            if (batchCount + i < batchSize) {
                outputs[2 + batchCount + i] = cigarElementNum[i];
                for (int j = 0; j < ((int )cigarElementNum[i]); j++) {
#pragma HLS pipeline
                    outputs[2 + batchSize + cigarElementPtr] = lengths[i][j];
                    cigarElementPtr++;
                    outputs[2 + batchSize + cigarElementPtr] = states[i][j];
                    cigarElementPtr++;
                }
                outputs[2 + batchSize + cigarElementPtr] = alignment_offset[i];
                cigarElementPtr++;
            }
        }
    }
    outputs[0] = (short)(cigarElementPtr & 0xffff);
    outputs[1] = (short)((cigarElementPtr >> 16) & 0xffff);
}
}
