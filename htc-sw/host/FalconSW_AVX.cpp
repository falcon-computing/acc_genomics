/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: Function definitions for Intel AVX Smith-Waterman implementation 
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <emmintrin.h>
#include <immintrin.h>

#include "common.h"

#define SW_BEST 0
#define SW_BASELINE 1
#define SW_ROWWISE 2
#define SW_ROWWISE_SIMD 3
#define SW_ROWWISE_SIMD_UNROLL 4


extern struct Cigar* retCigarBatch;
extern long SWPairwiseAlignment_C_time;
extern long isSWFailure_C_time;
extern long trimCigarByBases_C_time;
extern long getReferenceLength_C_time;
extern long leftAlignCigar_C_time;
extern long calMatrix_C_time;
extern long calCigar_C_time;
extern long malloc_time;
extern long SW_complexity;
extern long lastKernel_C_time;

int calculateMatrixRowWiseOpt(char*, char*, int**, int, int, int**, int, int, int, int, int, int);
int max_scan(int* sw_prime, int* best_gap_h, int* gap_size_h, int ncol, int w_open, int w_extend){
    int j;
    signed int lowInitValue = -1073741824;
    for(j = 0; j < ncol - 1; j++){
        best_gap_h[j] = sw_prime[j] + w_open + (ncol - j - 2) * w_extend;
    }
    for(j = 0; j < 1024; j += 2){
        if(best_gap_h[j + 1] > best_gap_h[j]){
            gap_size_h[j] = 1;
            gap_size_h[j + 1] = 1;
        }
        else{
            best_gap_h[j + 1] = best_gap_h[j];
            gap_size_h[j] = 1;
            gap_size_h[j + 1] = 2;
        }
    }
    for(j = 0; j < 1024; j += 4){
        if(best_gap_h[j + 3] <= best_gap_h[j + 1]){
            best_gap_h[j + 3] = best_gap_h[j + 1];
            gap_size_h[j + 3] = gap_size_h[j + 1] + 2;
        }
    }
    for(j = 0; j < 1024; j += 8){
        if(best_gap_h[j + 7] <= best_gap_h[j + 3]){
            best_gap_h[j + 7] = best_gap_h[j + 3];
            gap_size_h[j + 7] = gap_size_h[j + 3] + 4;
        }
    }
    for(j = 0; j < 1024; j += 16){
        if(best_gap_h[j + 15] <= best_gap_h[j + 7]){
            best_gap_h[j + 15] = best_gap_h[j + 7];
            gap_size_h[j + 15] = gap_size_h[j + 7] + 8;
        }
    }
    for(j = 0; j < 1024; j += 32){
        if(best_gap_h[j + 31] <= best_gap_h[j + 15]){
            best_gap_h[j + 31] = best_gap_h[j + 15];
            gap_size_h[j + 31] = gap_size_h[j + 15] + 16;
        }
    }
    for(j = 0; j < 1024; j += 64){
        if(best_gap_h[j + 63] <= best_gap_h[j + 31]){
            best_gap_h[j + 63] = best_gap_h[j + 31];
            gap_size_h[j + 63] = gap_size_h[j + 31] + 32;
        }
    }
    for(j = 0; j < 1024; j += 128){
        if(best_gap_h[j + 127] <= best_gap_h[j + 63]){
            best_gap_h[j + 127] = best_gap_h[j + 63];
            gap_size_h[j + 127] = gap_size_h[j + 63] + 64;
        }
    }
    for(j = 0; j < 1024; j += 256){
        if(best_gap_h[j + 255] <= best_gap_h[j + 127]){
            best_gap_h[j + 255] = best_gap_h[j + 127];
            gap_size_h[j + 255] = gap_size_h[j + 127] + 128;
        }
    }
    for(j = 0; j < 1024; j += 512){
        if(best_gap_h[j + 511] <= best_gap_h[j + 255]){
            best_gap_h[j + 511] = best_gap_h[j + 255];
            gap_size_h[j + 511] = gap_size_h[j + 255] + 256;
        }
    }
    if(best_gap_h[1023] <= best_gap_h[511]){
        best_gap_h[1024] = best_gap_h[511];
        gap_size_h[1024] = gap_size_h[511] + 512;
    }
    else{
        best_gap_h[1024] = best_gap_h[1023];
        gap_size_h[1024] = gap_size_h[1023];
    }
    best_gap_h[1023] = best_gap_h[511];
    gap_size_h[1023] = gap_size_h[511];
    best_gap_h[511] = lowInitValue;
    gap_size_h[511] = 0;
    //up sweep complete, start down sweep
    for(j = 0; j < 1024; j += 512){
        int tmp_best_gap_h = best_gap_h[j + 511];
        int tmp_gap_size_h = gap_size_h[j + 511];
        if(best_gap_h[j + 255] > best_gap_h[j + 511]){
            best_gap_h[j + 511] = best_gap_h[j + 255];
            gap_size_h[j + 511] = gap_size_h[j + 255];
        }
        else{
            gap_size_h[j + 511] = gap_size_h[j + 511] + 256;
        }
        best_gap_h[j + 255] = tmp_best_gap_h;
        gap_size_h[j + 255] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 256){
        int tmp_best_gap_h = best_gap_h[j + 255];
        int tmp_gap_size_h = gap_size_h[j + 255];
        if(best_gap_h[j + 127] > best_gap_h[j + 255]){
            best_gap_h[j + 255] = best_gap_h[j + 127];
            gap_size_h[j + 255] = gap_size_h[j + 127];
        }
        else{
            gap_size_h[j + 255] = gap_size_h[j + 255] + 128;
        }
        best_gap_h[j + 127] = tmp_best_gap_h;
        gap_size_h[j + 127] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 128){
        int tmp_best_gap_h = best_gap_h[j + 127];
        int tmp_gap_size_h = gap_size_h[j + 127];
        if(best_gap_h[j + 63] > best_gap_h[j + 127]){
            best_gap_h[j + 127] = best_gap_h[j + 63];
            gap_size_h[j + 127] = gap_size_h[j + 63];
        }
        else{
            gap_size_h[j + 127] = gap_size_h[j + 127] + 64;
        }
        best_gap_h[j + 63] = tmp_best_gap_h;
        gap_size_h[j + 63] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 64){
        int tmp_best_gap_h = best_gap_h[j + 63];
        int tmp_gap_size_h = gap_size_h[j + 63];
        if(best_gap_h[j + 31] > best_gap_h[j + 63]){
            best_gap_h[j + 63] = best_gap_h[j + 31];
            gap_size_h[j + 63] = gap_size_h[j + 31];
        }
        else{
            gap_size_h[j + 63] = gap_size_h[j + 63] + 32;
        }
        best_gap_h[j + 31] = tmp_best_gap_h;
        gap_size_h[j + 31] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 32){
        int tmp_best_gap_h = best_gap_h[j + 31];
        int tmp_gap_size_h = gap_size_h[j + 31];
        if(best_gap_h[j + 15] > best_gap_h[j + 31]){
            best_gap_h[j + 31] = best_gap_h[j + 15];
            gap_size_h[j + 31] = gap_size_h[j + 15];
        }
        else{
            gap_size_h[j + 31] = gap_size_h[j + 31] + 16;
        }
        best_gap_h[j + 15] = tmp_best_gap_h;
        gap_size_h[j + 15] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 16){
        int tmp_best_gap_h = best_gap_h[j + 15];
        int tmp_gap_size_h = gap_size_h[j + 15];
        if(best_gap_h[j + 7] > best_gap_h[j + 15]){
            best_gap_h[j + 15] = best_gap_h[j + 7];
            gap_size_h[j + 15] = gap_size_h[j + 7];
        }
        else{
            gap_size_h[j + 15] = gap_size_h[j + 15] + 8;
        }
        best_gap_h[j + 7] = tmp_best_gap_h;
        gap_size_h[j + 7] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 8){
        int tmp_best_gap_h = best_gap_h[j + 7];
        int tmp_gap_size_h = gap_size_h[j + 7];
        if(best_gap_h[j + 3] > best_gap_h[j + 7]){
            best_gap_h[j + 7] = best_gap_h[j + 3];
            gap_size_h[j + 7] = gap_size_h[j + 3];
        }
        else{
            gap_size_h[j + 7] = gap_size_h[j + 7] + 4;
        }
        best_gap_h[j + 3] = tmp_best_gap_h;
        gap_size_h[j + 3] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 4){
        int tmp_best_gap_h = best_gap_h[j + 3];
        int tmp_gap_size_h = gap_size_h[j + 3];
        if(best_gap_h[j + 1] > best_gap_h[j + 3]){
            best_gap_h[j + 3] = best_gap_h[j + 1];
            gap_size_h[j + 3] = gap_size_h[j + 1];
        }
        else{
            gap_size_h[j + 3] = gap_size_h[j + 3] + 2;
        }
        best_gap_h[j + 1] = tmp_best_gap_h;
        gap_size_h[j + 1] = tmp_gap_size_h;
    }
    for(j = 0; j < 1024; j += 2){
        int tmp_best_gap_h = best_gap_h[j + 1];
        int tmp_gap_size_h = gap_size_h[j + 1];
        if(best_gap_h[j] > best_gap_h[j + 1]){
            best_gap_h[j + 1] = best_gap_h[j];
            gap_size_h[j + 1] = gap_size_h[j];
        }
        else{
            gap_size_h[j + 1] = gap_size_h[j + 1] + 1;
        }
        best_gap_h[j] = tmp_best_gap_h;
        gap_size_h[j] = tmp_gap_size_h;
    }
    for(j = 1; j < ncol; j++){
        best_gap_h[j] = best_gap_h[j] - (ncol - j - 1) * w_extend;
    }

    int best_gap_h_golden[MAX_SEQ_LENGTH];
    int gap_size_h_golden[MAX_SEQ_LENGTH];
    int best_gap_h_target[MAX_SEQ_LENGTH];
    int sw_prime_after[MAX_SEQ_LENGTH];
    int i;
    int cur_best = lowInitValue;
    for(i = 0; i < ncol + 1; ++i){
        best_gap_h_golden[i] = lowInitValue;
        gap_size_h_golden[i] = 0;
    }
 
    for(j = 1; j < ncol; j++){
        best_gap_h_golden[j] = sw_prime[j - 1] + w_open + (ncol - j - 1) * w_extend;
    }
    for(j = 1; j < ncol; j++){
        if(best_gap_h_golden[j] > cur_best){
            cur_best = best_gap_h_golden[j];
            gap_size_h_golden[j] = 1;
        }
        else{
            best_gap_h_golden[j] = cur_best;
            gap_size_h_golden[j] = gap_size_h_golden[j - 1] + 1;
        }
    }

    for(j = 1; j < ncol; j++){
        sw_prime_after[j] = sw_prime[j - 1] + w_open + (ncol - j - 1) * w_extend;
    }
    for(j = 1; j < ncol; j++){
        best_gap_h_target[j] = best_gap_h[j] + (ncol - j - 1) * w_extend;
    }

    //now check the result
    for(j = 1; j < ncol; j++){
        if(best_gap_h_target[j] != best_gap_h_golden[j] || gap_size_h[j] != gap_size_h_golden[j]){
            printf("mismatch at %d, best_gap_h = %d, best_gap_h_golde = %d, gap_size_h = %d, gap_size_h_golden =%d\n", j, best_gap_h_target[j], best_gap_h_golden[j], gap_size_h[j], gap_size_h_golden[j]);
            int tmp;
            printf("sw_prime is \n");
            for(tmp = 0; tmp < ncol; tmp += 1){
                printf("%d, ", sw_prime_after[tmp]);
            }
            printf("\nbest_gap_h_golden is \n");
            for(tmp = 0; tmp < ncol; tmp += 1){
                printf("%d, ", best_gap_h_golden[tmp]);
            }
            printf("\ngap_size_h_golden is \n");
            for(tmp = 0; tmp < ncol; tmp += 1){
                printf("%d, ", gap_size_h_golden[tmp]);
            }
            printf("\nbest_gap_h is \n");
            for(tmp = 0; tmp < ncol; tmp += 1){
                printf("%d, ", best_gap_h_target[tmp]);
            }
            printf("\ngap_size_h is \n");
            for(tmp = 0; tmp < ncol; tmp += 1){
                printf("%d, ", gap_size_h[tmp]);
            }
            exit(-1);
        }
    }
    printf("max scan is right\n");

    return 0;
}



int SWPairwiseAlignmentMultiBatch(char* ref, int refLength, char alts[][MAX_SEQ_LENGTH], int batchSize, int* altLengths, struct Cigar* cigarResults, int* alignmentOffsets, int overhang_strategy, int option){
    int i = 0;
    for(i = 0; i < batchSize; ++i){
        if(SWPairwiseAlignmentOneBatch(ref, alts[i], refLength, altLengths[i], &cigarResults[i], &alignmentOffsets[i], overhang_strategy, option) == -1){
            printf("Exception happens in the SWPairwiseAlignmentOneBatch\n");
            return -1;
        }
    }
    return 0;
}

int SWPairwiseAlignmentOneBatch(char* reference, char* alternate, int refLength, int altLength, struct Cigar* cigarResult, int* alignment_offset_ptr, int overhang_strategy, int option){
    SW_complexity += refLength * altLength;
    int n = refLength + 1;
    int m = altLength + 1;
    int i = 0;
    int j = 0;
 
    //first dynamic allocate first
    struct timespec time1, time2, time_diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    
    int** sw;
    int* sw_oneD;
    int** btrack;
    int* btrack_oneD;

    if(!(sw = (int**)malloc(n * sizeof(int*)))){
        return -1;
    }
    if(!(sw_oneD = (int*)malloc(n * m * sizeof(int)))){
        return -1;
    }
    if(!(btrack = (int**)malloc(n * sizeof(int*)))){
        return -1;
    }
    if(!(btrack_oneD = (int*)malloc(n * m * sizeof(int)))){
        return -1;
    }
    for(i = 0; i < n ; ++i){
        sw[i] = &sw_oneD[i * m];
        for(j = 0; j < m; ++j){
            sw[i][j] = 0;
        }
    }
    
    for(i = 0; i < n ; ++i){
        btrack[i] = &btrack_oneD[i * m];
        for(j = 0; j < m; ++j){
            btrack[i][j] = 0;
        }
    }

    clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    malloc_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);

    int calculate_matrix_ret;
    if(option == SW_BEST){
         calculate_matrix_ret = calculateMatrixRowWiseSIMDUnroll4x(reference, alternate, sw, m, n, btrack, overhang_strategy, 0, W_MATCH, W_MISMATCH, W_OPEN, W_EXTEND);
         //calculate_matrix_ret = calculateMatrixRowWiseSIMD(reference, alternate, sw, m, n, btrack, overhang_strategy, 0);
    }
    else if(option == SW_BASELINE){
         calculate_matrix_ret = calculateMatrixOneBatch(reference, alternate, sw, m, n, btrack, overhang_strategy, 0);
    }
    else if(option == SW_ROWWISE){
         calculate_matrix_ret = calculateMatrixRowWise(reference, alternate, sw, m, n, btrack, overhang_strategy, 0);
    }
    else if(option == SW_ROWWISE_SIMD){
         calculate_matrix_ret = calculateMatrixRowWiseSIMD(reference, alternate, sw, m, n, btrack, overhang_strategy, 0);
    }
    else {
         calculate_matrix_ret = calculateMatrixRowWiseOpt(reference, alternate, sw, m, n, btrack, overhang_strategy, 0, W_MATCH, W_MISMATCH, W_OPEN, W_EXTEND);
        
    }

    if(calculate_matrix_ret < 0){
#ifndef NDEBUG
        printf("calculate_matrix has error\n");
#endif
        return -1;
    }
    
    clock_gettime(CLOCK_REALTIME, &time1);
    int calculate_cigar_ret = calculateCigarOneBatch(sw, btrack, n, m, overhang_strategy, cigarResult, alignment_offset_ptr);
    clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    calCigar_C_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);
    if(calculate_cigar_ret < 0){
#ifndef NDEBUG
        printf("calculate_cigar has error\n");
#endif
        return -1;
    }
    
    //free sw and btrack
    free(sw_oneD);
    free(btrack_oneD);
 //   for(i = 0; i < n; i++){
 //       free(sw[i]);
 //       free(btrack[i]);
 //   }
    free(sw);
    free(btrack);
    return 0;

    
}

int calculateMatrixRowWise(char* ref, char* alt, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff){
    signed int lowInitValue = -1073741824;
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < ncol; i++){
            currentValue += W_EXTEND;
            sw[0][i] = currentValue;
        }
        sw[1][0] = W_OPEN;
        currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw[i][0] = currentValue;
        }
        
    }
    int sw_prime[MAX_SEQ_LENGTH];
    int sw_prime_first_element[MAX_SEQ_LENGTH];
    int gap_size_h[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h[MAX_SEQ_LENGTH];
    int diagOrDown[MAX_SEQ_LENGTH]; //diag is 1, down is 0
    for(i = 0; i < ncol + 1; i++){
        best_gap_v[i] = lowInitValue; //pointing to F[i-1][j] when i is odd
        gap_size_v[i] = 0;
        sw_prime[i] = 0;
    }

    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw_prime_first_element[0] = 0;
        sw_prime_first_element[1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw_prime_first_element[i] = currentValue;
        }
    }
    else{
        for(i = 0; i < nrow; i++){
            sw_prime_first_element[i] = 0;
        }
    }
            

    for(i = 0; i < ncol + 1; ++i){
        best_gap_h[i] = lowInitValue;
        gap_size_h[i] = 0;
    }
    int j = 0;
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int step_down = 0;
    int step_right = 0;
    int prev_gap = 0;
    int kd = 0;
    int ki = 0;

    for(i = 1; i < nrow; i++){
        a_base = ref[i - 1];
        for(j = 1; j < ncol; j++){
            b_base = alt[j - 1];
            step_diag =  sw[i - 1][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }
            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            if(step_diag >= step_down){
                sw_prime[j] = max(step_diag, MATRIX_MIN_CUTOFF);
                btrack[i][j] = 0;
                diagOrDown[j] = 1;
            }
            else{
                sw_prime[j] = max(step_down, MATRIX_MIN_CUTOFF);
                btrack[i][j] = kd;
                diagOrDown[j] = 0;
            }
        }
        int cur_best = lowInitValue;
        /*
        for(j = 1; j < ncol; j++){
            cur_best  = lowInitValue;
            for(k = 0; k < j; k++){
                cur = sw_prime[k] + W_OPEN + (j - k - 1) * W_EXTEND;
                if(cur > cur_best){
                    cur_best = cur;
                    gap_size_h_best = j - k;
                }
            }
            
            gap_size_h[j] = gap_size_h_best;
            best_gap_h[j] = cur_best;
        }
        */
        sw_prime[0] = sw_prime_first_element[i];   
        for(j = 1; j < ncol; j++){
            best_gap_h[j] = sw_prime[j - 1] + W_OPEN + (ncol - j - 1) * W_EXTEND;
        }
        for(j = 1; j < ncol; j++){
            if(best_gap_h[j] > cur_best){
                cur_best = best_gap_h[j];
                gap_size_h[j] = 1;
            }
            else{
                best_gap_h[j] = cur_best;
                gap_size_h[j] = gap_size_h[j - 1] + 1;
            }
        }
        for(j = 1; j < ncol; j++){
            best_gap_h[j] = best_gap_h[j] - (ncol - j - 1) * W_EXTEND;
        }



        for(j = 1; j < ncol; j++){
            step_right = best_gap_h[j];
            ki = gap_size_h[j];
            
            if(diagOrDown[j]){
                if(step_right > sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }
            else{
                if(step_right >= sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }         
        }
    }
    return 0;
}

int calculateMatrixRowWiseOpt(char* ref, char* alt, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff, int w_match, int w_mismatch, int w_open, int w_extend){
    signed int lowInitValue = -1073741824;
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = w_open;
        int currentValue = w_open;
        for(i = 2; i < ncol; i++){
            currentValue += w_extend;
            sw[0][i] = currentValue;
        }
        sw[1][0] = w_open;
        currentValue = w_open;
        for(i = 2; i < nrow; i++){
            currentValue += w_extend;
            sw[i][0] = currentValue;
        }
        
    }
    int sw_prime[MAX_SEQ_LENGTH];
    int sw_prime_first_element[MAX_SEQ_LENGTH];
    int gap_size_h[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h[MAX_SEQ_LENGTH];
    int diagOrDown[MAX_SEQ_LENGTH]; //diag is 1, down is 0
    for(i = 0; i < ncol + 1; i++){
        best_gap_v[i] = lowInitValue; //pointing to F[i-1][j] when i is odd
        gap_size_v[i] = 0;
        sw_prime[i] = 0;
    }

    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw_prime_first_element[0] = 0;
        sw_prime_first_element[1] = w_open;
        int currentValue = w_open;
        for(i = 2; i < nrow; i++){
            currentValue += w_extend;
            sw_prime_first_element[i] = currentValue;
        }
    }
    else{
        for(i = 0; i < nrow; i++){
            sw_prime_first_element[i] = 0;
        }
    }
            

    for(i = 0; i < ncol + 1; ++i){
        best_gap_h[i] = lowInitValue;
        gap_size_h[i] = 0;
    }
    int j = 0;
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int step_down = 0;
    int step_right = 0;
    int prev_gap = 0;
    int kd = 0;
    int ki = 0;

    for(i = 1; i < nrow; i++){
        a_base = ref[i - 1];
        for(j = 1; j < ncol; j++){
            b_base = alt[j - 1];
            step_diag =  sw[i - 1][j - 1] + wd(a_base, b_base, w_match, w_mismatch);
            prev_gap = sw[i - 1][j] + w_open;
            best_gap_v[j] += w_extend;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }
            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            if(step_diag >= step_down){
                sw_prime[j] = max(step_diag, MATRIX_MIN_CUTOFF);
                btrack[i][j] = 0;
                diagOrDown[j] = 1;
            }
            else{
                sw_prime[j] = max(step_down, MATRIX_MIN_CUTOFF);
                btrack[i][j] = kd;
                //btrack[i][j] = 1;
                
                diagOrDown[j] = 0;
            }
        }
        //int cur_best = lowInitValue;
        sw_prime[0] = sw_prime_first_element[i];   
        /*for(j = 1; j < ncol; j++){
            best_gap_h[j] = sw_prime[j - 1] + w_open + (ncol - j - 1) * w_extend;
        }
        for(j = 1; j < ncol; j++){
            if(best_gap_h[j] > cur_best){
                cur_best = best_gap_h[j];
                gap_size_h[j] = 1;
            }
            else{
                best_gap_h[j] = cur_best;
                gap_size_h[j] = gap_size_h[j - 1] + 1;
            }
        }
        for(j = 1; j < ncol; j++){
            best_gap_h[j] = best_gap_h[j] - (ncol - j - 1) * w_extend;
        }*/
        max_scan(sw_prime, best_gap_h, gap_size_h, ncol, w_open, w_extend);



        for(j = 1; j < ncol; j++){
            step_right = best_gap_h[j];
            ki = gap_size_h[j];
            
            if(diagOrDown[j]){
                if(step_right > sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                    //btrack[i][j] = -1;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }
            else{
                if(step_right >= sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                    //btrack[i][j] = -1;
                    /*if(ki > 1){
                        int tmp;
                        printf("btrack[%d][%d] is supposed to be %d\n",i, j, kd);
                        for(tmp = 0; tmp < kd; tmp++){
                            printf("btrack[%d][%d] = %d\n", i - tmp , j , btrack[i-tmp][j]);
                            if(btrack[i-tmp][j] != 1){
                                printf("error\n");
                                exit(-1);
                            }
                        }
                    }*/
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }         
        }
    }
    return 0;
}
int calculateMatrixRowWiseSIMD(char* ref, char* alt, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff){
    struct timespec time1, time2, time_diff;
    clock_gettime(CLOCK_REALTIME, &time1);
  
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    int32_t lowInitValue = -1073741824;
    __m256i lowInitValue256 = _mm256_set1_epi32(lowInitValue);
    __m256i zero256         = _mm256_setzero_si256();

    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < ncol; i++){
            currentValue += W_EXTEND;
            sw[0][i] = currentValue;
        }
        sw[1][0] = W_OPEN;
        currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw[i][0] = currentValue;
        }
        
    }
    int sw_prime[MAX_SEQ_LENGTH];
    int sw_prime_first_element[MAX_SEQ_LENGTH];
    int gap_size_h[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h[MAX_SEQ_LENGTH];
    int diagOrDown[MAX_SEQ_LENGTH]; //diag is 1, down is 0
   
    int j = 0;
    for(i = 0; i + 8 < ncol + 1; i += 8){
        _mm256_storeu_si256((__m256i*)(best_gap_v + i), lowInitValue256);
        _mm256_storeu_si256((__m256i*)(gap_size_v + i), zero256);
        _mm256_storeu_si256((__m256i*)(best_gap_h + i), lowInitValue256);
        _mm256_storeu_si256((__m256i*)(gap_size_h + i), zero256);
        _mm256_storeu_si256((__m256i*)(sw_prime + i), zero256);
    }
    j = i;
    while(j < ncol + 1){
        best_gap_v[j] = lowInitValue;
        gap_size_v[j] = 0;
        best_gap_h[j] = lowInitValue;
        gap_size_v[j] = 0;
        sw_prime[j] = 0;
        j++;
    }
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw_prime_first_element[0] = 0;
        sw_prime_first_element[1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw_prime_first_element[i] = currentValue;
        }
    }
    else{
        for(i = 0; i < nrow; i++){
            sw_prime_first_element[i] = 0;
        }
    }
 
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int step_down = 0;
    int step_right = 0;
    int prev_gap = 0;
    int kd = 0;
    int ki = 0;
    __m256i w_open_256 = _mm256_set1_epi32((int)W_OPEN);
    __m256i w_extend_256 = _mm256_set1_epi32((int)W_EXTEND);
    __m256i w_match_256 = _mm256_set1_epi32((int)W_MATCH);
    __m256i w_mismatch_256 = _mm256_set1_epi32((int)W_MISMATCH);
    __m256i one_256 = _mm256_set1_epi32(1);
    int* altInt = (int*)malloc(ncol * sizeof(int));
    for(i = 0; i < ncol; ++i){
        altInt[i] = alt[i];
    }

    for(i = 1; i < nrow; i++){
        a_base = ref[i - 1];
        __m256i a_base_256 = _mm256_set1_epi32(a_base);
        for(j = 1; j + 8 < ncol; j += 8){
                
            __m256i tmp0 = _mm256_loadu_si256((__m256i*)(altInt + j - 1)); //b_base
            __m256i tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            __m256i tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(&(sw[i-1][j-1])));

            __m256i step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(&(sw[i-1][j])));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + j));
            __m256i best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            __m256i gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + j)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + j), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + j), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + j), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + j), tmp0);
            _mm256_storeu_si256((__m256i*)(&(btrack[i][j])), tmp1);

        }
        while(j < ncol){
            b_base = alt[j - 1];
            step_diag =  sw[i - 1][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }
            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            if(step_diag >= step_down){
                sw_prime[j] = max(step_diag, MATRIX_MIN_CUTOFF);
                btrack[i][j] = 0;
                diagOrDown[j] = 1;
            }
            else{
                sw_prime[j] = max(step_down, MATRIX_MIN_CUTOFF);
                btrack[i][j] = kd;
                diagOrDown[j] = 0;
            }
            j++;
        }

        sw_prime[0] = sw_prime_first_element[i];   
        int offset = W_OPEN + (ncol - 2) * W_EXTEND;
        __m256i offset_256 = _mm256_set1_epi32(offset);
        __m256i w_extend_ladder_256 = _mm256_setr_epi32(0, W_EXTEND, 2 * W_EXTEND, 3 * W_EXTEND, 4 * W_EXTEND, 5 * W_EXTEND, 6 * W_EXTEND, 7 * W_EXTEND);
        __m256i w_extend_8fold_256 = _mm256_set1_epi32(8 * W_EXTEND);
        for(j = 1; j + 8 < ncol; j += 8){
            __m256i sw_prime_lastCol_256 = _mm256_loadu_si256((__m256i*)(sw_prime + j - 1));
            sw_prime_lastCol_256 = _mm256_add_epi32(sw_prime_lastCol_256, offset_256);
            sw_prime_lastCol_256 = _mm256_sub_epi32(sw_prime_lastCol_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + j), sw_prime_lastCol_256);
        }
        int offset_local = W_OPEN + (ncol - j - 1) * W_EXTEND;
        while(j < ncol){
            best_gap_h[j] = sw_prime[j - 1] + offset_local;
            offset_local -= W_EXTEND;
            j++;
        }
        int cur_best = lowInitValue;
        int prev_gap_size_h = gap_size_h[0];
        for(j = 1; j < ncol; j++){
            if(best_gap_h[j] > cur_best){
                cur_best = best_gap_h[j];
                gap_size_h[j] = 1;
            }
            else{
                best_gap_h[j] = cur_best;
                gap_size_h[j] = prev_gap_size_h + 1;
            }
            prev_gap_size_h = gap_size_h[j];
        }

        offset = (ncol - 2) * W_EXTEND;
        offset_256 = _mm256_set1_epi32(offset);
        w_extend_ladder_256 = _mm256_setr_epi32(0, W_EXTEND, 2 * W_EXTEND, 3 * W_EXTEND, 4 * W_EXTEND, 5 * W_EXTEND, 6 * W_EXTEND, 7 * W_EXTEND);
        for(j = 1; j + 8 < ncol; j += 8){
            __m256i best_gap_h_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + j));
            best_gap_h_256 = _mm256_sub_epi32(best_gap_h_256, offset_256);
            best_gap_h_256 = _mm256_add_epi32(best_gap_h_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + j), best_gap_h_256);
            
        }
        offset_local = (ncol - j - 1) * W_EXTEND;
        while(j < ncol){
            best_gap_h[j] = best_gap_h[j] - offset_local;
            offset_local -= W_EXTEND;
            j++;
        }
            

        for(j = 1; j + 8 < ncol; j += 8){
            __m256i step_right_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + j));
            __m256i ki_256 = _mm256_loadu_si256((__m256i*)(gap_size_h + j));
            __m256i ki_256_neg = _mm256_sign_epi32(ki_256, lowInitValue256);
            __m256i diagOrDown_256 = _mm256_loadu_si256((__m256i*)(diagOrDown + j));
            __m256i sw_prime_256 = _mm256_loadu_si256((__m256i*)(sw_prime + j));
            __m256i right_gt_prime = _mm256_cmpgt_epi32(step_right_256, sw_prime_256);
            __m256i right_eq_prime = _mm256_cmpeq_epi32(step_right_256, sw_prime_256);
            __m256i right_ge_prime = _mm256_or_si256(right_gt_prime, right_eq_prime);
            __m256i right_cmp_prime = _mm256_blendv_epi8(right_ge_prime, right_gt_prime, diagOrDown_256);
            __m256i sw_curRow_curCol_256 = _mm256_blendv_epi8(sw_prime_256, step_right_256, right_cmp_prime);

            _mm256_maskstore_epi32(&btrack[i][j], right_cmp_prime, ki_256_neg);
            
            _mm256_storeu_si256((__m256i*)(&(sw[i][j])), sw_curRow_curCol_256);
        }
        while(j < ncol){
            step_right = best_gap_h[j];
            ki = gap_size_h[j];
            
            if(diagOrDown[j]){
                if(step_right > sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }
            else{
                if(step_right >= sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }
            j++;
        }
    }
     clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    calMatrix_C_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);

   return 0;
}
int calculateMatrixRowWiseSIMDUnroll4x(char* ref, char* alt, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff, int w_match, int w_mismatch, int w_open, int w_extend){
    struct timespec time1, time2, time_diff;
    clock_gettime(CLOCK_REALTIME, &time1);
  
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    int32_t lowInitValue = -1073741824;
    __m256i lowInitValue256 = _mm256_set1_epi32(lowInitValue);
    __m256i zero256         = _mm256_setzero_si256();

    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = w_open;
        int currentValue = w_open;
        for(i = 2; i < ncol; i++){
            currentValue += w_extend;
            sw[0][i] = currentValue;
        }
        sw[1][0] = w_open;
        currentValue = w_open;
        for(i = 2; i < nrow; i++){
            currentValue += w_extend;
            sw[i][0] = currentValue;
        }
        
    }
    int sw_prime[MAX_SEQ_LENGTH];
    int sw_prime_first_element[MAX_SEQ_LENGTH];
    int gap_size_h[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h[MAX_SEQ_LENGTH];
    int diagOrDown[MAX_SEQ_LENGTH]; //diag is 1, down is 0
   
    int j = 0;
    for(i = 0; i + 8 < ncol + 1; i += 8){
        _mm256_storeu_si256((__m256i*)(best_gap_v + i), lowInitValue256);
        _mm256_storeu_si256((__m256i*)(gap_size_v + i), zero256);
        _mm256_storeu_si256((__m256i*)(best_gap_h + i), lowInitValue256);
        _mm256_storeu_si256((__m256i*)(gap_size_h + i), zero256);
        _mm256_storeu_si256((__m256i*)(sw_prime + i), zero256);
    }
    j = i;
    while(j < ncol + 1){
        best_gap_v[j] = lowInitValue;
        gap_size_v[j] = 0;
        best_gap_h[j] = lowInitValue;
        gap_size_v[j] = 0;
        sw_prime[j] = 0;
        j++;
    }

    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw_prime_first_element[0] = 0;
        sw_prime_first_element[1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw_prime_first_element[i] = currentValue;
        }
    }
    else{
        for(i = 0; i < nrow; i++){
            sw_prime_first_element[i] = 0;
        }
    }
     
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int step_down = 0;
    int step_right = 0;
    int prev_gap = 0;
    int kd = 0;
    int ki = 0;
    __m256i w_open_256 = _mm256_set1_epi32((int)w_open);
    __m256i w_extend_256 = _mm256_set1_epi32((int)w_extend);
    __m256i w_match_256 = _mm256_set1_epi32((int)w_match);
    __m256i w_mismatch_256 = _mm256_set1_epi32((int)w_mismatch);
    __m256i one_256 = _mm256_set1_epi32(1);
    int* altInt = (int*)malloc(ncol * sizeof(int));
    for(i = 0; i < ncol; ++i){
        altInt[i] = alt[i];
    }
    int i_prev = 0;
    for(i = 1; i < nrow; i++){
        a_base = ref[i_prev];
        __m256i a_base_256 = _mm256_set1_epi32(a_base);
        int* sw_lastRow = sw[i_prev];
        int* sw_curRow = sw[i];
        int* btrack_curRow = btrack[i];
           
        for(j = 1; j + 64 < ncol; j += 64){
                
            int jj = j + 8;
            int jjj = j + 16;
            int jjjj = j + 24;
            int j4 = j + 32;
            int jj4 = j +40;
            int jjj4 = j + 48;
            int jjjj4 = j + 56;
            __m256i tmp0 = _mm256_loadu_si256((__m256i*)(altInt + j - 1)); //b_base
            __m256i tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            __m256i tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + j - 1));

            __m256i step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + j));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + j));
            __m256i best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            __m256i gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + j)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + j), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + j), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + j), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + j), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + j), tmp1);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + jj - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
             tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jj - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jj));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + jj));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + jj)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + jj), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + jj), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[jj:jj+7]
            //now best_v_curCol_256 is step_down[jj:jj+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + jj), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + jj), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + jj), tmp1);

            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + jjj - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
             tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjj - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjj));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + jjj));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + jjj)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + jjj), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + jjj), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[jjj:jjj+7]
            //now best_v_curCol_256 is step_down[jjj:jjj+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + jjj), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + jjj), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + jjj), tmp1);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + jjjj - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
             tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjjj - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjjj));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + jjjj));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + jjjj)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + jjjj), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + jjjj), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[jjjj:jjjj+7]
            //now best_v_curCol_256 is step_down[jjjj:jjjj+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + jjjj), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + jjjj), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + jjjj), tmp1);

            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + j4 - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + j4 - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + j4));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + j4));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + j4)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + j4), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + j4), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + j4), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + j4), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + j4), tmp1);
  
            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + jj4 - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jj4 - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jj4));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + jj4));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + jj4)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + jj4), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + jj4), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + jj4), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + jj4), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + jj4), tmp1);
   
            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + jjj4 - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjj4 - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjj4));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + jjj4));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + jjj4)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + jjj4), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + jjj4), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + jjj4), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + jjj4), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + jjj4), tmp1);
 
            tmp0 = _mm256_loadu_si256((__m256i*)(altInt + jjjj4 - 1)); //b_base
            tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjjj4 - 1));

            step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + jjjj4));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + jjjj4));
            best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + jjjj4)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + jjjj4), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + jjjj4), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + jjjj4), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + jjjj4), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + jjjj4), tmp1);

        }
        while(j + 8 < ncol){
            __m256i tmp0 = _mm256_loadu_si256((__m256i*)(altInt + j - 1)); //b_base
            __m256i tmp1 = _mm256_cmpeq_epi32(a_base_256, tmp0);
            __m256i tmp2 = _mm256_blendv_epi8(w_mismatch_256, w_match_256, tmp1);
        
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + j - 1));

            __m256i step_diag_256 = _mm256_add_epi32(tmp0, tmp2);
            
            tmp0 = _mm256_loadu_si256((__m256i*)(sw_lastRow + j));
            tmp1 = _mm256_add_epi32(tmp0, w_open_256); //prev_gap_256
            tmp2 = _mm256_loadu_si256((__m256i*)(best_gap_v + j));
            __m256i best_v_curCol_256 = _mm256_add_epi32(tmp2, w_extend_256);
            tmp0 = _mm256_cmpgt_epi32(tmp1, best_v_curCol_256);
            __m256i gap_v_curCol_256 = _mm256_loadu_si256((__m256i*)(gap_size_v + j)); 
            best_v_curCol_256 = _mm256_blendv_epi8(best_v_curCol_256, tmp1, tmp0);
            gap_v_curCol_256 = _mm256_add_epi32(gap_v_curCol_256, one_256);
            gap_v_curCol_256 = _mm256_blendv_epi8(gap_v_curCol_256, one_256, tmp0);
            _mm256_storeu_si256((__m256i*)(gap_size_v + j), gap_v_curCol_256);
            _mm256_storeu_si256((__m256i*)(best_gap_v + j), best_v_curCol_256);
            //now gap_v_curCol_256 is kd[j:j+7]
            //now best_v_curCol_256 is step_down[j:j+7]
            tmp0 = _mm256_cmpgt_epi32(step_diag_256, best_v_curCol_256);
            tmp1 = _mm256_cmpeq_epi32(step_diag_256, best_v_curCol_256);
            tmp2 = _mm256_or_si256(tmp0, tmp1);
            tmp0 = _mm256_blendv_epi8(best_v_curCol_256, step_diag_256, tmp2);
            tmp1 = _mm256_blendv_epi8(gap_v_curCol_256, zero256, tmp2);
            _mm256_storeu_si256((__m256i*)(diagOrDown + j), tmp2);
            _mm256_storeu_si256((__m256i*)(sw_prime + j), tmp0);
            _mm256_storeu_si256((__m256i*)(btrack_curRow + j), tmp1);
            

            j += 8;
        }
        while(j < ncol){
            b_base = alt[j - 1];
            step_diag =  sw[i_prev][j - 1] + wd(a_base, b_base, w_match, w_mismatch);
            prev_gap = sw[i_prev][j] + w_open;
            best_gap_v[j] += w_extend;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }
            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            if(step_diag >= step_down){
                sw_prime[j] = max(step_diag, MATRIX_MIN_CUTOFF);
                btrack[i][j] = 0;
                diagOrDown[j] = 1;
            }
            else{
                sw_prime[j] = max(step_down, MATRIX_MIN_CUTOFF);
                btrack[i][j] = kd;
                diagOrDown[j] = 0;
            }
            j++;
        }
        sw_prime[0] = sw_prime_first_element[i];   
        int offset = w_open + (ncol - 2) * w_extend;
        __m256i offset_256 = _mm256_set1_epi32(offset);
        __m256i w_extend_ladder_256 = _mm256_setr_epi32(0, w_extend, 2 * w_extend, 3 * w_extend, 4 * w_extend, 5 * w_extend, 6 * w_extend, 7 * w_extend);
        __m256i w_extend_8fold_256 = _mm256_set1_epi32(8 * w_extend);
        for(j = 1; j + 32 < ncol; j += 32){
            int jj = j + 8;
            int jjj = j + 16;
            int jjjj = j + 24;
            __m256i sw_prime_lastCol_256 = _mm256_loadu_si256((__m256i*)(sw_prime + j - 1));
            sw_prime_lastCol_256 = _mm256_add_epi32(sw_prime_lastCol_256, offset_256);
            sw_prime_lastCol_256 = _mm256_sub_epi32(sw_prime_lastCol_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + j), sw_prime_lastCol_256);
            sw_prime_lastCol_256 = _mm256_loadu_si256((__m256i*)(sw_prime + jj - 1));
            sw_prime_lastCol_256 = _mm256_add_epi32(sw_prime_lastCol_256, offset_256);
            sw_prime_lastCol_256 = _mm256_sub_epi32(sw_prime_lastCol_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + jj), sw_prime_lastCol_256);
            sw_prime_lastCol_256 = _mm256_loadu_si256((__m256i*)(sw_prime + jjj - 1));
            sw_prime_lastCol_256 = _mm256_add_epi32(sw_prime_lastCol_256, offset_256);
            sw_prime_lastCol_256 = _mm256_sub_epi32(sw_prime_lastCol_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + jjj), sw_prime_lastCol_256);
            sw_prime_lastCol_256 = _mm256_loadu_si256((__m256i*)(sw_prime + jjjj - 1));
            sw_prime_lastCol_256 = _mm256_add_epi32(sw_prime_lastCol_256, offset_256);
            sw_prime_lastCol_256 = _mm256_sub_epi32(sw_prime_lastCol_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + jjjj), sw_prime_lastCol_256);
        


        }
        while(j + 8 < ncol){
             __m256i sw_prime_lastCol_256 = _mm256_loadu_si256((__m256i*)(sw_prime + j - 1));
            sw_prime_lastCol_256 = _mm256_add_epi32(sw_prime_lastCol_256, offset_256);
            sw_prime_lastCol_256 = _mm256_sub_epi32(sw_prime_lastCol_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + j), sw_prime_lastCol_256);
            j += 8;
        }
        int offset_local = w_open + (ncol - j - 1) * w_extend;
        while(j < ncol){
            best_gap_h[j] = sw_prime[j - 1] + offset_local;
            offset_local -= w_extend;
            j++;
        }
        int cur_best = lowInitValue;
        int prev_gap_size_h = gap_size_h[0];
        for(j = 1; j+4 < ncol; j+=4){
            int jj = j + 1;
            int jjj = j + 2;
            int jjjj = j + 3;
             
            /*         int best_gap_h_0 = best_gap_h[j];
            int best_gap_h_1 = best_gap_h[jj];
            int best_gap_h_2 = best_gap_h[jjj];
            int best_gap_h_3 = best_gap_h[jjjj];

            int gap_size_h_0 = 1;
            int gap_size_h_1 = 1;
            int gap_size_h_2 = 1;
            int gap_size_h_3 = 1;


            if(best_gap_h_0 <= cur_best){
                best_gap_h_0 = cur_best;
                gap_size_h_0 = prev_gap_size_h + 1;
            }

            if(best_gap_h_1 <= best_gap_h_0){
                best_gap_h_1 = best_gap_h_0;
                gap_size_h_1 = gap_size_h_0 + 1;
            }

            if(best_gap_h_3 <= best_gap_h_2){
                best_gap_h_3 = best_gap_h_2;
                gap_size_h_3 = 0; //gap_size_h_2 + 1, else depends the comparison between best_gap_h_23 and best_gap_h_01
            }
            
            if(best_gap_h_2 <= best_gap_h_1){
                best_gap_h_2 = best_gap_h_1;
                gap_size_h_2 = gap_size_h_1 + 1;
            }
            int tmp_gap_size_h_3 = 1;
            if(best_gap_h_3 <= best_gap_h_1){
                best_gap_h_3 = best_gap_h_1;
                tmp_gap_size_h_3 = gap_size_h_1 + 2;
            }
            if(gap_size_h_3 == 0){
                gap_size_h_3 = gap_size_h_2 + 1;
            }
            else{
                gap_size_h_3 = tmp_gap_size_h_3;
            }
            cur_best = best_gap_h_3;
            prev_gap_size_h = gap_size_h_3;
            
            best_gap_h[j] = best_gap_h_0;
            best_gap_h[jj] = best_gap_h_1;
            best_gap_h[jjj] = best_gap_h_2;
            best_gap_h[jjjj] = best_gap_h_3;

            gap_size_h[j] = gap_size_h_0;
            gap_size_h[jj] = gap_size_h_1;
            gap_size_h[jjj] = gap_size_h_2;
            gap_size_h[jjjj] = gap_size_h_3;
*/


            if(best_gap_h[j] > cur_best){
                cur_best = best_gap_h[j];
                gap_size_h[j] = 1;
            }
            else{
                best_gap_h[j] = cur_best;
                gap_size_h[j] = prev_gap_size_h + 1;
            }
            prev_gap_size_h = gap_size_h[j];
            if(best_gap_h[jj] > cur_best){
                cur_best = best_gap_h[jj];
                gap_size_h[jj] = 1;
            }
            else{
                best_gap_h[jj] = cur_best;
                gap_size_h[jj] = prev_gap_size_h + 1;
            }
            prev_gap_size_h = gap_size_h[jj];
            if(best_gap_h[jjj] > cur_best){
                cur_best = best_gap_h[jjj];
                gap_size_h[jjj] = 1;
            }
            else{
                best_gap_h[jjj] = cur_best;
                gap_size_h[jjj] = prev_gap_size_h + 1;
            }
            prev_gap_size_h = gap_size_h[jjj];
            if(best_gap_h[jjjj] > cur_best){
                cur_best = best_gap_h[jjjj];
                gap_size_h[jjjj] = 1;
            }
            else{
                best_gap_h[jjjj] = cur_best;
                gap_size_h[jjjj] = prev_gap_size_h + 1;
            }
            prev_gap_size_h = gap_size_h[jjjj];
        
       }
        while(j < ncol){
            if(best_gap_h[j] > cur_best){
                cur_best = best_gap_h[j];
                gap_size_h[j] = 1;
            }
            else{
                best_gap_h[j] = cur_best;
                gap_size_h[j] = prev_gap_size_h + 1;
            }
            prev_gap_size_h = gap_size_h[j];
            j++;
        }

        offset = (ncol - 2) * w_extend;
        offset_256 = _mm256_set1_epi32(offset);
        w_extend_ladder_256 = _mm256_setr_epi32(0, w_extend, 2 * w_extend, 3 * w_extend, 4 * w_extend, 5 * w_extend, 6 * w_extend, 7 * w_extend);
        for(j = 1; j + 32 < ncol; j += 32){
            int jj = j + 8;
            int jjj = j + 16;
            int jjjj = j + 24;
            __m256i best_gap_h_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + j));
            best_gap_h_256 = _mm256_sub_epi32(best_gap_h_256, offset_256);
            best_gap_h_256 = _mm256_add_epi32(best_gap_h_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + j), best_gap_h_256);
            best_gap_h_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + jj));
            best_gap_h_256 = _mm256_sub_epi32(best_gap_h_256, offset_256);
            best_gap_h_256 = _mm256_add_epi32(best_gap_h_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + jj), best_gap_h_256);
            best_gap_h_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + jjj));
            best_gap_h_256 = _mm256_sub_epi32(best_gap_h_256, offset_256);
            best_gap_h_256 = _mm256_add_epi32(best_gap_h_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + jjj), best_gap_h_256);
            best_gap_h_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + jjjj));
            best_gap_h_256 = _mm256_sub_epi32(best_gap_h_256, offset_256);
            best_gap_h_256 = _mm256_add_epi32(best_gap_h_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + jjjj), best_gap_h_256);
        }
        while(j + 8 < ncol){
             __m256i best_gap_h_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + j));
            best_gap_h_256 = _mm256_sub_epi32(best_gap_h_256, offset_256);
            best_gap_h_256 = _mm256_add_epi32(best_gap_h_256, w_extend_ladder_256);
            w_extend_ladder_256 = _mm256_add_epi32(w_extend_ladder_256, w_extend_8fold_256);
            _mm256_storeu_si256((__m256i*)(best_gap_h + j), best_gap_h_256);
            j += 8;
        }
        offset_local = (ncol - j - 1) * w_extend;
        while(j < ncol){
            best_gap_h[j] = best_gap_h[j] - offset_local;
            offset_local -= w_extend;
            j++;
        }
            

        for(j = 1; j + 32 < ncol; j += 32){
            int jj = j + 8;
            int jjj = j + 16;
            int jjjj = j + 24;
            __m256i step_right_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + j));
            __m256i ki_256 = _mm256_loadu_si256((__m256i*)(gap_size_h + j));
            __m256i ki_256_neg = _mm256_sign_epi32(ki_256, lowInitValue256);
            __m256i diagOrDown_256 = _mm256_loadu_si256((__m256i*)(diagOrDown + j));
            __m256i sw_prime_256 = _mm256_loadu_si256((__m256i*)(sw_prime + j));
            __m256i right_gt_prime = _mm256_cmpgt_epi32(step_right_256, sw_prime_256);
            __m256i right_eq_prime = _mm256_cmpeq_epi32(step_right_256, sw_prime_256);
            __m256i right_ge_prime = _mm256_or_si256(right_gt_prime, right_eq_prime);
            __m256i right_cmp_prime = _mm256_blendv_epi8(right_ge_prime, right_gt_prime, diagOrDown_256);
            __m256i sw_curRow_curCol_256 = _mm256_blendv_epi8(sw_prime_256, step_right_256, right_cmp_prime);

            _mm256_maskstore_epi32(btrack_curRow + j, right_cmp_prime, ki_256_neg);
            
            _mm256_storeu_si256((__m256i*)(sw_curRow + j), sw_curRow_curCol_256);
        
            step_right_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + jj));
            ki_256 = _mm256_loadu_si256((__m256i*)(gap_size_h + jj));
            ki_256_neg = _mm256_sign_epi32(ki_256, lowInitValue256);
            diagOrDown_256 = _mm256_loadu_si256((__m256i*)(diagOrDown + jj));
            sw_prime_256 = _mm256_loadu_si256((__m256i*)(sw_prime + jj));
            right_gt_prime = _mm256_cmpgt_epi32(step_right_256, sw_prime_256);
            right_eq_prime = _mm256_cmpeq_epi32(step_right_256, sw_prime_256);
            right_ge_prime = _mm256_or_si256(right_gt_prime, right_eq_prime);
            right_cmp_prime = _mm256_blendv_epi8(right_ge_prime, right_gt_prime, diagOrDown_256);
            sw_curRow_curCol_256 = _mm256_blendv_epi8(sw_prime_256, step_right_256, right_cmp_prime);

            _mm256_maskstore_epi32(btrack_curRow + jj, right_cmp_prime, ki_256_neg);
            
            _mm256_storeu_si256((__m256i*)(sw_curRow + jj), sw_curRow_curCol_256);
            
            step_right_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + jjj));
            ki_256 = _mm256_loadu_si256((__m256i*)(gap_size_h + jjj));
            ki_256_neg = _mm256_sign_epi32(ki_256, lowInitValue256);
            diagOrDown_256 = _mm256_loadu_si256((__m256i*)(diagOrDown + jjj));
            sw_prime_256 = _mm256_loadu_si256((__m256i*)(sw_prime + jjj));
            right_gt_prime = _mm256_cmpgt_epi32(step_right_256, sw_prime_256);
            right_eq_prime = _mm256_cmpeq_epi32(step_right_256, sw_prime_256);
            right_ge_prime = _mm256_or_si256(right_gt_prime, right_eq_prime);
            right_cmp_prime = _mm256_blendv_epi8(right_ge_prime, right_gt_prime, diagOrDown_256);
            sw_curRow_curCol_256 = _mm256_blendv_epi8(sw_prime_256, step_right_256, right_cmp_prime);

            _mm256_maskstore_epi32(btrack_curRow + jjj, right_cmp_prime, ki_256_neg);
            
            _mm256_storeu_si256((__m256i*)(sw_curRow + jjj), sw_curRow_curCol_256);
            
            step_right_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + jjjj));
            ki_256 = _mm256_loadu_si256((__m256i*)(gap_size_h + jjjj));
            ki_256_neg = _mm256_sign_epi32(ki_256, lowInitValue256);
            diagOrDown_256 = _mm256_loadu_si256((__m256i*)(diagOrDown + jjjj));
            sw_prime_256 = _mm256_loadu_si256((__m256i*)(sw_prime + jjjj));
            right_gt_prime = _mm256_cmpgt_epi32(step_right_256, sw_prime_256);
            right_eq_prime = _mm256_cmpeq_epi32(step_right_256, sw_prime_256);
            right_ge_prime = _mm256_or_si256(right_gt_prime, right_eq_prime);
            right_cmp_prime = _mm256_blendv_epi8(right_ge_prime, right_gt_prime, diagOrDown_256);
            sw_curRow_curCol_256 = _mm256_blendv_epi8(sw_prime_256, step_right_256, right_cmp_prime);

            _mm256_maskstore_epi32(btrack_curRow + jjjj, right_cmp_prime, ki_256_neg);
            
            _mm256_storeu_si256((__m256i*)(sw_curRow + jjjj), sw_curRow_curCol_256);
        


        }
        while(j + 8 < ncol){
            __m256i step_right_256 = _mm256_loadu_si256((__m256i*)(best_gap_h + j));
            __m256i ki_256 = _mm256_loadu_si256((__m256i*)(gap_size_h + j));
            __m256i ki_256_neg = _mm256_sign_epi32(ki_256, lowInitValue256);
            __m256i diagOrDown_256 = _mm256_loadu_si256((__m256i*)(diagOrDown + j));
            __m256i sw_prime_256 = _mm256_loadu_si256((__m256i*)(sw_prime + j));
            __m256i right_gt_prime = _mm256_cmpgt_epi32(step_right_256, sw_prime_256);
            __m256i right_eq_prime = _mm256_cmpeq_epi32(step_right_256, sw_prime_256);
            __m256i right_ge_prime = _mm256_or_si256(right_gt_prime, right_eq_prime);
            __m256i right_cmp_prime = _mm256_blendv_epi8(right_ge_prime, right_gt_prime, diagOrDown_256);
            __m256i sw_curRow_curCol_256 = _mm256_blendv_epi8(sw_prime_256, step_right_256, right_cmp_prime);

            _mm256_maskstore_epi32(btrack_curRow + j, right_cmp_prime, ki_256_neg);
            
            _mm256_storeu_si256((__m256i*)(sw_curRow + j), sw_curRow_curCol_256);
        

            j += 8;
        }
        while(j < ncol){
            step_right = best_gap_h[j];
            ki = gap_size_h[j];
            
            if(diagOrDown[j]){
                if(step_right > sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }
            else{
                if(step_right >= sw_prime[j]){
                    sw[i][j] = max(step_right, MATRIX_MIN_CUTOFF);
                    btrack[i][j] = -ki;
                }
                else{
                    sw[i][j] = sw_prime[j];
                }
            }
            j++;
        }
        i_prev = i;
    }
     clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    calMatrix_C_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);

   return 0;
}
int calculateMatrixOneBatch(char* reference, char* alternate, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff){
    struct timespec time1, time2, time_diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    signed int lowInitValue = -1073741824;
    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int gap_size_h[MAX_SEQ_LENGTH];
   /* int* best_gap_v;
    int* best_gap_h;
    int* gap_size_v;
    int* gap_size_h;
    if(!(best_gap_v = (int*)malloc((ncol + 1) * sizeof(int)))){
        return -1;
    }
    if(!(best_gap_h = (int*)malloc((nrow + 1) * sizeof(int)))){
        return -1;
    }
    if(!(gap_size_v = (int*)malloc((ncol + 1) * sizeof(int)))){
        return -1;
    }
    if(!(gap_size_h = (int*)malloc((nrow + 1) * sizeof(int)))){
        return -1;
    }*/


    for(i = 0; i < ncol + 1; ++i){
        best_gap_v[i] = lowInitValue;
        gap_size_v[i] = 0;
    }
    for(i = 0; i < nrow + 1; ++i){
        best_gap_h[i] = lowInitValue;
        gap_size_h[i] = 0;
    }
    
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < ncol; i++){
            currentValue += W_EXTEND;
            sw[0][i] = currentValue;
        }
        sw[1][0] = W_OPEN;
        currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw[i][0] = currentValue;
        }
        
    }


    //build smith-waterman matrix and keep backtrack info
    int curRow_id = 0;
    int lastRow_id = 0;
    int curBackTrackRow_id = 0;
    int j = 0;
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int prev_gap = 0;
    int step_down = 0;
    int kd = 0;
    int step_right = 0;
    int ki = 0;
    int diagHighestOrEqual = 0; 
    for(i = 1; i < nrow; ++i){
        a_base = reference[i - 1];
        lastRow_id = curRow_id;
        curRow_id = i;
        curBackTrackRow_id = i;
        for(j = 1; j< ncol; ++j){
            b_base = alternate[j - 1];
            
            step_diag = sw[lastRow_id][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[lastRow_id][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }

            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            
            prev_gap = sw[curRow_id][j - 1] + W_OPEN;
            best_gap_h[i] += W_EXTEND;
            if(prev_gap > best_gap_h[i]){
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            }
            else{
                gap_size_h[i]++;
            }
            step_right = best_gap_h[i];
            ki = gap_size_h[i];
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            
            if(diagHighestOrEqual){
                sw[curRow_id][j] = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[curBackTrackRow_id][j] = 0;
            }
            else if(step_right >= step_down){
                sw[curRow_id][j] = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[curBackTrackRow_id][j] = -ki;
            }
            else{
                sw[curRow_id][j] = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[curBackTrackRow_id][j] = kd;
            }
            //printf("golden sw[%d][%d] = %d, btrack[%d][%d] = %d\n", curRow_id,j,sw[curRow_id][j],curBackTrackRow_id,j, btrack[curBackTrackRow_id][j]);
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    calMatrix_C_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);

//free(best_gap_h);
    //free(best_gap_v);
    //free(gap_size_h);
    //free(gap_size_v);
    return 0; 
}

int calculateMatrixOneBatchUnroll2x(char* reference, char* alternate, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff){
    struct timespec time1, time2, time_diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    signed int lowInitValue = -1073741824;
    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int gap_size_h[MAX_SEQ_LENGTH];
   
    for(i = 0; i < ncol + 1; ++i){
        best_gap_v[i] = lowInitValue;
        gap_size_v[i] = 0;
    }
    for(i = 0; i < nrow + 1; ++i){
        best_gap_h[i] = lowInitValue;
        gap_size_h[i] = 0;
    }
    
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < ncol; i++){
            currentValue += W_EXTEND;
            sw[0][i] = currentValue;
        }
        sw[1][0] = W_OPEN;
        currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw[i][0] = currentValue;
        }
        
    }


    //build smith-waterman matrix and keep backtrack info
    int j = 0;
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int prev_gap = 0;
    int step_down = 0;
    int kd = 0;
    int step_right = 0;
    int ki = 0;
    int diagHighestOrEqual = 0; 
    for(i = 1; i < nrow; ++i){
        a_base = reference[i - 1];
        for(j = 1; j + 1 < ncol; j += 2){
            //inner loop
            //read:
            //  alternate[j - 1 : j], sw[i - 1][j - 1 : j], sw[i][j - 1 : j]
            //  best_gap_v[j : j + 1], best_gap_h[i]
            //  gap_size_v[j : j + 1], gap_size_h[i]
            //write
            //  btrack[i][j : j + 1]
            //  sw[i][j : j + 1]

            b_base = alternate[j - 1];
            
            step_diag = sw[i - 1][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }

            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            
            prev_gap = sw[i][j - 1] + W_OPEN;
            best_gap_h[i] += W_EXTEND;
            if(prev_gap > best_gap_h[i]){
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            }
            else{
                gap_size_h[i]++;
            }
            step_right = best_gap_h[i];
            ki = gap_size_h[i];
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            
            int sw_i_j;
            if(diagHighestOrEqual){
                sw_i_j = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j] = 0;
            }
            else if(step_right >= step_down){
                sw_i_j = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j] = -ki;
            }
            else{
                sw_i_j = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j] = kd;
            }
            sw[i][j] = sw_i_j;
            b_base = alternate[j];
            step_diag = sw[i - 1][j] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j + 1] + W_OPEN;
            best_gap_v[j + 1] += W_EXTEND;
            if(prev_gap > best_gap_v[j + 1]){
                best_gap_v[j + 1] = prev_gap;
                gap_size_v[j + 1] = 1;
            }
            else
                gap_size_v[j + 1]++;
            step_down = best_gap_v[j + 1];
            kd = gap_size_v[j + 1];
            
            prev_gap = sw_i_j + W_OPEN;
            best_gap_h[i] += W_EXTEND;
            if(prev_gap > best_gap_h[i]){
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            }
            else{
                gap_size_h[i]++;
            }
            step_right = best_gap_h[i];
            ki = gap_size_h[i];
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            if(diagHighestOrEqual){
                sw[i][j + 1] = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j + 1] = 0;
            }
            else if(step_right >= step_down){
                sw[i][j + 1] = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j + 1] = -ki;
            }
            else{
                sw[i][j + 1] = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j + 1] = kd;
            }
        
        }
        while(j < ncol){
            b_base = alternate[j - 1];
        
            step_diag = sw[i - 1][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap; gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }

            step_down = best_gap_v[j];
            kd = gap_size_v[j];
        
            prev_gap = sw[i][j - 1] + W_OPEN;
            best_gap_h[i] += W_EXTEND;
            if(prev_gap > best_gap_h[i]){
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            }
            else{
                gap_size_h[i]++;
            }
            step_right = best_gap_h[i];
            ki = gap_size_h[i];
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
        
            if(diagHighestOrEqual){
                sw[i][j] = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j] = 0;
            }
            else if(step_right >= step_down){
                sw[i][j] = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j] = -ki;
            }
            else{
                sw[i][j] = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j] = kd;
            }
            j++;
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    calMatrix_C_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);

//free(best_gap_h);
    //free(best_gap_v);
    //free(gap_size_h);
    //free(gap_size_v);
    return 0; 
}

int calculateMatrixOneBatchUnroll4x(char* reference, char* alternate, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy, int cutoff){
    struct timespec time1, time2, time_diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int i = 0;
    signed int lowInitValue = -1073741824;
    int best_gap_v[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
   
    for(i = 0; i < ncol + 1; ++i){
        best_gap_v[i] = lowInitValue;
        gap_size_v[i] = 0;
    }
   
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = W_OPEN;
        int currentValue = W_OPEN;
        for(i = 2; i < ncol; i++){
            currentValue += W_EXTEND;
            sw[0][i] = currentValue;
        }
        sw[1][0] = W_OPEN;
        currentValue = W_OPEN;
        for(i = 2; i < nrow; i++){
            currentValue += W_EXTEND;
            sw[i][0] = currentValue;
        }
        
    }


    //build smith-waterman matrix and keep backtrack info
    int j = 0;
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int prev_gap = 0;
    int step_down = 0;
    int kd = 0;
    int step_right = 0;
    int ki = 0;
    int diagHighestOrEqual = 0; 
    for(i = 1; i < nrow; ++i){
        a_base = reference[i - 1];
        int best_gap_h_i = lowInitValue;
        int gap_size_h_i = 0;
        for(j = 1; j + 3 < ncol; j += 4){
            b_base = alternate[j - 1];
            
            step_diag = sw[i - 1][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }

            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            
            prev_gap = sw[i][j - 1] + W_OPEN;
            best_gap_h_i += W_EXTEND;
            if(prev_gap > best_gap_h_i){
                best_gap_h_i = prev_gap;
                gap_size_h_i = 1;
            }
            else{
                gap_size_h_i++;
            }
            step_right = best_gap_h_i;
            ki = gap_size_h_i;
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            
            int sw_i_j;
            if(diagHighestOrEqual){
                sw_i_j = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j] = 0;
            }
            else if(step_right >= step_down){
                sw_i_j = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j] = -ki;
            }
            else{
                sw_i_j = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j] = kd;
            }
            sw[i][j] = sw_i_j;
            b_base = alternate[j];
            step_diag = sw[i - 1][j] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j + 1] + W_OPEN;
            best_gap_v[j + 1] += W_EXTEND;
            if(prev_gap > best_gap_v[j + 1]){
                best_gap_v[j + 1] = prev_gap;
                gap_size_v[j + 1] = 1;
            }
            else
                gap_size_v[j + 1]++;
            step_down = best_gap_v[j + 1];
            kd = gap_size_v[j + 1];
            
            prev_gap = sw_i_j + W_OPEN;
            best_gap_h_i += W_EXTEND;
            if(prev_gap > best_gap_h_i){
                best_gap_h_i = prev_gap;
                gap_size_h_i= 1;
            }
            else{
                gap_size_h_i++;
            }
            step_right = best_gap_h_i;
            ki = gap_size_h_i;
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            int sw_i_j1;
            if(diagHighestOrEqual){
                sw_i_j1 = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j + 1] = 0;
            }
            else if(step_right >= step_down){
                sw_i_j1 = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j + 1] = -ki;
            }
            else{
                sw_i_j1 = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j + 1] = kd;
            }

            sw[i][j + 1] = sw_i_j1; 

            b_base = alternate[j + 1];
            step_diag = sw[i - 1][j + 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j + 2] + W_OPEN;
            best_gap_v[j + 2] += W_EXTEND;
            if(prev_gap > best_gap_v[j + 2]){
                best_gap_v[j + 2] = prev_gap;
                gap_size_v[j + 2] = 1;
            }
            else
                gap_size_v[j + 2]++;
            step_down = best_gap_v[j + 2];
            kd = gap_size_v[j + 2];
            
            prev_gap = sw_i_j1 + W_OPEN;
            best_gap_h_i += W_EXTEND;
            if(prev_gap > best_gap_h_i){
                best_gap_h_i = prev_gap;
                gap_size_h_i = 1;
            }
            else{
                gap_size_h_i++;
            }
            step_right = best_gap_h_i;
            ki = gap_size_h_i;
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            int sw_i_j2;
            if(diagHighestOrEqual){
                sw_i_j2 = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j + 2] = 0;
            }
            else if(step_right >= step_down){
                sw_i_j2 = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j + 2] = -ki;
            }
            else{
                sw_i_j2 = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j + 2] = kd;
            }
            sw[i][j + 2] = sw_i_j2;

            b_base = alternate[j + 2];
            step_diag = sw[i - 1][j + 2] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j + 3] + W_OPEN;
            best_gap_v[j + 3] += W_EXTEND;
            if(prev_gap > best_gap_v[j + 3]){
                best_gap_v[j + 3] = prev_gap;
                gap_size_v[j + 3] = 1;
            }
            else
                gap_size_v[j + 3]++;
            step_down = best_gap_v[j + 3];
            kd = gap_size_v[j + 3];
            
            prev_gap = sw_i_j2 + W_OPEN;
            best_gap_h_i+= W_EXTEND;
            if(prev_gap > best_gap_h_i){
                best_gap_h_i = prev_gap;
                gap_size_h_i = 1;
            }
            else{
                gap_size_h_i++;
            }
            step_right = best_gap_h_i;
            ki = gap_size_h_i;
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            if(diagHighestOrEqual){
                sw[i][j + 3] = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j + 3] = 0;
            }
            else if(step_right >= step_down){
                sw[i][j + 3] = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j + 3] = -ki;
            }
            else{
                sw[i][j + 3] = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j + 3] = kd;
            }
        }
        while(j < ncol){
            b_base = alternate[j - 1];
        
            step_diag = sw[i - 1][j - 1] + wd(a_base, b_base, W_MATCH, W_MISMATCH);
            prev_gap = sw[i - 1][j] + W_OPEN;
            best_gap_v[j] += W_EXTEND;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap; gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }

            step_down = best_gap_v[j];
            kd = gap_size_v[j];
        
            prev_gap = sw[i][j - 1] + W_OPEN;
            best_gap_h_i += W_EXTEND;
            if(prev_gap > best_gap_h_i){
                best_gap_h_i = prev_gap;
                gap_size_h_i = 1;
            }
            else{
                gap_size_h_i++;
            }
            step_right = best_gap_h_i;
            ki = gap_size_h_i;
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
        
            if(diagHighestOrEqual){
                sw[i][j] = max(MATRIX_MIN_CUTOFF, step_diag);
                btrack[i][j] = 0;
            }
            else if(step_right >= step_down){
                sw[i][j] = max(MATRIX_MIN_CUTOFF, step_right);
                btrack[i][j] = -ki;
            }
            else{
                sw[i][j] = max(MATRIX_MIN_CUTOFF, step_down);
                btrack[i][j] = kd;
            }
            j++;
        }
      //  best_gap_h[i] = best_gap_h_i;
      //  gap_size_h[i] = gap_size_h_i;
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    time_diff = diff_time(time1, time2);
    calMatrix_C_time += (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);

    return 0; 
}

inline int32_t MoveMask_256_8(__m256i x){
    int32_t x_32 = _mm256_movemask_epi8(x);
    int32_t x_8 = (x_32 | (x_32 >> 3)) & 0x03030303;
    x_8 = (x_8 | (x_8 >> 6)) & 0x000f000f;
    return (x_8 | (x_8 >> 12)) & 0x000000ff;
}

int calculateCigarOneBatch(int** sw, int** btrack, int nrow, int ncol, int overhang_strategy, struct Cigar* cigarResult, int* alignment_offset_ptr){
    int p1 = 0;
    int p2 = 0;
    int refLength = nrow - 1;
    int altLength = ncol - 1;
    int maxscore = -2147483648;
    int segment_length = 0;
    int i = 0;
    int j = 0;
    int curScore = 0;

    if(overhang_strategy == OVERHANG_STRATEGY_INDEL){
        p1 = refLength;
        p2 = altLength;
    }
    else{
        p2 = altLength;
        for(i = 1; i < nrow; ++i){
            curScore = sw[i][altLength];
            if(curScore >= maxscore){
                p1 = i;
                maxscore = curScore;
            }
        }
        //now look for a larger score on the bottom-most row
        if(overhang_strategy != OVERHANG_STRATEGY_LEADING_INDEL){
            for(j = 1; j < ncol; ++j){
                curScore = sw[refLength][j];
                if(curScore > maxscore || (curScore == maxscore && (abs(refLength - j) < abs(p1 - p2)))){
                    p1 = refLength;
                    p2 = j;
                    maxscore = curScore;
                    segment_length = altLength - j;
                }
            }
        }
    }
    
    cigarResult->CigarElementNum = 0;    
    if(segment_length > 0 && overhang_strategy == OVERHANG_STRATEGY_SOFTCLIP){
        addCigarElement(cigarResult, segment_length, STATE_CLIP);
        segment_length = 0;
    }

    int state = STATE_MATCH;
    int btr = 0;
    int new_state =0;
    int step_length = 0;
    do{
        btr = btrack[p1][p2];
        step_length = 1;
        if(btr > 0){
            new_state = STATE_DELETION;
            step_length = btr;
        }
        else if(btr < 0){
            new_state = STATE_INSERTION;
            step_length = (-btr);
        }
        else
            new_state = STATE_MATCH;

        switch(new_state){
            case STATE_MATCH: p1--; p2--; break;
            case STATE_INSERTION: p2 -= step_length; break; 
            case STATE_DELETION: p1 -= step_length; break;
        }

        if(new_state == state) segment_length += step_length;
        else{
            addCigarElement(cigarResult, segment_length, state);
            segment_length = step_length;
            state = new_state;
        }
    }while(p1 > 0 && p2 > 0);

    if(overhang_strategy == OVERHANG_STRATEGY_SOFTCLIP){
        addCigarElement(cigarResult, segment_length, state);

        if(p2 > 0){
            addCigarElement(cigarResult, p2, STATE_CLIP);
        }
        *alignment_offset_ptr = p1;
    }
    else if(overhang_strategy == OVERHANG_STRATEGY_IGNORE){
        addCigarElement(cigarResult, segment_length + p2, state);
        *alignment_offset_ptr = p1 - p2;
    }
    else{
        addCigarElement(cigarResult, segment_length, state);
 
        if(p1 > 0){
            addCigarElement(cigarResult, p1, STATE_DELETION);
        }
        else if(p2 > 0){
            addCigarElement(cigarResult, p2, STATE_INSERTION);
        }
        *alignment_offset_ptr = 0;
    }
    //reverse the cigarResult->CigarElementNum
    if(cigarResult->CigarElementNum <= 0){
        //printf("cigarResult has non pos elements %d \n", cigarResult->CigarElementNum);
        return -1;
    }

    i = 0;
    j = cigarResult->CigarElementNum - 1;
    while(j > i){
        struct CigarElement tmp;
        tmp = cigarResult->cigarElements[j];
        cigarResult->cigarElements[j] = cigarResult->cigarElements[i];
        cigarResult->cigarElements[i] = tmp;
        j--;
        i++;
    }
    return 0;
}

