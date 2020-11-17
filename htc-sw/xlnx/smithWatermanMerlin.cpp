/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: Xilinx Smith-Waterman Merlin compiler implementation 
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#define MAX_SEQ_LENGTH 512
#define MAX_BATCH_SIZE 260
#define OVERHANG_STRATEGY_SOFTCLIP 0
#define OVERHANG_STRATEGY_INDEL 1
#define OVERHANG_STRATEGY_LEADING_INDEL 2
#define OVERHANG_STRATEGY_IGNORE 3
#define W_MATCH 200
#define W_MISMATCH -150
#define W_OPEN -260
#define W_EXTEND -11
#define STATE_MATCH 0       //equal to CigarOperator.M
#define STATE_INSERTION 1   //equal to CigarOperator.I
#define STATE_DELETION 2    //equal to CigarOperator.D
#define STATE_CLIP 4        //equal to CigarOperator.S

#define PARALLEL_FACTOR 8

int smithWatermanKernel(char*, char*, int, int, int, int, int, int, int, int, short*, short*, short*, short*);

inline int smithWatermanCellWiseMiniKernel(int, int, int, int, int, int, int*, int*, int, int , int, char, char, int, int, int, int, int*, int*, int*, int*, int*, int*, int, int, int*, int*, int, short*, int, int, int, int, int, int, int, int, int, int, int, char, int);

#pragma ACCEL kernel
void smithWatermanMerlin(char* inputs, int refLength, int batchSize, int overhang_strategy, int w_match, int w_mismatch, int w_open, int w_extend, short* outputs, int* cigarInfoLength){
/*
 * in inputs, the first 2*batchSize bytes are altLengths, each 2 bytes is a short
 * The rest data is ref and alts
 * in outputs, the first batchSize of short value is cigarSize, the rest value is cigarInfo.
 * */
#pragma ACCEL interface variable=inputs depth=134152 max_depth=134152
#pragma ACCEL interface variable=outputs depth=266762 max_depth=266762
#pragma ACCEL interface variable=cigarInfoLength depth=1 max_depth=1
    int batchCount = 0;
    char ref[MAX_SEQ_LENGTH];
    short altLengths[MAX_BATCH_SIZE];
    for(int i = 0; i < batchSize; i++){
        altLengths[i] = (short)(inputs[2 * i + 1] << 8) | (signed short)(inputs[2 * i] & 0xff);
    }
    memcpy(ref, inputs + 2 * batchSize, MAX_SEQ_LENGTH * sizeof(char));
    
    int cigarElementPtr = 0;
#define KERNEL_NUM 3
    assert(batchSize >= 0 && batchSize < MAX_BATCH_SIZE);
    for(batchCount = 0; batchCount < batchSize; batchCount += KERNEL_NUM){
        int altLengthsLocal[KERNEL_NUM];
        char alt[KERNEL_NUM][MAX_SEQ_LENGTH];
        short lengths[KERNEL_NUM][MAX_SEQ_LENGTH];
        short states[KERNEL_NUM][MAX_SEQ_LENGTH];
        short cigarElementNum[KERNEL_NUM];
        short alignment_offset[KERNEL_NUM];
        for(int i = 0; i < KERNEL_NUM; ++i){
            if(batchCount + i < batchSize){
                memcpy(&alt[i][0], inputs + 2 * batchSize + (batchCount + i + 1) * MAX_SEQ_LENGTH, MAX_SEQ_LENGTH * sizeof(char));
                altLengthsLocal[i] = altLengths[batchCount + i];
            }
        }
#pragma ACCEL parallel factor=3
        for(int i = 0; i < KERNEL_NUM; ++i){
            if(batchCount + i < batchSize){
                smithWatermanKernel(ref, &alt[i][0], refLength, altLengthsLocal[i], overhang_strategy, 0, w_match, w_mismatch, w_open, w_extend, &lengths[i][0], &states[i][0], &cigarElementNum[i], &alignment_offset[i]);
            }
        }
        for(int i = 0; i < KERNEL_NUM; ++i){
            if(batchCount + i < batchSize){
                outputs[batchCount + i] = cigarElementNum[i];
                for(int j = 0; j < cigarElementNum[i]; j++){
                    outputs[batchSize + cigarElementPtr] = lengths[i][j];
                    cigarElementPtr++;
                    outputs[batchSize + cigarElementPtr] = states[i][j];
                    cigarElementPtr++;
                }
                outputs[batchSize + cigarElementPtr] = alignment_offset[i];
                cigarElementPtr++;
            }
        }
    }

    *cigarInfoLength = cigarElementPtr;
}

int smithWatermanKernel(char* ref, char* alt, int refLength, int altLength, int overhang_strategy, int cutoff, int w_match, int w_mismatch, int w_open, int w_extend, short* lengths, short* states, short* cigarElementNum, short* alignment_offset_ptr){
    signed int lowInitValue = -1073741824;
    int MATRIX_MIN_CUTOFF; 
    if(cutoff) MATRIX_MIN_CUTOFF = 0;
    else MATRIX_MIN_CUTOFF = (int)(-1e8);
    int ncol = altLength + 1;
    int nrow = refLength + 1;
    assert(nrow > 0 && nrow < MAX_SEQ_LENGTH);
    assert(ncol > 0 && ncol < MAX_SEQ_LENGTH);
 
    short btrack[MAX_SEQ_LENGTH * MAX_SEQ_LENGTH];
    
    int i = 0;
    int j = 0;
    int lastLastDiag[MAX_SEQ_LENGTH];
    int lastDiag[MAX_SEQ_LENGTH];
    int curDiag[MAX_SEQ_LENGTH];

    int best_gap_v[MAX_SEQ_LENGTH];
    int best_gap_h_last[MAX_SEQ_LENGTH];
    int best_gap_h_curr[MAX_SEQ_LENGTH];
    int gap_size_v[MAX_SEQ_LENGTH];
    int gap_size_h_last[MAX_SEQ_LENGTH];
    int gap_size_h_curr[MAX_SEQ_LENGTH];
    int sw_firstRow[MAX_SEQ_LENGTH + 1];
    int sw_firstCol[MAX_SEQ_LENGTH + 1];
    int sw_finalRow[MAX_SEQ_LENGTH + 1];
    int sw_finalCol[MAX_SEQ_LENGTH + 1];
 
//#pragma ACCEL parallel factor=2
    for(i = 0; i < MAX_SEQ_LENGTH; ++i){
        best_gap_v[i] = lowInitValue;
        gap_size_v[i] = 0;
        best_gap_h_last[i] = lowInitValue;
        gap_size_h_last[i] = 0;
        best_gap_h_curr[i] = lowInitValue;
        gap_size_h_curr[i] = 0;

    }

    sw_firstRow[0] = 0;
    sw_firstCol[0] = 0;
//#pragma ACCEL parallel factor=2
    for(i = 0; i < MAX_SEQ_LENGTH; i++){
        if(i < ncol){
            sw_firstRow[i] = 0;
            sw_finalRow[i] = 0;
        }
    }
//#pragma ACCEL parallel factor=2
    for(i = 0; i < MAX_SEQ_LENGTH; i++){
        if(i < nrow){
            sw_firstCol[i] = 0;
            sw_finalCol[i] = 0;
        }
    }
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL || overhang_strategy == OVERHANG_STRATEGY_LEADING_INDEL){
        sw_firstRow[1] = w_open;
        int currentValue = w_open;
//#pragma ACCEL parallel factor=2
        for(i = 2; i < MAX_SEQ_LENGTH; i++){
            if(i < ncol){
#ifdef MCC_ACC
                sw_firstRow[i] = w_open + (i - 1) * w_extend;
#else
                currentValue += w_extend;
                sw_firstRow[i] = currentValue;
#endif
            }
        }
        sw_firstCol[1] = w_open;
        currentValue = w_open;
//#pragma ACCEL parallel factor=2
        for(i = 2; i < MAX_SEQ_LENGTH; i++){
            if(i < nrow){
#ifdef MCC_ACC
                sw_firstCol[i] = w_open + (i - 1) * w_extend;
#else
                currentValue += w_extend;
                sw_firstCol[i] = currentValue;
#endif
            }
        }
    }
    sw_finalRow[0] = sw_firstCol[nrow - 1];
    sw_finalCol[0] = sw_firstRow[ncol - 1];
 
    char a_base;
    char b_base;
    int k;
    i = 0;
    j = 0;
    char imod3 = 0;

#pragma ACCEL false_dependence variable=lastLastDiag
#pragma ACCEL false_dependence variable=lastDiag
#pragma ACCEL false_dependence variable=curDiag
#pragma ACCEL false_dependence variable=best_gap_v
#pragma ACCEL false_dependence variable=gap_size_v
#pragma ACCEL false_dependence variable=best_gap_h_last
#pragma ACCEL false_dependence variable=best_gap_h_curr
#pragma ACCEL false_dependence variable=gap_size_h_last
#pragma ACCEL false_dependence variable=gap_size_h_curr
    
    int j_PE = 0;
    for(k = 0; k < (2 * MAX_SEQ_LENGTH - 1) * MAX_SEQ_LENGTH / 2; k++){
        j = j_PE * 2;
        //if(i < refLength + altLength - 1 && i >= j && i <= j + refLength - 1 && j < altLength){
        {   int step_diag;
            int step_down;
            int step_right;
            int kd;
            int ki;
            int sw_lastRow_lastCol;
            int sw_lastRow_curCol;
            int sw_curRow_lastCol;
            int prev_gap_v = 0;
            int prev_gap_h = 0;
            int diagHighestOrEqual = 0;
            char imod2 = i % 2;
            //working on j + 1, i - j + 1 element in matrix
            //b_base = alt[j];
            //a_base = ref[i_j_gap];
            
            int lastDiag_pre[2];
            int lastDiag_cur[2];
            int lastLastDiag_pre[2];
            int lastLastDiag_cur[2];
            int curDiag_pre[2];
            int curDiag_cur[2];
            char j_PEmod2 = j_PE & 0x1;
            int jdiv4 = j >> 2;
            if(j_PEmod2 == 0){
                if(j != 0)
                    lastDiag_pre[0] = lastDiag[4 * jdiv4 - 1];
                lastDiag_cur[0] = lastDiag[4 * jdiv4];
                if(j != 0)
                    lastLastDiag_pre[0] = lastLastDiag[4 * jdiv4 - 1];
                lastLastDiag_cur[0] = lastLastDiag[4 * jdiv4];
                if(j != 0 )
                    curDiag_pre[0] = curDiag[4 * jdiv4 - 1];
                curDiag_cur[0] = curDiag[4 * jdiv4]; 

                lastDiag_pre[1] = lastDiag_cur[0];
                lastDiag_cur[1] = lastDiag[4 * jdiv4 + 1];
                lastLastDiag_pre[1] = lastLastDiag_cur[0];
                lastLastDiag_cur[1] = lastLastDiag[4 * jdiv4 + 1];
                curDiag_pre[1] = curDiag_cur[0];
                curDiag_cur[1] = curDiag[4 * jdiv4 + 1];
                
                int curDiagVal[2];
                int i_j_gap_local[2];
                int j_local[2];
                i_j_gap_local[0] = i - j;
                i_j_gap_local[1] = i - j - 1;
                j_local[0] = j;
                j_local[1] = j + 1;
                char a_base_local[2];
                char b_base_local[2];
                a_base_local[0] = ref[i - j];
                a_base_local[1] = ref[i - j - 1];
                b_base_local[0] = alt[j];
                b_base_local[1] = alt[j + 1];
            

#pragma ACCEL parallel factor=2
                for(int index = 0; index < 2; index++){
                    if(i < refLength + altLength - 1 && i >= j_local[index] && i <= j_local[index] + refLength - 1 && j_local[index] < altLength){
                        curDiagVal[index] = smithWatermanCellWiseMiniKernel(lastLastDiag_pre[index], lastLastDiag_cur[index], lastDiag_pre[index], lastDiag_cur[index], curDiag_pre[index], curDiag_cur[index], sw_firstCol, sw_firstRow, i_j_gap_local[index], j_local[index], imod3, a_base_local[index], b_base_local[index], w_match, w_mismatch, w_open, w_extend, best_gap_v, gap_size_v, best_gap_h_curr, best_gap_h_last, gap_size_h_curr, gap_size_h_last, refLength, altLength, sw_finalCol, sw_finalRow, MATRIX_MIN_CUTOFF, btrack, step_diag, step_down, step_right, kd, ki, sw_lastRow_lastCol, sw_lastRow_curCol, sw_curRow_lastCol, prev_gap_v, prev_gap_h, diagHighestOrEqual, imod2, lowInitValue);
                    }
                }
                if(imod3 == 0){
                    curDiag[4 * jdiv4] = curDiagVal[0];
                    curDiag[4 * jdiv4 + 1] = curDiagVal[1];
                }
                else if(imod3 == 1){
                    lastLastDiag[4 * jdiv4] = curDiagVal[0];
                    lastLastDiag[4 * jdiv4 + 1] = curDiagVal[1];
                }
                else{
                    lastDiag[4 * jdiv4] = curDiagVal[0];    
                    lastDiag[4 * jdiv4 + 1] = curDiagVal[1];
                }
            }
            else{
                if(j != 0)
                    lastDiag_pre[0] = lastDiag[4 * jdiv4 + 1];
                lastDiag_cur[0] = lastDiag[4 * jdiv4 + 2];
                if(j != 0)
                    lastLastDiag_pre[0] = lastLastDiag[4 * jdiv4 + 1];
                lastLastDiag_cur[0] = lastLastDiag[4 * jdiv4 + 2];
                if(j != 0 )
                    curDiag_pre[0] = curDiag[4 * jdiv4 + 1];
                curDiag_cur[0] = curDiag[4 * jdiv4 + 2]; 

                lastDiag_pre[1] = lastDiag_cur[0];
                lastDiag_cur[1] = lastDiag[4 * jdiv4 + 3];
                lastLastDiag_pre[1] = lastLastDiag_cur[0];
                lastLastDiag_cur[1] = lastLastDiag[4 * jdiv4 + 3];
                curDiag_pre[1] = curDiag_cur[0];
                curDiag_cur[1] = curDiag[4 * jdiv4 + 3];
                
                int curDiagVal[2];
                int i_j_gap_local[2];
                int j_local[2];
                i_j_gap_local[0] = i - j;
                i_j_gap_local[1] = i - j - 1;
                j_local[0] = j;
                j_local[1] = j + 1;
                char a_base_local[2];
                char b_base_local[2];
                a_base_local[0] = ref[i - j];
                a_base_local[1] = ref[i - j - 1];
                b_base_local[0] = alt[j];
                b_base_local[1] = alt[j + 1];
            

#pragma ACCEL parallel factor=2
                for(int index = 0; index < 2; index++){
                    if(i < refLength + altLength - 1 && i >= j_local[index] && i <= j_local[index] + refLength - 1 && j_local[index] < altLength){
                        curDiagVal[index] = smithWatermanCellWiseMiniKernel(lastLastDiag_pre[index], lastLastDiag_cur[index], lastDiag_pre[index], lastDiag_cur[index], curDiag_pre[index], curDiag_cur[index], sw_firstCol, sw_firstRow, i_j_gap_local[index], j_local[index], imod3, a_base_local[index], b_base_local[index], w_match, w_mismatch, w_open, w_extend, best_gap_v, gap_size_v, best_gap_h_curr, best_gap_h_last, gap_size_h_curr, gap_size_h_last, refLength, altLength, sw_finalCol, sw_finalRow, MATRIX_MIN_CUTOFF, btrack, step_diag, step_down, step_right, kd, ki, sw_lastRow_lastCol, sw_lastRow_curCol, sw_curRow_lastCol, prev_gap_v, prev_gap_h, diagHighestOrEqual, imod2, lowInitValue);
                    }
                }
                if(imod3 == 0){
                    curDiag[4 * jdiv4 + 2] = curDiagVal[0];
                    curDiag[4 * jdiv4 + 3] = curDiagVal[1];
                }
                else if(imod3 == 1){
                    lastLastDiag[4 * jdiv4 + 2] = curDiagVal[0];
                    lastLastDiag[4 * jdiv4 + 3] = curDiagVal[1];
                }
                else{
                    lastDiag[4 * jdiv4 + 2] = curDiagVal[0];    
                    lastDiag[4 * jdiv4 + 3] = curDiagVal[1];
                }
            }
        }
        j_PE = (j_PE == MAX_SEQ_LENGTH / 2 - 1) ? 0 : j_PE + 1;
        if(j_PE == MAX_SEQ_LENGTH / 2 - 1){
            i = i + 1;
        }
        if(j_PE == MAX_SEQ_LENGTH / 2 - 1){
            imod3 = (imod3 == 2) ? 0 : imod3 + 1;
        }
    }

    int p1 = 0;
    int p2 = 0;
    int maxscore = -2147483648;
    int segment_length = 0;
    int curScore = 0;
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL){
        p1 = refLength;
        p2 = altLength;
    }
    else{
        p2 = altLength;
#pragma ACCEL pipeline 
        for(i = 1; i < nrow; ++i){
            curScore = sw_finalCol[i];
            if(curScore >= maxscore){
                p1 = i;
                maxscore = curScore;
            }
        }
        //now look for a larger score on the bottom-most row
        if(overhang_strategy != OVERHANG_STRATEGY_LEADING_INDEL){
#pragma ACCEL pipeline 
            for(j = 1; j < ncol; ++j){
                curScore = sw_finalRow[j];
                int abs1 = refLength - j;
                if(abs1 < 0)
                    abs1 = -abs1;
                int abs2 = p1 - p2;
                if(abs2 < 0)
                    abs2 = -abs2;
                if(curScore > maxscore || (curScore == maxscore && (abs1 < abs2))){
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

    short cigarSize = 0;
    if(segment_length > 0 && overhang_strategy == OVERHANG_STRATEGY_SOFTCLIP){
        lengths[cigarSize] = segment_length;
        states[cigarSize] = STATE_CLIP;
        cigarSize++;
        segment_length = 0;
    }

    short state = STATE_MATCH;
    int btr = 0;
    int new_state =0;
    int step_length = 0;
    do{
        btr = btrack[p1 * MAX_SEQ_LENGTH + p2];
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
            lengths[cigarSize] = segment_length;
            states[cigarSize] = state;
            cigarSize++;
            segment_length = step_length;
            state = new_state;
        }
    }while(p1 >= 0 && p2 >= 0);
    
    p1++;
    p2++;
    if(overhang_strategy == OVERHANG_STRATEGY_SOFTCLIP){
        lengths[cigarSize] = segment_length;
        states[cigarSize] = state;
        cigarSize++;

        if(p2 > 0){
            lengths[cigarSize] = p2;
            states[cigarSize] = STATE_CLIP;
            cigarSize++;
        }
        *alignment_offset_ptr = p1;
    }
    else if(overhang_strategy == OVERHANG_STRATEGY_IGNORE){
        lengths[cigarSize] = segment_length + p2;
        states[cigarSize] = state;
        cigarSize++;
        *alignment_offset_ptr = p1 - p2;
    }
    else{
        lengths[cigarSize] = segment_length;
        states[cigarSize] = state;
        cigarSize++;
 
        if(p1 > 0){
            lengths[cigarSize] = p1;
            states[cigarSize] = STATE_DELETION;
            cigarSize++;
        }
        else if(p2 > 0){
            lengths[cigarSize] = p2;
            states[cigarSize] = STATE_INSERTION;
            cigarSize++;
        }
        *alignment_offset_ptr = 0;
    }

    *cigarElementNum = cigarSize;
    return 0;
}
inline int smithWatermanCellWiseMiniKernel(int lastLastDiag_pre, int lastLastDiag_cur, int lastDiag_pre, int lastDiag_cur, int curDiag_pre, int curDiag_cur, int* sw_firstCol, int* sw_firstRow, int i_j_gap, int j, int imod3, char a_base, char b_base, int w_match, int w_mismatch, int w_open, int w_extend, int* best_gap_v, int* gap_size_v, int* best_gap_h_curr, int* best_gap_h_last, int* gap_size_h_curr, int* gap_size_h_last, int refLength, int altLength, int* sw_finalCol, int* sw_finalRow, int MATRIX_MIN_CUTOFF, short* btrack, int step_diag, int step_down, int step_right, int kd, int ki, int sw_lastRow_lastCol, int sw_lastRow_curCol, int sw_curRow_lastCol, int prev_gap_v, int prev_gap_h, int diagHighestOrEqual, char imod2, int lowInitValue){
    int sw_firstCol_cur = sw_firstCol[i_j_gap];
    int sw_firstCol_nxt = sw_firstCol[i_j_gap + 1];
    int sw_firstRow_cur = sw_firstRow[j];
    int sw_firstRow_nxt = sw_firstRow[j + 1];

    if(j == 0){
        sw_lastRow_lastCol = sw_firstCol_cur;   
    }
    else if(i_j_gap == 0){
        sw_lastRow_lastCol = sw_firstRow_cur;
    }
    else{
        if(imod3 == 0)
            sw_lastRow_lastCol = lastLastDiag_pre;
        else if(imod3 == 1)
            sw_lastRow_lastCol = lastDiag_pre;
        else
            sw_lastRow_lastCol = curDiag_pre;
    }
    if(j == 0){
        sw_curRow_lastCol = sw_firstCol_nxt;
    }
    else{
        if(imod3 == 0)
            sw_curRow_lastCol = lastDiag_pre;
        else if(imod3 == 1)
            sw_curRow_lastCol = curDiag_pre;
        else
            sw_curRow_lastCol = lastLastDiag_pre;
    }
    if(i_j_gap == 0){
        sw_lastRow_curCol = sw_firstRow_nxt;   
    }
    else{
        if(imod3 == 0)
            sw_lastRow_curCol = lastDiag_cur;
        else if(imod3 == 1)
            sw_lastRow_curCol = curDiag_cur;
        else
            sw_lastRow_curCol = lastLastDiag_cur;
    }
    step_diag = sw_lastRow_lastCol +((a_base == b_base) ? w_match : w_mismatch);
    int cur_best_gap_v;
    int cur_gap_size_v;
    prev_gap_v = sw_lastRow_curCol + w_open;
    cur_best_gap_v = best_gap_v[j] + w_extend;
    if(prev_gap_v > cur_best_gap_v){
        cur_best_gap_v = prev_gap_v;
        cur_gap_size_v = 1;
    }
    else{
        cur_gap_size_v = gap_size_v[j] + 1;
    }
    best_gap_v[j] = cur_best_gap_v;
    gap_size_v[j] = cur_gap_size_v;
    step_down = cur_best_gap_v;
    kd = cur_gap_size_v;
    int cur_best_gap_h;
    int cur_gap_size_h;
    prev_gap_h = sw_curRow_lastCol + w_open;
    if(j == 0){
        cur_best_gap_h = lowInitValue + w_extend;
    }
    else{
        if(imod2 == 0){
            cur_best_gap_h = best_gap_h_last[j - 1] + w_extend;
        }
        else{
            cur_best_gap_h = best_gap_h_curr[j - 1] + w_extend;
        }
    }
    if(prev_gap_h > cur_best_gap_h){
        cur_best_gap_h = prev_gap_h;
        cur_gap_size_h = 1;
    }
    else{
        if(j == 0)
            cur_gap_size_h = 1;
        else{
            if(imod2 == 0)
                cur_gap_size_h = gap_size_h_last[j - 1] + 1;
            else
                cur_gap_size_h = gap_size_h_curr[j - 1] + 1;
        }

    }
    if(imod2 == 0){
        best_gap_h_curr[j] = cur_best_gap_h;
        gap_size_h_curr[j] = cur_gap_size_h;
    }
    else{
        best_gap_h_last[j] = cur_best_gap_h;
        gap_size_h_last[j] = cur_gap_size_h;
    }
    step_right = cur_best_gap_h;
    ki = cur_gap_size_h;
       
    diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
    
    int curDiagVal;
    int curBtrack;
    if(diagHighestOrEqual){
        curDiagVal = (MATRIX_MIN_CUTOFF > step_diag ? MATRIX_MIN_CUTOFF : step_diag);
        curBtrack = 0;
    }
    else if(step_right >= step_down){
        curDiagVal = (MATRIX_MIN_CUTOFF > step_right ? MATRIX_MIN_CUTOFF : step_right);
        curBtrack = -ki;
    }
    else{
        curDiagVal = (MATRIX_MIN_CUTOFF > step_down ? MATRIX_MIN_CUTOFF : step_down);
        curBtrack = kd;
    }
    btrack[i_j_gap * MAX_SEQ_LENGTH + j] = curBtrack;

    if(j == altLength - 1){
        sw_finalCol[i_j_gap + 1] = curDiagVal;
    }
    if(i_j_gap == refLength - 1){
        sw_finalRow[j + 1] = curDiagVal;
    }

    return curDiagVal;   
}


