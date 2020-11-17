/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: AVX function declarations and macros
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/


#ifndef COMMON_H
#define COMMON_H

#define MAX_SEQ_LENGTH 1536
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

#define max(a,b) \
    ({__typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
            _a > _b ? _a : _b; })

#define min(a,b) \
    ({__typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
            _a < _b ? _a : _b; })

//#define abs(a) ({__typeof__ (a) _a = (a); _a > 0  ? _a : -_a; })

#define wd(x, y, w1, w2) \
    ({__typeof__ (x) _x = (x); \
      __typeof__ (y) _y = (y); \
      __typeof__ (w1) _w1 = (w1); \
      __typeof__ (w2) _w2 = (w2); \
            _x == _y ? _w1 : _w2 ; })


struct CigarElement{
    int length;
    int state;
};

struct Cigar{
    struct CigarElement cigarElements[MAX_SEQ_LENGTH];
    int CigarElementNum;
};


int CigarUtilsCalculateCigar(char*, int, char[][MAX_SEQ_LENGTH], int, int*);
int CigarUtilsCalculateCigarBatch(char*, int, char[][MAX_SEQ_LENGTH], int, int*, int);
int createIndelString(struct Cigar*, int, char*, int, char*, int, int, int, char*, int*);
int moveCigarLeft(struct Cigar*, int, struct Cigar*);
int cigarHasZeroSizeElement(struct Cigar*);
int arraysEqual(char*, int, char*, int);
int copyCigar(struct Cigar*, struct Cigar*);
int leftAlignSingleIndel(struct Cigar*, char*, int, char*, int, int, int, struct Cigar*);
int leftAlignCigarSequentially(struct Cigar*, char*, int, char*, int, int, int, struct Cigar*);
int getReadLength(struct Cigar*);
int getReferenceLength(struct Cigar*);
int trimCigarByBases(struct Cigar*, int, int, struct Cigar*);
int addCigarElement(struct Cigar*, int, int);
int needsConsolidation(struct Cigar*);
int consolidate(struct Cigar*, struct Cigar*);
int isSWFailure(struct Cigar*, int);
int SWPairwiseAlignmentOneBatch(char*, char*, int, int, struct Cigar*, int*, int, int);
int calculateMatrixOneBatch(char*, char*, int**, int, int, int**, int, int);
int calculateMatrixRowWise(char*, char*, int**, int, int, int**, int, int);
int calculateMatrixRowWiseSIMD(char*, char*, int**, int, int, int**, int, int);
int calculateMatrixRowWiseSIMDUnroll4x(char*, char*, int**, int, int, int**, int, int, int, int, int, int);
int calculateMatrixOneBatchUnroll2x(char*, char*, int**, int, int, int**, int, int);
int calculateMatrixOneBatchUnroll4x(char*, char*, int**, int, int, int**, int, int);
int calculateCigarOneBatch(int**, int**, int, int, int, struct Cigar*, int*);
int SWPairwiseAlignmentMultiBatch(char*, int, char[][MAX_SEQ_LENGTH], int, int*, struct Cigar*, int*, int, int);

struct timespec diff_time(struct timespec, struct timespec);
int printCigar(struct Cigar*);

#endif
