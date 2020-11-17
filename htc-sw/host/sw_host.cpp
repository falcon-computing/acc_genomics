#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "common.h"
#include "intel_avx/avx2_impl.h"

long calMatrix_C_time = 0;
long calCigar_C_time = 0;
long malloc_time = 0;
long SW_complexity = 0;

void FalconSWFPGA_init(char*);
double FalconSWFPGA_run(char*, int, char[][MAX_SEQ_LENGTH], int*, int, int, int, int, int, int, struct Cigar*, int*, bool);
void FalconSWFPGA_release();
int addCigarElement(struct Cigar* cigar, int length, int state){
    if(cigar->CigarElementNum < 0)
        return -1;
    if(length > 0){
        (cigar->cigarElements[cigar->CigarElementNum]).length = length;
        (cigar->cigarElements[cigar->CigarElementNum]).state = state;
        cigar->CigarElementNum++;
    }
    return 0;
}


int GetInputs(int batch_id, char* ref, int* refLength, char altArray[][MAX_SEQ_LENGTH], int* altArraySize, int* altArrayLengths){
   // printf("start  getInputs function\n");
    char common_path[100] = "/curr/jysheng/gatk3_tmp/input";
    char new_common_path[100];
    sprintf(new_common_path, "%s%d", common_path, batch_id);
    FILE *fp;
    //printf("input file is %s\n", new_common_path);
    fp = fopen(new_common_path, "r");
    if(!fp){
        printf("file %s open failed\n", new_common_path);
        exit(-1);
    }
    fgets(ref, MAX_SEQ_LENGTH, fp);
    //printf("ref is %s\n", ref);
    int i = 0;
    while(ref[i])
        i++;
    *refLength = i - 1;
    char altArraySizeString[10];
    fgets(altArraySizeString, 10, fp);
    *altArraySize = atoi(altArraySizeString);
    for(i = 0; i < *altArraySize; ++i){
        fgets(altArray[i], MAX_SEQ_LENGTH, fp);
        int j = 0;
        while(altArray[i][j])
            j++;
        altArrayLengths[i] = j - 1;
    }
    fclose(fp);
    return 0;

}
int printCigar(struct Cigar* cigar){
    int i = 0;
    if(cigar->CigarElementNum == 0){
        printf("this cigar has not element");
    }
    else{

        printf("this cigar %p  has %d elements\n", cigar, cigar->CigarElementNum);
    }
    for(i = 0; i < cigar->CigarElementNum; ++i){
        printf("%d",cigar->cigarElements[i].length);
        if(cigar->cigarElements[i].state == STATE_MATCH){
            printf("M");
        }
        else if(cigar->cigarElements[i].state == STATE_INSERTION){
            printf("I");
        }
        else if(cigar->cigarElements[i].state == STATE_DELETION){
            printf("D");
        }
        else if(cigar->cigarElements[i].state == STATE_CLIP){
            printf("S");
        }
        else{
            printf("Illegal cigar Op in printCigar(), which is %d\n", cigar->cigarElements[i].state);
            return -1;
        }
    }
    printf("\n");
    return 0;
}
struct timespec diff_time(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}


int cmpCigarResults(struct Cigar* golden, struct Cigar* target, int batchSize, int& error_batch_id){
    int i = 0;
    int j = 0;
    for(i = 0; i < batchSize; ++i){
        if(golden[i].CigarElementNum != target[i].CigarElementNum){
            printf("batch id = %d, golden[%d] num is %d, target[%d], num is %d\n", i, i, golden[i].CigarElementNum, i,  target[i].CigarElementNum);
            printf("golden is ");
            printCigar(&golden[i]);
            printf("target address is %p\n", &target[i]);
            printf("target is ");
            printCigar(&target[i]);
            error_batch_id = i;
            return 0;
        }
        for(j = 0; j < golden[i].CigarElementNum; ++j){
            if(golden[i].cigarElements[j].length != target[i].cigarElements[j].length || golden[i].cigarElements[j].state != target[i].cigarElements[j].state){
                printf("i = %d, j = %d, golden length is %d, ret length is %d, golden state is %d, ret state is %d\n", i, j, golden[i].cigarElements[j].length, target[i].cigarElements[j].length, golden[i].cigarElements[j].state, target[i].cigarElements[j].state);
                printf("golden is ");
                printCigar(&golden[i]);
                printf("target is ");
                printCigar(&target[i]);
                error_batch_id = i;
                return 0;
            }
        }
    }
            
    return 1;
}

int sum_array(int* array, int size){
    int i = 0;
    int ret = 0;
    for(i = 0; i < size; i++){
        ret += array[i];
    }
    return ret;
}

int GenInputs(char* ref, int* refLength, char alts[][MAX_SEQ_LENGTH], int* batchSize, int* altLengths){
    int avgLength = 286;
    float diffRatio = 0.05;
#define MAX_FPGA_SEQ_LENGTH 512
    int refRange = 2 * (MAX_FPGA_SEQ_LENGTH - avgLength - 1);
    int refBase = 2 * avgLength - MAX_FPGA_SEQ_LENGTH;
    *refLength = refBase + rand() % refRange;
    int altRange = refRange * diffRatio;
    for(int i = 0; i < *batchSize; i++){
        altLengths[i] = *refLength - altRange / 2 + rand() % altRange;
        if(altLengths[i] > MAX_FPGA_SEQ_LENGTH - 2){
            altLengths[i] = MAX_FPGA_SEQ_LENGTH - 2;
        }
        if(altLengths[i] == 0){
            altLengths[i] = 1;
        }
    }
    char GENE[4] = {'A', 'T', 'C', 'G'};
    for(int i = 0; i < *refLength; i++){
        ref[i] = GENE[rand() % 4];
    }
    for(int i = 0; i < *batchSize; i++){
        for(int j = 0; j < altLengths[i]; j++){
            if(rand() % 10 == 0){
                alts[i][j] = GENE[rand() % 4];
            }
            else{
                if(j < *refLength){
                    alts[i][j] = ref[j];
                }
                else{
                    alts[i][j] = GENE[rand() % 4];
                }
            }
        }
    }
    return 0;
}

int main(int argc, char*argv[]){
    srand(time(NULL));
    //test CigarUtilsCalculateCigar
    //gen rand inputs from input 

    printf("load bitstream %s\n", argv[1]);
    FalconSWFPGA_init(argv[1]);
    char ref[MAX_SEQ_LENGTH];
    char alt[MAX_BATCH_SIZE][MAX_SEQ_LENGTH];
    int refLength;
    int batchSize;
    int altArrayLengths[MAX_BATCH_SIZE];

    double total_golden_time = 0;
    double total_target_time = 0;
    double total_target_kernel_time = 0;
    double total_SW_complexity = 0;
    struct Cigar golden[MAX_BATCH_SIZE];
    struct Cigar target[MAX_BATCH_SIZE];
    int alignment_offset_golden[MAX_BATCH_SIZE];
    int alignment_offset_target[MAX_BATCH_SIZE];
    int i = 0;
    int minBatchSize = 1000000;
    int maxBatchSize = 0;
    int minRefLength = 1000000;
    int maxRefLength = 0;
    float avgRefLength = 0;
    float avgBatchSize = 0;
    int testNum = 1;
    //int slow_count = 0;
    //int fast_singleBatchCount = 0;
 
    struct timespec time1, time2, time_diff;
    
    for(batchSize = 1; batchSize < 256; batchSize *= 2){ 
        total_golden_time = 0;
        total_target_time = 0;
        total_target_kernel_time = 0;
        total_SW_complexity = 0;
        minBatchSize = 1000000;
        maxBatchSize = 0;
        minRefLength = 1000000;
        maxRefLength = 0;
        avgRefLength = 0;
        avgBatchSize = 0;
        
        for(i = 0; i < testNum; ++i){
            if(GenInputs(ref, &refLength, alt, &batchSize, altArrayLengths) < 0){
                exit(-1);
            }
            
            if(batchSize > maxBatchSize)
                maxBatchSize = batchSize;
            if(batchSize < minBatchSize)
                minBatchSize = batchSize;
            if(refLength > maxRefLength)
                maxRefLength = refLength;
            if(refLength < minRefLength)
                minRefLength = refLength;
            avgRefLength = (avgRefLength * i + refLength) / (i + 1);
            avgBatchSize = (avgBatchSize * i + batchSize) / (i + 1);
            int overhang_strategy = rand() % 4;
            overhang_strategy = 3;

            clock_gettime(CLOCK_REALTIME, &time1);
            SWPairwiseAlignmentMultiBatch(ref, refLength, alt, batchSize, altArrayLengths, golden, alignment_offset_golden, overhang_strategy, 0);
            clock_gettime(CLOCK_REALTIME, &time2);
            time_diff = diff_time(time1, time2);
            long cur_golden_time = 0;
            cur_golden_time = (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);
            total_golden_time += cur_golden_time;
            clock_gettime(CLOCK_REALTIME, &time1);
/*#ifdef FPGA
            total_target_kernel_time += FalconSWFPGA_run(ref, refLength, alt, altArrayLengths, batchSize, overhang_strategy, W_MATCH, W_MISMATCH, W_OPEN, W_EXTEND, target, alignment_offset_target, true);
#else 
            total_target_kernel_time += FalconSWFPGA_run(ref, refLength, alt, altArrayLengths, batchSize, overhang_strategy, W_MATCH, W_MISMATCH, W_OPEN, W_EXTEND, target, alignment_offset_target, false);
#endif*/
            for(int ii = 0; ii < batchSize; ii++){
                uint8_t* ref_ptr = (uint8_t*)ref;
                uint8_t* alt_ptr = (uint8_t*)alt[ii];
                alignment_offset_target[ii] = runSWOnePairBT_fp_avx2(W_MATCH, W_MISMATCH, W_OPEN, W_EXTEND, ref_ptr, alt_ptr, refLength, altArrayLengths[ii], overhang_strategy, &target[ii]);
            }


            clock_gettime(CLOCK_REALTIME, &time2);
            time_diff = diff_time(time1, time2);
            long cur_target_time = 0;
            cur_target_time = (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);
            total_target_time += cur_target_time;
     /*       if(cur_golden_time < cur_target_time){
                slow_count++;
                printf("batch id %d, target is slower than golden, batch size is %d, total %d out of %d is slower\n", i, batchSize, slow_count, i);
            }
            else{
                if(batchSize == 1){
                    fast_singleBatchCount++;
                    printf("batch id %d size is 1, targt is faster than golden, total %d out %d is faster\n", i, fast_singleBatchCount, i);
                }
            }*/

            
            for(int index = 0; index < batchSize; index++){
                if(alignment_offset_golden[index] != alignment_offset_target[index]){
                    printf("alignment offset error happens in test%d, batch id %d\n", i, index);
                    printf("target %d, golden %d\n", alignment_offset_target[index], alignment_offset_golden[index]);
                    exit(-1);
                }
            }
            int error_id;
            if(cmpCigarResults(golden, target, batchSize, error_id) == 0){
                printf("error happens in test %d, batch id = %d \n", i, error_id);
                printf("overhang strategy is %d\n", overhang_strategy);
                printf("refLength = %d, batchSize = %d\n", refLength, batchSize);
                printf("altLength[%d] = %d \n", error_id, altArrayLengths[error_id]);
                int index;
                printf("ref is \n");
                for(index = 0; index < refLength; index++){
                    printf("%d, ", ref[index]);
                }
                printf("\nalt is \n");
                for(index = 0; index < altArrayLengths[error_id]; index++){
                    printf("%d, ", alt[error_id][index]);
                }
                //return -1;
            }   
            else{
                if(i % 1000 == 0)
                    printf("this test %d passed\n", i);
            }
            printf("this test %d passed\n", i);
            total_SW_complexity += refLength * sum_array(altArrayLengths, batchSize);

        }
        
        printf("================This is for batch size %d===================\n", batchSize);
        printf("min ref length is %d, max ref length is %d, avg ref length is %f\n", minRefLength, maxRefLength, avgRefLength);
        printf("min batch size is %d, max batch size is %d, avg batch size is %f\n", minBatchSize, maxBatchSize, avgBatchSize);
        printf("total golden time is %lf secs\n", total_golden_time*(1e-9)); 
        printf("total target time is %lf secs\n", total_target_time*(1e-9)); 
        printf("total pure target time is %lf secs\n", total_target_kernel_time*(1e-9)); 
        printf("golden GCUPs is %lf \n", total_SW_complexity/total_golden_time);
        printf("target GCUPs is %lf \n", total_SW_complexity/total_target_time);

    }

    return 0;

}

