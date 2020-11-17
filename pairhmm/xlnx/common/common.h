#ifndef COMMON_H
#define COMMON_H
#define MAX_READ_LEN 192
#define MAX_HAP_LEN 1024
#define MAX_RSDATA_NUM 2048   //default 2048, minimum 32, maximum 4096
#define MAX_HAPDATA_NUM 128 //default 128, minimum 8, maximum 256
#define DEP_DIST 42
#define READ_BLOCK_SIZE 2
#define HAP_BLOCK_SIZE 4
#define FPGA_PERF 16.5
#define AVX_PERF 0.3
//this is the pipeline stages of the computeEngine pipeline, modify this based on that shown in HLS report
#endif

