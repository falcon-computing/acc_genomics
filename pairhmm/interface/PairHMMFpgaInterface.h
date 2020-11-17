#ifndef PAIRHMM_FPGA_INTERFACE_H
#define PAIRHMM_FPGA_INTERFACE_H
#include <cstdint>
#include <string>

#include "gkl-pairhmm/Context.h"
#include "PairHMMHostInterface.h"
#include "m2m.h"

extern double peak_est1_gcups;
extern double peak_est2_gcups;
extern double curr_est1_gcups;
extern double curr_est2_gcups;

// Xilinx FPGA constants
#define MAX_READ_LEN      192
#define MAX_HAP_LEN       1024
#define MAX_RSDATA_NUM    2048   //default 2048, minimum 32, maximum 4096
#define MAX_HAPDATA_NUM   128    //default 128, minimum 8, maximum 256
#define READ_BLOCK_SIZE   2
#define HAP_BLOCK_SIZE    4
#define DEP_DIST          42

// FPGA input data structures
typedef struct {
  int   haplen;
  float one_over_haplen; // 1 / haplen
} HaplenBundle;

typedef struct {
  // 31:28 readBases
  // 27:21 readQuals
  // 20:14 insertionGOP
  // 13:7 deletionGOP
  // 6:0 overallGCP
  uint32_t scores; 
  float    m2m;
} ReadBundle;

typedef struct {
  int  idx;
  int  num_readinfo;
  int* readinfo_idx;
  uint64_t  trip_count;
} PU;

typedef struct {
  uint8_t  read_len[READ_BLOCK_SIZE]; // length of the reads
  int      output_idx;                // idx to output likelihood array
  int      trip_count;                // loop trip count of a PU
  bool     is_full;
} ReadInfo;

typedef struct {
  uint64_t     read_info[MAX_RSDATA_NUM / READ_BLOCK_SIZE];
  ReadBundle   read_data[MAX_RSDATA_NUM][MAX_READ_LEN];
  HaplenBundle hap_len[MAX_HAPDATA_NUM];
  uint16_t     hap_data[MAX_HAPDATA_NUM / HAP_BLOCK_SIZE][MAX_HAP_LEN];
  uint16_t     num_reads_per_pu[64];
  uint64_t     trip_count_per_pu[64]; 
} FpgaInputBundle;

// helper functions for data conversions
static inline uint8_t read2char(char c) {

  uint8_t ret = 0;
  switch (c) {
    case 'A': ret = 0; break;
    case 'C': ret = 1; break;
    case 'T': ret = 2; break;
    case 'G': ret = 3; break;
    case 'N': ret = 4; break;
  }
  return ret;
}

// serialize and distribute fpga data
// assuming all input data is fpga ready
void pack_fpga_input(
    const int num_pu,
    const int num_read,
    const int num_hap,
    const read_t* reads,
    const hap_t*  haps,
    FpgaInputBundle *bundle);

#endif
