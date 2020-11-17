#ifndef PAIRHMM_FPGA_H
#define PAIRHMM_FPGA_H

#ifdef XILINX
#include "acc_lib/XCLEnv.h" // load xilinx opencl framework
#else
#include "acc_lib/AOCLEnv.h" // load intel opencl framework
#endif

#include "PairHMMHostInterface.h"
#include "PairHMMFpgaInterface.h"

extern double peak_kernel_gcups;
extern double curr_kernel_gcups;

float* compute_fpga(
    const char* bit_path,
    std::string read_data,
    std::string hap_data,
    uint64_t num_cell);

void __attribute__((destructor)) cleanup() ;

#endif

