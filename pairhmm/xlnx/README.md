# PairHMM

## Description
This design contains the PairHMM acceleration kernel for haplotypecaller in GATK built by Falcon Computing Solutions.

## Major Optimization
This design accelerates both on Intel-AVX and Xilinx-FPGA, when the task size is small, it uses the AVX to accelerate. Otherwise it uses Xilinx FPGA
1. For Intel-AVX solution, we port it from Intel GKL library
2. For Xilinx-FPGA solution, we explore the parallelism in the anti-diagonal direction.
   Instead of traditional Systolic-Array fashion, we explore parallelism in coarse-grain level. On each FPGA die, we deploy a certain number of processing units (PUs).
   Differnt PUs work on different reads. In each PU, we deploy 8 processing elements (PEs), different PEs work on different haplotypes and reads. 
   Each PE works on a dynamic programming scoreing matrix along the anti-diagnoal direction. This hierarchical architecture blance the resource utilization of LUTS, BRAMs and DSPs, and achieve good timing.

## Usage:
    make host:              make the host executable
    make help:              get the help information from Makefile
	make test:              run the execution. Change the FPGA and BORAD_SETTING in the Makefile. If FPGA=1 and BOARD_SETTING=1, will run the bitstream generation and on-board execution, takes many hours! Be careful about this. If FPGA=1 and BOARD_SETTING=0, will run xocc software emulation flow. If FPGA = 0 and BOARD_SETTING=0, will run everything on AVX
	make hls:               run the HLS only, use this to do hardware resource utilization estimation
    make bitgen:            run the bitstream generation, which takes many hours
	make clean:             remove all the temp files and DST files

 
## QoR results on Xilinx vcu1525 FPGA and Intel(R) Xeon(R) CPU E5-2687W v4 @ 3.00GHz 
tools:  xocc version: v2017.4_sdx (64-bit) DSA versoin: 5.0
        gcc/g++ version: 4.8.5 20150623 (Red Hat 4.8.5-4)

Dataset: 
        NA12878-Garvan-Vial1

FPGA board:
        Xilinx VCU1525

|PE#| LUTs % | FFs % | BRAMs % | DSPs % |freq MHz| Avg Throughput on FPGA GCUPs | Peak Throughput on AVX GCUPs |
|:-:|:------:|:-----:|:-------:|:------:|:------:|:----------------------------:|:----------------------------:|
|152|47      |32     |48       |56      |200     |15                            |30                            |               


