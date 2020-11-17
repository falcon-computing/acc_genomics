# Smith-Waterman

## Description
This design contains the Smith-Waterman acceleration kernel for haplotypecaller in GATK built by Falcon Computing Solutions.

## Major Optimization
This design accelerates both on Intel-AVX and Xilinx-FPGA, when the task size is small, it uses the AVX to accelerate. Otherwise it uses Xilinx FPGA
1. For Intel-AVX solution, we explore the parallelism in the row direction, which is called striped smith-waterman. 
   Basic idea is that we first ingore the dependency on the column direction and do the calculation in a SIMD fashion.
   Then we do a max scan on column direction
2. For Xilinx-FPGA solution, we explore the parallelism in the anti-diagonal direction. we do 8-PE at same time on each anti-diagonal. 

## Usage:
    make help: display all the helping message
	make all: make the falcon .so file only.
    before doing any FPGA runs, please change the platform in the Makefile
	make test_sw_emu: run the software emulation, if you want to turn off the xilinx software emulation and run AVX only, You can comment out the FPGA_FLAG.
	make test_hw: run the hardware compilation, this will take many hours. Run it only once.
	make clean: remove all the temp files and binary files (except the FPGA bitstream files)

## QoR results on Xilinx ku115 FPGA and Intel(R) Xeon(R) CPU E5-2687W v4 @ 3.00GHz 
tools:  xocc version: v2017.1_sdx (64-bit) DSA versoin: 4.0 
        gcc/g++ version: 4.8.5 20150623 (Red Hat 4.8.5-4)

|PE#| LUTs % | FFs % | BRAMs % | DSPs % |freq MHz| Throughput on FPGA(W/O PCIe latency) GCUPs | Throughput on FPGA GCUPs |
|:-:|:------:|:-----:|:-------:|:------:|:------:|:------------------------------------------:|:-------------------:|
|96|35|19|45|0.0|200|1.45|0.3|


