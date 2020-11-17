# SMEM (supermaximal exact match)

## Description
This design contains the SMEM acceleration kernel for BWA-MEM built by Falcon Computing Solutions.

## Major Optimization
This design is the replacement of mem_collect_intv() function in BWA-MEM. This kernel originally takes 40% of overall BWA-MEM. 
SMEM is an DRAM random access bound application. The performance is bound by massive amount of random queries of BWT of reference (>3 GB), which can only be fit on DRAM.
The three major optimization techniques in this design are:
1. Reorganize the software SMEM algorithm into 7 stages so that they can be fully dataflowed. The 7 stages are input duplication, forward extension (fe), backward extension (be), long SMEM forward extension (lfe), long SMEM backward extension, all-like forward extension (afe)
2. Apply loop tiling on all the fe, be, lfe, lbe, and afe stages so that the DRAM access latency can be hidden
3. Isolate BWT query from SMEM dataflow engine, so that enough DRAM accesses can be batched together in order to hide DRAM access latency
 
## QoR results on Xilinx vcu1525 FPGA and Intel(R) Xeon(R) CPU E5-2687W v4 @ 3.00GHz 
tools:  xocc version: v2018.2_sdx (64-bit) DSA versoin: 5.0
        gcc/g++ version: 4.8.5 20150623 (Red Hat 4.8.5-4)

FPGA board:
        Xilinx VCU1525

|BANK#| LUTs % | LUTMem % | FFs % | BRAMs % | DSPs % | URAMs % |freq MHz| Avg DRAM Bandwidth GBps | Singe Core CPU baseline (GBps) |
|:---:|:------:|:--------:|:-----:|:-------:|:------:|:-------:|:------:|------------------------:|:------------------------------:|
|4    |19.5    |2.41      |12.6   |44.7     |0       |10.4     |286     |3                        |0.33                            |               


