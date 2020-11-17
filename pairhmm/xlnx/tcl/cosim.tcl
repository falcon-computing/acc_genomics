# Project settings
open_project pairhmm -reset
set XILINX_SDX /curr/software/Xilinx/SDx/2017.1
add_files device/pairhmm.cpp -cflags "-DHLS_SIM -D DIE_NUM=3 -g -I. -I./common -I./device -I./host"
add_files -tb host/avx_impl.cpp -cflags "-O3 -march=native -mavx -std=c++0x -DSIM -DHLS_SIM -g -I. -I./common -I./device -I./host"
add_files -tb host/baseline_impl.cpp -cflags "-O3 -march=native -mavx -std=c++0x -DSIM -DHLS_SIM -g -I. -I./common -I./device -I./host"
add_files -tb host/FalconPairHMM.cpp -cflags "-O3 -march=native -mavx -std=c++0x -g -DSIM -DHLS_SIM -D DIE_NUM=3 -D SLR0_PE_NUM=8 -D SLR1_PE_NUM=8 -D SLR2_PE_NUM=8 -I. -I./common -I./device -I./host -I$XILINX_SDX/runtime/include/1_2 -I$XILINX_SDX/Vivado_HLS/include"
add_files -tb pairhmm_test.cpp -cflags "-O3 -march=native -mavx -std=c++0x -g -DSIM -DHLS_SIM -D DIE_NUM=3 -D SLR0_PE_NUM=8 -D SLR1_PE_NUM=8 -D SLR2_PE_NUM=8 -I. -I./common -I./device -I./host -I$XILINX_SDX/runtime/include/1_2 -I$XILINX_SDX/Vivado_HLS/include"
set_top computePairhmmFpgaMerlin

# Solution settings
open_solution -reset solution1
set_part xcvu9p-flgb2104-2-i
#set_part xcku115-flvb2104-2-e
create_clock -period 200MHz -name default

config_rtl -register_reset

csim_design -ldflags "-L$XILINX_SDX/runtime/lib/x86_64 -lxilinxopencl"
csynth_design 
#cosim_design -ldflags "-L$XILINX_SDX/runtime/lib/x86_64 -lxilinxopencl"
cosim_design -bc -ldflags "-L$XILINX_SDX/runtime/lib/x86_64 -lxilinxopencl"
#cosim_design -trace_level all -ldflags "-L$XILINX_SDX/runtime/lib/x86_64 -lxilinxopencl"

