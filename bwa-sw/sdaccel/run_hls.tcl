set top sw_top

# Project settings
open_project proj_new
add_file smithwaterman.cpp -cflags "-DHLS_"
add_file -tb main_cl.cpp  -cflags "-DHLS_";
set_top $top

# Solution settings
open_solution -reset solution3
#set_part xc7vx690tffg1157-2
set_part xcku060-ffva1156-2-e
create_clock -period 250MHz
set_clock_uncertainty 1.080000
config_interface -m_axi_addr64
config_compile -pipeline_loops 64

csim_design -argv "0 [pwd]/testdata/input [pwd]/testdata/golden"
#exit
config_rtl -register_reset
#config_array_partition -throughput_driven
csynth_design
#export_design -evaluate verilog
#cosim_design -bc
#cosim_design 
#cosim_design -argv "49 9 3" -bc
#cosim_design -argv "49 9 3"

