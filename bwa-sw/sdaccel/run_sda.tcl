#*******************************************************************************
# Define the solution for SDAccel
set bit [lindex $argv 0]
create_solution -name $bit -dir . -force
add_device -vbnv xilinx:adm-pcie-ku3:2ddr:3.1
#add_device -vbnv xilinx:adm-pcie-ku3:2ddr-xpr:3.1
#add_device -vbnv xilinx:adm-pcie-7v3:1ddr:2.0
#add_device -vbnv baidu:leopard-pcie-7v2:fpga_card:1.0

set PWD [pwd]

# Host Compiler Flags
set_property -name host_cflags -value "-g -Wall -D FPGA_DEVICE -D C_KERNEL" -objects [current_solution]

# Host Source Files
#add_files "main_cl.c"

# set to 0 for building host only
if {1} {
# Kernel Definition
create_kernel sw_top -type c
add_files -kernel [get_kernels sw_top] "smithwaterman.cpp"

# Define Binary Containers

create_opencl_binary smithwaterman
set_property region "OCL_REGION_0" [get_opencl_binary smithwaterman]
create_compute_unit -opencl_binary [get_opencl_binary smithwaterman] -kernel [get_kernels sw_top] -name k1

if {0} {
# Compile the design for CPU based emulation
compile_emulation -flow cpu -opencl_binary [get_opencl_binary smithwaterman]

# Run the compiled application in CPU based emulation mode
run_emulation -flow cpu -args "smithwaterman.xclbin $PWD/data/total_input.dat"

compile_emulation -flow hardware -opencl_binary [get_opencl_binary smithwaterman]

# Run the compiled application in CPU based emulation mode
run_emulation -flow hardware -args "smithwaterman.xclbin $PWD/data/total_input.dat"
}
}

#compile_emulation -flow cpu -opencl_binary [get_opencl_binary smithwaterman]
# Compile the application to run on the accelerator card
#set_param compiler.worstNegativeSlack -1.0
build_system

# Package the application binaries
#package_system

