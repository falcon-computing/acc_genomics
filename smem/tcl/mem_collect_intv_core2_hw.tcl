catch {::common::set_param -quiet hls.xocc.mode csynth};
catch {::common::set_param -quiet hls.xocc.enable_rule_check_client 1};
# 
# Hls run script generated by SDAccel
# 

set xocc_optimize_level 0
open_project mem_collect_intv_core2
set_top mem_collect_intv_core2
add_files "/curr/jysheng/acc_lib_smem/acc_lib/smem/device/smem.cpp" -cflags " -D CPP -I /curr/jysheng/acc_lib_smem/acc_lib/smem "
open_solution solution
set_part xcvu9p-fsgd2104-2-i
create_clock -period 300MHz -name default
config_sdx -target xocc -optimization_level $xocc_optimize_level
config_dataflow -strict_mode warning
set_clock_uncertainty 27.000000%
config_interface -m_axi_addr64
csynth_design
export_design -format ipxact -kernel_drc -sdaccel -ipname mem_collect_intv_core2
close_project
puts "Vivado HLS completed successfully"
exit
