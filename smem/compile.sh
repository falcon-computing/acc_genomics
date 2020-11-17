PLATFORM=xilinx_vcu1525_dynamic_5_0
HLS_FREQ=300
IMP_FREQ=300
PWD=`pwd`
KEN_NME=mem_collect_intv
KEN_SRC=./device/smem.cpp
DATE=`date +%m%d`
XCLBIN_NAME=${KEN_NME}_${PLATFORM}_${DATE}.xclbin
IMP_FLG="--xp vivado_param:bd.ForceAppCoreUpgrade=1 \
--xp param:compiler.userPostSysLinkTcl=${PWD}/tcl/syslink.tcl \
--xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.TCL.POST=${PWD}/tcl/postopt.tcl \
--xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.ARGS.DIRECTIVE=Explore \
--xp vivado_prop:run.impl_1.STEPS.PLACE_DESIGN.ARGS.DIRECTIVE=Explore \
--xp vivado_prop:run.impl_1.STEPS.PHYS_OPT_DESIGN.IS_ENABLED=true \
--xp vivado_prop:run.impl_1.STEPS.PHYS_OPT_DESIGN.ARGS.DIRECTIVE=AggressiveExplore \
--xp vivado_prop:run.impl_1.STEPS.ROUTE_DESIGN.ARGS.DIRECTIVE=Explore \
--sp ${KEN_NME}_core0_1.m_axi_bwt_para:bank0 \
--sp ${KEN_NME}_core0_1.m_axi_seq:bank0 \
--sp ${KEN_NME}_core0_1.m_axi_seq_len:bank0 \
--sp ${KEN_NME}_core0_1.m_axi_mem:bank0 \
--sp ${KEN_NME}_core0_1.m_axi_mem_num:bank0 \
--sp ${KEN_NME}_core0_1.m_axi_bwt:bank0 \
--sp ${KEN_NME}_core1_1.m_axi_bwt_para:bank1 \
--sp ${KEN_NME}_core1_1.m_axi_seq:bank1 \
--sp ${KEN_NME}_core1_1.m_axi_seq_len:bank1 \
--sp ${KEN_NME}_core1_1.m_axi_mem:bank1 \
--sp ${KEN_NME}_core1_1.m_axi_mem_num:bank1 \
--sp ${KEN_NME}_core1_1.m_axi_bwt:bank1 \
--sp ${KEN_NME}_core2_1.m_axi_bwt_para:bank2 \
--sp ${KEN_NME}_core2_1.m_axi_seq:bank2 \
--sp ${KEN_NME}_core2_1.m_axi_seq_len:bank2 \
--sp ${KEN_NME}_core2_1.m_axi_mem:bank2 \
--sp ${KEN_NME}_core2_1.m_axi_mem_num:bank2 \
--sp ${KEN_NME}_core2_1.m_axi_bwt:bank2 \
--sp ${KEN_NME}_core3_1.m_axi_bwt_para:bank3 \
--sp ${KEN_NME}_core3_1.m_axi_seq:bank3 \
--sp ${KEN_NME}_core3_1.m_axi_seq_len:bank3 \
--sp ${KEN_NME}_core3_1.m_axi_mem:bank3 \
--sp ${KEN_NME}_core3_1.m_axi_mem_num:bank3 \
--sp ${KEN_NME}_core3_1.m_axi_bwt:bank3"


xocc -t hw --kernel_frequency $HLS_FREQ  --platform $PLATFORM -s -o ${KEN_NME}_core0.xo -c ${KEN_SRC} -DCPP -k ${KEN_NME}_core0 
xocc -t hw --kernel_frequency $HLS_FREQ  --platform $PLATFORM -s -o ${KEN_NME}_core1.xo -c ${KEN_SRC} -DCPP -k ${KEN_NME}_core1 
xocc -t hw --kernel_frequency $HLS_FREQ  --platform $PLATFORM -s -o ${KEN_NME}_core2.xo -c ${KEN_SRC} -DCPP -k ${KEN_NME}_core2 
xocc -t hw --kernel_frequency $HLS_FREQ  --platform $PLATFORM -s -o ${KEN_NME}_core3.xo -c ${KEN_SRC} -DCPP -k ${KEN_NME}_core3

xocc -t hw --kernel_frequency $IMP_FREQ --platform $PLATFORM ${IMP_FLG} -s -o ${XCLBIN_NAME} --link ${KEN_NME}_core0.xo ${KEN_NME}_core1.xo ${KEN_NME}_core2.xo ${KEN_NME}_core3.xo


