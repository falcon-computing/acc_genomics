SHELL = /bin/bash
KERNEL_DEBUG=1

ifndef XILINX_SDX
$(error Environment variable XILINX_SDX is required and should point to SDAccel install area)
endif

ifndef XILINX_VIVADO
$(error Environment variable XILINX_VIVADO is required and should point to SDAccel install area)
endif

BWAFLOW_DIR=/curr/diwu/prog/bwa-flow
HTSLIB_DIR=/curr/diwu/tools/htslib-1.3.1

HOST_SRCS := main.cpp
HOST_HDRS := $(wildcard ./host/*.h)
CMM_HDRS := $(wildcard ./common/*.h)
OBJECTS := main.o host/baseline.o host/ocl.o

BANK_NUM=3

HOST = ./host_top
HOST_CFLAGS= -g -O3 -fPIC -std=c++11 \
	     -I$(BWAFLOW_DIR) \
	     -I$(HTSLIB_DIR) \
	     -I$(XILINX_SDX)/runtime/include/1_2 \
	     -I$(XILINX_VIVADO)/include \
	     -I. -D BANK_NUM=$(BANK_NUM)


ifeq (${EMU_TARGET},sw_emu)
	HOST_CFLAGS += -DCL
else 
	HOST_CFLAGS += -DCPP
endif

LINK= -L$(XILINX_SDX)/runtime/lib/x86_64 -lxilinxopencl\
      -L$(BWAFLOW_DIR)/bwa -lbwa \
      -L$(HTSLIB_DIR) -lhts \
      -lpthread -lm -ldl -lrt -lstdc++ -lz

CL_KERNEL_SRCS = device/smem.cl
CPP_KERNEL_SRCS = device/smem.cpp
KERNEL_SWEMU_OBJ = smem_sw_emu.xo
KERNEL_HW_OBJ0 = smem_hw0.xo
KERNEL_HW_OBJ1 = smem_hw1.xo
KERNEL_HW_OBJ2 = smem_hw2.xo
KERNEL_HW_OBJ3 = smem_hw3.xo

KERNEL_HWEMU_OBJ0 = smem_hw_emu0.xo
KERNEL_HWEMU_OBJ1 = smem_hw_emu1.xo
KERNEL_HWEMU_OBJ2 = smem_hw_emu2.xo
KERNEL_HWEMU_OBJ3 = smem_hw_emu3.xo

KERNEL_NAME = mem_collect_intv_core
KERNEL_DEFS = 
KERNEL_INCS = 
KERNEL_FREQ = 300
CLCC_OPT_LEVEL=
#set target device for XCLBIN
#aws or huawei or vcu1525
XDEVICE=aws 
DSA_MAJOR_VER=5
DSA_MINOR_VER=0

ifeq (${XDEVICE},vcu1525)
	PLATFORM=xilinx_$(XDEVICE)_dynamic_$(DSA_MAJOR_VER)_$(DSA_MINOR_VER)
else ifeq (${XDEVICE},aws)
	PLATFORM=xilinx_$(XDEVICE)-vu9p-f1-04261818_dynamic_$(DSA_MAJOR_VER)_$(DSA_MINOR_VER)
else ifeq (${XDEVICE},huawei)
	PLATFORM=$(XDEVICE)_vu9p_dynamic_fp1_$(DSA_MAJOR_VER)_$(DSA_MINOR_VER)
endif

KEEP_TEMP=1
KERNEL_DEBUG=
XCLBIN_NAME=smem
DATE=$(shell date +%m%d)
XCLBIN_HW = ${XCLBIN_NAME}_${PLATFORM}_${DATE}.xclbin
XCLBIN_SWEMU = ${XCLBIN_NAME}_sw_emu.xclbin
XCLBIN_HWEMU = ${XCLBIN_NAME}_hw_emu.xclbin

TEST_FILES=/tmp/smem_testcases_WGS
#REF_FILES=/genome/ref/human_g1k_v37.fasta
REF_FILES=/genome/ref/v38/Homo_sapiens_assembly38.fasta

HOST_ARGS_SW_EMU = ${TEST_FILES}/ 1024 $(shell echo $(BANK_NUM)\* 16 | bc)
HOST_ARGS_HW_EMU = ${TEST_FILES}/ 4 $(shell echo $(BANK_NUM)\* 9 | bc)
HOST_ARGS_HW = ${TEST_FILES}/ 1024 $(shell echo $(BANK_NUM)\* 256 | bc)

CC = gcc
CLCC = xocc
CLCC_OPT += $(CLCC_OPT_LEVEL) --platform ${PLATFORM} ${KERNEL_DEFS} ${KERNEL_INCS} 

IMP_OPT = --xp vivado_param:bd.ForceAppCoreUpgrade=1 \
		--xp param:compiler.userPostSysLinkTcl=${PWD}/tcl/syslink.tcl \
		--xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.TCL.POST=${PWD}/tcl/postopt.tcl \
		--xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.ARGS.DIRECTIVE=Explore \
		--xp vivado_prop:run.impl_1.STEPS.PLACE_DESIGN.ARGS.DIRECTIVE=Explore \
		--xp vivado_prop:run.impl_1.STEPS.PHYS_OPT_DESIGN.IS_ENABLED=true \
		--xp vivado_prop:run.impl_1.STEPS.PHYS_OPT_DESIGN.ARGS.DIRECTIVE=AggressiveExplore \
		--xp vivado_prop:run.impl_1.STEPS.ROUTE_DESIGN.ARGS.DIRECTIVE=Explore

ifneq (${EMU_TARGET},sw_emu)
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_bwt:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_bwt:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_bwt:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_bwt:bank3
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_seq:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_seq:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_seq:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_seq:bank3
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_mem:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_mem:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_mem:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_mem:bank3
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_bwt_para:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_bwt_para:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_bwt_para:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_bwt_para:bank3
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_mem_num:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_mem_num:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_mem_num:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_mem_num:bank3
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_seq_len:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_seq_len:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_seq_len:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_seq_len:bank3
else 
	CLCC_OPT += --sp $(KERNEL_NAME)0_1.m_axi_gmem:bank0 \
				--sp $(KERNEL_NAME)1_1.m_axi_gmem:bank1 \
				--sp $(KERNEL_NAME)2_1.m_axi_gmem:bank2 \
				--sp $(KERNEL_NAME)3_1.m_axi_gmem:bank3
endif

CLCC_OPT += -I. -D BANK_NUM=$(BANK_NUM)

ifeq (${KEEP_TEMP},1)
    CLCC_OPT += -s
endif

ifeq (${KERNEL_DEBUG},1)
    CLCC_OPT += -g
endif


.PHONY: all clean cleanall bit sim hw_sim run hls run_emu run_sw_emu run_hw_emu

data: ${TEST_FILES} ${REF_FILES}

all: ${HOST} ${XCLBIN_HW}

sim: 	${XCLBIN_SWEMU}	

hw_sim: ${XCLBIN_HWEMU}

bit: 	${XCLBIN_HW}

run: ${HOST} ${XCLBIN_HW}
	${HOST} ${XCLBIN_HW} ${HOST_ARGS_HW}

${TEST_FILES}:
ifeq ("$(wildcard ${TEST_FILES})","")
	cd /tmp; pwd; aws s3 cp s3://fcs-genome-data/data-suite/smem_tb/smem_testcases_WGS.tgz .; tar xvf smem_testcases_WGS.tgz
endif


ifeq (${EMU_TARGET},sw_emu)
run_emu: run_sw_emu host
else
run_emu: run_hw_emu host
endif

run_sw_emu: emconfig.json ${HOST} ${XCLBIN_SWEMU} ${TEST_FILES}
	XCL_EMULATION_MODE=sw_emu ${HOST} ${XCLBIN_SWEMU} ${REF_FILES} ${HOST_ARGS_SW_EMU}

run_hw_emu: emconfig.json ${HOST} ${XCLBIN_HWEMU} ${TEST_FILES}
	echo hw_emu requires sdx 18.2 and above, please check vivado_hls/vivado version
	XCL_EMULATION_MODE=hw_emu ${HOST} ${XCLBIN_HWEMU} ${REF_FILES} ${HOST_ARGS_HW_EMU}

hls: ${KERNEL_HWEMU_OBJ0} 
	echo hw_emu requires sdx 18.2 and above, please check vivado_hls/vivado version

host: ${HOST}

${HOST}: $(OBJECTS)
	$(CC) -o $@ $^ $(HOST_CFLAGS) $(LINK)  \
			-L./ \
			-Wl,-rpath=./

%.o: %.cpp
	${CC} ${HOST_CFLAGS} -c $< -o $@

host/%.o: host/%.cpp $(HOST_HDRS) $(CMM_HDRS) .FORCE
	${CC} ${HOST_CFLAGS} -c $< -o $@

.FORCE:

device/%.o: device/%.cpp
	${CC} ${HOST_CFLAGS} -c $< -o $@

emconfig.json:
	emconfigutil --platform ${PLATFORM} ${DEVICE_REPO_OPT} --od .

${XCLBIN_SWEMU}: ${KERNEL_SWEMU_OBJ}
	${CLCC} -o $@ ${CLCC_OPT} -t sw_emu --link $< --kernel_frequency ${KERNEL_FREQ} --xp param:compiler.version=31

${KERNEL_SWEMU_OBJ}: ${CL_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t sw_emu -c $< --kernel_frequency ${KERNEL_FREQ} -DCL --xp param:compiler.version=31

${XCLBIN_HWEMU}: ${KERNEL_HWEMU_OBJ0} ${KERNEL_HWEMU_OBJ1} ${KERNEL_HWEMU_OBJ2} ${KERNEL_HWEMU_OBJ3}
	${CLCC} -o $@ ${CLCC_OPT} -t hw_emu --link ${KERNEL_HWEMU_OBJ0} ${KERNEL_HWEMU_OBJ1} ${KERNEL_HWEMU_OBJ2} ${KERNEL_HWEMU_OBJ3} --kernel_frequency ${KERNEL_FREQ}

${KERNEL_HWEMU_OBJ0}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -g -t hw_emu -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}0 -DCPP --custom_script ${KERNEL_NAME}0:./tcl/${KERNEL_NAME}0_hw.tcl

${KERNEL_HWEMU_OBJ1}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -g -t hw_emu -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}1 -DCPP --custom_script ${KERNEL_NAME}1:./tcl/${KERNEL_NAME}1_hw.tcl

${KERNEL_HWEMU_OBJ2}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -g -t hw_emu -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}2 -DCPP --custom_script ${KERNEL_NAME}2:./tcl/${KERNEL_NAME}2_hw.tcl

${KERNEL_HWEMU_OBJ3}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -g -t hw_emu -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}3 -DCPP --custom_script ${KERNEL_NAME}3:./tcl/${KERNEL_NAME}3_hw.tcl


${KERNEL_HW_OBJ0}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t hw -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}0 -DCPP

${KERNEL_HW_OBJ1}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t hw -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}1 -DCPP

${KERNEL_HW_OBJ2}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t hw -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}2 -DCPP

${KERNEL_HW_OBJ3}: ${CPP_KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t hw -c $< --kernel_frequency ${KERNEL_FREQ} -k ${KERNEL_NAME}3 -DCPP

${XCLBIN_HW}: ${KERNEL_HW_OBJ0} ${KERNEL_HW_OBJ1} ${KERNEL_HW_OBJ2} ${KERNEL_HW_OBJ3}
	${CLCC} -o $@ ${CLCC_OPT} ${IMP_OPT} -t hw --link ${KERNEL_HW_OBJ0} ${KERNEL_HW_OBJ1} ${KERNEL_HW_OBJ2} ${KERNEL_HW_OBJ3} --kernel_frequency ${KERNEL_FREQ}


# by default let's only clean host compile files
.PHONY: clean cleanall

clean:
	${RM} -rf ${HOST} ${OBJECTS} ${XCLBIN_HWEMU} ${XCLBIN_SWEMU} ${KERNEL_SWEMU_OBJ} ${KERNEL_HWEMU_OBJ0} ${KERNEL_HWEMU_OBJ1} ${KERNEL_HWEMU_OBJ2} ${KERNEL_HWEMU_OBJ3}

cleanall: clean
	${RM} -rf _x .run/ sdaccel_profile* _xocc_compile* _xocc_link* _xocc_krnl* TempConfig emconfig.json _xocc_*_${XCLBIN_NAME}_*.dir .Xil
