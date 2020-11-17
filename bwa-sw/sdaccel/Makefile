SHELL = /bin/bash

ifdef XILINX_SDX
    XILINX_SDACCEL=${XILINX_SDX}
endif

ifndef XILINX_SDACCEL
$(error Environment variable XILINX_SDACCEL or XILINX_SDX is required and should point to SDAccel install area)
endif

BWAFLOW_DIR=/curr/diwu/prog/bwa-flow
HTSLIB_DIR=/curr/diwu/tools/htslib-1.3.1

HOST_SRCS := main_cl.cpp
OBJECTS := main_cl.o smithwaterman.o

HOST = ./host_top
HOST_CFLAGS= -g -O2 -fPIC \
             -DTARGET_DEVICE=\"${XDEVICE}\" -DFPGA_DEVICE -DC_KERNEL \
	     -I../../common/include \
	     -I$(BWAFLOW_DIR) \
	     -I$(HTSLIB_DIR) \
	     -I$(XILINX_SDACCEL)/runtime/include/1_2 \
	     -I$(XILINX_VIVADO)/include
#-DDUAL_BANK \

ifneq ($(DO_CSIM),)
	HOST_CFLAGS := $(HOST_CFLAGS) -DHLS_
endif

LINK= -L$(XILINX_SDACCEL)/runtime/lib/x86_64 -lxilinxopencl \
      -L$(BWAFLOW_DIR)/bwa -lbwa \
      -L$(HTSLIB_DIR) -lhts \
      -lpthread -lm -lz

KERNEL_SRCS = smithwaterman.cpp
KERNEL_OBJS = smithwaterman.xo
KERNEL_NAME = sw_top
KERNEL_DEFS = 
KERNEL_INCS = 
KERNEL_FREQ = 300
CLCC_OPT_LEVEL=
#set target device for XCLBIN

XDEVICE=huawei_vu9p
DSA_MAJOR_VERSION=5
DSA_MINOR_VERSION=0
ifeq (${XDEVICE},ku115)
    PLATFORM=xilinx:xil-accel-rd-$(XDEVICE):4ddr-xpr:$(DSA_MAJOR_VERSION).$(DSA_MINOR_VERSION)
else ifeq (${XDEVICE},vcu1525)
    PLATFORM=xilinx_$(XDEVICE)_dynamic_$(DSA_MAJOR_VERSION)_$(DSA_MINOR_VERSION)
else ifeq (${XDEVICE},aws-vu9p)
	ifeq ($(DSA_MAJOR_VERSION),4)
    	PLATFORM=xilinx:$(XDEVICE)-f1:4ddr-xpr-2pr:$(DSA_MAJOR_VERSION).$(DSA_MINOR_VERSION)
	else 
    	PLATFORM=xilinx_$(XDEVICE)-f1-04261818_dynamic_$(DSA_MAJOR_VERSION)_$(DSA_MINOR_VERSION)
	endif
else ifeq (${XDEVICE},huawei_vu9p)
	ifeq (${DSA_MAJOR_VERSION},4)
		PLATFORM=xilinx:$(XDEVICE)-fp1:4ddr-xpr:$(DSA_MAJOR_VERSION).$(DSA_MINOR_VERSION)
	else 
		PLATFORM=$(XDEVICE)_dynamic_fp1_$(DSA_MAJOR_VERSION)_$(DSA_MINOR_VERSION)
	endif
endif

KEEP_TEMP=1
KERNEL_DEBUG=
XCLBIN_NAME=sw_3x20_dup_$(KERNEL_FREQ)MHz_$(XDEVICE)_$(DSA_MAJOR_VERSION)_$(DSA_MINOR_VERSION)
XCLBIN_HW = ${XCLBIN_NAME}.xclbin
XCLBIN_SWEMU = ${XCLBIN_NAME}_sw_emu.xclbin
XCLBIN_HWEMU = ${XCLBIN_NAME}_hw_emu.xclbin

HOST_ARGS = $(PWD)/testdata/small_gcat_m $(PWD)/testdata/small_gcat_golden

CC = gcc
PP = g++
CLCC = xocc
CLCC_OPT += $(CLCC_OPT_LEVEL) --platform ${PLATFORM} ${KERNEL_DEFS} ${KERNEL_INCS}

ifeq (${KEEP_TEMP},1)
    CLCC_OPT += -s
endif

ifeq (${KERNEL_DEBUG},1)
    CLCC_OPT += -g
endif

CLCC_OPT += --kernel ${KERNEL_NAME}

ifeq (${DSA_MAJOR_VERSION},4)
	CLCC_OPT += --xp misc:map_connect=add.kernel.$(KERNEL_NAME)_1.M_AXI_GMEM_A.core.OCL_REGION_0.M00_AXI \
			--xp misc:map_connect=add.kernel.$(KERNEL_NAME)_1.M_AXI_GMEM_B.core.OCL_REGION_0.M02_AXI \
			--xp misc:map_connect=add.kernel.$(KERNEL_NAME)_1.M_AXI_GMEM_C.core.OCL_REGION_0.M01_AXI \
			--xp misc:map_connect=add.kernel.$(KERNEL_NAME)_1.M_AXI_GMEM_D.core.OCL_REGION_0.M03_AXI
else 
	CLCC_OPT += --sp $(KERNEL_NAME)_1.m_axi_gmem_a:bank0 \
				--sp $(KERNEL_NAME)_1.m_axi_gmem_b:bank2 \
				--sp $(KERNEL_NAME)_1.m_axi_gmem_c:bank1 \
				--sp $(KERNEL_NAME)_1.m_axi_gmem_d:bank3
endif

.PHONY: all clean cleanall bit sim hw_sim csim run run_emu run_sw_emu run_hw_emu

all: ${HOST} ${XCLBIN_HW}

sim: 	${XCLBIN_SWEMU}	

hw_sim: ${XCLBIN_HWEMU}

# for quick verfication using GCC
csim: 
	make host DO_CSIM=1

bit: 	${XCLBIN_HW}

run: ${HOST} ${XCLBIN_HW}
	${HOST} ${XCLBIN_HW} ${HOST_ARGS}

run_emu: run_sw_emu

run_sw_emu: emconfig.json ${HOST} ${XCLBIN_SWEMU}
	XCL_EMULATION_MODE=true ${HOST} ${XCLBIN_SWEMU} ${HOST_ARGS}

run_hw_emu: emconfig.json ${HOST} ${XCLBIN_HWEMU}
	XCL_EMULATION_MODE=true ${HOST} ${XCLBIN_HWEMU} ${HOST_ARGS}
	
host: ${HOST} 

${HOST}: $(OBJECTS)
	$(PP) -o $@ $^ $(HOST_CFLAGS) $(LINK)  \
			-L./ \
			-Wl,-rpath=./


%.o: %.cpp
	${CC} ${HOST_CFLAGS} -c $< -o $@

hls:
	${CLCC} -c -s -t hw_emu -o kernel_top.xo --platform ${PLATFORM} --report estimate --kernel_frequency ${KERNEL_FREQ} --kernel ${KERNEL_NAME} ${KERNEL_SRCS}

emconfig.json:
	emconfigutil --platform ${PLATFORM} ${DEVICE_REPO_OPT} --od .

${XCLBIN_SWEMU}: ${KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t sw_emu ${KERNEL_SRCS} --kernel_frequency ${KERNEL_FREQ}

${XCLBIN_HWEMU}: ${KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -g -t hw_emu ${KERNEL_SRCS} --kernel_frequency ${KERNEL_FREQ}

${XCLBIN_HW}: ${KERNEL_OBJS}
	${CLCC} -o $@ ${CLCC_OPT} -t hw --link ${KERNEL_OBJS} --kernel_frequency ${KERNEL_FREQ}

${KERNEL_OBJS}: ${KERNEL_SRCS}
	${CLCC} -o $@ ${CLCC_OPT} -t hw -c ${KERNEL_SRCS} --kernel_frequency ${KERNEL_FREQ}

# by default let's only clean host compile files
clean:
	${RM} -rf ${HOST} ${OBJECTS} ${XCLBIN_HWEMU} ${XCLBIN_SWEMU}

cleanall: clean
	${RM} -rf *.xclbin sdaccel_profile_summary.* _xocc_compile* _xocc_link* _xocc_krnl* TempConfig emconfig.json _xocc_${XCLBIN_NAME}_*.dir .Xil
