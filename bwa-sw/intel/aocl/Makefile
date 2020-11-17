SHELL = /bin/bash

DATE := $(shell date +'%Y%m%d')
AOCL_VER := $(shell echo "$(ALTERAOCLSDKROOT)" | grep -Eo '[0-9]+([.][0-9]+)+[^/]+')

# Where is the Intel(R) FPGA SDK for OpenCL(TM) software?
ifeq ($(wildcard $(ALTERAOCLSDKROOT)),)
$(error Set ALTERAOCLSDKROOT to the root directory of the Intel(R) FPGA SDK for OpenCL(TM) software installation)
endif
ifeq ($(wildcard $(ALTERAOCLSDKROOT)/host/include/CL/opencl.h),)
$(error Set ALTERAOCLSDKROOT to the root directory of the Intel(R) FPGA SDK for OpenCL(TM) software installation.)
endif

# OpenCL compile and link flags.
AOCL_COMPILE_CONFIG := $(shell aocl compile-config )
AOCL_LINK_CONFIG := $(shell aocl link-config )

# Host compile related stuff
HOST_DEBUG=
HOST_SRCS := ./src/main_cl.cpp
OBJECTS := main_cl.o

HOST = ./bin/host_top
# Compilation flags
ifeq ($(HOST_DEBUG),1)
    HOST_CFLAGS += -g
else
    HOST_CFLAGS += -O2
endif
HOST_CFLAGS += -std=c++0x -fPIC -I../../../common/include -D K_NUM=${K_NUM}

LINK= -lpthread -lm -ldl -lrt -lstdc++ \
	   	-lboost_system \
		-lboost_thread \
		-lboost_iostreams \

# Kernel compile related stuff
K_NUM=2
BLOCK_PE_NUM=10
READ_PE_NUM=2
SEED_PE_NUM=10
KERNEL_SRCS = ./src/smithwaterman.cl
KERNEL_DEFS = -D K_NUM=${K_NUM}
KERNEL_DEFS += -D READ_PE_NUM=${READ_PE_NUM}
KERNEL_DEFS += -D BLOCK_PE_NUM=${BLOCK_PE_NUM}
KERNEL_DEFS += -D SEED_PE_NUM=${SEED_PE_NUM}

PLATFORM=a10gx
KERNEL_DEBUG=
#AOCXBIN_NAME := sw_${K_NUM}x${READ_PE_NUM}x${SEED_PE_NUM}_${AOCL_VER}_${DATE}
AOCXBIN_NAME := sw_${K_NUM}x${READ_PE_NUM}x${SEED_PE_NUM}_${DATE}
AOCX_HW := ./bin/${AOCXBIN_NAME}.aocx
AOCX_PROF_HW := ./bin/${AOCXBIN_NAME}_prof.aocx
AOCO_HW := ./bin/${AOCXBIN_NAME}.aoco
AOCX_SWEMU := ./bin/${AOCXBIN_NAME}_sw_emu.aocx

HOST_ARGS = $(PWD)/testdata/input $(PWD)/testdata/golden

CC = gcc
CLCC = aoc
CLCC_INC = -I $$ALTERAOCLSDKROOT/include/kernel_headers/

.PHONY: default all clean cleanall bit sim run run_sw_emu host est

default:
	@echo "Usage"
	@echo "$$ make <target>"
	@echo "Available targets:"
	@echo "          host - Compile Host Code"
	@echo "        sw_emu - Generate kernel binary for simulation on CPU"
	@echo "    run_sw_emu - Run Kernel simulation on host"
	@echo "           est - Get resource and performance estimates"
	@echo "           bit - Generate kernel binary for Acceleration Platform"
	@echo "      bit_prof - Generate kernel binary for Acceleration Platform with Profiling"
	@echo "           run - Run Accelerated Kernel on Platform"
	@echo "         clean - Remove all output products (except for bitgen outputs)"
	
all: ${HOST} ${AOCX_HW}

host: ${HOST}

sw_emu: ${AOCX_SWEMU}	

est: ${AOCO_HW}

bit_prof: ${AOCX_PROF_HW}

bit: ${AOCX_HW}

run_sw_emu: ${HOST} ${AOCX_SWEMU}
	CL_CONTEXT_EMULATOR_DEVICE_ALTERA=1 ${HOST} ${AOCX_SWEMU} ${HOST_ARGS}

run: ${HOST} ${AOCX_HW}
	${HOST} ${AOCX_HW} ${HOST_ARGS}

${HOST}: $(HOST_SRCS)
	$(CC) -o $@ $^ $(HOST_CFLAGS) $(LINK) $(AOCL_COMPILE_CONFIG) $(AOCL_LINK_CONFIG) \
			-L./ \
			-Wl,-rpath=./

%.o: %.cpp
	${CC} ${HOST_CFLAGS} -c $< -o $@

${AOCX_SWEMU}: ${KERNEL_SRCS}
	${CLCC} -march=emulator ${KERNEL_SRCS} $(KERNEL_DEFS) ${CLCC_INC} -o $@ --board ${PLATFORM} -v

${AOCO_HW}: ${KERNEL_SRCS}
	${CLCC} --high-effort -c ${KERNEL_SRCS} ${CLCC_INC} $(KERNEL_DEFS) -o $@ --board ${PLATFORM} -v -v -v --report

${AOCX_PROF_HW}: ${KERNEL_SRCS}
	${CLCC} --high-effort --profile -g0 ${KERNEL_SRCS} ${CLCC_INC} $(KERNEL_DEFS) -o $@ --board ${PLATFORM} -v -v -v --report

${AOCX_HW}: ${KERNEL_SRCS}
	${CLCC} --high-effort ${KERNEL_SRCS} ${CLCC_INC} $(KERNEL_DEFS) -o $@ --board ${PLATFORM} -v -v -v --report

# by default let's only clean host compile files
clean:
	${RM} -rvf ${HOST} ${OBJECTS} bin/*sw_emu.aoc* bin/*sw_emu

cleanall: clean
	${RM} -rf *.xclbin sdaccel_profile_summary.* _xocc_compile* _xocc_link* _xocc_krnl* TempConfig emconfig.json _xocc_*_${XCLBIN_NAME}_*.dir .Xil
