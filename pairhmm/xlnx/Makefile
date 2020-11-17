CC=gcc
PP=g++
XCC=xocc

CPPFLAGS  := -g -fPIC -Wall -O3 -march=native -m64 -std=c++11
CFLAGS 	  := -g -fPIC -Wall -O3 -march=native -m64

CLCOMPILE := -I$(XILINX_XRT)/include -I$(XILINX_SDX)/runtime/include/1_2 -I$(XILINX_SDX)/Vivado_HLS/include 
CLLINK    := -lpthread -lm -ldl -lrt -lstdc++ -L$(XILINX_XRT)/lib -L$(XILINX_SDX)/runtime/lib/x86_64 -lxilinxopencl 

INCLUDES  := -I.

SRC       := $(wildcard ./host/*.cpp)
HDR       := $(wildcard ./host/*.h)
CMM_HDR       := $(wildcard ./common/*.h)
OBJ       := $(patsubst ./host/%.cpp,./host/%.o,$(SRC))

TEST_SRC := pairhmm_test.cpp
TEST_OBJ := pairhmm_test.o
TEST_DST  := test_bin

KEN_DIR   := ./device
KEN_NME   := pmm_core_top
KEN_SRC    = $(KEN_DIR)/pairhmm.cpp
SW_EMU_KEN0_OBJ:= $(KEN_NME)0_sw_emu.xo
SW_EMU_KEN1_OBJ:= $(KEN_NME)1_sw_emu.xo
SW_EMU_KEN2_OBJ:= $(KEN_NME)2_sw_emu.xo
HLS_OBJ:= $(KEN_NME)_hls.xo
HW_KEN0_OBJ:= $(KEN_NME)0_hw.xo
HW_KEN1_OBJ:= $(KEN_NME)1_hw.xo
HW_KEN2_OBJ:= $(KEN_NME)2_hw.xo
SW_EMU_DST:= pairhmm_kernel_sw_emu.xclbin
HW_EMU_DST:= pairhmm_kernel_hw_emu.xclbin

IMP_FLAG= --xp vivado_param:bd.ForceAppCoreUpgrade=1 \
		--xp vivado_prop:run.pfm_dynamic_${KEN_NME}0_1_0_synth_1.STEPS.SYNTH_DESIGN.TCL.PRE=$(shell pwd)/tcl/presynth.tcl \
		--xp vivado_prop:run.pfm_dynamic_${KEN_NME}1_1_0_synth_1.STEPS.SYNTH_DESIGN.TCL.PRE=$(shell pwd)/tcl/presynth.tcl \
		--xp vivado_prop:run.pfm_dynamic_${KEN_NME}2_1_0_synth_1.STEPS.SYNTH_DESIGN.TCL.PRE=$(shell pwd)/tcl/presynth.tcl \
		--xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.TCL.POST=$(shell pwd)/tcl/postopt.tcl \
		--xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.ARGS.DIRECTIVE=Explore \
		--xp vivado_prop:run.impl_1.STEPS.PLACE_DESIGN.ARGS.DIRECTIVE=Explore \
		--xp vivado_prop:run.impl_1.STEPS.PHYS_OPT_DESIGN.IS_ENABLED=true \
		--xp vivado_prop:run.impl_1.STEPS.PHYS_OPT_DESIGN.ARGS.DIRECTIVE=AggressiveExplore \
		--xp vivado_prop:run.impl_1.STEPS.ROUTE_DESIGN.ARGS.DIRECTIVE=Explore

HW_NAME=vcu1525
DSA_MAJOR_VER=5
DSA_MINOR_VER=1

ifeq ($(HW_NAME),aws)
	ifeq ($(DSA_MAJOR_VER),4)
		DIE_NUM = 3
		PLATFORM := xilinx:$(HW_NAME)-vu9p-f1:4ddr-xpr-2pr:$(DSA_MAJOR_VER).$(DSA_MINOR_VER)
		KERNEL0_BANKID = 3
		KERNEL1_BANKID = 1
		KERNEL2_BANKID = 2
		KERNEL_DRAM_MAP=--xp misc:map_connect=add.kernel.$(KEN_NME)0_1.M_AXI_INPUT_DATA.core.OCL_REGION_0.M0$(KERNEL0_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)0_1.M_AXI_OUTPUT_DATA.core.OCL_REGION_0.M0$(KERNEL0_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)1_1.M_AXI_INPUT_DATA.core.OCL_REGION_0.M0$(KERNEL1_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)1_1.M_AXI_OUTPUT_DATA.core.OCL_REGION_0.M0$(KERNEL1_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)2_1.M_AXI_INPUT_DATA.core.OCL_REGION_0.M0$(KERNEL2_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)2_1.M_AXI_OUTPUT_DATA.core.OCL_REGION_0.M0$(KERNEL2_BANKID)_AXI
		IMP_FLAG+= --xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.TCL.PRE=$(shell pwd)/tcl/pblock_3die_aws.tcl $(KERNEL_DRAM_MAP)
		SLR0_PE_NUM = 48
		SLR1_PE_NUM = 40
		SLR2_PE_NUM = 64
	else
		DIE_NUM = 3
		PLATFORM := xilinx_$(HW_NAME)-vu9p-f1-04261818_dynamic_$(DSA_MAJOR_VER)_$(DSA_MINOR_VER)
		KERNEL0_BANKID = 3
		KERNEL1_BANKID = 1
		KERNEL2_BANKID = 0
		KERNEL_DRAM_MAP=--sp $(KEN_NME)0_1.m_axi_input_data:bank$(KERNEL0_BANKID) \
						--sp $(KEN_NME)0_1.m_axi_output_data:bank$(KERNEL0_BANKID) \
						--sp $(KEN_NME)1_1.m_axi_input_data:bank$(KERNEL1_BANKID) \
						--sp $(KEN_NME)1_1.m_axi_output_data:bank$(KERNEL1_BANKID) \
						--sp $(KEN_NME)2_1.m_axi_input_data:bank$(KERNEL2_BANKID) \
						--sp $(KEN_NME)2_1.m_axi_output_data:bank$(KERNEL2_BANKID)
		IMP_FLAG+= --xp param:compiler.userPostSysLinkTcl=$(shell pwd)/tcl/syslink.tcl $(KERNEL_DRAM_MAP)
		SLR0_PE_NUM = 48
		SLR1_PE_NUM = 56
		SLR2_PE_NUM = 80
	endif
else ifeq ($(HW_NAME),huawei)
	KERNEL0_BANKID = 0
	KERNEL1_BANKID = 1
	KERNEL2_BANKID = 3
	ifeq ($(DSA_MAJOR_VER),4)
		DIE_NUM = 3
		PLATFORM := xilinx:$(HW_NAME)-vu9p-fp1:4ddr-xpr:$(DSA_MAJOR_VER).$(DSA_MINOR_VER)
		KERNEL_DRAM_MAP=--xp misc:map_connect=add.kernel.$(KEN_NME)0_1.M_AXI_INPUT_DATA.core.OCL_REGION_0.M0$(KERNEL0_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)0_1.M_AXI_OUTPUT_DATA.core.OCL_REGION_0.M0$(KERNEL0_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)1_1.M_AXI_INPUT_DATA.core.OCL_REGION_0.M0$(KERNEL1_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)1_1.M_AXI_OUTPUT_DATA.core.OCL_REGION_0.M0$(KERNEL1_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)2_1.M_AXI_INPUT_DATA.core.OCL_REGION_0.M0$(KERNEL2_BANKID)_AXI \
		--xp misc:map_connect=add.kernel.$(KEN_NME)2_1.M_AXI_OUTPUT_DATA.core.OCL_REGION_0.M0$(KERNEL2_BANKID)_AXI
		IMP_FLAG+= --xp vivado_prop:run.impl_1.STEPS.OPT_DESIGN.TCL.PRE=$(shell pwd)/tcl/pblock_3die_hwc.tcl $(KERNEL_DRAM_MAP)
		SLR0_PE_NUM = 64
		SLR1_PE_NUM = 24
		SLR2_PE_NUM = 64
	else
		DIE_NUM = 3
		PLATFORM := $(HW_NAME)_vu9p_dynamic_fp1_$(DSA_MAJOR_VER)_$(DSA_MINOR_VER)
		KERNEL_DRAM_MAP=--sp $(KEN_NME)0_1.m_axi_input_data:bank$(KERNEL0_BANKID) \
						--sp $(KEN_NME)0_1.m_axi_output_data:bank$(KERNEL0_BANKID) \
						--sp $(KEN_NME)1_1.m_axi_input_data:bank$(KERNEL1_BANKID) \
						--sp $(KEN_NME)1_1.m_axi_output_data:bank$(KERNEL1_BANKID) \
						--sp $(KEN_NME)2_1.m_axi_input_data:bank$(KERNEL2_BANKID) \
						--sp $(KEN_NME)2_1.m_axi_output_data:bank$(KERNEL2_BANKID)
		IMP_FLAG+= --xp param:compiler.userPostSysLinkTcl=$(shell pwd)/tcl/syslink.tcl $(KERNEL_DRAM_MAP)
		SLR0_PE_NUM = 80
		SLR1_PE_NUM = 40
		SLR2_PE_NUM = 80
	endif
else ifeq ($(HW_NAME),vcu1525)
	KERNEL0_BANKID = 0
	KERNEL1_BANKID = 1
	KERNEL2_BANKID = 3
	DIE_NUM = 3
	PLATFORM := xilinx_$(HW_NAME)_dynamic_$(DSA_MAJOR_VER)_$(DSA_MINOR_VER) 
	KERNEL_DRAM_MAP=--sp $(KEN_NME)0_1.m_axi_input_data:bank$(KERNEL0_BANKID) \
						--sp $(KEN_NME)0_1.m_axi_output_data:bank$(KERNEL0_BANKID) \
						--sp $(KEN_NME)1_1.m_axi_input_data:bank$(KERNEL1_BANKID) \
						--sp $(KEN_NME)1_1.m_axi_output_data:bank$(KERNEL1_BANKID) \
						--sp $(KEN_NME)2_1.m_axi_input_data:bank$(KERNEL2_BANKID) \
						--sp $(KEN_NME)2_1.m_axi_output_data:bank$(KERNEL2_BANKID)
	IMP_FLAG+= --xp param:compiler.userPostSysLinkTcl=$(shell pwd)/tcl/syslink.tcl $(KERNEL_DRAM_MAP)
	SLR0_PE_NUM = 80
	SLR1_PE_NUM = 40
	SLR2_PE_NUM = 80
endif

HW_DST = bin/$(HW_NAME).xclbin
TEST_FILES = /tmp/pairHmm_testcases_WES

FPGA = 1
BOARD_SETTING = 1

ifeq ($(FPGA), 1)
	FPGA_FLAG = -DFPGA -D DIE_NUM=$(DIE_NUM) -D SLR0_PE_NUM=$(SLR0_PE_NUM) -D SLR1_PE_NUM=$(SLR1_PE_NUM) \
				-D SLR2_PE_NUM=$(SLR2_PE_NUM) -D KERNEL0_BANKID=$(KERNEL0_BANKID) \
				-D KERNEL1_BANKID=$(KERNEL1_BANKID) -D KERNEL2_BANKID=$(KERNEL2_BANKID) \
				-DPLATFORM='"$(HW_NAME)"' \
				-DDSA='"$(DSA_MAJOR_VER).$(DSA_MINOR_VER)"'
else
	FPGA_FLAG = -D DIE_NUM=$(DIE_NUM) -D SLR0_PE_NUM=$(SLR0_PE_NUM) -D SLR1_PE_NUM=$(SLR1_PE_NUM) \
				-D SLR2_PE_NUM=$(SLR2_PE_NUM) -D KERNEL0_BANKID=$(KERNEL0_BANKID) \
				-D KERNEL1_BANKID=$(KERNEL1_BANKID) -D KERNEL2_BANKID=$(KERNEL2_BANKID) \
				-DPLATFORM='"$(HW_NAME)"' \
				-DDSA='"$(DSA_MAJOR_VER).$(DSA_MINOR_VER)"'
endif

ifeq ($(BOARD_SETTING), 0)
	SIM_FLAG = -DSIM
else
	SIM_FLAG =
endif

ifneq ("$(wildcard $(HW_DST))","")
	HW_DST_EXISTS = 1
else
	HW_DST_EXISTS = 0
endif

ifneq ("$(wildcard $(TEST_FILES))","")
	TEST_FILE_EXISTS = 1
else
	TEST_FILE_EXISTS = 0
endif


TEST_DST  := test_bin

all: $(DST)

help:
	@echo "make all:     make the falcon .so file only."; 
	@echo "make HW_NAME="FPGA device name" test:    run the execution. If FPGA=1 and BOARD_SETTING=1, will run the bitstream generation and on-board execution, takes many hours! Be careful about this. If FPGA=1 and BOARD_SETTING=0, will run xocc software emulation flow. If FPGA = 0 and BOARD_SETTING=0, will run everything on AVX";
	@echo "make hw_emu"  run hardware emulation
	@echo "make hls:     run the HLS only, use this to do hardware resource utilization estimation";
	@echo "make bitgen: run the bitstream generation, which takes many hours";
	@echo "make clean:  remove all the temp files and destination files.";
	@echo "make HW_NAME="FPGA device name" host:    make host";
	@echo "Availabe FPGA devices aws: xilinx:aws-vu9p-f1:4ddr-xpr-2pr:4.0 and xilinx_aws-vu9p-f1-04261818_dynamic_5_0"
	@echo "                      hwc: xilinx:huawei-vu9p-fp1:4ddr-xpr:4.1 and huawei_vu9p_dynamic_fp1_5_0"
	@echo "                  vcu1525: xilinx_vcu1525_dynamic_5_0"

emconfig.json:
	emconfigutil --platform $(PLATFORM) --od .

ifeq ($(FPGA), 1)
ifeq ($(BOARD_SETTING), 1)
test: $(TEST_DST) $(DST) $(HW_DST) $(TEST_FILES)
	./test_bin $(HW_DST) -h
else
test: emconfig.json $(TEST_DST) $(DST) $(SW_EMU_DST) $(TEST_FILES)
	XCL_EMULATION_MODE=sw_emu ./test_bin -h
endif
else
test: $(TEST_DST) $(DST) $(TEST_FILES)
	./test_bin
endif

hw_emu: emconfig.json $(HW_EMU_DST) $(DST) $(TEST_DST) $(TEST_FILES)
	XCL_EMULATION_MODE=hw_emu ./test_bin -h

$(HW_EMU_DST): $(HW_KEN0_OBJ) $(HW_KEN1_OBJ) $(HW_KEN2_OBJ)
	$(XCC) -t hw_emu -g --kernel_frequency 200 --platform $(PLATFORM) $(SIM_FLAG) $(KERNEL_DRAM_MAP) -o $(HW_EMU_DST) --link $(HW_KEN0_OBJ) $(HW_KEN1_OBJ) $(HW_KEN2_OBJ) -I.

hls: $(HLS_OBJ)
	@echo "Run HLS"

bitgen: $(HW_DST)
	@echo "Run bit gen"

host: $(TEST_DST) $(DST)

$(TEST_FILES):
ifeq ($(TEST_FILE_EXISTS), 0)
	cd /tmp; pwd; aws s3 cp s3://fcs-genome-data/data-suite/pmm_tb/pmm_testcases_WES.tgz .; tar xvf pmm_testcases_WES.tgz
endif


$(TEST_DST): $(TEST_OBJ) $(OBJ)
	$(PP) -o $@ $(TEST_OBJ) $(OBJ) $(CLLINK)

$(TEST_OBJ): $(TEST_SRC) $(CMM_HDR) $(HDR) 
	$(CC) $(CPPFLAGS) $(FPGA_FLAG) $(SIM_FLAG) $(INCLUDES) $(CLCOMPILE) -c $< -o $@



$(SW_EMU_DST): $(SW_EMU_KEN0_OBJ) $(SW_EMU_KEN1_OBJ) $(SW_EMU_KEN2_OBJ) 
	$(XCC) -t sw_emu --kernel_frequency 250 --platform $(PLATFORM) $(KERNEL_DRAM_MAP) -o $(SW_EMU_DST) --link $(SW_EMU_KEN0_OBJ) $(SW_EMU_KEN1_OBJ) $(SW_EMU_KEN2_OBJ) -I.

$(SW_EMU_KEN0_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t sw_emu --kernel_frequency 250 --platform $(PLATFORM) -D SLR_PE_NUM=$(SLR0_PE_NUM) $(SIM_FLAG) -D SLRID=0 -o $(SW_EMU_KEN0_OBJ) -c $(KEN_SRC) -k $(KEN_NME)0 -I.

$(SW_EMU_KEN1_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t sw_emu --kernel_frequency 250 --platform $(PLATFORM) -D SLR_PE_NUM=$(SLR1_PE_NUM) $(SIM_FLAG) -D SLRID=1 -o $(SW_EMU_KEN1_OBJ) -c $(KEN_SRC) -k $(KEN_NME)1 -I.

$(SW_EMU_KEN2_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t sw_emu --kernel_frequency 250 --platform $(PLATFORM) -D SLR_PE_NUM=$(SLR2_PE_NUM) $(SIM_FLAG) -D SLRID=2 -o $(SW_EMU_KEN2_OBJ) -c $(KEN_SRC) -k $(KEN_NME)2 -I.

ifeq ($(HW_DST_EXISTS), 0)
$(HW_DST): $(HW_KEN0_OBJ) $(HW_KEN1_OBJ) $(HW_KEN2_OBJ) 
	echo "hw .xclbin file does not exist, start the full hw compilation flow, which takes multiple hours!!"
	$(XCC) -t hw --kernel_frequency 250 --platform $(PLATFORM) $(IMP_FLAG) -s -o $(HW_NAME).xclbin --link $(HW_KEN0_OBJ) $(HW_KEN1_OBJ) $(HW_KEN2_OBJ)
else
$(HW_DST): 
	echo "hw .xclbin file exists"
endif

$(HLS_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t hw -s -g --kernel_frequency 250 --platform $(PLATFORM) -D SLR_PE_NUM=8 -D SLRID=0 --optimize 3 -o $(HLS_OBJ) -c $(KEN_SRC) -k $(KEN_NME)0 -I.

$(HW_KEN0_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t hw -s -g --kernel_frequency 200 --platform $(PLATFORM) -D SLR_PE_NUM=$(SLR0_PE_NUM) -D SLRID=0 --optimize 3 -o $(HW_KEN0_OBJ) -c $(KEN_SRC) -k $(KEN_NME)0 -I.

$(HW_KEN1_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t hw -s -g --kernel_frequency 200 --platform $(PLATFORM) -D SLR_PE_NUM=$(SLR1_PE_NUM) -D SLRID=1 --optimize 3 -o $(HW_KEN1_OBJ) -c $(KEN_SRC) -k $(KEN_NME)1 -I.

$(HW_KEN2_OBJ): $(KEN_SRC) $(CMM_HDR)
	$(XCC) -t hw -s -g --kernel_frequency 200 --platform $(PLATFORM) -D SLR_PE_NUM=$(SLR2_PE_NUM) -D SLRID=2 --optimize 3 -o $(HW_KEN2_OBJ) -c $(KEN_SRC) -k $(KEN_NME)2 -I.


host/%.o: host/%.cpp $(HDR) $(CMM_HDR)
	$(CC) -c $(CPPFLAGS) $(FPGA_FLAG) $(SIM_FLAG) $(INCLUDES) $(CLCOMPILE)$< -o $@

clean: 
	rm -rf $(OBJ) $(SW_EMU_KEN0_OBJ) $(SW_EMU_KEN1_OBJ) $(SW_EMU_KEN2_OBJ) $(TEST_OBJ) $(HW_KEN0_OBJ) $(HW_KEN1_OBJ) $(HW_KEN2_OBJ)
	rm -rf $(TEST_DST) $(SW_EMU_DST) $(HW_EMU_DST)
	rm -rf .Xil/ sdaccel_* *.dir/
	rm -rf _x/
	rm -rf .run/
	rm -rf *.json
	rm -rf *.log

