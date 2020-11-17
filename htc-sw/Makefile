PP=g++
CC=gcc
XCC=xocc

CPPFLAGS  := -fPIC -Wall -O3 -march=native -m64 -std=c++11
CFLAGS 	  := -fPIC -Wall -O3 -march=native -m64

CLCOMPILE := -fPIC -Wall -I$(XILINX_SDX)/runtime/include/1_2 -I$(XILINX_SDX)/Vivado_HLS/include
CLLINK    := -lpthread -lm -ldl -lrt -lstdc++ -L$(XILINX_SDX)/runtime/lib/x86_64 -lxilinxopencl 



INCLUDES  := -I.

LIBS	  := -lpthread -lm

SRC 	  := $(wildcard host/*.cpp)
OBJ       := $(patsubst host/%.cpp,host/%.o,$(SRC))

INTEL_AVX_DIR = intel_avx
INTEL_AVX_SRC := $(wildcard intel_avx/*.cc)
INTEL_AVX_HDR := $(wildcard intel_avx/*.h)
INTEL_AVX_OBJ := $(patsubst intel_avx/%.cc,intel_avx/%.o,$(INTEL_AVX_SRC)) 

KEN_DIR   := ./xlnx
KEN_NME   := smithWatermanMerlin
KEN_SRC   := $(KEN_DIR)/kernel_$(KEN_NME).cpp
SW_EMU_KEN_OBJ:= $(KEN_NME)_sw_emu.xo
HW_KEN_OBJ:= $(KEN_NME)_hw.xo
SW_EMU_DST:= sw_kernel_sw_emu.xclbin
HW_DST:= sw_kernel_hw.xclbin
PLATFORM  := xilinx:xil-accel-rd-ku115:4ddr-xpr:4.0

#FPGA_FLAG := -DFPGA #remove this flag if only AVX


TEST_DST := test_bin
TEST_SW_OBJ := sw_host.o
TEST_SW_DST := sw_host
TEST_SW_SRC := host/sw_host.cpp

ifneq ($(shell uname -s ), Darwin)
DST := lib_htc_sw.so
else # assume it's mac
DST := lib_htc_sw.dylib
endif

all: $(DST)

help:
	@echo "make all: make the falcon .so file only."; 
	@echo "make test_sw_emu: run the software emulation, if you want to turn off the xilinx software emulation and run AVX only, You can comment out the FPGA_FLAG."; 
	@echo "make test_hw: run the hardware compilation, this will take many hours. Run it only once.";
	@echo "make clean: remove all the temp files and DST files.";

sw_emu: $(TEST_SW_DST) $(SW_EMU_DST)

hw: $(TEST_SW_DST) $(HW_DST)

test_sw_emu: $(TEST_SW_DST) $(SW_EMU_DST) $(DST)
	XCL_EMULATION_MODE=sw_emu ./sw_host $(SW_EMU_DST)

$(SW_EMU_DST): $(SW_EMU_KEN_OBJ) 
	$(XCC) -t sw_emu --platform $(PLATFORM) -DCUSTOM -o $(SW_EMU_DST) --link $(SW_EMU_KEN_OBJ)

$(SW_EMU_KEN_OBJ): $(KEN_SRC)
	$(XCC) -t sw_emu --platform $(PLATFORM) -DCUSTOM -o $(SW_EMU_KEN_OBJ) -c $(KEN_SRC) -k $(KEN_NME)

#Warning: this will take many hours!!
test_hw_emu: $(TEST_SW_DST) $(HW_DST) $(DST)
	./sw_host $(HW_DST)

$(HW_DST): $(HW_KEN_OBJ) 
	$(XCC) -t hw --platform $(PLATFORM) -DCUSTOM -o $(HW_DST) --link $(HW_KEN_OBJ)

$(HW_KEN_OBJ): $(KEN_SRC)
	$(XCC) -t hw --platform $(PLATFORM) -DCUSTOM -o $(SW_EMU_KEN_OBJ) -c $(KEN_SRC) -k $(KEN_NME)



$(TEST_SW_DST): $(TEST_SW_OBJ) $(DST)
	$(CC) -o $(TEST_SW_DST) $(TEST_SW_OBJ) $(DST)

$(TEST_SW_OBJ): $(TEST_SW_SRC)
	$(CC) -c $(FPGA_FLAG) $(INCLUDES) -o $@ $<

$(TEST_DST): $(TEST_OBJ) $(OBJ) $(DST)
	$(PP) -o $@ $(OBJ) $(TEST_OBJ) $(LIBS) 

$(DST):	$(OBJ) $(HDR) $(INTEL_AVX_OBJ) 
	$(PP) -shared $(CLLINK) $(OBJ) $(INTEL_AVX_OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) -c $(CLCOMPILE) $(CPPFLAGS) $(INCLUDES) $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(INTEL_AVX_OBJ): $(INTEL_AVX_SRC) $(INTEL_AVX_HDR)
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@



.PHONY: all clean

clean:
	rm -rf $(DST) $(TEST_DST) $(TEST_SW_DST) $(SW_EMU_DST)
	rm -rf $(OBJ) $(TEST_OBJ) $(TEST_SW_OBJ) $(SW_EMU_KEN_OBJ) $(HW_KEN_OBJ)
	rm -rf *.dir/ sdaccel_profile_summary.* .Xil/
