#ifndef BWAOCLENV_H
#define BWAOCLENV_H

#include <string>
#include <stdexcept>

#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "acc_lib/XCLEnv.h" // load xilinx opencl framework

class BWAOCLEnv : public OpenCLEnv {
 public:
  BWAOCLEnv(const char* bin_path, 
      const char* kernel_name,
      int64_t pac_size,
      char* pac): OpenCLEnv(bin_path, kernel_name) 
  {
    int err = 0;
    cl_mem_ext_ptr_t ext_c, ext_d;
#ifdef DUAL_BANK
    ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
    ext_d.flags = XCL_MEM_DDR_BANK3; ext_d.obj = 0; ext_d.param = 0;
#else
    ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
    ext_d.flags = XCL_MEM_DDR_BANK1; ext_d.obj = 0; ext_d.param = 0;
#endif

    cl_context context = this->getContext();
    cl_command_queue commands = this->getCmdQueue();

    pac_input_a_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
        pac_size, &ext_c, &err);
    pac_input_b_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
        pac_size, &ext_d, &err);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("Failed to create reference OpenCL buffer!");
    }
    cl_event event[2];
    err = clEnqueueWriteBuffer(commands, pac_input_a_, CL_TRUE, 0, pac_size, pac, 0, NULL, &event[0]);
    err = clEnqueueWriteBuffer(commands, pac_input_b_, CL_TRUE, 0, pac_size, pac, 0, NULL, &event[1]);
    clWaitForEvents(2, event);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("Failed to write reference to DDR!");
    }
    clReleaseEvent(event[0]);
    clReleaseEvent(event[1]);
  }

  ~BWAOCLEnv(){
    clReleaseMemObject(pac_input_a_);
    clReleaseMemObject(pac_input_b_);
  }

  cl_mem pac_input_a_;
  cl_mem pac_input_b_;
};
#endif
