#ifndef PAIRHMM_TASK_H
#define PAIRHMM_TASK_H

#include <stdlib.h>
#include <stdexcept>
#include <time.h>

#include <glog/logging.h>

#include "blaze/Block.h" 
#include "blaze/Task.h" 

#include "blaze/xlnx_opencl/OpenCLEnv.h" 

class PairHMMInput {
 public:
  PairHMMInput(blaze::OpenCLEnv* env, int bank) {
    //PLACE_TIMER;

    bundle = (FpgaInputBundle*)aligned_alloc(4096, sizeof(FpgaInputBundle));
    memset(bundle, 0, sizeof(bundle));

    // allocate input buffer
    static unsigned bankID[4] = {
      XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1,
      XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};

    cl_mem_ext_ptr_t input_ext;

    input_ext.flags  = bankID[bank];
    input_ext.obj    = bundle;
    input_ext.param  = 0;

    cl_int err = 0;
    cl_context context = env->getContext();
    buf = clCreateBuffer(context, 
        CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
        sizeof(FpgaInputBundle), &input_ext, &err);

    if (!buf || err != CL_SUCCESS)
      throw blaze::invalidParam("failed to allocate CL buffer for input");

    DVLOG(1) << "create one pairhmm input for bank: " << bank;
  }

  ~PairHMMInput() {
    clReleaseMemObject(buf);
    free(bundle);  
    DVLOG(1) << "free one pairhmm input";
  }
  
  cl_mem           buf;
  FpgaInputBundle* bundle;
};

typedef boost::shared_ptr<PairHMMInput> PairHMMInput_ptr;

class PairHMM : public blaze::Task {
 public:
  // extends the base class constructor
  // to indicate how many input blocks
  // are required
  PairHMM();
  virtual ~PairHMM();

  virtual uint64_t estimateClientTime() { return 0; }
  virtual uint64_t estimateTaskTime() { return 0; }

  virtual void compute();
  virtual void prepare();

 private:
  int conf_str2int(std::string key) {
    std::string conf;
    get_conf(key, conf);
    return std::stoi(conf);
  }

  ksight::IntvTimer timer_;
  blaze::OpenCLEnv* env;

  uint64_t num_cell;
  int      num_read;
  int      num_hap;
  PairHMMInput_ptr     input_;
  blaze::DataBlock_ptr output_;
};

// define the constructor and destructor for dlopen()
extern "C" blaze::Task* create() {
  return new PairHMM();
}

extern "C" void destroy(blaze::Task* p) {
  delete p;
}

#endif
