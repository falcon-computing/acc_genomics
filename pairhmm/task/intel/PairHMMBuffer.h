#ifndef PAIRHMMBUFFER_H
#define PAIRHMMBUFFER_H
#include <CL/opencl.h>
#include <glog/logging.h>
#include <stdlib.h>

#include "blaze/altr_opencl/OpenCLEnv.h" 

class PairHMMBuffer {
 public:
  PairHMMBuffer(blaze::OpenCLEnv *env, 
      size_t max_size, 
      cl_int flag = CL_MEM_READ_ONLY): env_(env), max_size_(max_size)
  {
    cl_int err = 0;
    cl_context context = env->getContext();

    //h_buf_ = aligned_alloc(4096, max_size);
    int rc = posix_memalign(&h_buf_, 64, max_size);
    d_buf_ = clCreateBuffer(context, flag, 
                  max_size, NULL, &err);

    if (!d_buf_ || err != CL_SUCCESS) {
      throw blaze::invalidParam("failed to allocate CL buffer for tile results data");
    }

    DVLOG(1) << "create one set of pairhmm tile buffer";
  }

  ~PairHMMBuffer() {
    clReleaseMemObject(d_buf_);
    free(h_buf_);

    DVLOG(1) << "free one set of pairhmm tile buffer";
  }

  void write_buffer(size_t size = 0, int q = 0) {
    if (size == 0) {
      size = max_size_;
    }
    if (size > max_size_) {
      DVLOG(1) << "size: " << size << ", max_size: " << max_size_;
      throw blaze::invalidParam("maximum buffer size exceeded");
    }

    cl_command_queue cmd = env_->getCmdQueue(q);
    cl_int status = clEnqueueWriteBuffer(cmd,
                        d_buf_,
                        CL_FALSE,
                        0,
                        size,
                        h_buf_,
                        0, NULL, NULL);

    if (status != CL_SUCCESS) {
      throw blaze::invalidParam("fails to write buffer");
    }
  }

  void read_buffer(size_t size, int q) {
    if (size == 0) {
      size = max_size_;
    }
    if (size > max_size_) {
      DVLOG(1) << "size: " << size << ", max_size: " << max_size_;
      throw blaze::invalidParam("maximum buffer size exceeded");
    }

    cl_command_queue cmd = env_->getCmdQueue(q);
    cl_int status = clEnqueueReadBuffer(cmd,
                        d_buf_,
                        CL_FALSE,
                        0,
                        size,
                        h_buf_,
                        0, NULL, NULL);

    if (status != CL_SUCCESS) {
      throw blaze::invalidParam("fails to write buffer");
    }
  }

  void*  get_host_buffer() { return h_buf_; }
  cl_mem get_device_buffer() { return d_buf_; }

 private:
  blaze::OpenCLEnv *env_;

  void*   h_buf_;
  cl_mem  d_buf_;
  size_t  max_size_;
};
#endif
