#ifndef XCLENV_H
#define XCLENV_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

class OpenCLEnv {
 public:
  OpenCLEnv(
      const char* bin_path,
      const char* kernel_name)
  {
    // start platform setting up
    int err;

    cl_platform_id platform_id[2];
    cl_device_id device_id;

    char cl_platform_vendor[1001];

    cl_uint num_platforms = 0;

    // Connect to first platform
    err = clGetPlatformIDs(2, platform_id, &num_platforms);

    if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "Failed to find an OpenCL platform!");
    }
    //printf("Found %d OpenCL Platform\n", num_platforms);

    int platform_idx;
    for (int i = 0; i < num_platforms; i++) {
      char cl_platform_name[1001];

      err = clGetPlatformInfo(
          platform_id[i], 
          CL_PLATFORM_NAME, 
          1000, 
          (void *)cl_platform_name,NULL);
      if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!");
      }

      if (strstr(cl_platform_name, "Xilinx") != NULL) {
        // found platform
        //printf("Found Xilinx platform\n");
        platform_idx = i;
        break;
      }
    }

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id[platform_idx], CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id, NULL);

    if (err != CL_SUCCESS) {
        fprintf(stderr, "clGetDeviceIDs return %d\n", err);
        throw std::runtime_error("Failed to create a device group!");
    }

    // Create a compute context 
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);

    if (!context) {
        throw std::runtime_error("Failed to create a compute context!");
    }

    // Create a command commands
    commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);

    if (!commands) {
        throw std::runtime_error("Failed to create a command queue context!");
    }

    // Load binary from disk
    unsigned char *kernelbinary;
    
    int n_i = load_file(bin_path, (char **) &kernelbinary);

    if (n_i < 0) {
        throw std::runtime_error(
            "failed to load kernel from xclbin");
    }
    size_t n_t = n_i;

    int status = 0;

    // Create the compute program from offline
    program = clCreateProgramWithBinary(context, 1, &device_id, &n_t,
            (const unsigned char **) &kernelbinary, &status, &err);

    if ((!program) || (err!=CL_SUCCESS)) {
        throw std::runtime_error(
            "Failed to create compute program from binary");
    }

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to build program executable!");
    }

    // Create the compute kernel in the program we wish to run
    kernel = clCreateKernel(program, kernel_name, &err);

    if (!kernel || err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create compute kernel!");
    }
  }

  ~OpenCLEnv() {
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
  }

  cl_context& getContext() {
    return context;
  }

  cl_command_queue& getCmdQueue() {
    return commands;
  }

  cl_kernel& getKernel() {
    return kernel;
  }

 private:
  int load_file(
      const char *filename, 
      char **result) { 
    int size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) 
    { 
      *result = NULL;
      return -1; // -1 means file opening fail 
    } 
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if (size != fread(*result, sizeof(char), size, f)) 
    { 
      free(*result);
      return -2; // -2 means file reading fail 
    } 
    fclose(f);
    (*result)[size] = 0;
    return size;
  }

  cl_context context;                 // compute context
  cl_command_queue commands;          // compute command queue
  cl_program program;                 // compute program
  cl_kernel kernel;                   // compute kernel
};
#endif
