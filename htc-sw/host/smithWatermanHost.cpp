/*
 * Copyright (C) 2015-2018 Falcon Computing Solutions, Inc. All Rights Reserved.
 * Description: OpenCL function definitions for Xilinx FPGA implementation 
 * Author: Jiayi Sheng
 * Date: Apr 30th 2018
 * Version: 1.0
*/

#include <string.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <assert.h>
#include <string.h>
#include "smithWatermanHost.h"

#define KERNEL_NUM 1

static int count_init=0;
cl_platform_id platform_id;         // platform id
cl_device_id device_id;             // compute device id 
cl_context context;                 // compute context
cl_command_queue commands[KERNEL_NUM];          // compute command queue
cl_program program;                 // compute program

cl_mem _inputs_buffer;
cl_event __event__inputs_buffer;
cl_mem _outputs_buffer;
cl_event __event__outputs_buffer;
cl_kernel __smithWaterman_kernel;
cl_event __event_smithWaterman_kernel;

int _create_context() {
    int err;
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    if (!context) {
      printf("Error: Failed to create a compute context!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    return 0;
}

int _get_device_id() {
    int err;
    int fpga = 1;
    err = clGetDeviceIDs(platform_id, fpga ? CL_DEVICE_TYPE_ACCELERATOR : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    if (err != CL_SUCCESS) {
      printf("Error: Failed to create a device group!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    return 0;
}

int _get_platform_info() {
    char cl_platform_vendor[1001];
    char cl_platform_name[1001];
    cl_platform_vendor[0] = 0;
    cl_platform_name[0] = 0;
    int err;
    
    err = clGetPlatformIDs(1,&platform_id,NULL);
    if (err != CL_SUCCESS) {
      printf("Error: Failed to find an OpenCL platform!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    err = clGetPlatformInfo(platform_id,CL_PLATFORM_VENDOR,1000,(void *)cl_platform_vendor,NULL);
    if (err != CL_SUCCESS) {
      printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    printf("CL_PLATFORM_VENDOR %s\n",cl_platform_vendor);
    err = clGetPlatformInfo(platform_id,CL_PLATFORM_NAME,1000,(void *)cl_platform_name,NULL);
    if (err != CL_SUCCESS) {
      printf("Error: clGetPlatformInfo(CL_PLATFORM_NAME) failed!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    printf("CL_PLATFORM_NAME %s\n",cl_platform_name);
    return 0;
}

int _create_command_queue() {
    int err;
    int i;
    for(i=0; i < 1; i++) {
        commands[i] = clCreateCommandQueue(context, device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
        if (!commands[i]) {
          printf("Error: Failed to create a command queue commands[%d]!\n",i);
          printf("Error: code %i\n",err);
          exit(EXIT_FAILURE);
        }
    }
    return 0;
}

int load_file_to_memory(const char *filename, char **result) { 
  size_t size = 0;
  FILE *f = fopen(filename, "rb");
  if (f == NULL) { 
    printf("ERROR : Kernel binary %s not exist!\n", filename);
    *result = NULL;
    return -1; // -1 means file opening fail 
  } 
  fseek(f, 0, SEEK_END);
  size = ftell(f);
  fseek(f, 0, SEEK_SET);
  *result = (char *)malloc(size+1);
  if (size != (size_t)fread(*result, sizeof(char), size, f)) { 
    free(*result);
    return -2; // -2 means file reading fail 
  } 
  fclose(f);
  (*result)[size] = 0;
  return size;
}

int _create_program(const char * bitstream) {
    int err;
    int n_i = 0;
    unsigned char *kernelbinary;
    char * bit_file;
    if(bitstream == NULL) {
        bit_file=(char *)"kernel_top.xclbin";
    } else { 
        bit_file=(char *)bitstream; 
    }
    printf("loading %s\n", bit_file);
    n_i = load_file_to_memory(bit_file, (char **) &kernelbinary);
    if (n_i < 0) {
      printf("ERROR : failed to load kernel from binary: %s\n", bit_file);
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    
    int status;
    size_t n = n_i;
    program = clCreateProgramWithBinary(context, 1, &device_id, &n,
                                        (const unsigned char **) &kernelbinary, &status, &err);
    if ((!program) || (err!=CL_SUCCESS)) {
      printf("Error: Failed to create compute program from binary %d!\n", err);
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
      size_t len;
      char buffer[2048];

      printf("Error: Failed to build program executable!\n");
      clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
      printf("%s\n", buffer);
      exit(EXIT_FAILURE);
    }
    return 0;
}

int _init_opencl(const char * bitstream) {
    if(count_init == 0) {
        _get_platform_info();
        _get_device_id();
        _create_context();
        _create_command_queue();
        _create_program(bitstream);
        count_init++;
        return 1;
    } else {
        return 0;
    }
}

int opencl_create_kernel(cl_kernel * kernel, char * kernel_name) {
    int err;
    *kernel = clCreateKernel(program, kernel_name, &err);
    if (!(*kernel) || err != CL_SUCCESS) {
      printf("Error: Failed to create compute kernel!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    return CL_SUCCESS;
}

int opencl_create_buffer(cl_mem * cl_buffer, long long size, int type) {
    if (size >= (1LL<<32)) {
        printf("ERROR: Invalid Size: size=%lld\n", size);
        exit(EXIT_FAILURE);
    }
    size_t flag = 0;
    if (type == 0) flag = CL_MEM_READ_ONLY;
    if (type == 1) flag = CL_MEM_WRITE_ONLY;
    if (type == 2) flag = CL_MEM_READ_WRITE;
    *cl_buffer = clCreateBuffer(context,  flag, size, NULL, NULL);
    if (!(*cl_buffer)) {
      printf("Error: Failed to allocate device memory!\n");
      exit(EXIT_FAILURE);
    }    
    return CL_SUCCESS;
}

int _init_kernel_buffer() {
    opencl_create_kernel(&__smithWaterman_kernel, (char *)"smithWatermanMerlin");
    opencl_create_buffer(&_inputs_buffer, 1*134152, 2);
    opencl_create_buffer(&_outputs_buffer, 2*266764, 2);
    return 0;
}

int opencl_write_buffer(cl_mem cl_buffer, cl_command_queue commands, long long offset, void * host_buffer, long long size) {
    if (offset < 0 || offset + size - 1 >= (1LL<<32)) {
        printf("ERROR: offset or size overflow: offset=%lld size=%lld\n", offset, size);
        exit(EXIT_FAILURE);
    }
    cl_event event;
    int err = clEnqueueWriteBuffer(commands, cl_buffer, CL_FALSE, offset, size, host_buffer, 0, NULL, &event);
    clWaitForEvents(1, &event);
    clReleaseEvent(event);
    if (err != CL_SUCCESS) {
      printf("Error: Failed to write to device buffer!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    return err;
}

int opencl_read_buffer(cl_mem cl_buffer, cl_command_queue commands, long long offset, void * host_buffer, long long size) {
    if (offset < 0 || offset + size - 1 >= (1LL<<32)) {
        printf("ERROR: offset or size overflow: offset=%lld size=%lld\n", offset, size);
        exit(EXIT_FAILURE);
    }
    cl_event readevent;
    int err = clEnqueueReadBuffer(commands, cl_buffer, CL_FALSE, offset, size, host_buffer, 0, NULL, &readevent);
    clWaitForEvents(1, &readevent);
    clReleaseEvent(readevent);
    if (err != CL_SUCCESS) {
      printf("Error: Failed to read from device buffer!\n");
      printf("Error: code %i\n",err);
      exit(EXIT_FAILURE);
    }
    return err;
}

int opencl_set_kernel_arg(cl_kernel kernel, int index, size_t size, const void * content) {
    return clSetKernelArg(kernel, index, size, content);
}

int opencl_enqueue_kernel(cl_kernel kernel, cl_command_queue commands, int dim, size_t global_in[], cl_event *event_out) {
    int i;
    size_t global[100], local[100];
    assert(dim < 100);
    for (i = 0; i < dim; i++) {
        global[i] = global_in[i];
        local[i] = 1;
    }
//    cl_event event;
//    int err = clEnqueueNDRangeKernel(commands, kernel, dim, NULL, (size_t*)&global, (size_t*)&local, 0, NULL, &event);
//    clWaitForEvents(1, &event);
    int err = clEnqueueNDRangeKernel(commands, kernel, dim, NULL, (size_t*)&global, (size_t*)&local, 0, NULL, event_out);
    return err;
}

void _smithWatermanRun(char *inputs,int refLength,int batchSize,int overhang_strategy,int w_match,int w_mismatch,int w_open,int w_extend,short *outputs)
{
    if (inputs != 0) {
        opencl_write_buffer(_inputs_buffer,commands[0],0 * sizeof(char ),((char *)inputs),((unsigned long )134152) * sizeof(char ));
    }
    else {
        printf("Error! The external memory pointed by 'inputs' is invalid.\n");
        exit(1);
    }
  
	opencl_set_kernel_arg(__smithWaterman_kernel, 0, sizeof(cl_mem), &_inputs_buffer);
	opencl_set_kernel_arg(__smithWaterman_kernel, 1, sizeof(int), &refLength);
	opencl_set_kernel_arg(__smithWaterman_kernel, 2, sizeof(int), &batchSize);
	opencl_set_kernel_arg(__smithWaterman_kernel, 3, sizeof(int), &overhang_strategy);
	opencl_set_kernel_arg(__smithWaterman_kernel, 4, sizeof(int), &w_match);
	opencl_set_kernel_arg(__smithWaterman_kernel, 5, sizeof(int), &w_mismatch);
	opencl_set_kernel_arg(__smithWaterman_kernel, 6, sizeof(int), &w_open);
	opencl_set_kernel_arg(__smithWaterman_kernel, 7, sizeof(int), &w_extend);
	opencl_set_kernel_arg(__smithWaterman_kernel, 8, sizeof(cl_mem), &_outputs_buffer);
	size_t __gid[1];
	__gid[0] = 1;
   
	opencl_enqueue_kernel(__smithWaterman_kernel, commands[0], 1, __gid, &__event_smithWaterman_kernel);
	clWaitForEvents(1, &__event_smithWaterman_kernel);
   
      
            
    if (outputs != 0) {
        int cigarInfoLength = 0;
        opencl_read_buffer(_outputs_buffer,commands[0],0 * sizeof(short ),((short *)outputs),((unsigned long )32) * sizeof(short ));
        cigarInfoLength = outputs[1] << 16 | outputs[0];

        if(batchSize + cigarInfoLength > 30){
            opencl_read_buffer(_outputs_buffer,commands[0],32 * sizeof(short ),((short *)outputs + 32),((unsigned long )(batchSize + cigarInfoLength - 30)) * sizeof(short ));
        }
    }
    else {
        printf("Error! The external memory pointed by 'outputs' is invalid.\n");
        exit(1);
    }
}

int _release_smithWaterman(){
	if(__event_smithWaterman_kernel) {
		clReleaseEvent(__event_smithWaterman_kernel);
	}
	if(__smithWaterman_kernel) {
		clReleaseKernel(__smithWaterman_kernel);
	}
	if(_inputs_buffer) {
		clReleaseMemObject(_inputs_buffer);
	}
	if(_outputs_buffer) {
		clReleaseMemObject(_outputs_buffer);
	}
    return 0;
}


