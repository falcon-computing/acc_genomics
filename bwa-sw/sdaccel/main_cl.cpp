#include <dirent.h> 
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <stdexcept>

#include "bwa/bwa.h"
#ifndef HLS_
#include "BWAOCLEnv.h" // load xilinx opencl framework
#endif

timespec diff(timespec start, timespec end) {   
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

timespec tic( ) {
	timespec start_time;
	clock_gettime(CLOCK_REALTIME, &start_time);
	return start_time;
}

char* toc( timespec* start_time ) {
	timespec current_time;
	//printTimeSpec(diff(*start_time, current_time));
	*start_time = current_time;
}

#ifdef HLS_
void sw_top(int *a, int *b,
    int *output_a, int *output_b,
    int *pac_input_a, int *pac_input_b,
    int size_a, int size_b);

int single_run(char* fname_input, char* fname_golden, int* pac, int pac_size) {
#else
int single_run(BWAOCLEnv *env, char* fname_input, char* fname_golden, int* pac, int pac_size) {
#endif

  static bool first_set = true;
  FILE* f_inp = fopen(fname_input, "rb");
  if (f_inp == NULL) { 
    printf("ERROR: Cannot open the input data file %s\n", fname_input);
    return 1;
  }

  int size = 0;
  fseek(f_inp, 0, SEEK_END);
  size = ftell(f_inp);
  fseek(f_inp, 0, SEEK_SET);
  size = (size+3)/4;

  int *a = (int *) malloc(sizeof(int)*size);
  fread(a, sizeof(int), size, f_inp);
  fclose(f_inp);

  // Deserialize data into bwa-sw kernel input pattern
  int task_num = 0;
  int k = 0;
  while (k < size) {
    int next_idx = a[k++];
    int read_len = a[k++];
    k += (read_len + 7) / 8; // skip read seq

    int chain_num = a[k++];
    for (int chain = 0; chain < chain_num; chain++) {
      int seq_len = a[k+2] - a[k];
      k += 4; // skip rmax

      int seed_num = a[k++];
      task_num += seed_num;
      k += 5*seed_num; // skip seed data
    }
  }

  int* results_a = (int *) malloc(sizeof(int)*task_num*5);
  int* results_b = (int *) malloc(sizeof(int)*task_num*5);
  for(int i = 0; i < task_num*5; i++) {
    results_a[i] = 0;
    results_b[i] = 0;
  }

#ifndef HLS_
  // Setup OpenCL environment
  cl_int           err = 0;
  cl_context       context = env->getContext();
  cl_kernel        kernel  = env->getKernel();
  cl_command_queue command = env->getCmdQueue();
  
  cl_mem_ext_ptr_t ext_a, ext_b;
#ifdef DUAL_BANK
  ext_a.flags = XCL_MEM_DDR_BANK0; ext_a.obj =0; ext_a.param=0;
  ext_b.flags = XCL_MEM_DDR_BANK2; ext_b.obj =0; ext_b.param=0;
#else
  ext_a.flags = XCL_MEM_DDR_BANK1; ext_a.obj =0; ext_a.param=0;
  ext_b.flags = XCL_MEM_DDR_BANK1; ext_b.obj =0; ext_b.param=0;
#endif

  // Create the input and output arrays in device memory for our calculation
  cl_mem input_a = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
      sizeof(int) * size, &ext_a, NULL);
  cl_mem input_b = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
      sizeof(int) * size, &ext_b, NULL);
  cl_mem output_a = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX,
      sizeof(int) * (task_num*5), &ext_a, NULL);
  cl_mem output_b = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX,
      sizeof(int) * (task_num*5), &ext_b, NULL);
  
  if (!input_a || !output_a) {
    fprintf(stderr, "Error: Failed to allocate device memory\n");
    return EXIT_FAILURE;
  }    

  // Write our data set into the input array in device memory 
  timespec timer = tic();
  err = clEnqueueWriteBuffer(command, input_a, CL_TRUE, 0, sizeof(int) * size, a, 0, NULL, NULL);
  err = clEnqueueWriteBuffer(command, input_b, CL_TRUE, 0, sizeof(int) * size, a, 0, NULL, NULL);
  
  if (err != CL_SUCCESS) {
    fprintf(stderr, "Error: Failed to write to source array\n");
    return EXIT_FAILURE;
  }
  timespec e_time = diff(timer, tic());
  printf("input used %d.%06d sec, ", (int)e_time.tv_sec, (int)e_time.tv_nsec/1e3);

  // Set the arguments to our compute kernel
  int i_arg = 0;
  err = 0;
  err  = clSetKernelArg(kernel, i_arg++, sizeof(cl_mem), &input_a);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(cl_mem), &input_b);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(cl_mem), &output_a);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(cl_mem), &output_b);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(cl_mem), &env->pac_input_a_);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(cl_mem), &env->pac_input_b_);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(int), &size);
  err |= clSetKernelArg(kernel, i_arg++, sizeof(int), &size);

  if (err != CL_SUCCESS) {
    fprintf(stderr, "Error: Failed to set kernel arguments %d\n", err);
    return EXIT_FAILURE;
  }
  cl_event taskevent;

  timer = tic();

  err = clEnqueueTask(command, kernel, 0, NULL, &taskevent);
  clWaitForEvents(1, &taskevent);

  if (err) {
    fprintf(stderr, "Error: Failed to execute kernel! %d\n", err);
    return EXIT_FAILURE;
  }

  e_time = diff(timer, tic());
  printf("kernel %d.%06d sec ", (int)e_time.tv_sec, (int)e_time.tv_nsec/1e3);

  err = clEnqueueReadBuffer(command, output_a, CL_TRUE, 0, sizeof(int) * (task_num*5), results_a, 0, NULL, NULL);
  err = clEnqueueReadBuffer(command, output_b, CL_TRUE, 0, sizeof(int) * (task_num*5), results_b, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    fprintf(stderr, "Error: Failed to read output_a array! %d\n", err);
    return EXIT_FAILURE;
  }

  clReleaseMemObject(input_a);
  clReleaseMemObject(input_b);
  clReleaseMemObject(output_a);
  clReleaseMemObject(output_b);
#else
  sw_top(a, a, results_a, results_b, pac, pac, size, size);
#endif
  free(a);

  // Validate results
  if (fname_golden) {
    // assuming the results data is the same dumped from FPGA
    int* results_raw  = (int *)malloc(sizeof(int)*task_num*5);
    int* results_base = (int *)malloc(sizeof(int)*task_num*4);
    FILE* fout = fopen(fname_golden, "rb");
    if (!fout) {
      fprintf(stderr, "cannot find baseline file %s\n", fname_golden);
      return 1;
    }
    fread(results_raw, sizeof(int), 5*task_num, fout);
    fclose(fout);

    // re-order results-base
    for (int i = 0; i < task_num; i++) {
      int seed_idx = results_raw[i*5];
      memcpy(&results_base[seed_idx*4], &results_raw[i*5+1], 4*sizeof(int));
    }
    free(results_raw);

    int num_diffs = 0;
    for (int i = 0; i < task_num; i++) {
      int seed_index_a = results_a[i*5];
      int seed_index_b = results_b[i*5];
      for (int j = 0; j < 4; j++) {
        if (results_base[seed_index_a*4 + j] != results_a[i*5 + j + 1] ||
            results_base[seed_index_b*4 + j] != results_b[i*5 + j + 1]) {

          if (num_diffs < 10) {
            printf("%d != %d\n", 
                results_base[seed_index_b*4 +j],
                results_a[i*5 + j + 1]);
          }
          num_diffs ++;
        }
      }
    }
    free(results_base);
    free(results_a);
    free(results_b);

    if (num_diffs == 0) {
      printf("Pass\n");
      return 0;
    }
    else {
      printf("Failed\n");
      fprintf(stderr, "%d out of %d results are different.\n", num_diffs, 4*task_num);
      return 1;
    }
  }
  else {
    free(results_a);
    free(results_b);
    printf("\n");
    return 0;
  }
}

// NOTE: the FPGA kernel requires a different layout for the pac
// buffer, hence the _set_pac() and _get_pac() macros are different
// from bwa/bntseq.c
// _get_pac_2() is introduced to get the origin pac array from
// aux->idx, reorganizing it for FPGA
#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((l&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((l&3)<<1)&3)
#define _get_pac_2(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
int64_t get_full_pac(char* &pac, bwaidx_t* bwaidx) { 
  int64_t l_pac = bwaidx->bns->l_pac;
  int64_t pac_size = (bwaidx->bns->l_pac * 2 + 3) / 4;

  pac = (char*)calloc(pac_size, 1);

  int64_t k = l_pac;
  
  // forward, using different set/get_pac schemes
  for (int64_t l = 0; l < l_pac; l++) {
    _set_pac(pac, l, _get_pac_2(bwaidx->pac, l));
  }
  // backward
  for (int64_t l = l_pac - 1; l >= 0; --l) {
    _set_pac(pac, k, 3-_get_pac(pac, l));
    k++;
  }

  return pac_size;
}

int main(int argc, char** argv) {
  
  if (argc < 4) {
    fprintf(stderr, "Usage: %s xclbin ref.fa input_data [output_data] \n", argv[0]);
    return 1;
  }
  // load bwa index
  bwaidx_t* bwaidx = bwa_idx_load(argv[2], BWA_IDX_ALL);
  if (bwaidx == 0) {
    printf("failed to load reference from %s\n", argv[2]);
  }

  char* pac;
  int64_t pac_size = get_full_pac(pac, bwaidx);

  try { 
#ifndef HLS_
    BWAOCLEnv env(argv[1], "sw_top", pac_size, pac);
#endif
    int case_count = 0;
    int case_pass = 0;

    // loop through directory
    DIR           *d = opendir(argv[3]);
    struct dirent *dir;
    if (d) {
      while ((dir = readdir(d)) != NULL) {
        if (dir->d_name[0] == '.') {
          continue;
        }
        printf("%s ", dir->d_name);
        char fname_in[4096];
        sprintf(fname_in, "%s/%s", argv[3], dir->d_name);

        int ret = 0;

        if (argc == 5) { // has output file means comparing results
          char fname_out[4096];
          sprintf(fname_out, "%s/%s", argv[4], dir->d_name);
          
#ifdef HLS_
          ret = single_run(fname_in, fname_out, pac, pac_size);
#else
          ret = single_run(&env, fname_in, fname_out, (int*)pac, pac_size);
#endif
        }
        else {
#ifdef HLS_
          ret = single_run(fname_in, NULL, pac, pac_size);
#else
          ret = single_run(&env, fname_in, NULL, (int*)pac, pac_size);
#endif
        }
        case_count ++;
        if (ret == 0) {
          case_pass ++;
        }
      }
      printf("%d out of %d cases passed the test.\n", case_pass, case_count);
    }
    else {
      throw std::runtime_error("input folder is not valid");
    }
  }
  catch (std::runtime_error &e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }

  free(pac);
  bwa_idx_destroy(bwaidx);

  return 0;
}
