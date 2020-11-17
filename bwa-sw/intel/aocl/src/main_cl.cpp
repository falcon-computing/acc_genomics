#include <boost/atomic.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>
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
#include <iostream>
#include <queue>
#include <string>
#include <stdexcept>

#include "acc_lib/CLEnv.h" // load opencl framework

#define AOCL_ALIGNMENT 64

class KernelWorker;

struct KernelTask {
  int      i_size[2];
  int      o_size[2];
  cl_mem   i_buf[2];
  cl_mem   o_buf[2];
  cl_event write_events[2];
  cl_event kernel_events[2];
  KernelWorker* worker;
};

void *alignedMalloc(size_t size) {
  void *result = NULL;
  posix_memalign (&result, AOCL_ALIGNMENT, size);
  return result;
}

inline uint64_t getNs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);

  return (uint64_t)tr.tv_sec*1e9 + tr.tv_nsec;
}

class KernelWorker {
  public:
    KernelWorker(OpenCLEnv* env) {
      context_ = env->getContext();

      cl_int        err       = 0;
      cl_device_id  device_id = env->getDeviceID();
      cl_program    program   = env->getProgram();

      for (int i = 0; i < 4; i++) {
        cmd_[i] = clCreateCommandQueue(context_, device_id, 0, &err);
        if (err != CL_SUCCESS) {
          throw std::runtime_error("Failed to create a command queue context!");
        }
      }
      for (int i = 0; i < 2; i++) {
        char kernel_in_name[100];
        char kernel_out_name[100];

        sprintf(kernel_in_name,"data_parse%d", i);
        sprintf(kernel_out_name,"upload_results%d", i);

        kernels_[2*i+0] = clCreateKernel(program, kernel_in_name, &err);
        kernels_[2*i+1] = clCreateKernel(program, kernel_out_name, &err);
      }
    }

    ~KernelWorker() {
      for (int i=0; i<4; i++) {
        clReleaseCommandQueue(cmd_[i]);
        clReleaseKernel(kernels_[i]);
      }
    }

    void write_buffer(cl_mem buf, void* data, size_t size, int i) {
      uint64_t start_ts = getNs();
      cl_int err = clEnqueueWriteBuffer(cmd_[2*i], buf, CL_FALSE, 0, size, 
          data, 0, NULL, NULL);
      if (err != CL_SUCCESS) {
        throw std::runtime_error("failed to write data");
      }
      std::cerr << "write_buffer-" << i << ": " << getNs() - start_ts << " ns\n";
    }

    void read_buffer(cl_mem buf, void* data, size_t size, int i) {
      cl_int err = clEnqueueReadBuffer(cmd_[2*i+1], buf, CL_TRUE, 0, size, 
          data, 0, NULL, NULL);
      if (err != CL_SUCCESS) {
        throw std::runtime_error("failed to read data");
      }
    }

    void post_task(KernelTask* task, KernelWorker* prev_worker) {
      // kernel execution
      for (int k = 0; k < 2; k++) {
        cl_int err = 0;
        err  = clSetKernelArg(kernels_[2*k+0], 0, sizeof(cl_mem), &task->i_buf[k]);
        err |= clSetKernelArg(kernels_[2*k+0], 1, sizeof(int), &task->i_size[k]);
        err |= clSetKernelArg(kernels_[2*k+1], 0, sizeof(cl_mem), &task->o_buf[k]);
        err |= clSetKernelArg(kernels_[2*k+1], 1, sizeof(int), &task->o_size[k]);

        if (err) {
          std::cerr << "failed to set kernel args" << std::endl;
        }
      }

      uint64_t start_ts = getNs();
      for (int k = 0; k < 2; k++) {
        cl_int err = 0;
        if (!prev_worker) {
          err  = clEnqueueTask(cmd_[2*k+0], kernels_[2*k+0], 0, NULL, &events_[2*k+0]);
          err |= clEnqueueTask(cmd_[2*k+1], kernels_[2*k+1], 0, NULL, &events_[2*k+1]);
        }
        else {
          err  = clEnqueueTask(cmd_[2*k+0], kernels_[2*k+0], 4, prev_worker->events_, &events_[2*k+0]);
          err |= clEnqueueTask(cmd_[2*k+1], kernels_[2*k+1], 4, prev_worker->events_, &events_[2*k+1]);
        }

        if (err) {
          std::cerr << "failed to execute kernels: " << err << std::endl;
        }
      }
    }

    void finish_task() {
      for (int k = 0; k < 4; k++) {
        clReleaseEvent(events_[k]);
      }
    }

  private:
    cl_context       context_;
    cl_kernel        kernels_[4];
    cl_command_queue cmd_[4];
    cl_event         events_[4];
};

int run_kernel(OpenCLEnv* env,
    int max_iter, 
    int i_size, int o_size,
    int* input,
    int** results) 
{

  cl_int       err = 0;
  cl_context   context = env->getContext();
  cl_device_id device_id = env->getDeviceID();

  cl_command_queue cmd = clCreateCommandQueue(context, device_id, 0, &err);
  if (err != CL_SUCCESS) {
    std::cerr << "failed to create a command queue context" << std::endl;
  }

  // for ping-pong buffering
  std::queue<KernelTask*>   task_queue; 
  std::queue<KernelWorker*> worker_queue;

  for (int k = 0; k < 2; k++) {
    KernelWorker* worker = new KernelWorker(env);
    worker_queue.push(worker);
  }

  for (int iter = 0; iter < max_iter + 1; iter++) {
    if (iter < max_iter) {

      KernelWorker* worker = worker_queue.front();
      KernelTask* task = new KernelTask;
      task->worker = worker;

      // switch kernel to the back of the queue
      worker_queue.push(worker);
      worker_queue.pop();

      for (int k = 0; k < 2; k++) {
        task->i_size[k] = i_size; 
        task->o_size[k] = o_size;
        task->i_buf[k] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int)*i_size, NULL, NULL);
        task->o_buf[k] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(int)*o_size, NULL, NULL);

        if (!task->i_buf[k] || !task->o_buf[k]) {
          fprintf(stderr, "Failed to allocate device memory\n");
          return 1;
        }
      }

      uint64_t start_ts = getNs();
      for (int k = 0; k < 2; k++) {
        worker->write_buffer(task->i_buf[k], input, i_size*sizeof(int), k);
      }
      std::cerr << "write buf: " << getNs()-start_ts <<  " ns\n";

      if (task_queue.empty()) {
        worker->post_task(task, NULL);
      }
      else {
        // get prev worker
        worker->post_task(task, task_queue.front()->worker);
      }
      task_queue.push(task);
    }

    if (iter > 0) {
      KernelTask*   task = task_queue.front();
      KernelWorker* worker = task->worker;
      task_queue.pop();

      uint64_t start_ts = getNs();

      // copy data
      for (int k = 0; k < 2; k++) {
        worker->read_buffer(task->o_buf[k], results[k], o_size*sizeof(int), k);
      }
      worker->finish_task();
      std::cerr << "read buf: " << getNs()-start_ts << " ns\n";

      // release events
      for (int k = 0; k < 2; k++) {
        clReleaseMemObject(task->i_buf[k]);
        clReleaseMemObject(task->o_buf[k]);
      }

      delete task;
    }
  }
}

int single_run(OpenCLEnv *env, char* fname_input, char* fname_golden) {

  FILE* f_inp = fopen(fname_input, "rb");
  if (f_inp == NULL) { 
    printf("ERROR: Cannot open the input data file %s\n", fname_input);
    return 1;
  }

  int size = 0;
  int out_size = 0;
  fseek(f_inp, 0, SEEK_END);
  size = ftell(f_inp);
  fseek(f_inp, 0, SEEK_SET);
  size = (size+3)/4;

  int *a = NULL;
  // allocate size that is a multiple of AOCL_ALIGNMENT
  int a_size = sizeof(int)*((size/AOCL_ALIGNMENT) + 1)*AOCL_ALIGNMENT;
  a = (int *)alignedMalloc(a_size); 
  
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
      k += (seq_len + 7) / 8; // skip chain seq

      int seed_num = a[k++];
      task_num += seed_num;
      k += 5*seed_num; // skip seed data
    }
  }
  //printf("task num: %d\n", task_num);
  out_size = task_num*5;

  int* results[2];
  for (int i = 0; i < 2; i++) {
    results[i] = (int*)alignedMalloc(sizeof(int)*out_size);
    for (int j = 0; j < out_size; j++) {
      results[i][j] = 0;
    }
  }

  uint64_t start_ts = getNs();
  run_kernel(env, 1, size, out_size, a, results);
  fprintf(stderr, "total: %d ns\n", getNs()-start_ts);

  free(a);
    
  // Validate results
  if (fname_golden) {
    int* results_base = (int *)malloc(sizeof(int)*task_num*4);
    FILE* fout = fopen(fname_golden, "rb");
    FILE* fout_txt = fopen("out.txt", "w");
    FILE* fres = fopen("res.txt", "w");
    if (!fout) {
      fprintf(stderr, "cannot find baseline file %s\n", fname_golden);
      return 1;
    }
    fread(results_base, sizeof(int), 4*task_num, fout);
    fclose(fout);

    int num_diffs = 0;
    int seed_index[K_NUM];
    for (int i=0; i<K_NUM; i++) {
      for (int j=0; j<task_num; j++) {
        seed_index[i] = results[i][j*5];
        for (int k=0; k<4; k++) {
          fprintf(fout_txt,"%d\n",results_base[seed_index[i]*4 + k]);
          //fprintf(fres,"%d\n",results[i][j*5 + k + 1]);
          if (results_base[seed_index[i]*4 + k] != results[i][j*5 + k + 1]) {
            num_diffs ++;
          }
        }
      }
    }

    for (int i=0; i<K_NUM; i++) {
      for (int j=0; j<out_size; j++) {
         fprintf(fres,"%d\n",results[i][j]);
      }
    }
    fclose(fout_txt);
    fclose(fres);
    free(results_base);
    for (int i=0; i<K_NUM; i++) {
      free(results[i]);
    }

    if (num_diffs == 0) {
      printf("Pass\n");
      return 0;
    }
    else {
      printf("Failed\n");
      fprintf(stderr, "%d out of %d results are different.\n", num_diffs, K_NUM*4*task_num);
      return 1;
    }
  }
  else {
    for (int i=0; i<K_NUM; i++) {
      free(results[i]);
    }
    printf("\n");
    return 0;
  }
}

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s xclbin/aocx input_data [output_data] \n", argv[0]);
    return 1;
  }

  try { 
#ifndef HLS_
    OpenCLEnv env(argv[1], "sw_top");
#endif

    int case_count = 0;
    int case_pass = 0;
#ifdef SINGLE_TEST
    int ret = 0;
    char fname_in[4096] = "testdata/input/sample_3.dat";
    char fname_out[4096] = "testdata/golden/sample_3.dat";
    ret = single_run(&env, fname_in, fname_out);
    if (ret == 0) 
        printf("Passed!\n");
    else
        printf("Failed!\n");
#else
    // loop through directory
    DIR           *d = opendir(argv[2]);
    struct dirent *dir;
    if (d) {
      while ((dir = readdir(d)) != NULL) {
        if (dir->d_name[0] == '.') {
          continue;
        }
        printf("%s ", dir->d_name);
        char fname_in[4096];
        sprintf(fname_in, "%s/%s", argv[2], dir->d_name);

        int ret = 0;

        if (argc == 4) { // has output file means comparing results
          char fname_out[4096];
          sprintf(fname_out, "%s/%s", argv[3], dir->d_name);
#ifdef HLS_
          ret = single_run(fname_in, fname_out);
#else
          ret = single_run(&env, fname_in, fname_out);
#endif
        }
        else {
#ifdef HLS_
          ret = single_run(fname_in, NULL);
#else
          ret = single_run(&env, fname_in, NULL);
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
#endif
  }
  catch (std::runtime_error &e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }

  return 0;
}
