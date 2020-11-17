#include <cstdint>
#include <stdexcept>

#include "ksight/tools.h"
#include "PairHMMFpga.h"
#include "PairHMMFpgaInterface.h"

OpenCLEnv *env = NULL;
float *ret_buf = NULL;
FpgaInputBundle* bundle = NULL;

double peak_kernel_gcups = 0;
double curr_kernel_gcups = 0;

static Context<float>  g_ctxf;
static Context<double> g_ctxd;

#ifdef XILINX
int compute(OpenCLEnv* env,
    int num_pu,
    int num_read,
    int num_hap,
    uint64_t num_cell,
    FpgaInputBundle *bundle) 
{
  PLACE_TIMER;
  cl_int err = 0;

  static unsigned bankID[4] = {
          XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1,
          XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};

  cl_context       context = env->getContext();
  cl_kernel        kernel  = env->getKernel();
  cl_command_queue command = env->getCmdQueue();

  cl_mem_ext_ptr_t input_ext;
  cl_mem_ext_ptr_t output_ext;

  input_ext.flags  = bankID[KERNEL_BANKID];
  input_ext.obj    = bundle;
  input_ext.param  = 0;

  output_ext.flags = bankID[KERNEL_BANKID];
  output_ext.obj   = ret_buf;
  output_ext.param = 0; 

  cl_mem input_cl_buffer = clCreateBuffer(context, 
		  CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
		  sizeof(FpgaInputBundle), &input_ext, &err);

  if (!input_cl_buffer || err != CL_SUCCESS)
    throw std::runtime_error("failed to allocate CL buffer for input");

  cl_mem output_cl_buffer = clCreateBuffer(context, 
		  CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
		  sizeof(float) * MAX_RSDATA_NUM * MAX_HAPDATA_NUM, &output_ext, &err);

  if (!output_cl_buffer || err != CL_SUCCESS)
	  throw std::runtime_error("failed to allocate CL buffer for output");

  // start tranfering the data
  cl_event event_write;
  cl_event event_kernel;
  cl_event event_read;

  {
    PLACE_TIMER1("write_buffer");
    err = clEnqueueMigrateMemObjects(command, 1, &input_cl_buffer, 0, 0, NULL, &event_write);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to migrate input buffer");
    }
    clWaitForEvents(1, &event_write);
  }

  clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_cl_buffer);
  clSetKernelArg(kernel, 1, sizeof(int), &num_read);
  clSetKernelArg(kernel, 2, sizeof(int), &num_hap);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), &output_cl_buffer);

  {
    PLACE_TIMER1("kernel");
    err = clEnqueueTask(command, kernel, 0, NULL, &event_kernel);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to enqueue kernel");
    }
    clWaitForEvents(1, &event_kernel);
  }

  cl_ulong k_start, k_end;
  clGetEventProfilingInfo(event_kernel, CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
  clGetEventProfilingInfo(event_kernel, CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);

  double gcups = (double)num_cell / (k_end - k_start);
  curr_kernel_gcups = gcups;
  peak_kernel_gcups = std::max(gcups, peak_kernel_gcups);

#if 0
  printf("Kernel timer from profiler: %d ns\n", k_end - k_start);
  printf("Current Kernel GCUPS: %.3f #read %d #hap %d #cell %d\n", 
      gcups,
      num_read, num_hap, num_cell);
#endif

  {
    PLACE_TIMER1("read_buffer");
    err = clEnqueueMigrateMemObjects(command, 1, &output_cl_buffer, CL_MIGRATE_MEM_OBJECT_HOST,
        0, NULL, &event_read);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to enqueue read output");
    }

    clWaitForEvents(1, &event_read);
  }

  clReleaseEvent(event_write);
  clReleaseEvent(event_kernel);
  clReleaseEvent(event_read);
  
  clReleaseMemObject(input_cl_buffer);
  clReleaseMemObject(output_cl_buffer);
  return 0;
}
    
float* compute_fpga(
    const char* bit_path,
    std::string read_data,
    std::string hap_data,
    uint64_t num_cell) 
{
  if (!env) {
    env = new OpenCLEnv(bit_path, KERNEL_NAME);
  }
  int num_pu = KERNEL_NUM_PE / (READ_BLOCK_SIZE * HAP_BLOCK_SIZE);

  if (!bundle) {
    bundle = (FpgaInputBundle*)aligned_alloc(4096, sizeof(FpgaInputBundle));
  }
  memset(bundle, 0, sizeof(bundle));

  // deserialize data
  read_t* reads = NULL;
  hap_t*  haps  = NULL;
  int num_read = deserialize(read_data, reads);
  int num_hap  = deserialize(hap_data, haps);

  // prepare input
  pack_fpga_input(num_pu, num_read, num_hap, reads, haps, bundle);

  free_reads(reads, num_read);
  free_haps(haps, num_hap);

  if (!ret_buf) {
    // allocate output buffer
    ret_buf = (float*)aligned_alloc(4096, 
        sizeof(float) * MAX_RSDATA_NUM * MAX_HAPDATA_NUM);
  }

  compute(env, num_pu, num_read, num_hap, num_cell, bundle);

  return ret_buf;
}
#else
int compute(OpenCLEnv* env,
    int num_read,
    int num_hap,
    int num_results,
    int num_rows,
    int seg_num_cell,
    cl_mem   &read_cl_buffer,
    ReadData *read_host_buf,
    int       total_read_data,
    cl_mem   &hap_cl_buffer,
    HapData  *hap_host_buf,
    int       total_hap_data,
    cl_mem   &output_cl_buffer,
    ResultData *results) 
{
  PLACE_TIMER;
  cl_int err = 0;

  cl_context       context  = env->getContext();
  cl_kernel        kernel0  = env->getKernel(0);
  cl_kernel        kernel1  = env->getKernel(1);
  cl_command_queue command0 = env->getCmdQueue(0);
  cl_command_queue command1 = env->getCmdQueue(1);

  cl_event event_kernel0;
  cl_event event_kernel1;

  // start tranfering the data
  {
    PLACE_TIMER1("write_buffer");
    err = clEnqueueWriteBuffer( command0, read_cl_buffer, CL_FALSE, 0,
                                total_read_data*sizeof(ReadData), read_host_buf, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to write input buffer for read");
    }
    err = clEnqueueWriteBuffer( command0, hap_cl_buffer, CL_FALSE, 0,
                                total_hap_data*sizeof(HapData), hap_host_buf, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to write input buffer for read");
    }
    //clFinish( command0 );
  }

  clSetKernelArg(kernel1, 0, sizeof(cl_mem), &output_cl_buffer);
  clSetKernelArg(kernel1, 1, sizeof(cl_int), &num_results);

  clSetKernelArg(kernel0, 0, sizeof(cl_mem), &read_cl_buffer);
  clSetKernelArg(kernel0, 1, sizeof(cl_mem), &hap_cl_buffer);
  clSetKernelArg(kernel0, 2, sizeof(cl_int), &total_read_data);
  clSetKernelArg(kernel0, 3, sizeof(cl_int), &total_hap_data);
  clSetKernelArg(kernel0, 4, sizeof(cl_int), &num_rows);

  //uint64_t start_ts = ::getUs();
  {
    PLACE_TIMER1("kernel");
    err  = clEnqueueTask(command0, kernel0, 0, NULL, NULL);
    err |= clEnqueueTask(command1, kernel1, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to enqueue kernel");
    }
    //clWaitForEvents(1, &event_kernel0);
    //clWaitForEvents(1, &event_kernel1);
    //clFinish( command1 );
  }
  //uint64_t end_ts = ::getUs();


#if 0
  double gcups = (double)seg_num_cell / (double)(end_ts - start_ts) / 1.e3 ;
  curr_kernel_gcups = gcups;
  peak_kernel_gcups = std::max(gcups, peak_kernel_gcups);
  printf("Kernel timer from profiler: %d ns\n", (end_ts - start_ts) * 1000 );
  printf("Current Kernel GCUPS: %.3f #read %d #hap %d #cell %d\n", 
      gcups,
      num_read, num_hap, seg_num_cell);
#endif

  {
    PLACE_TIMER1("read_buffer");
    err = clEnqueueReadBuffer( command1, output_cl_buffer, CL_TRUE, 0,
                               num_results * sizeof(ResultData), results, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to enqueue read output");
    }
  }
  //clFinish( command1 );

  return 0;
}

float* compute_fpga(
    const char* bit_path,
    std::string read_data,
    std::string hap_data,
    uint64_t num_cell) 
{
  if (!env) {
    std::vector<std::string> kernel_names_str;
    kernel_names_str.push_back("pairhmm_driver");
    kernel_names_str.push_back("pairhmm_output");

    std::vector<const char *> kernel_names;
    kernel_names.push_back(kernel_names_str[0].c_str());
    kernel_names.push_back(kernel_names_str[1].c_str());

    env = new OpenCLEnv(bit_path, kernel_names);
  }


  if (!ret_buf) {
    // allocate output buffer
    ret_buf = (float*)aligned_alloc(4096, 
        sizeof(float) * MAX_RSDATA_NUM * MAX_HAPDATA_NUM);
  }

  std::vector<ReadINTL> reads;
  std::vector<HaplotypeINTL> haps;
  int total_num_reads = deserialize(read_data, reads);
  int total_num_haps  = deserialize(hap_data, haps);

  int read_segment = 0;
  int hap_segment = 0;
  std::vector<std::vector<ReadINTL>*> concat_reads_array;
  std::vector<std::vector<HaplotypeINTL>*> haps_array;
  std::vector<int> max_read_lengths;
  int read_count = 0;
  int total_read_size = 0;
  int max_read_length = 0;

  int hap_count = 0;
  int total_hap_size = 0;
 
  concat_reads_array.push_back(new std::vector<ReadINTL>());
  for(int i = 0; i < total_num_reads; i++) {
      if( concat_reads_array[read_segment]->size() > MAX_NUM_RESULTS ||
          total_read_size + padded_length( reads[i].length ) >= MAX_READ_LENGTH ) {
          concat_reads_array.push_back(new std::vector<ReadINTL>());
          max_read_lengths.push_back( max_read_length );
          read_segment++;
          total_read_size = 0;
      }

      concat_reads_array[read_segment]->push_back( reads[i] );
      total_read_size += padded_length( reads[i].length );
      max_read_length = std::max( max_read_length, (int)reads[i].length );
      read_count++;
  }
  max_read_lengths.push_back( max_read_length );

  // Sort each read segment
  for( int segment = 0; segment <= read_segment; segment++ ) {
    std::sort( concat_reads_array[segment]->begin(),
               concat_reads_array[segment]->end() );
  }

  // Add the haps
  haps_array.push_back(new std::vector<HaplotypeINTL>());
  for(int i = 0; i < total_num_haps; i++) {
      haps[i].y_init = g_ctxf.INITIAL_CONSTANT / (float) haps[i].length;
      hap_count++;
      if( total_hap_size + (int) haps[i].length > MAX_HAP_LENGTH ) {
          haps_array.push_back(new std::vector<HaplotypeINTL>());
          hap_segment++;
          total_hap_size = 0;
      }

      haps_array[hap_segment]->push_back( haps[i] );
      total_hap_size += (int)haps[i].length;
  }

#if 0
  if(!useFPGA(num_cell, (read_segment + 1), (hap_segment + 1))){
      //INFO("use AVX this batch, num reads = %d, num haps = %d, read_segment = %d, hap_segment = %d", reads.size(), haps.size(), read_segment + 1, hap_segment + 1);
      //run_batch_avx(batch);
      return NULL;
  }
  INFO("use FPGA this batch, num reads = %d, num haps = %d, read_segment = %d, hap_segment = %d", reads.size(), haps.size(), read_segment + 1, hap_segment + 1);
#endif

  ReadData *read_host_buf = (ReadData *)aligned_alloc(4096, MAX_READ_LENGTH * sizeof(ReadData));
  HapData  *hap_host_buf  = (HapData  *)aligned_alloc(4096, MAX_HAP_LENGTH * sizeof(HapData));
  ResultData *results     = (ResultData *)aligned_alloc(4096, MAX_READS * MAX_HAPS * sizeof(cl_float) );

  cl_int err;
  cl_context context = env->getContext();

  cl_mem read_cl_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, MAX_READ_LENGTH * sizeof(ReadData), NULL, &err);
  if (!read_cl_buffer || err != CL_SUCCESS)
    throw std::runtime_error("failed to allocate CL buffer for read");

  cl_mem hap_cl_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, MAX_HAP_LENGTH * sizeof(HapData), NULL, &err);
  if (!hap_cl_buffer || err != CL_SUCCESS)
    throw std::runtime_error("failed to allocate CL buffer for hap");

  cl_mem output_cl_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, MAX_READS * MAX_HAPS * sizeof(cl_float), NULL, &err);
  if (!output_cl_buffer || err != CL_SUCCESS)
    throw std::runtime_error("failed to allocate CL buffer for output");

  for( int cur_read_segment = 0; cur_read_segment <= read_segment; cur_read_segment++ ) {
      ReadData *read_host_data = read_host_buf;
      int       seg_read_size = 0; 
      int       num_reads = concat_reads_array[cur_read_segment]->size();

      for( auto read : *(concat_reads_array[cur_read_segment]) ) {
          add_read(read_host_data, read, seg_read_size);
      }

      for( int cur_hap_segment = 0; cur_hap_segment <= hap_segment; cur_hap_segment++ ) {
          HapData *hap_host_data = hap_host_buf;
          int       seg_hap_size = 0;
          int       num_haps = haps_array[cur_hap_segment]->size();

          for( auto hap : *(haps_array[cur_hap_segment]) ) {
              add_hap(hap_host_data, hap, seg_hap_size);
          }

          int num_rows = max_read_length + seg_hap_size - 1;
          int num_results = num_reads * num_haps;
          int seg_num_cell = seg_read_size * seg_hap_size;

#if 0
          DBG( "num_results = %d num_reads %d num_haps %d",
               num_results,
               num_reads,
               num_haps );
#endif

          //compute(env, num_pu, num_read, num_hap, num_cell, bundle);
          compute(env, num_reads, num_haps, num_results, num_rows, seg_num_cell,
                  read_cl_buffer, read_host_buf, seg_read_size,
                  hap_cl_buffer, hap_host_buf, seg_hap_size,
                  output_cl_buffer, results);

          for( int i = 0; i < num_results; i++ ) {
              int result_index = results[i].result_read_num * total_num_haps + results[i].result_hap_num;
              ret_buf[result_index] = results[i].result;
          }
      } // Hap segment
  } // Read segment

  clReleaseMemObject(read_cl_buffer);
  clReleaseMemObject(hap_cl_buffer);
  clReleaseMemObject(output_cl_buffer);

  return ret_buf;
}

#endif

void cleanup() {
  free(bundle);
  free(ret_buf);
  free(env);
}
