#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <time.h>

#include <glog/logging.h>

#include "blaze/altr_opencl/OpenCLEnv.h" 
#include "ksight/tools.h"
#include "PairHMMBuffer.h"
#include "PairHMMHostInterface.h"
#include "PairHMMFpgaInterface.h"
#include "PairHMMTask.h"

PairHMM::PairHMM(): blaze::Task(3),
  ROWS(40), COLS(8),
  max_read_length(0),
  read_segment(0),
  hap_segment(0),
  concat_reads_array(2),
  haps_array(2),
  max_read_lengths(2)
{ ; }

PairHMM::~PairHMM() {
  if (env) {
    // push back the scratch blocks to TaskEnv
    for (auto m : buf_table_) {
      // NOTE: this should execute before member variables' destruction
      env->putScratch(m.first, m.second);
      DVLOG(1) << "Putting " << m.first << " back to scratch";
    }
    env->putScratch("output", output_);
    DVLOG(1) << "Putting output back to scratch";
  }
  DVLOG(1) << "Destroying PairHMM Task";
}

void PairHMM::prepare() {
  PLACE_TIMER;
  
  env = (blaze::OpenCLEnv*)getEnv();
  kernel_table_["pairhmm_driver"]  = env->getKernel(0);
  kernel_table_["pairhmm_output"]  = env->getKernel(1);

  arg_index_["pairhmm_driver"] = 0;
  arg_index_["pairhmm_output"] = 0;

  // initialize ROW/COL from config
  try {
    ROWS = conf_str2int("ROWS");
    COLS = conf_str2int("COLS");
  }
  catch (std::exception &e) {
    DLOG(ERROR) << "Failed to obtain kernel paramter, use default";
  }
  DVLOG(1) << "ROWS = " << ROWS << ", COLS = " << COLS;

  num_cell_ = *((int*)getInput(0));
  num_read_ = deserialize(getInput(1), reads_);
  num_hap_  = deserialize(getInput(2), haps_);

  // get buffers from scratch, create if it does not exists
  create_buffer("read_data", 
                CL_MEM_READ_ONLY,
                MAX_READ_LENGTH*sizeof(ReadData));
  create_buffer("hap_data",
                CL_MEM_READ_ONLY,
                MAX_HAP_LENGTH * sizeof(HapData));
  create_buffer("results",
                CL_MEM_WRITE_ONLY,
                (MAX_READS * MAX_HAPS * sizeof(ResultData)));

  // create output block, cache it in scratch
  if (!env->getScratch("output", output_)) {
    uint64_t output_size = MAX_RSDATA_NUM * MAX_HAPDATA_NUM;
    blaze::DataBlock_ptr b(new blaze::DataBlock(
          1, output_size, output_size*sizeof(float),
          0, blaze::DataBlock::OWNED));

    output_ = b;
    DLOG(INFO) << "Creating a new PairHMM output buffer";
  }
  setOutput(0, output_);

  // start preparing host data
  reset_read_buffers();
  reset_hap_buffers();
  uint readCnt = 0;

  size_t total_read_size = 0;

  ushort read_count = 0;

  for (int i = 0; i < num_read_; i++) {
    if (concat_reads_array[read_segment].size() > MAX_NUM_RESULTS ||
        total_read_size + padded_read_length( reads_[i].len ) >= MAX_READ_LENGTH) {

      max_read_lengths[read_segment] = max_read_length;

      // realloc
      read_segment++;
      if (read_segment >= concat_reads_array.size()) {
        concat_reads_array.resize(concat_reads_array.size() * 2);
        max_read_lengths.resize(max_read_lengths.size() * 2);
      }

      total_read_size = 0;
    }

    concat_reads_array[read_segment].push_back( std::make_pair(reads_[i], i) );
    total_read_size += padded_read_length( reads_[i].len );
    max_read_length = std::max( max_read_length, (uint) reads_[i].len );
    read_count++;
  }

  max_read_lengths[read_segment] = max_read_length;

  // Sort each read segment
  for (uint segment = 0; segment < concat_reads_array.size(); segment++ ) {
    std::sort( 
        concat_reads_array[segment].begin(),
        concat_reads_array[segment].end(),
        [](const Read & a, const Read & b) { return a.first.len < b.first.len; }
        );
  }

  // Add the haps
  ushort hap_count = 0;
  size_t total_hap_size = 0;

  for (int i = 0; i < num_hap_; i++) {

    hap_count++;

    // consider padded hap_len rather than hap_len
    int hap_length = haps_[i].len;
    int pad_hap_data = hap_length;
    if (hap_length <= pad_hap_data_cap) {
      pad_hap_data_cap -= hap_length;
    } 
    else {
      pad_hap_data += pad_hap_data_cap;
      pad_hap_data_cap = (hap_length % ROWS) ? (ROWS - (hap_length % ROWS)) : 0;
    }

    if (total_hap_size + pad_hap_data > MAX_HAP_LENGTH) {
      hap_segment++;

      // realloc
      if (hap_segment >= haps_array.size()) {
        haps_array.resize(haps_array.size() * 2);
      }

      // similar to reset_hap_buffers
      pad_hap_data_cap = 0;
      total_hap_size = 0;
    }

    haps_array[hap_segment].push_back( std::make_pair(haps_[i], i) );
    total_hap_size += pad_hap_data;
  }
}

void PairHMM::compute() {
  PLACE_TIMER;

  for (int cur_read_segment = 0; cur_read_segment <= read_segment; cur_read_segment++ ) {
    reset_read_buffers();
    for( auto read : concat_reads_array[cur_read_segment] ) {
      add_read( read );
    }

    // Compute adjusted value for read length
    uint adj_total_read_data = comp_adj_param( total_read_data, COLS );

    for( uint cur_hap_segment = 0; cur_hap_segment <= hap_segment; cur_hap_segment++ ) {
      reset_hap_buffers();
      for( auto hap : haps_array[cur_hap_segment] ) {
        add_hap( hap );
      }

      uint num_rows = max_read_length + total_hap_data - 1;
      uint num_results = num_reads * num_haps;

      DVLOG(1) << "total_read_data = " << total_read_data;;
      DVLOG(1) << "total_hap_data = " << total_hap_data;;
      DVLOG(1) << "num_rows = " << num_rows;
      DVLOG(1) << "num_results = " << num_results 
               << " num_reads " << num_reads
               << " num_haps " << num_haps;

      { PLACE_TIMER1("fpga_run");
      // Copy data from CPU buffers to FPGA global memory
      write_buffer( "read_data", total_read_data * sizeof(ReadData) );
      write_buffer( "hap_data", total_hap_data * sizeof(HapData) );
      
      // Start pairhmm_output kernel in resultq
      add_buffer_arg( "pairhmm_output", "results" );
      add_scalar_arg( "pairhmm_output", &num_results, sizeof(uint) );

      // Compute adjusted value for number of rows
      uint adj_num_rows = comp_adj_param( num_rows, ROWS );

      // Start pairhmm_driver kernel
      add_buffer_arg( "pairhmm_driver", "read_data" );
      add_buffer_arg( "pairhmm_driver", "hap_data" );
      add_scalar_arg( "pairhmm_driver", &total_read_data, sizeof(uint) );
      add_scalar_arg( "pairhmm_driver", &adj_total_read_data, sizeof(uint) );
      add_scalar_arg( "pairhmm_driver", &total_hap_data, sizeof(uint) );
      add_scalar_arg( "pairhmm_driver", &num_rows, sizeof(uint) );
      add_scalar_arg( "pairhmm_driver", &adj_num_rows, sizeof(uint) );

      run_task("pairhmm_driver" );
      run_task("pairhmm_output", 1);

      // Read results
      read_buffer( "results",
          num_results * sizeof(ResultData),
          1);

      finish(1);
      }

      float *output_ptr = static_cast<float*>(output_->getData());
      ResultData *results_ptr = static_cast<ResultData *>(buf_table_["results"]->get_host_buffer());
      for( int i = 0; i < num_results; i++ ) {
        uint result_index = (results_ptr[i].result_read_num * num_hap_) + results_ptr[i].result_hap_num;
        output_ptr[result_index] = results_ptr[i].result;

        DVLOG(2) << i 
                 << " results " << results_ptr[i].result
                 << " result_num_read " << results_ptr[i].result_read_num
                 << " result_num_hap " << results_ptr[i].result_hap_num
                 << " result_index " << result_index;
      }
    }
  }
}
