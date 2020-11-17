#include <glog/logging.h>
#include <unordered_set>
#include <vector>

#include "blaze/Client.h"
#include "gkl-pairhmm/avx_impl.h"
#include "gkl-pairhmm/Context.h"
#ifndef NO_PROFILE
#include "ksight/tools.h"
#endif
#include "PairHMMClient.h"
#include "PairHMMFpgaInterface.h"
#include "PairHMMHostInterface.h"
#include "PairHMMWorker.h"

using namespace blaze;

static Context<float>  g_ctxf;
static Context<double> g_ctxd;

PairHMMWorker::PairHMMWorker(
    PairHMMClient* client,
    int num_read, 
    int num_hap, 
    read_t* reads, 
    hap_t* haps): 
  client_(client),
  use_cpu_(false),
  num_read_(num_read),
  num_hap_(num_hap),
  host_reads_(reads),
  host_haps_(haps)
{
#ifndef NO_PROFILE
  PLACE_TIMER;
#endif

  output_.resize(num_read * num_hap);

  // calculate total cells
  uint64_t ttl_rl     = 0;
  uint64_t ttl_hl     = 0;
  uint64_t ttl_rl_cpu = 0;
  uint64_t ttl_hl_cpu = 0;

  uint64_t num_cell      = 0;
  uint64_t num_cell_cpu  = 0;
  uint64_t num_cell_fpga = 0;

  for (int i = 0; i < num_read; ++i) ttl_rl += reads[i].len;
  for (int i = 0; i < num_hap; ++i)  ttl_hl += haps[i].len;
  num_cell = ttl_rl * ttl_hl;
#ifndef NO_PROFILE
  ksight::ksight.add("num cells", num_cell);
#endif

  // empirical threshold, if data too small, don't
  // bother running fpga
  if (num_read < 32 || num_hap < 2 || num_cell < 5e6) {
    use_cpu_ = true;
    DLOG(INFO) << "Use AVX for all computation";
#ifndef NO_PROFILE
    ksight::ksight.add("num cells on cpu", num_cell);
#endif
    return;
  }

  // scan all reads to make sure it does not violate FPGA constraint
  int max_haplen = 0;
  for (int i = 0; i < num_hap; ++i) {
    if (haps[i].len >= MAX_HAP_LEN) {
      cpu_hap_idx_.push_back(i);
      ttl_hl_cpu += haps[i].len;
    }
    else {
      max_haplen = std::max(max_haplen, haps[i].len);
      fpga_haps_.push_back(haps[i]);
      fpga_hap_idx_.push_back(i);
    }
  }

  for (int i = 0; i < num_read; ++i) {
    if (reads[i].len >= MAX_READ_LEN || 
        reads[i].len * (reads[i].len + max_haplen) <= 
          1 + READ_BLOCK_SIZE * MAX_READ_LEN) {
      cpu_read_idx_.push_back(i);
      ttl_rl_cpu += reads[i].len;
    }
    else {
      fpga_reads_.push_back(reads[i]);
      fpga_read_idx_.push_back(i);
    }
  }

  num_cell_cpu = ttl_rl_cpu * ttl_hl + 
                 ttl_rl * ttl_hl_cpu - 
                 ttl_hl_cpu * ttl_rl_cpu;

#ifndef NO_PROFILE
  ksight::ksight.add("num cells on cpu", num_cell_cpu);
#endif
  
  DLOG(INFO) << "cpu workload: num_read = " << cpu_read_idx_.size()
             << ", num_hap = " << cpu_hap_idx_.size();

  DLOG(INFO) << "fpga workload: num_read = " << fpga_read_idx_.size()
             << ", num_hap = " << fpga_hap_idx_.size();
}

PairHMMWorker:: ~PairHMMWorker() {
  DLOG(INFO) << "Destroying client";
}

inline testcase convert_avx_input(
    read_t & read,
    hap_t  & hap) 
{
  testcase ret;
  ret.rslen  = read.len; 
  ret.rs     = read._b; 
  ret.i      = read._i;
  ret.d      = read._d;
  ret.c      = read._c;
  ret.q      = read._q;
  ret.haplen = hap.len;
  ret.hap    = hap._b;
  return ret;
}

// avx compute routine
void PairHMMWorker::compute() {
#ifndef NO_PROFILE
  ksight::AutoTimer __timer("compute on cpu");
#endif

  if (use_cpu_) {
    for (int i = 0; i < num_read_; ++i) {
      for (int j = 0; j < num_hap_; ++j) {
        // use GKL routines to perform computation
        testcase avx_v = convert_avx_input(host_reads_[i], host_haps_[j]);
        output_[i*num_hap_ + j] = compute_fp_avxs(&avx_v);
      }
    }
  }
  else {
    for (auto i : cpu_read_idx_) {
      for (auto j : cpu_hap_idx_) {
        // use GKL routines to perform computation
        testcase avx_v = convert_avx_input(host_reads_[i], host_haps_[j]);
        output_[i * num_hap_ + j] = compute_fp_avxs(&avx_v);
      } 
    }
  }
}

// NOTE: need to guarantee output has size num_hap * num_read
void PairHMMWorker::getOutput(double* output) {
#ifndef NO_PROFILE
  PLACE_TIMER;
#endif

  // indexes to cpu output and fpga output
  int cpu_row  = 0;
  int cpu_col  = 0;
  int fpga_row = 0;
  int fpga_col = 0;

  int num_recal = 0;
  uint64_t recal_etime = 0;

  for (int i = 0; i < num_read_; ++i) {
    for (int j = 0; j < num_hap_; ++j) {
      float v = output_[i * num_hap_ + j];
        
      // set output values
      if (v < MIN_ACCEPTED) {
        //Timer timer("AVX Double Recalculate");

        //uint64_t start_ts = blaze::getUs();
        // recalculate using double precision 
        testcase avx_v = convert_avx_input(host_reads_[i], host_haps_[j]);
        double d = compute_fp_avxd(&avx_v);

        output[i * num_hap_ + j] = log10(d) - g_ctxd.LOG10_INITIAL_CONSTANT;

        //recal_etime += blaze::getUs() - start_ts;
        ++num_recal;
      }
      else {
        output[i * num_hap_ + j] = (double)(log10f(v) - g_ctxf.LOG10_INITIAL_CONSTANT);
      }
    }
  }
  DLOG_IF(INFO, num_recal > 0) << "Recalculate in double precision " << num_recal
             << "/" << num_read_ * num_hap_ << " times";
             //<< " times, taking " << recal_etime << " us";
} 

// in this version, use a separate thread to run compute()
// for all reads/haps that is impossible for fpga to compute
void PairHMMWorker::run() {
  DLOG(INFO) << "Computing mode: " << (use_cpu_ ? "CPU" : "CFX");
  if (use_cpu_) {
    compute();
  }
  else {
    //PLACE_TIMER1("compute on fpga");
#ifndef NO_PROFILE
    ksight::Timer timer;
    timer.start();
#endif

    // start a thread to run cpu
    boost::thread t(boost::bind(&PairHMMWorker::compute, this));

    // each fpga client invokation takes a maximum data size
    for (int row = 0; row < fpga_reads_.size(); row += MAX_RSDATA_NUM) {
      for (int col = 0; col < fpga_haps_.size(); col += MAX_HAPDATA_NUM) {

        int read_bound = std::min((int)fpga_reads_.size(), row + MAX_RSDATA_NUM);
        int hap_bound  = std::min((int)fpga_haps_.size(), col + MAX_HAPDATA_NUM);

        // start performing fpga computation
        int num_read = read_bound - row;
        int num_hap  = hap_bound - col;
        {
#ifndef NO_PROFILE
          PLACE_TIMER1("serialize input");
#endif

          client_->setup(
              &fpga_reads_[row], num_read,
              &fpga_haps_[col],  num_hap);
        }

        // start fpga run
        {
#ifndef NO_PROFILE
          PLACE_TIMER1("start()");
#endif
          client_->start();
        }

        // if cpu fallback is called, skip output copy
        {
#ifndef NO_PROFILE
          PLACE_TIMER1("copy output");
#endif
          
          // process output with data copy
          float* results = (float*)client_->getOutputPtr(0); 

          for (int i = row; i < read_bound; ++i) {
            for (int j = col; j < hap_bound; ++j) {
              // write output to global output buffer     
              int ii = fpga_read_idx_[i];
              int jj = fpga_hap_idx_[j];
              output_[ii * num_hap_ + jj] = 
                results[(i - row) * num_hap + j - col];
            }
          }
        }
      }
    }
    // don't count cpu time
#ifndef NO_PROFILE
    ksight::ksight.add("compute on client", timer.stop());
#endif
    t.join();
  }
}
