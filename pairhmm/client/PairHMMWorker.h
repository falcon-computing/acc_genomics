#ifndef PAIRHMMWORKER_H
#define PAIRHMMWORKER_H

#include <stdexcept>

#include "blaze/Client.h"
#include "PairHMMClient.h"

class PairHMMWorker {
 public:
  PairHMMWorker(PairHMMClient* client,
      int num_read, int num_hap, 
      read_t* reads, hap_t* haps);

  ~PairHMMWorker();

  // perform all computation
  void run();

  void getOutput(double* output);

  void compute();     // fallback compute() when client has problem

 private:
  void distribute() {;} // TODO: distribute workloads using wait time

  PairHMMClient* client_;

  bool use_cpu_;      // if true, use cpu compute everything
  int num_read_;
  int num_hap_;
  read_t* host_reads_;
  hap_t*  host_haps_;

  std::vector<int> cpu_read_idx_;
  std::vector<int> cpu_hap_idx_;

  std::vector<read_t> fpga_reads_;
  std::vector<hap_t>  fpga_haps_;
  // TODO: see if there are any better solutions than this
  std::vector<int>    fpga_read_idx_;
  std::vector<int>    fpga_hap_idx_;

  std::vector<float>  output_;

  uint64_t cpu_end_ts_;  // finish timestamp of cpu compute
  uint64_t fpga_end_ts_; // finish timestamp of fpga compute
};


#endif
