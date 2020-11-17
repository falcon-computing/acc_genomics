#ifndef PairHMMCLIENT_H
#define PairHMMCLIENT_H
#include <stdexcept>

#include "blaze/Client.h"
#include "ksight/tools.h"
#include "PairHMMHostInterface.h"

class PairHMMClient : public blaze::Client {

 public:
  PairHMMClient();

  void setup(read_t* reads, int num_read, 
             hap_t*  haps,  int num_hap);

  void compute();

  int getMaxInputNumItems(int idx);
  int getMaxInputItemLength(int idx);
  int getMaxInputDataWidth(int idx);

 private:
  // used to perform easy calculation on cpu
  int      num_read_;
  int      num_hap_;
  uint64_t num_cell_;
  read_t*  reads_;
  hap_t*   haps_;
};

#endif
