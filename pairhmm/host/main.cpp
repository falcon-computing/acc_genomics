#include <cassert>
#include <cmath>
#include <dirent.h> 
#include <fcntl.h>
#include <fstream>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include <glog/logging.h>
#include <iterator>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <vector>

#include "ksight/tools.h"

#ifdef LOCAL_BLAZE
#include "blaze/PlatformManager.h"
#include "blaze/AppCommManager.h"
#endif

#include "gkl-pairhmm/avx_impl.h"
#include "gkl-pairhmm/Context.h"

#ifdef CLIENT_MODE
#include "PairHMMWorker.h"
#include "PairHMMClient.h"
#else

extern double peak_kernel_gcups;
extern double curr_kernel_gcups;

static Context<float>  g_ctxf;
static Context<double> g_ctxd;

float* compute_fpga(
    const char* bit_path,
    std::string read_data,
    std::string hap_data,
    uint64_t num_cell);
#endif

#include "PairHMMHostInterface.h"
#include "PairHMMFpgaInterface.h"

int fail_num = 0;

inline std::vector<std::string> tok_nline(std::ifstream &fin) {
  if (!fin.good()) return std::vector<std::string>();
  std::string line;
  std::getline(fin, line);
  std::stringstream ss(line);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> ret(begin, end);

  return ret;
}

void get_input(
    int &num_read,
    int &num_hap,
    read_t * & reads,
    hap_t * & haps,
    const char* filename)
{
  std::ifstream ifs(filename, std::ifstream::in);
  std::string line;

  // first get num_reads and num_haps
  auto tokens = tok_nline(ifs);
  assert(tokens.size() == 4);
  num_read = std::stoi(tokens[1]);
  num_hap  = std::stoi(tokens[3]);

  // allocate read_t and hap_t
  reads = (read_t*)malloc(num_read * sizeof(read_t));
  haps = (hap_t*)malloc(num_hap * sizeof(hap_t));

  for (int i = 0; i < num_read; ++i) {

    auto tokens = tok_nline(ifs);
    assert(tokens.size() == 1);
    int read_len = std::stoi(tokens[0]);

    alloc_data(&reads[i], read_len);
    
    // skip line, then do the following order: _b, _r, _i, _d, _c
    for (int j = 0; j < 5; ++j) {
      std::getline(ifs, line); // to be ignored

      auto tokens = tok_nline(ifs);
      assert(tokens.size() == read_len);

      for  (int k = 0; k < read_len; ++k) {
        switch (j) {
          case 0: reads[i]._b[k] = (char)std::stoi(tokens[k]); break;
          case 1: reads[i]._q[k] = (char)std::stoi(tokens[k]); break;
          case 2: reads[i]._i[k] = (char)std::stoi(tokens[k]); break;
          case 3: reads[i]._d[k] = (char)std::stoi(tokens[k]); break;
          case 4: reads[i]._c[k] = (char)std::stoi(tokens[k]); break;
          default: ;
        }
      }
    }
  }

  std::getline(ifs, line); // skip an empty line

  for (int i = 0; i < num_hap; ++i) {
    auto tokens = tok_nline(ifs);
    assert(tokens.size() == 1);
    int hap_len = std::stoi(tokens[0]);

    haps[i].len = hap_len;

    std::getline(ifs, line); // skip the header line

    std::getline(ifs, line);
    assert(line.size() == hap_len);

    // hap bases are literal chars, so use strdup
    haps[i]._b = strdup(line.c_str());
  }
  ifs.close();
}

int get_output(
    double* likelihood, 
    int size, 
    const char* filename)
{
  std::ifstream ifs(filename, std::ifstream::in);

  if (!ifs.good()) {
    printf("bad file name %s\n", filename);
    return 1;
  }

  for (int i = 0; i < size; ++i) {
    double ref;
    ifs >> ref;
    union {
      long long i;
      double d;
    } value;
    ifs >> value.i;
    likelihood[i] = value.d;
  }

  return 0;
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

bool use_fpga(
    int num_read,
    int num_hap,
    read_t* reads, 
    hap_t*  haps) 
{
  if (num_read < 64 || num_hap < 4) return false;

  if (num_read > MAX_RSDATA_NUM ||
      num_hap > MAX_HAPDATA_NUM) {
    return false;
  }

  int max_haplen = 0;
  for (int i = 0; i < num_hap; i ++) {
    if (haps[i].len >= MAX_HAP_LEN) {
      return false;
    }
    max_haplen = std::max(max_haplen, haps[i].len);
  }
  for (int i = 0; i < num_read; i ++) {
    if (reads[i].len >= MAX_READ_LEN) {
      return false;
    }

    if (reads[i].len * (reads[i].len + max_haplen) <= 
        1 + READ_BLOCK_SIZE * MAX_READ_LEN) {
      return false;
    }
  }  
  return true;
}

typedef struct {
  uint64_t ts;
  int id;
  int num_read;
  int num_hap;
  read_t* reads;
  hap_t*  haps;
} pairhmm_input_t;

typedef struct {
  int id;
  uint64_t num_cell;
  int size;
  double* v;
} pairhmm_output_t;

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("USAGE: %s bit_path input_dir\n", argv[0]);
    return 1;
  }

  DIR           *d = opendir(argv[2]);
  struct dirent *dir;

  int test_num = 0;
  if (d) {
    while ((dir = readdir(d)) != NULL) {
      test_num ++;
    }
    closedir(d);
  }
  else {
    throw std::runtime_error("input folder is not valid");
  }

  // init for avx
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  ConvertChar::init();

  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
#ifdef LOCAL_BLAZE
  int file_handle = open(argv[1], O_RDONLY);
  if (file_handle < 0) {
    printf("cannot find configure file: %s\n",
        argv[3]);
    return -1;
  }

  google::protobuf::io::FileInputStream fin(file_handle);

  // config manager
  blaze::ManagerConf conf;
  if (!google::protobuf::TextFormat::Parse(&fin, &conf)) {
    throw std::runtime_error("cannot parse protobuf message");
  }
  FLAGS_v = conf.verbose();
  std::cout << conf.DebugString();

  // start manager
  blaze::PlatformManager platform_manager(&conf);
  blaze::AppCommManager comm(&platform_manager, "127.0.0.1", 1027); 
#endif

  test_num = (test_num - 2) / 2;
  uint64_t total_num_cell = 0;

  printf("test id, num_read, num_haps, est gcups, est best gcups, real gcups\n");

#ifdef CLIENT_MODE
  PairHMMClient client;
#endif

  //for (int i = 2346; i < 2347; i++) {
  int test_run_num = 0;
  for (int i = 0; i < test_num; i++) {
    
    char fname_in[4096];
    char fname_out[4096];
    sprintf(fname_in,  "%s/input%d",  argv[2], i);
    sprintf(fname_out, "%s/output%d", argv[2], i);
 
    read_t* reads = NULL;
    hap_t*  haps  = NULL;

    int num_read = 0;
    int num_hap  = 0;

    get_input(num_read, num_hap, reads, haps, fname_in);

    int output_size = num_read * num_hap;;
    //printf("Test #%d: %d, %d\n", i, num_read, num_hap);

    uint64_t total_rl = 0;
    uint64_t total_hl = 0;
    for (int i = 0; i < num_read; ++i) {
      total_rl += reads[i].len;
    }
    for (int i = 0; i < num_hap; ++i) {
      total_hl += haps[i].len;
    }
    uint64_t num_cell = total_rl * total_hl;

    if (!use_fpga(num_read, num_hap, reads, haps)) {
      continue;
    }
    DLOG(INFO) << "Run batch #" << i;

#ifndef CLIENT_MODE
    std::string read_data = serialize(reads, num_read);
    std::string hap_data  = serialize(haps, num_hap);

    float* output = NULL;
    {
      ksight::AutoTimer timer("compute::host");
      output = compute_fpga(argv[1], read_data, hap_data, num_cell);
    }

    if (output == NULL) {
      printf("Test #%d: Skipped\n", i);
      continue;
    }
    total_num_cell += num_cell;
#else
    total_num_cell += num_cell;
    double* output = (double*)malloc(output_size*sizeof(double));
    PairHMMWorker worker(&client, 
        num_read, num_hap,
        reads, haps);
    {
      ksight::AutoTimer timer("compute::host");
      worker.run();
      worker.getOutput(output);
    }
#endif
    DLOG(INFO) << "Finish batch #" << i;

#ifndef CLIENT_MODE
    printf("%d, %d, %d, %.3f, %.3f, %.3f\n", i, num_read, num_hap,
        curr_est2_gcups, curr_est1_gcups, curr_kernel_gcups);
#endif

#ifndef NO_COMPARISON
    double* baseline = new double[output_size];
    get_output(baseline, output_size, fname_out);

    // comparison
    int recal_count = 0;
    int error_count = 0;
    for (int k = 0; k < output_size; k++) {
#ifndef CLIENT_MODE
      double target;
      if (output[k] < MIN_ACCEPTED) {
        testcase avx_v = convert_avx_input(reads[k/num_hap], haps[k%num_hap]);
        double d = compute_fp_avxd(&avx_v);
        target = log10(d) - g_ctxd.LOG10_INITIAL_CONSTANT;
        recal_count++;
      }
      else
        target = log10(output[k]) - g_ctxf.LOG10_INITIAL_CONSTANT;
#else
      double target = output[k];
#endif
      
      if (std::isnan(target)) {
        printf("target is nan\n");
        error_count++;
      }
      double error = fabs((target - baseline[k]) / baseline[k]);
      if (error > 5e-3){ 
        printf("result has significant error, golden=%f, target=%f\n", baseline[k], target);
        error_count++;
      }
    }

    free_reads(reads, num_read);
    free_haps(haps, num_hap);

    if (error_count > 0) {
      fail_num ++;
      printf("Test #%d: %d errors\n", i, error_count);
    }
    else {
#ifndef CLIENT_MODE
      if (recal_count == 0)
        printf("Test #%d: Pass (no recalc)\n", i);
      else
        printf("Test #%d: Pass (recalc %d/%d)\n", i, recal_count, output_size);
#else
      printf("Test #%d: Pass\n", i);
#endif
    }
    test_run_num ++;

    delete [] baseline;
#ifdef CLIENT_MODE
    free(output);
#endif
#endif
  }

  printf("%d out of %d failed test\n", fail_num, test_run_num);
  if (fail_num != 0) {
    return 1;
  }
  ksight::ksight.print_total();

  //printf("Average estimated GCUPS: %.3f\n", (double)total_num_cell / ksight::ksight.get_total<uint64_t>("Estimated"));
  //printf("Average best GCUPS: %.3f\n", (double)total_num_cell  / ksight::ksight.get_total<uint64_t>("Estimated best"));
  //printf("Average kernel GCUPS: %.3f\n", (double)total_num_cell / ksight::ksight.get_total<uint64_t>("compute::kernel"));
  //printf("Average host GCUPS: %.3f\n", (double)total_num_cell  / ksight::ksight.get_total<uint64_t>("compute::host"));

  return 0;
}
