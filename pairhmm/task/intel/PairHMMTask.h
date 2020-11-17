#ifndef PAIRHMM_TASK_H
#define PAIRHMM_TASK_H

#include <CL/opencl.h>
#include <fstream>
#include <glog/logging.h>
#include <stdlib.h>
#include <stdexcept>
#include <string>
#include <time.h>
#include <vector>

#include "blaze/altr_opencl/OpenCLEnv.h" 
#include "blaze/Block.h" 
#include "blaze/Task.h" 
#include "PairHMMBuffer.h"
#include "PairHMMHostInterface.h"

#define MAX_READS          1024
#define MAX_READ_LENGTH    (32 * 1024)
#define MAX_HAPS           1024
#define MAX_HAP_LENGTH     (32 * 1024)

#define MAX_NUM_RESULTS    (MAX_READS * MAX_HAPS)

// assign position with OR, test position with AND
// example: position |= FIRST; if (position & LAST) ...
#define FIRST    (1 << 0)
#define LAST     (1 << 1)
#define OTHER    0

// 2016/09/15: GSP: make 0 INVALID base
// enum values for hap encoding
#define INVALID '-'

static Context<float> ctx;
typedef boost::shared_ptr<PairHMMBuffer>   PairHMMBuffer_ptr;

// read/hap data and origin index pair
typedef std::pair<read_t, int> Read;
typedef std::pair<hap_t,  int> Haplotype;

#define PRIM_TYPE_UCHAR cl_uchar
#define PRIM_TYPE_CHAR cl_char
#define PRIM_TYPE_SHORT cl_short
#define PRIM_TYPE_USHORT cl_ushort
#define PRIM_TYPE_FLOAT cl_float
#define PRIM_TYPE_INT cl_int
#define PRIM_TYPE_BOOL bool
#define PRIM_TYPE_UINT cl_uint

typedef struct {
    PRIM_TYPE_CHAR base;
    PRIM_TYPE_UCHAR position;
    PRIM_TYPE_USHORT hap_num;
    PRIM_TYPE_FLOAT y_init;
} HapData;

typedef struct {
    PRIM_TYPE_CHAR base;
    PRIM_TYPE_UCHAR position;
    PRIM_TYPE_USHORT read_num;
    PRIM_TYPE_FLOAT mx; // Pr(match to insert gap)
    PRIM_TYPE_FLOAT my; // Pr(match to delete gap)
    PRIM_TYPE_FLOAT gg; // Pr(gap to gap)
    PRIM_TYPE_FLOAT mm_1m_qual; // Pr(match to match) * (1 - read_qual)
    PRIM_TYPE_FLOAT mm_qual_div3; // Pr(match to match) * read_qual / 3
    PRIM_TYPE_FLOAT gm_1m_qual; // Pr(gap to match) * (1 - read_qual)
    PRIM_TYPE_FLOAT gm_qual_div3; // Pr(gap to match) * read_qual / 3
} ReadData;

typedef struct {
    PRIM_TYPE_USHORT result_read_num;
    PRIM_TYPE_USHORT result_hap_num;
    PRIM_TYPE_FLOAT result;
} ResultData;


class PairHMM : public blaze::Task {
 public:
  PairHMM();
  virtual ~PairHMM();

  virtual uint64_t estimatedClientTime() { return 0; }
  virtual uint64_t estimatedTaskTime() { return 0; }

  virtual void compute();
  virtual void prepare();

 private:
  blaze::OpenCLEnv* env;

  // kernel configurations
  int ROWS;
  int COLS;

  // raw inputs from client
  uint64_t num_cell_;
  int      num_read_;
  int      num_hap_;
  read_t*  reads_;
  hap_t*   haps_;

  // static pointers that stores the actual data
  std::map<std::string, cl_kernel> kernel_table_;
  std::map<std::string, uint>      arg_index_;
  std::map<std::string, PairHMMBuffer_ptr> buf_table_;

  blaze::DataBlock_ptr  output_;

  // dynamic variables used by computation
  ReadData   *read_data;
  HapData    *hap_data;
  ResultData *results;
  ushort     num_reads;
  ushort     num_haps;
  uint       total_read_data;
  uint       total_hap_data;
  size_t     pad_hap_data_cap;
  uint       max_read_length;
  uint       read_segment;
  uint       hap_segment;

  std::vector<std::vector<Read> >      concat_reads_array;
  std::vector<std::vector<Haplotype> > haps_array;
  std::vector<uint>                    max_read_lengths;

  int conf_str2int(std::string key) {
    std::string conf;
    get_conf(key, conf);
    return std::stoi(conf);
  }

  uint comp_adj_param( uint orig_val, uint incr_val )
  {
    uint mod_val = orig_val % incr_val;
    uint offset_val = (mod_val == 0) ? 0 : 1;
    uint adj_val = ((orig_val - mod_val) / incr_val) + offset_val;

    return adj_val;
  }

  size_t padded_read_length( size_t size )
  {
    return size + (size % COLS ? COLS - (size % COLS) : 0);
  }

  void add_buffer_arg(std::string kernel_name, std::string buffer_name) {
    cl_kernel kernel = kernel_table_[kernel_name];
    uint arg_index   = arg_index_[kernel_name];

    arg_index_[kernel_name]++;
    DVLOG(2) << "Add buf arg: " << kernel_name 
             << "(" << arg_index << ") = " << buffer_name;

    cl_mem buffer = buf_table_[buffer_name]->get_device_buffer();

    cl_int status = clSetKernelArg( kernel, (cl_uint) arg_index, sizeof(cl_mem), &buffer );
    if (status != CL_SUCCESS) {
      throw blaze::invalidParam("fails to set kernel args");
    }
  }

  void add_scalar_arg(std::string kernel_name, void* arg, size_t size) {
    cl_kernel kernel = kernel_table_[kernel_name];
    uint arg_index   = arg_index_[kernel_name];

    arg_index_[kernel_name]++;
    DVLOG(2) << "Add scalar arg: " << kernel_name << "(" << arg_index << ")";

    cl_int status = clSetKernelArg(kernel, (cl_uint) arg_index, size, arg);
    if (status != CL_SUCCESS) {
      throw blaze::invalidParam("fails to set kernel args");
    }
  }

  void run_task(std::string kernel_name, int q = 0) {
    cl_int status = clEnqueueTask(env->getCmdQueue(q), 
        kernel_table_[kernel_name], 0, NULL, NULL );

    arg_index_[kernel_name] = 0;

    if (status != CL_SUCCESS) {
      throw blaze::invalidParam(std::string("fails to enqueue kernel: ") + kernel_name);
    }
  }

  void finish(int q) {
    DVLOG(2) << "finish queue " << q;
    clFinish(env->getCmdQueue(q));
  }

  void create_buffer(std::string name, cl_int flag, size_t size) {
    PLACE_TIMER1(std::string("create buffer: ") + name);

    PairHMMBuffer_ptr p;
    if (!env->getScratch(name, p)) {
      PairHMMBuffer_ptr np(new PairHMMBuffer(env, size, flag));
      DVLOG(1) << "Creating " << name << " add add it to scratch";
      p = np;
    }
    else {
      DVLOG(1) << "Getting " << name << " from scratch";
    }
    buf_table_[name] = p; 
  }

  void write_buffer(std::string name, size_t size, int q = 0) {
    if (!buf_table_.count(name)) 
      throw blaze::invalidParam(std::string("cannot find buf: ") + name);

    DVLOG(2) << "write_buffer " << name << ", size = " << size << " using queue " << q;
    PairHMMBuffer_ptr p = buf_table_[name]; 
    p->write_buffer(size, q); 
  }

  void read_buffer(std::string name, size_t size, int q = 0) {
    if (!buf_table_.count(name)) 
      throw blaze::invalidParam(std::string("cannot find buf: ") + name);

    DVLOG(2) << "read_buffer " << name << ", size = " << size << " using queue " << q;
    PairHMMBuffer_ptr p = buf_table_[name]; 
    p->read_buffer(size, q); 
  }

  void reset_read_buffers( void )
  {
    read_data = static_cast<ReadData *>(buf_table_["read_data"]->get_host_buffer());
    total_read_data = 0;
    num_reads = 0;
  }

  void reset_hap_buffers( void )
  {
    hap_data = static_cast<HapData *>(buf_table_["hap_data"]->get_host_buffer());
    total_hap_data = 0;
    num_haps = 0;
    pad_hap_data_cap = 0;
  }

  void add_read( Read &read )
  {
    read_t r           = read.first;
    size_t read_length = r.len;
    for( size_t i = 0; i < read_length; i++ ) {
      read_data->read_num = read.second;
      read_data->base = r._b[i];
      read_data->position = OTHER;
      read_data->position |= (i == 0) ? FIRST : OTHER;
      read_data->position |= (i == read_length - 1) ? LAST : OTHER;

      float mm = ctx.set_mm_prob( r._i[i], r._d[i] );
      float gm = 1.0f - ctx.ph2pr[r._c[i]];
      read_data->mx = ctx.ph2pr[r._i[i]];
      read_data->my = ctx.ph2pr[r._d[i]];
      read_data->gg = ctx.ph2pr[r._c[i]];
      read_data->mm_1m_qual = mm * (1.0f - ctx.ph2pr[r._q[i]]);
      read_data->mm_qual_div3 = mm * ctx.ph2pr[r._q[i]] / 3.0f;
      read_data->gm_1m_qual = gm * (1.0f - ctx.ph2pr[r._q[i]]);
      read_data->gm_qual_div3 = gm * ctx.ph2pr[r._q[i]] / 3.0f;
      read_data++;
      total_read_data++;
    }

    size_t pad_length = (read_length % COLS) ? (COLS - (read_length % COLS)) : 0;
    DVLOG(2) << "read pad length = " << pad_length;

    for ( size_t i = 0; i < pad_length; i++ ) {
      read_data->base = INVALID;
      read_data->position = OTHER;
      read_data++;
      total_read_data++;
    }

    num_reads++;
  }

  void add_hap( Haplotype &hap )
  {
    hap_t  h = hap.first;
    size_t hap_length = h.len;

    if( hap_length <= pad_hap_data_cap) {
      // Entire hap sequence can fit in remaining space
      pad_hap_data_cap -= hap_length;
      DVLOG(2) <<  "updated hap pad cap = " << pad_hap_data_cap;
    } 
    else {
      // Need to pad the remaining space and align hap sequence to next row block
      for( size_t i = 0; i < pad_hap_data_cap; i++ ) {
        hap_data->base = (char) 0;
        hap_data->y_init = 0.0f;
        hap_data->position = OTHER;
        hap_data++;
        total_hap_data++;
      }

      pad_hap_data_cap = (hap_length % ROWS) ? (ROWS - (hap_length % ROWS)) : 0;
      DVLOG(2) << "new hap pad cap = " << pad_hap_data_cap;
    }
    float y_init = ctx.INITIAL_CONSTANT / hap_length; 
    DVLOG(2) << "add_hap #" << hap.second;
    for( size_t i = 0; i < hap_length; i++ ) {
      hap_data->base = h._b[i];
      hap_data->y_init = y_init;
      hap_data->hap_num = hap.second;
      hap_data->position = OTHER;
      hap_data->position |= (i == 0) ? FIRST : OTHER;
      hap_data->position |= (i == hap_length - 1) ? LAST : OTHER;

      hap_data++;
      total_hap_data++;
    }

    num_haps++;
  }
};

// define the constructor and destructor for dlopen()
extern "C" blaze::Task* create() {
  //return new PairHMM();
  return new PairHMM();
}

extern "C" void destroy(blaze::Task* p) {
  delete p;
}

#endif
