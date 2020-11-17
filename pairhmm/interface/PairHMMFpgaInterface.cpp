#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <functional>
#include <queue>
#include <math.h>
#include <vector>

#include "gkl-pairhmm/Context.h"
#include "PairHMMFpgaInterface.h"

#include "ksight/tools.h"

double curr_est1_gcups = 0;
double curr_est2_gcups = 0;
double peak_est1_gcups = 0;
double peak_est2_gcups = 0;

static Context<float>  g_ctxf;
static Context<double> g_ctxd;

static inline HaplenBundle pack_haplen(int len) {

  HaplenBundle ret;
  ret.haplen = len;
  ret.one_over_haplen = g_ctxf.INITIAL_CONSTANT / (float)len;

  return ret;
}

static float m2m_table[] = M2M_INIT;

static inline void pack_read_bundle(const read_t &r, ReadBundle* b) {
  int len = r.len;

  for (int i = 0; i < len; i++) {
    uint8_t _rs = read2char(r._b[i]);
    uint8_t _q  = r._q[i];
    uint8_t _i  = r._i[i];
    uint8_t _d  = r._d[i];
    uint8_t _c  = r._c[i];

    uint32_t _score = ((uint32_t)(_rs & 0x7) << 28) | 
                     ((uint32_t)(_q & 127) << 21) | 
                     ((uint32_t)(_i & 127) << 14) | 
                     ((uint32_t)(_d & 127) << 7) |
                     (uint32_t)(_c & 127);

    float _m2m = m2m_table[128 * _i + _d];

    b[i].scores = _score;
    b[i].m2m    = _m2m;
  }
}

// return a bundled 64-bit int
static inline uint64_t serialize_read_info(ReadInfo &v) { 
  uint64_t ret = 0;
  ret = v.read_len[0] + (v.read_len[1] << 8) + 
        ((uint64_t)v.output_idx << 16) + 
        (((uint64_t)v.trip_count - 1) << 37) + 
        ((uint64_t)v.is_full << 58);
  return ret;
}

std::vector<ReadInfo> prepare_readinfo(
    const int   num_read,
    const int   num_hap,
    const int   max_haplen,
    const read_t* reads,
    const hap_t*  haps)
{
  PLACE_TIMER;
  std::vector<ReadInfo> readinfo((num_read+READ_BLOCK_SIZE-1) / READ_BLOCK_SIZE);

  // prepare read info
  for (int i = 0; i < num_read; i += READ_BLOCK_SIZE) {
    ReadInfo* r = &readinfo[i / READ_BLOCK_SIZE];  // idx to read info list

    // until proven otherwize
    r->is_full = true;

    // should work for READ_BLOCK_SIZE=2, but not sure for larger
    uint8_t max_read_len = 0;
    for (int j = 0; j < READ_BLOCK_SIZE; j++) {
      if (i + j < num_read) {
        r->read_len[j] = reads[i+j].len; 
        //printf("%d: %d\n", i+j, reads[i+j].bases.size());
      }
      else {
        r->read_len[j] = 0;
        r->is_full = false;
      }
      max_read_len = std::max(max_read_len, r->read_len[j]);
    }
    r->trip_count = (std::max(max_read_len + 1, DEP_DIST) + 1) * (max_read_len + 1 + max_haplen);
    r->output_idx = num_hap * i;
  }

  return readinfo;
}

uint64_t total_tripcount = 0;

std::vector<PU> dist_reads_to_pu(
    const int num_pu,
    const int max_haplen,
    std::vector<ReadInfo> &readinfo)
{
  PLACE_TIMER; 

  // initialize a heap to storage pu workloads
  // always assign reads to the pu with the least
  // amount of work
  auto cmp = [](PU &l, PU &r) { 
    return l.trip_count > r.trip_count;
  };

  // priority queue must use an actual object rather than pointers
  //std::priority_queue<PU, std::vector<PU>, decltype(cmp)> q(cmp);
  std::queue<PU> q;

  std::vector<PU> pu_list(num_pu);

  // initialize 
  for (int i = 0; i < num_pu; i++) {
    pu_list[i].idx = i;
    pu_list[i].num_readinfo = 0;
    pu_list[i].readinfo_idx = (int*)malloc(MAX_RSDATA_NUM*sizeof(int));
    pu_list[i].trip_count = 0;

    q.push(pu_list[i]);
  }

  // distribute reads to pu
  int max_tripcount = 0;
  total_tripcount = 0;
  for (int i = 0; i < readinfo.size(); i ++) {

    // get the PU with the smallest tripcount
    //PU p = q.top();
    PU p = q.front();
    q.pop();

    // update tripcount and assignment
    p.trip_count += readinfo[i].trip_count;
    p.readinfo_idx[p.num_readinfo++] = i;
    
    pu_list[p.idx] = p;

    //max_tripcount = std::max(max_tripcount, p.idx*384 + readinfo[i].trip_count);
    max_tripcount = std::max(max_tripcount, readinfo[i].trip_count);

    if (i % num_pu == 0) {
      //printf("max_tripcount: %d\n", max_tripcount);
      total_tripcount += max_tripcount;
      max_tripcount = 0;
    }

    // this will heap-sort the PU according to the tripcount
    q.push(p); 
  }
  if (max_tripcount != 0) {
    //printf("max_tripcount: %d\n", max_tripcount);
    total_tripcount += max_tripcount;
  }

  return pu_list;
}

void pack_hap_data(
    int num_hap,
    const hap_t* haps,
    FpgaInputBundle *bundle) 
{
  PLACE_TIMER;
  assert(num_hap <= MAX_HAPDATA_NUM);
  assert(HAP_BLOCK_SIZE == 4);

  for (int i = 0; i < (num_hap + HAP_BLOCK_SIZE - 1) / HAP_BLOCK_SIZE; i ++) {
    for (int j = 0; j < MAX_HAP_LEN; j++) {
      int finished_haps = 0;
      // packing 4 hap (1 byte each) into a 16-bit uint
      uint16_t hap_pack = 0;
      for (int k = 0; k < HAP_BLOCK_SIZE; k ++) {
        uint8_t b = 0;
        int hap_idx = i * HAP_BLOCK_SIZE + k;
        if (hap_idx < num_hap && j < haps[hap_idx].len) {
          b = read2char(haps[hap_idx]._b[j]);
        }
        else {
          finished_haps += 1;
          b = read2char('N');
        }
        hap_pack = hap_pack | ((uint16_t)b << (k * 4));
      }
      bundle->hap_data[i][j] = hap_pack;
      if (finished_haps >= HAP_BLOCK_SIZE) break;
    } 
  }
}

void pack_read_data(
    const int   num_pu,
    const int   num_read,
    const int   num_hap,
    const int   max_haplen,
    const read_t* reads,
    const hap_t*  haps,
    FpgaInputBundle *bundle) 
{
  PLACE_TIMER;
  std::vector<ReadInfo> readinfo = prepare_readinfo(
            num_read, num_hap, max_haplen, 
            reads, haps);

  std::vector<PU> pu_list = dist_reads_to_pu(
            num_pu, max_haplen, readinfo);

  int hap_batch_size = (num_hap + HAP_BLOCK_SIZE - 1) / HAP_BLOCK_SIZE;

  uint64_t total_trip_count = 0;
  uint64_t max_trip_count = 0;
  for (int i = 0; i < num_pu; i ++) {
    PU p = pu_list[i];

    // write PU info to bundle
    bundle->trip_count_per_pu[i] = p.trip_count * hap_batch_size;
    bundle->num_reads_per_pu[i] = 0;

    total_trip_count += bundle->trip_count_per_pu[i];
    max_trip_count = std::max(max_trip_count, bundle->trip_count_per_pu[i]);
  }

  uint64_t total_rl = 0;
  uint64_t total_hl = 0;
  for (int i = 0; i < num_read; i++) {
    total_rl += reads[i].len;
  }
  for (int i = 0; i < num_hap; i++) {
    total_hl += haps[i].len;
  }

  uint64_t avg_trip_count = total_trip_count / num_pu;
  uint64_t estimate_time1 = avg_trip_count * 1000 / 185;
  uint64_t estimate_time2 = total_tripcount*hap_batch_size * 1000 / 185;
  //uint64_t estimate_time2 = max_trip_count * 1000 / 185;

  ksight::ksight.add("Estimated", estimate_time2);
  ksight::ksight.add("Estimated best", estimate_time1);

  uint64_t num_cell = total_rl * total_hl;

  curr_est1_gcups = (double)num_cell / estimate_time1;
  peak_est1_gcups = std::max(curr_est1_gcups, peak_est1_gcups);

  curr_est2_gcups = (double)num_cell / estimate_time2;
  peak_est2_gcups = std::max(curr_est2_gcups, peak_est2_gcups);

  // global index to the bundle data
  int read_info_idx = 0;
  while (read_info_idx < (num_read + READ_BLOCK_SIZE - 1) / READ_BLOCK_SIZE) {
    for (int i = 0; i < num_pu; i ++) {
      if (pu_list[i].num_readinfo == 0) {
        continue;
      }
      // read_info idx to PU assignments
      int idx = pu_list[i].readinfo_idx[--pu_list[i].num_readinfo];

      // pack readinfo and reads to bundle according to
      // PU assignments
      bundle->num_reads_per_pu[i] += readinfo[idx].is_full ? 2 : 1;

      // write read_info
      bundle->read_info[read_info_idx] = serialize_read_info(readinfo[idx]);

      // write read_data
      for (int j = 0; j < READ_BLOCK_SIZE; j++) {
        int read_idx = idx * READ_BLOCK_SIZE + j;
        if (j < 1 || readinfo[idx].is_full) {
          pack_read_bundle(reads[read_idx], 
              bundle->read_data[read_info_idx * READ_BLOCK_SIZE + j]);
        }
      }
      read_info_idx++;
    }
  }
  // free PU info
  for (int i = 0; i < num_pu; i ++) {
    //printf("PU%d: %016lx %d %d\n", i,
    //    bundle->read_info[i],
    //    bundle->trip_count_per_pu[i],
    //    bundle->num_reads_per_pu[i]);
    PU p = pu_list[i];
    free(p.readinfo_idx);
  }
}

void pack_fpga_input(
    const int   num_pu,
    const int   num_read,
    const int   num_hap,
    const read_t* reads,
    const hap_t*  haps,
    FpgaInputBundle *bundle) 
{
  PLACE_TIMER;
  assert(num_read <= MAX_RSDATA_NUM);
  assert(num_hap <= MAX_HAPDATA_NUM);

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  // pack haplen and calculate the max hap length
  int max_haplen = 0;
  for (int i = 0; i < num_hap; i++) {
    int haplen = haps[i].len;
    max_haplen = std::max(max_haplen, haplen + 1);

    bundle->hap_len[i] = pack_haplen(haplen);
  }

  // pack hap data
  pack_hap_data(num_hap, haps, bundle);

  // pack read data
  pack_read_data(num_pu, num_read, num_hap, max_haplen,
    reads, haps, bundle);

  //FILE* fout = fopen("test.dat", "wb+");
  //fwrite(bundle, sizeof(FpgaInputBundle), 1, fout);
  //fclose(fout);
}


