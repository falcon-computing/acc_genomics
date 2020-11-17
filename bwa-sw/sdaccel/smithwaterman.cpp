#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <ap_int.h>
#include <string.h>
#include <ap_utils.h>
#include <hls_stream.h>
using namespace hls;

typedef ap_uint<2> uint2_t;
typedef ap_uint<3> uint3_t;
typedef ap_uint<1> uint1_t;
typedef ap_uint<4> uint4_t;
typedef ap_uint<10> uint10_t;
typedef ap_uint<11> uint11_t;
typedef ap_int<11> int11_t;
typedef ap_uint<33> uint33_t;
typedef ap_uint<29> uint29_t;
typedef unsigned char uint8_t;
typedef short int16_t;
typedef unsigned short uint16_t;

#define READ_PE_NUM 3
#define SEED_PE_NUM 20
#define BLOCK_NUM 6
#define READ_BLOCK_NUM 2
#define BLOCK_PE_NUM 10
#define o_del 6
#define e_del 1
#define o_ins 6
#define e_ins 1
#define pen_clip 5
#define w_in 100
#define zero (int11_t)(0)

template <typename T>
inline void increment(T& cnt, T max) {
  cnt ++;
  if (cnt >= max) cnt = 0;
}

template <typename T>
inline T readStream(stream<T>& input) {
  T tmp;
  input.read(tmp);
  return tmp;
}

inline static void nextPE(char& pe_idx, char total_pe = READ_PE_NUM) {
  if ( pe_idx == total_pe -1 ) {
    pe_idx = 0;
  }
  else {
    pe_idx ++;
  }
}

void sw_extend(uint8_t qs_baddr, uint11_t ts_baddr, uint3_t *seed_seq, 
    uint8_t qlen, uint11_t tlen, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins,
    int16_t max_del, int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret,
    int16_t *gscore_ret, int16_t *maxoff_ret);

void data_parse(int *totalinp, int total_size, 
    stream<int> readTask[READ_PE_NUM],
    stream<uint29_t> &chain_rseq_span,
    stream<uint1_t> readTask_ctrl[READ_PE_NUM],
    stream<uint1_t> &chain_rseq_span_ctrl);

void seed_proc(stream<uint1_t>& seed_ctrl,
    uint3_t* seed_seq,
    stream<uint16_t>& seed_param,
    stream<int>& match);

void sw_extend(uint8_t qs_baddr, uint11_t ts_baddr, uint3_t *seed_seq, 
    uint8_t qlen, uint11_t tlen, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins,
    int16_t max_del, int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret,
    int16_t *gscore_ret, int16_t *maxoff_ret)
{
  uint11_t i;
  uint8_t j;
  ap_uint<2> k;
  ap_int<12> max_i, max_ie, max_off;
  ap_int<12> gscore;
  uint8_t qs_baddr_t;
  uint11_t ts_baddr_t;
  short max_j;
  char oe_del = o_del + e_del;
  char oe_ins = o_ins + e_ins;
  uint8_t beg, end;
  uint8_t backw_tmp;
  char forw_update;
  uint8_t forw_tmp;
  ap_int<10> abs_mj_m_i;
  char tmp_ehh_m_eins;
  int11_t tmp_eme;
  int11_t h1_init_val;
  int11_t max;
  int11_t h, e, M;
  int11_t e_tmp;
  int11_t h_tmp;
  int11_t h1_reg;
  int11_t t, f, h1, m;
  short mj;
  uint3_t q_i, q_j;
  ap_int<10> prev;
  char isBreak;
  uint8_t aw1;
  uint8_t aw_tmp;

  char h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0

  const char my_mat[5][5]={{1, -4, -4, -4, -1}, {-4, 1, -4, -4, -1}, {-4, -4, 1, -4, -1}, {-4, -4, -4, 1, -1}, {-1, -1, -1, -1, -1}};
#pragma HLS ARRAY_PARTITION variable=my_mat complete dim=0
  int11_t eh_h [512];
#pragma HLS ARRAY_MAP variable=eh_h instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_h core=RAM_2P_BRAM
  int11_t eh_e [512];
#pragma HLS ARRAY_MAP variable=eh_e instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_e core=RAM_2P_BRAM

  max = h0;
  max_i = max_j = -1;
  max_ie = -1;
  gscore = -1;
  max_off = 0;

  k = 0;
  isBreak = 0;
  for(j=0;j<=qlen;j++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=125 max=125
    eh_e[j]= 0;
    eh_h[j]= 0;
  }

ext_while_loop : 
  while ((k < 2) && (!isBreak)) {
#pragma HLS LOOP_TRIPCOUNT min=2 max=2
    prev = *regScore;
    aw_tmp = w_in << k;
    aw1 = aw_tmp < max_ins ? aw_tmp : max_ins;
    aw1 = aw1 < max_del ? aw1 : max_del;
    beg = 0;
    end = qlen;
    tmp_eme = h0 - oe_ins;
    tmp_eme = (tmp_eme > 0) ? tmp_eme : (int11_t)(0);
    h1_init_val = h0 - o_del;
target_loop : 
    for (i = 0; i < tlen; i++) {
#pragma HLS LOOP_TRIPCOUNT min=50 max=50
      f = 0; m = 0; mj = -1;
      // ts_baddr_t = ts_baddr + i;
      ts_baddr_t = ts_baddr + i;
      q_i = seed_seq[ts_baddr_t];

      if (beg < i - aw1) beg = i - aw1;
      //#pragma HLS resource variable=beg core=AddSub_DSP
      if (end > i + aw1 + 1) end = i + aw1 + 1;
      //#pragma HLS resource variable=end core=AddSub_DSP
      if (end > qlen) end = qlen;
      if(beg ==0){
        h1_init_val -= e_del;
        h1 = h1_init_val;
        if (h1 < 0) h1 = 0;
      }
      else h1 = 0;
      backw_tmp = 0; 
      forw_tmp = 0; 
      forw_update = 0;
query_loop : 
      for (j = beg; j < end; ++j) {
#pragma HLS LOOP_TRIPCOUNT min=50 max=50
#pragma HLS pipeline II=1
        #pragma AP dependence variable=eh_e array inter false
        #pragma AP dependence variable=eh_h array inter false
        qs_baddr_t = qs_baddr + j;
        q_j = seed_seq[qs_baddr_t];
        h_tmp = eh_h[j];// get H(i-1,j-1) and E(i-1,j)
        e_tmp = eh_e[j];
        if (i == 0) {
          e = 0;
          if (j == 0) {
            h = h0;
            M = h0;
          }
          else if (j == 1) {
            h = tmp_eme;
            M = tmp_eme;
          }
          else {
            tmp_eme -= e_ins;
            h = (tmp_eme > 0) ? tmp_eme : zero;
            M = (tmp_eme > 0) ? tmp_eme : zero;
          }
        }
        else {
          e = e_tmp;
          h = h_tmp;
          M = h_tmp;
        }
        h1_reg = h1;
        M = M?(int11_t)(M+my_mat[q_i][q_j]):zero;
        h = M > e? M : e;
        h = h > f? h : f;
        h1 = h;             // save H(i,j) to h1 for the next column

        t = M - oe_del;
        t = (t > 0) ? t : zero;
        e -= e_del;
        e = (e > t) ? e : t;   // computed E(i+1,j)
        t = M - oe_ins;
        t = t > 0? t : zero;
        f -= e_ins;
        f = (f > t) ? f : t;   // computed F(i,j+1)
        eh_e[j] = e; // save E(i+1,j) for the next row
        eh_h[j] = h1_reg;          // set H(i,j-1) for the next row
        if (m <= h)
        {
          mj = j;
          m = h;
        }
        if (forw_update == 0) { //((h1_reg == 0) &&
          if (h1_reg == 0 && e ==0) {
            forw_tmp++;
          }
          else {
            forw_update = 1;
          }
        }
        if(h1_reg ==0 && e ==0){
          backw_tmp++;
        } 
        else backw_tmp = 0;


      }
      eh_h[end] = h1;
      eh_e[end] = 0; 
      if( h1 == 0) 
        backw_tmp++; 
      else 
        backw_tmp = 0; 
      if (j == qlen) {
        if (gscore <= h1) {
          max_ie = i;
          gscore = h1;
        }
      }
      if (m == 0) break;
      if (m > max) {
        max = m;
        max_i = i;
        max_j = mj;
        abs_mj_m_i = abs(mj - i);
        //#pragma HLS resource variable=abs_mj_m_i core=AddSub_DSP
        if (max_off < abs_mj_m_i) max_off = abs_mj_m_i;
      }
      beg = beg + forw_tmp;
      end = end - backw_tmp + 2 <qlen ? end-backw_tmp +2:qlen;
    }
    *qle_ret = max_j + 1;
    *tle_ret = max_i + 1;
    *gtle_ret = max_ie + 1;
    *gscore_ret = gscore;
    *maxoff_ret = max_off;
    *regScore = max;
    if (max == prev || ( max_off < (aw_tmp >> 1) + (aw_tmp >> 2))) isBreak = 1;
    k++;
  }
  *w_ret = aw_tmp;
}

// distribute the data by read
void data_parse(int *totalinp, int total_size, 
    stream<int> read_task[READ_PE_NUM],
    stream<uint29_t> &chain_rseq_span,
    stream<uint1_t> read_ctrl[READ_PE_NUM],
    stream<uint1_t> &chain_rseq_span_ctrl)
{
#pragma HLS inline off
  int readEndIndex = 0;
  char pe_idx = 0;

  int read_seq_length;
  int read_seq_length_div8;
  int rseq_length;
  int rseq_length_div8;

  uint33_t rmax_0;
  uint33_t rmax_1;

  uint11_t chain_num;
  uint11_t seed_num;

  int seed_index;
  int seed_qbeg;
  uint33_t seed_rbeg;
  int seed_len;

  uint11_t chain_count = 0;
  uint11_t seed_count = 0;
  uint11_t read_count = 0;

  uint4_t stage = 0;
  uint32_t tmpval_hi = 0;
  uint32_t tmpval_lo = 0;
  int stage_cnt = 0; // counter for substate

data_parse_loop :
  for (int idx = 0; idx < total_size; idx ++) {
#pragma HLS LOOP_TRIPCOUNT min=20000 max=20000
#pragma HLS pipeline
    // get an input from dram
    int data_input = totalinp[idx];

    if (idx < readEndIndex) {
      if (stage == 0) { // parse read param
        read_seq_length = data_input;
        read_seq_length_div8 = (read_seq_length + 7) >> 3;

        // write read seq length
        read_task[pe_idx].write(read_seq_length);

        read_count = 0;
        stage = 1;
      }
      else if (stage == 1) { // parse read seqs
        if (read_count < read_seq_length_div8) {
          // write read seq data
          read_task[pe_idx].write(data_input);
          read_count += 1 ;
        }
        else {
          chain_num = data_input;

          // write chain_num
          read_task[pe_idx].write(chain_num);

          read_count = 0;
          chain_count = 0;
          if (chain_num > 0) {
            stage = 2;
          }
          else {
            stage = 0;
          }
        } 
      }
      else if (stage == 2) { // parse the chain param
        const int state_cnt_max = 4;
        switch (stage_cnt) {
          case 0: // read in tmpval_lo
            tmpval_lo = data_input;
            read_task[pe_idx].write(data_input);
            chain_rseq_span.write(pe_idx);
            increment(stage_cnt, state_cnt_max);
            break;
          case 1: // read in tmpval_hi, set rmax_0
            tmpval_hi = (uint32_t)data_input;
            rmax_0 = ((uint33_t)tmpval_hi << 32) | tmpval_lo;
            chain_rseq_span.write(rmax_0>>4);
            //read_task[pe_idx].write(data_input);
            increment(stage_cnt, state_cnt_max);
            break;
          case 2: // read in tmpval_lo
            tmpval_lo = (uint32_t)data_input;

            //read_task[pe_idx].write(data_input);
            increment(stage_cnt, state_cnt_max);
            break;
          case 3: // read in tmpval_hi, set rmax_1
            tmpval_hi = (uint32_t)data_input;
            rmax_1 = ((uint33_t)tmpval_hi << 32) | tmpval_lo;

            read_task[pe_idx].write((rmax_1>>4)-(rmax_0>>4)+1);
            chain_rseq_span.write(rmax_1>>4);
            //read_task[pe_idx].write(data_input);
            //  this chain_rseq part should be independent with the parameters part
            chain_count += 1;
            stage = 3;
            increment(stage_cnt, state_cnt_max);
            break;
          default: ;
        }
      }
      else if (stage == 3) { // parse the chain seq
        seed_num = data_input;
        read_task[pe_idx].write(seed_num);
        read_count = 0;
        if (seed_num == 0) {
          if (chain_count < chain_num) {
            stage = 2;
          }
          else {
            stage = 0;
            chain_count = 0;
          }
        }
        else{
          stage = 4;
        }
      }
      else if (stage == 4){ // parse seed param
        const int state_cnt_max = 5;
        int leftRlen;
        int rightRlen;
        switch (stage_cnt) {
          case 0:
            seed_index = data_input;
            read_task[pe_idx].write(seed_index);
            increment(stage_cnt, state_cnt_max);
            break;
          case 1:
            tmpval_lo = (uint32_t)data_input;
            increment(stage_cnt, state_cnt_max);
            break;
          case 2:
            tmpval_hi = (uint32_t)data_input;
            seed_rbeg = ((uint33_t)tmpval_hi << 32) | tmpval_lo;

            leftRlen = seed_rbeg - rmax_0;
            read_task[pe_idx].write(leftRlen);

            increment(stage_cnt, state_cnt_max);
            break;
          case 3:
            seed_qbeg = data_input;
            read_task[pe_idx].write(seed_qbeg);
            increment(stage_cnt, state_cnt_max);
            break;
          case 4:
            seed_len = data_input;
            rightRlen = rmax_1 - seed_rbeg - seed_len;

            if (seed_len > 65535 || rightRlen > 65535) {
              fprintf(stderr, "seed_len or rightRlen overflow\n");
            }

            read_task[pe_idx].write(seed_len << 16 | rightRlen);

            seed_count += 1;
            if (seed_count >= seed_num) {
              seed_count = 0;
              if (chain_count < chain_num) {
                stage = 2;
              }
              else {
                stage = 0 ;
                chain_count = 0;
              }
            }
            increment(stage_cnt, state_cnt_max);
            break;
          default: ;
        }
      }
    }
    else {
      // update readEndIndx
      readEndIndex = data_input; 
      stage = 0;

      // switch to a new PE
      if (idx) {
        nextPE(pe_idx);
      }
    }
  }

  // write the end signal
  for (int i = 0; i < READ_PE_NUM; i++){
#pragma HLS UNROLL
    read_ctrl[i].write(1);
  }
  chain_rseq_span_ctrl.write(1);
}

#define MEMORY_BATCH_SIZE 100
void chain_rseq_proc ( 
    int* pac_input,
    stream<uint29_t> &chain_rseq_span,
    stream<uint1_t> &chain_rseq_span_ctrl,
    stream<int> chain_rseq_buffer[READ_PE_NUM]) {
  int pac_buffer[100];
  while (1) {
#pragma HLS LOOP_TRIPCOUNT min=1300 max=1300
    if (!chain_rseq_span_ctrl.empty() && chain_rseq_span.empty()) {
      uint1_t flag =0;
      chain_rseq_span_ctrl.read_nb(flag);
      break;
    }
    
    uint29_t lower_bound = 0;
    uint29_t upper_bound = 0;
    uint29_t pe_idx = 0;
    chain_rseq_span.read(pe_idx);
    chain_rseq_span.read(lower_bound);
    chain_rseq_span.read(upper_bound);
  
    memcpy(pac_buffer, &pac_input[lower_bound], sizeof(int)*(upper_bound - lower_bound + 1));
    for (uint29_t i=0; i<=upper_bound-lower_bound; i++) {
#pragma HLS LOOP_TRIPCOUNT min=15 max=30
#pragma HLS PIPELINE
      chain_rseq_buffer[pe_idx].write(pac_buffer[i]);
    }
  }
}

void seed_proc( stream<uint1_t>& seed_ctrl,
    stream<uint3_t>& seed_seq,
    stream<uint16_t>& seed_param,
    stream<uint1_t>& match_ctrl,
    stream<int>& match)
{
#pragma HLS inline off

  // for now I dont use flexible opts 
  uint8_t qlen[2];
#pragma HLS ARRAY_PARTITION variable=qlen complete dim=0
  uint11_t tlen[2];
#pragma HLS ARRAY_PARTITION variable=tlen complete dim=0
  int16_t max_ins[2];
#pragma HLS ARRAY_PARTITION variable=max_ins complete dim=0
  int16_t max_del[2];
#pragma HLS ARRAY_PARTITION variable=max_del complete dim=0
  int16_t h0 = 0;
  int16_t regScore = 0;

  int16_t qle =0;
  int16_t tle =0;
  int16_t gtle =0;
  int16_t gscore =0;
  int16_t maxoff =0;
  int16_t qBeg =0;
  int16_t rBeg=0;
  int16_t qEnd=0;
  int16_t rEnd=0;
  int16_t score=0;
  int16_t trueScore=0;
  int16_t width=0;
  int16_t aw[2];
#pragma HLS ARRAY_PARTITION variable=aw complete dim=0
  int16_t sc0 =0;
  int16_t h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0
  uint8_t qs_baddr;
  uint11_t ts_baddr;

  uint2_t i = 0;
  uint16_t leftQlen = 0;
  uint16_t leftRlen = 0;
  uint16_t rightQlen = 0;
  uint16_t rightRlen = 0;
  uint16_t seed_len = 0;
  uint16_t seed_qbeg = 0;
  uint16_t seed_index = 0;
  uint3_t seq_mem[2048];

seed_proc_outer:
  while (1) {
#pragma HLS LOOP_TRIPCOUNT min=80 max=80
    if (!seed_ctrl.empty() && seed_seq.empty()) {
      uint1_t flag = 0;
      seed_ctrl.read_nb(flag);
      match_ctrl.write(1);
      break;
    }

    seed_param.read(leftQlen);
    seed_param.read(leftRlen);
    seed_param.read(rightQlen);
    seed_param.read(rightRlen);
    seed_param.read(seed_len);
    seed_param.read(seed_qbeg);
    seed_param.read(seed_index);

seed_proc_seq:
    for (uint16_t k = 0; k<leftQlen + leftRlen + rightQlen + rightRlen; k++){
#pragma HLS LOOP_TRIPCOUNT min=250 max=250
#pragma HLS pipeline
      seed_seq.read(seq_mem[k]);
    }

    qlen[0] = leftQlen;
    tlen[0] = leftRlen;
    qlen[1] = rightQlen;
    tlen[1] = rightRlen;
    h0 = seed_len;
    regScore = seed_len;
    max_ins[0] = qlen[0]; // for static opt
    max_del[0] = qlen[0];
    max_ins[1] = qlen[1];
    max_del[1] = qlen[1];
    aw[0] = w_in; 
    aw[1] = w_in;
    qBeg = 0;
    qEnd = qlen[1];
    rBeg = 0;
    rEnd = 0;
    trueScore = regScore;
    qle = -1;
    tle = -1;
    gtle = -1;
    gscore = -1;
    maxoff = -1;

    qs_baddr = 0;
    ts_baddr = qlen[0] + qlen[1];

seed_proc_sw_lr:
    for (i=0; i<2; i++){
      sc0 = regScore;
      h0_arr[0] = h0;
      h0_arr[1] = sc0;
      sw_extend(
          qs_baddr,
          ts_baddr,
          seq_mem,
          qlen[i],
          tlen[i],
          h0_arr[i],
          &regScore,
          seed_qbeg,
          max_ins[i],
          max_del[i],
          &aw[i],
          &qle,
          &tle,
          &gtle,
          &gscore,
          &maxoff
          );
      score = regScore;
      if (gscore <= 0 || gscore <= (regScore - pen_clip)) {
        if (i == 0) {
          qBeg = seed_qbeg - qle;
          rBeg = -tle;
          trueScore = regScore;
        }
        else {
          qEnd = qle;
          rEnd = tle;
          trueScore += regScore - sc0;
        }
      }
      else {
        if (i == 0) {
          qBeg = 0;
          rBeg = -gtle;
          trueScore = gscore;
        }
        else {
          qEnd = qlen[1];
          rEnd = gtle;
          trueScore += gscore - sc0;
        }
      }
      qs_baddr += qlen[i] ;
      ts_baddr += tlen[i] ;
    }
    if (aw[0] > aw[1]) width = aw[0];
    else width = aw[1];
    // get the results
    match.write(seed_index);
    match.write((qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000));
    match.write((rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000));
    match.write((score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000));
    match.write(width & 0xFFFF);
  }
}

void read_proc(
    // input data fifo
    stream<uint1_t> &read_ctrl,
    stream<int>     &read_data,
    stream<int>     &chain_rseq_buffer,
    // output data fifo
    stream<uint1_t> seed_ctrl[SEED_PE_NUM],
    stream<uint3_t> seed_seq[SEED_PE_NUM],
    stream<uint16_t> seed_param[SEED_PE_NUM]
    ) {
#pragma HLS inline off

  uint3_t read_seq[512];
  uint2_t chain_rseqs[2048];
#pragma HLS ARRAY_PARTITION variable=read_seq cyclic factor=8
#pragma HLS ARRAY_PARTITION variable=chain_rseqs cyclic factor=16
#pragma HLS resource variable=read_seq core=RAM_S2P_BRAM latency=4
//#pragma HLS resource variable=chain_rseqs core=RAM_S2P_BRAM latency=4

  char pe_idx = 0;

read_proc_outer:
  while (1) {
#pragma HLS LOOP_TRIPCOUNT min=200 max=200
    if (!read_ctrl.empty() && read_data.empty()) {
      uint1_t flag = 0;
      read_ctrl.read_nb(flag); // received the end signal
      break;
    }
   // int chain_num;


    uint16_t read_seq_len = 0;
    read_seq_len = readStream(read_data);
read_proc_read_seq:
    for (int i = 0; i < ((read_seq_len + 7) >> 3); i++) {
#pragma HLS LOOP_TRIPCOUNT min=20 max=20
#pragma HLS pipeline
      int data_input = readStream(read_data);
      //read_data.read(data_input);
      for (int k = 0; k < 8; k++) {
        read_seq[i*8 + k] = (data_input & 0xF0000000) >> 28 ; 
        data_input <<= 4;
      }
    }
    uint11_t chain_num;
    //read_data.read(chain_num);
    chain_num = readStream(read_data);
read_proc_chain_outer:
    for (uint11_t chain = 0; chain < chain_num; chain ++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=2
      int addr = 0;
      int length = 0;
      read_data.read(addr);
      read_data.read(length);
      uint4_t addr_t4 = (uint4_t)(addr);
      int first_int_size = 0;
      uint2_t m = 0;
      for (uint11_t i=0; i< length; i++) {
#pragma HLS LOOP_TRIPCOUNT min=12 max=12
#pragma HLS PIPELINE
        int ref_input;
        chain_rseq_buffer.read(ref_input);
        if (i==0 && (addr_t4 !=0)){
          for ( int j=0; j<16; j++) {
#pragma HLS UNROLL
            if (addr_t4 + j <16) {
              chain_rseqs[j] = ref_input >> ((addr_t4+j)<<1);
            }
            else if(addr_t4 + j == 16) {
              first_int_size = j; 
            }
          }
        }
        else if (first_int_size) {
          for (int k=0; k<16; k++) {
#pragma HLS UNROLL
            chain_rseqs[(i-1)*16+k+first_int_size] = ref_input >>(k<<1);
          }          
        }
        else {
          for (int k=0; k<16; k++) {
#pragma HLS UNROLL
            chain_rseqs[i*16+k] = ref_input>>(k<<1);
          }
        }
      }

      uint11_t seed_num = 0;
      //read_data.read(seed_num);
      seed_num = readStream(read_data);
read_proc_seed_outer:
      for (uint11_t seed = 0; seed < seed_num; seed ++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=2
        uint16_t seed_index = readStream(read_data);
        uint16_t leftRlen = readStream(read_data);
        uint16_t seed_qbeg = readStream(read_data);
        int tmp_data = readStream(read_data);
        uint16_t rightRlen = tmp_data & 0x0000FFFF;
        uint16_t seed_len = (ap_uint<32>)tmp_data >> 16;

        uint16_t leftQlen = seed_qbeg;
        uint16_t rightQlen = read_seq_len - seed_qbeg - seed_len;

        // send the params
        seed_param[pe_idx].write(leftQlen);
        seed_param[pe_idx].write(leftRlen);
        seed_param[pe_idx].write(rightQlen);
        seed_param[pe_idx].write(rightRlen);
        seed_param[pe_idx].write(seed_len);
        seed_param[pe_idx].write(seed_qbeg);
        seed_param[pe_idx].write(seed_index);

        uint8_t leftQlen_u8 = leftQlen;
        uint11_t leftRlen_u11 = leftRlen;
        uint8_t rightQlen_u8 = rightQlen;
        uint11_t rightRlen_u11 = rightRlen;
        // send the seqs
        ap_uint<12> seed_seq_length = leftQlen_u8 + leftRlen_u11 + rightQlen_u8 + rightRlen_u11 ;
read_proc_seed_seq:
        for (ap_uint<12> i = 0; i < seed_seq_length; i++){
#pragma HLS LOOP_TRIPCOUNT min=230 max=250
#pragma HLS pipeline II=1
          if ( i< leftQlen_u8){
            seed_seq[pe_idx].write(read_seq[leftQlen_u8-1-i]);
          }
          else if (i < leftQlen_u8 + rightQlen_u8){
            // seed_qbeg + seed_len +i - leftQlen
            seed_seq[pe_idx].write(read_seq[seed_len+i]);
          }
          else if (i < leftQlen_u8 + leftRlen_u11 + rightQlen_u8){
            // leftRlen -1 -(i- leftQlen - rightQlen) 
            seed_seq[pe_idx].write(chain_rseqs[leftRlen_u11 +leftQlen_u8 +rightQlen_u8-1-i]);
          }
          else {
            // seed_rbeg + seed_len - rmax0 + i - leftQlen - rightQlen - leftRlen
            seed_seq[pe_idx].write(chain_rseqs[2*seed_len - read_seq_len + i]);
          }
        }
        nextPE(pe_idx, SEED_PE_NUM);
      }
    }
  }
  for (int i = 0; i < SEED_PE_NUM; i++){
#pragma HLS unroll
    seed_ctrl[i].write(1);
  }
}

void result_proc(
    stream<uint1_t> match_ctrl[READ_PE_NUM][SEED_PE_NUM],
    stream<int>     match_data[READ_PE_NUM][SEED_PE_NUM],
    int* output_buf) 
{
#pragma HLS inline off

  int pe_finish_cnt = 0;
  bool pe_finish[READ_PE_NUM][SEED_PE_NUM];
  bool pe_final[READ_PE_NUM][SEED_PE_NUM];

#pragma HLS resource variable=pe_finish core=RAM_2P_LUTRAM latency=4
#pragma HLS resource variable=pe_final core=RAM_2P_LUTRAM latency=4

  // initialize pe status array
  for (int i = 0; i < READ_PE_NUM; i++) {
    for (int j = 0; j < SEED_PE_NUM; j++) {
#pragma HLS pipeline
      pe_finish[i][j] = false;
    }
  }

  int output_offset = 0;
  int output_idx = 0;

  //const int buf_size = READ_PE_NUM*SEED_PE_NUM*5;
  //int buf_idx = 0;
  //int result_buf[buf_size];

  // scan all PE and put results to buffer
  while (pe_finish_cnt < READ_PE_NUM*SEED_PE_NUM) {
    for (int i = 0; i < READ_PE_NUM; i ++) {
      for (int j = 0; j < SEED_PE_NUM; j ++) {
        if (!pe_finish[i][j]) {
          if (!match_ctrl[i][j].empty() && match_data[i][j].empty()) {
            uint1_t flag;
            match_ctrl[i][j].read_nb(flag);
            pe_finish[i][j] = true;
            pe_finish_cnt ++;
          }
          else if (!match_data[i][j].empty()){
            for (int k = 0; k < 5; k ++) {
#pragma HLS pipeline
              int tmp;
              match_data[i][j].read(tmp);
              output_buf[output_offset+k] = tmp;
            }
            output_offset += 5;
          }
        }
      }
    }     
  }
}

void match_assemble(stream<int> match_data[BLOCK_PE_NUM], 
    stream<uint1_t> match_ctrl[BLOCK_PE_NUM],
    stream<int>& block_data,
    stream<uint1_t>& block_ctrl) 
{
  ap_uint<4> i = 0;
  ap_uint<4> pe_finish_cnt = 0;
  int one_match_buffer[5];
match_assemble_outer:
  while(pe_finish_cnt < BLOCK_PE_NUM) {
match_assemble_inner:
    for ( i=0; i<BLOCK_PE_NUM; i++ ) {
      if (!match_data[i].empty()){
        for (int k = 0; k < 5; k++) {
#pragma HLS PIPELINE II=1
          match_data[i].read(one_match_buffer[k]);
        }
        for (int k =0; k < 5; k++) {
#pragma HLS PIPELINE II=1
          block_data.write(one_match_buffer[k]);
        }
      }  
      else if(!match_ctrl[i].empty()) {
        uint1_t flag;
        match_ctrl[i].read_nb(flag);
        pe_finish_cnt++;
      }
    }
  }
  block_ctrl.write(1);
}

void block_assemble(stream<int> block_data[BLOCK_NUM],
    stream<uint1_t> block_ctrl[BLOCK_NUM],
    stream<int> &results,
    stream<uint1_t> &results_ctrl
    )
{
  ap_uint<4> i =0;
  ap_uint<4> block_finish_cnt = 0;
  int one_block_buffer[5];
  while(block_finish_cnt < BLOCK_NUM) {
    for ( i=0; i<BLOCK_NUM; i++) {
      if(!block_data[i].empty()) {
        for (ap_uint<3> k =0; k<5; k++){
#pragma HLS PIPELINE II=1
          block_data[i].read(one_block_buffer[k]);
        }
        for (ap_uint<3> k =0; k<5; k++){
#pragma HLS PIPELINE II=1
          results.write(one_block_buffer[k]);
        }
      }
      else if (!block_ctrl[i].empty()) {
        uint1_t flag;
        block_ctrl[i].read_nb(flag);
        block_finish_cnt++;
      }
    }
  }
  results_ctrl.write(1);
}

void upload_results(stream<int> &results, stream<uint1_t> &results_ctrl, int *output_a)
{
  int output_offset = 0;
  while(1)
  {
    if(!results.empty()) {
      int tmp;
      results.read(tmp);
      output_a[output_offset] = tmp;
      output_offset++;
    }
    else if(!results_ctrl.empty()) {
      uint1_t flag;
      results_ctrl.read_nb(flag);
      break;
    }
  }
}



void sw_core(int *a, int *output_a, int *pac_input_a, int total_size) {
#pragma HLS inline off
#pragma HLS dataflow

  const int EXTRA = 0;
  stream<int> readTask[READ_PE_NUM];
#pragma HLS STREAM variable=readTask depth=2048+EXTRA
  stream<uint29_t> chain_rseq_span;
#pragma HLS STREAM variable=chain_rseq_span depth=128+EXTRA
  stream<uint1_t> chain_rseq_span_ctrl;
#pragma HLS STREAM variable=chain_rseq_span_ctrl depth=1+EXTRA
  stream<int> chain_rseq_buffer[READ_PE_NUM];
#pragma HLS STREAM variable=chain_rseq_buffer depth=2048+EXTRA
  stream<uint1_t> readTask_ctrl[READ_PE_NUM];
#pragma HLS STREAM variable=readTask_ctrl depth=1+EXTRA
  stream<uint1_t> seed_ctrl[READ_PE_NUM][SEED_PE_NUM];
#pragma HLS STREAM variable=seed_ctrl depth=1+EXTRA
  stream<uint3_t> seed_seq[READ_PE_NUM][SEED_PE_NUM];
#pragma HLS STREAM variable=seed_seq depth=1024+EXTRA
  stream<uint16_t> seed_param[READ_PE_NUM][SEED_PE_NUM];
#pragma HLS STREAM variable=seed_param depth=1024+EXTRA
  stream<int> match_data[READ_PE_NUM][SEED_PE_NUM];
#pragma HLS STREAM variable=match_data depth=128+EXTRA
  stream<uint1_t> match_ctrl[READ_PE_NUM][SEED_PE_NUM];
#pragma HLS STREAM variable=match_ctrl depth=1+EXTRA
  stream<int> block_data[BLOCK_NUM];
#pragma HLS STREAM variable=block_data depth=512+EXTRA
  stream<uint1_t> block_ctrl[BLOCK_NUM];
#pragma HLS STREAM variable=block_ctrl depth=1+EXTRA
  stream<int> results;
#pragma HLS STREAM variable=results depth=1024+EXTRA
  stream<uint1_t> results_ctrl;
#pragma HLS STREAM variable=results_ctrl depth=1+EXTRA

  data_parse(a, total_size, readTask, chain_rseq_span, readTask_ctrl, chain_rseq_span_ctrl);
  chain_rseq_proc(pac_input_a, chain_rseq_span, chain_rseq_span_ctrl, chain_rseq_buffer);

  for (int i = 0; i < READ_PE_NUM; i++) {
#pragma HLS UNROLL
    read_proc(readTask_ctrl[i], readTask[i], chain_rseq_buffer[i],
        seed_ctrl[i], seed_seq[i], seed_param[i]);

    for (int j = 0; j < SEED_PE_NUM; j++) {
#pragma HLS UNROLL
      seed_proc(
          seed_ctrl[i][j],
          seed_seq[i][j],
          seed_param[i][j],
          match_ctrl[i][j],
          match_data[i][j]);
    }
    for (int k =0; k < READ_BLOCK_NUM; k++) {
#pragma HLS UNROLL
      match_assemble(
          &match_data[i][BLOCK_PE_NUM*k],
          &match_ctrl[i][BLOCK_PE_NUM*k],
          block_data[i*READ_BLOCK_NUM + k],
          block_ctrl[i*READ_BLOCK_NUM + k]
          );
    }
  }
//  match_assemble(&match_data[0][0],  &match_ctrl[0][0],  block_data[0], block_ctrl[0]);
//  match_assemble(&match_data[0][10], &match_ctrl[0][10], block_data[1], block_ctrl[1]);
//  match_assemble(&match_data[0][20], &match_ctrl[0][20], block_data[2], block_ctrl[2]);
//
//  match_assemble(&match_data[1][0],  &match_ctrl[1][0],  block_data[3], block_ctrl[3]);
//  match_assemble(&match_data[1][10], &match_ctrl[1][10], block_data[4], block_ctrl[4]);
//  match_assemble(&match_data[1][20], &match_ctrl[1][20], block_data[5], block_ctrl[5]);

  block_assemble(block_data, block_ctrl, results, results_ctrl);
  upload_results(results, results_ctrl, output_a);
}

void sw_core_1(int *a, int *output_a, int *pac_input_a, int total_size) {
  sw_core(a, output_a, pac_input_a, total_size);
}

void sw_core_2(int *a, int *output_a, int *pac_input_a, int total_size) {
  sw_core(a, output_a, pac_input_a, total_size);
}

#ifndef HLS_
extern "C" {
#endif
  void sw_top(
      int *input_a, 
      int *input_b, 
      int *output_a, 
      int *output_b, 
      int *pac_input_a,
      int *pac_input_b,
      int size_a,
      int size_b) 
  {
#ifndef HLS_
#pragma HLS INTERFACE m_axi port=input_a offset=slave bundle=gmem_a
#pragma HLS INTERFACE m_axi port=output_a offset=slave bundle=gmem_a
#pragma HLS INTERFACE m_axi port=pac_input_a offset=slave bundle=gmem_c
#pragma HLS INTERFACE m_axi port=input_b offset=slave bundle=gmem_b
#pragma HLS INTERFACE m_axi port=output_b offset=slave bundle=gmem_b
#pragma HLS INTERFACE m_axi port=pac_input_b offset=slave bundle=gmem_d
#else
#pragma HLS INTERFACE m_axi port=input_a offset=slave bundle=gmem depth=252916
#pragma HLS INTERFACE m_axi port=input_b offset=slave bundle=gmem depth=252916
#pragma HLS INTERFACE m_axi port=output_a offset=slave bundle=gmem depth=8037*5
#pragma HLS INTERFACE m_axi port=output_b offset=slave bundle=gmem depth=8037*5
#endif
#pragma HLS INTERFACE s_axilite port=input_a bundle=control
#pragma HLS INTERFACE s_axilite port=input_b bundle=control
#pragma HLS INTERFACE s_axilite port=pac_input_a bundle=control
#pragma HLS INTERFACE s_axilite port=pac_input_b bundle=control
#pragma HLS INTERFACE s_axilite port=output_a bundle=control
#pragma HLS INTERFACE s_axilite port=output_b bundle=control
#pragma HLS INTERFACE s_axilite port=size_a bundle=control
#pragma HLS INTERFACE s_axilite port=size_b bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    sw_core_1(input_a, output_a, pac_input_a, size_a);
    sw_core_2(input_b, output_b, pac_input_b, size_b);
  }
#ifndef HLS_
}
#endif
