//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
//#include <ap_int.h>
//#include <string.h>
//#include <ap_utils.h>
#include "ihc_apint.h"
typedef ap_uint<1> uint1_t;
//typedef bool uint1_t;

//typedef unsigned char uint8_t;
//typedef short int16_t;
//typedef unsigned short uint16_t;

#ifndef K_NUM
#define K_NUM 2 
#endif

#ifndef READ_PE_NUM
#define READ_PE_NUM 2
#endif

#ifndef BLOCK_PE_NUM
#define BLOCK_PE_NUM 10
#endif

#ifndef SEED_PE_NUM
#define SEED_PE_NUM 20 // must be multiple of BLOCK_PE_NUM 
#endif

#define READ_BLOCK_NUM SEED_PE_NUM/BLOCK_PE_NUM // (3) 
#define BLOCK_NUM READ_PE_NUM*READ_BLOCK_NUM // (6)
#define o_del 6
#define e_del 1
#define o_ins 6
#define e_ins 1
#define pen_clip 5
#define w_in 100
#define zero (int11_t)(0)

/*
 *
 *  for (int h = 0; h < K_NUM; h++) {
 *      data_parse(input[h], total_size, read_ch[h]);
 *  
 *      for (int i = 0; i < READ_PE_NUM; i++) {
 *          read_proc(read_ch[h][i],
 *                    seed_param_ctrl[h][i], seed_seq[h][i], seed_param[h][i]);
 *                                                                                                
 *          for (int j = 0; j < SEED_PE_NUM; j++) {
 *              seed_proc(
 *                seed_param_ctrl[h][i][j],
 *                seed_seq[h][i][j],
 *                seed_param[h][i][j],
 *                match_ctrl1[h][i][j],
 *                match_data[h][i][j]);
 *          }
 *          for (int k =0; k < READ_BLOCK_NUM; k++) {
 *              match_assemble(
 *                &match_data[h][i][BLOCK_PE_NUM*k],
 *                &match_ctrl1[h][i][BLOCK_PE_NUM*k],
 *                block_data[h][i*READ_BLOCK_NUM + k],
 *                block_ctrl1[h][i*READ_BLOCK_NUM + k]);
 *          }
 *      }
 *                                                                                              
 *      block_assemble(block_data[h], block_ctrl1[h], results[h]);
 *      upload_results(results[h], output[h], output_size);
 *  }
 *
 */

// define channels (FIFOs) between kernels
channel int      read_seqlen_ch[K_NUM][READ_PE_NUM] __attribute__((depth(8)));
channel int      read_seq_ch[K_NUM][READ_PE_NUM] __attribute__((depth(1024)));
channel int      read_cnum_ch[K_NUM][READ_PE_NUM] __attribute__((depth(8)));
channel int      read_cseqlen_ch[K_NUM][READ_PE_NUM] __attribute__((depth(32)));
channel int      read_cseq_ch[K_NUM][READ_PE_NUM] __attribute__((depth(1024)));
channel int      read_snum_ch[K_NUM][READ_PE_NUM] __attribute__((depth(32)));
channel int      read_sparams_ch[K_NUM][READ_PE_NUM] __attribute__((depth(256)));
channel uint16_t seed_param_ch[K_NUM][READ_PE_NUM][SEED_PE_NUM] __attribute__((depth(1024)));
channel uint8_t  seed_param_ctrl_ch[K_NUM][READ_PE_NUM][SEED_PE_NUM] __attribute__((depth(8)));
channel uint8_t  seed_seq_ch[K_NUM][READ_PE_NUM][SEED_PE_NUM] __attribute__((depth(1024)));
channel int      match_ch[K_NUM][READ_PE_NUM][SEED_PE_NUM] __attribute__((depth(128)));
channel uint8_t  match_ctrl1_ch[K_NUM][READ_PE_NUM][SEED_PE_NUM] __attribute__((depth(8)));
channel int      block_ch[K_NUM][BLOCK_NUM] __attribute__((depth(512)));
channel uint8_t  block_ctrl1_ch[K_NUM][BLOCK_NUM] __attribute__((depth(128)));
channel int      results_ch[K_NUM] __attribute__((depth(1024)));

inline void increment(int *cnt, const int max) {
  *cnt = *cnt + 1;
  if (*cnt >= max) *cnt = 0;
}

inline static void nextPE(char *pe_idx, char total_pe) {
    if (*pe_idx == total_pe-1) {
        *pe_idx = 0;
    }
    else {
        *pe_idx = *pe_idx + 1;
    }
}

void sw_extend(uint8_t qs_baddr, uint11_t ts_baddr, uint8_t *seed_seq, 
    uint8_t qlen, uint11_t tlen, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins,
    int16_t max_del, int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret,
    int16_t *gscore_ret, int16_t *maxoff_ret)
{
    uint2_t k = 0;
    uint8_t aw_tmp;
    char oe_del = o_del + e_del;
    char oe_ins = o_ins + e_ins;
    int11_t max = h0; // ********** assigning int16 to int11
    char isBreak = 0;

    const char my_mat[5][5]={{1, -4, -4, -4, -1}, {-4, 1, -4, -4, -1}, {-4, -4, 1, -4, -1}, {-4, -4, -4, 1, -1}, {-1, -1, -1, -1, -1}};
    #pragma HLS ARRAY_PARTITION variable=my_mat complete dim=0
    int11_t eh_h [512];
    #pragma HLS ARRAY_MAP variable=eh_h instance=eh_arr vertical
    #pragma HLS RESOURCE variable=eh_h core=RAM_2P_BRAM
    int11_t eh_e [512];
    #pragma HLS ARRAY_MAP variable=eh_e instance=eh_arr vertical
    #pragma HLS RESOURCE variable=eh_e core=RAM_2P_BRAM

    for(uint8_t j=0;j<=qlen;j++) {
    #pragma HLS PIPELINE II=1
    #pragma HLS LOOP_TRIPCOUNT min=125 max=125
        eh_e[j]= 0;
        eh_h[j]= 0;
    }

    ext_while_loop : 
    while ((k < 2) && (!isBreak)) {
    #pragma HLS LOOP_TRIPCOUNT min=2 max=2
        int12_t max_i = -1;
        int12_t max_ie = -1;
        int12_t max_off = 0;
        short max_j = -1;
        int12_t gscore = -1;
        uint8_t end;
        int11_t tmp_eme;
        int11_t h1_init_val;
        int10_t prev;
        uint8_t aw1;

        prev = *regScore;
        aw_tmp = w_in << k;
        aw1 = aw_tmp < max_ins ? aw_tmp : max_ins;
        aw1 = aw1 < max_del ? aw1 : max_del;
        uint8_t beg = 0;
        end = qlen;
        tmp_eme = h0 - oe_ins;
        tmp_eme = (tmp_eme > 0) ? tmp_eme : (int11_t)(0);
        h1_init_val = h0 - o_del;
        target_loop : 
        for (uint11_t i = 0; i < tlen; i++) {
        #pragma HLS LOOP_TRIPCOUNT min=10 max=10
            int11_t f = 0; 
            int11_t m = 0; 
            short mj = -1;
            uint11_t ts_baddr_t;
            uint8_t backw_tmp;
            char forw_update;
            uint8_t forw_tmp;
            int11_t h1;
            uint3_t q_i;
            uint8_t j;

            ts_baddr_t = ts_baddr + i;
            q_i = ((uint3_t)seed_seq[ts_baddr_t]);

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
            #pragma ivdep array(eh_e)
            #pragma ivdep array(eh_h)
            query_loop :
            for (j = beg; j < end; ++j) {
            #pragma HLS LOOP_TRIPCOUNT min=10 max=10
            #pragma HLS pipeline II=1
            #pragma AP dependence variable=eh_e array inter false
            #pragma AP dependence variable=eh_h array inter false
                uint8_t qs_baddr_t;
                int11_t h, e, M;
                int11_t e_tmp;
                int11_t h_tmp;
                int11_t h1_reg;
                int11_t t;
                uint3_t q_j;

                qs_baddr_t = qs_baddr + j;
                q_j = ((uint3_t)seed_seq[qs_baddr_t]);
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
                if (m <= h){
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
                int10_t abs_mj_m_i;
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

void data_parse(__global int * restrict totalinp, int total_size, int k_idx) {

    int readEndIndex = 0;
    char read_pe_idx = 0;

    uint11_t chain_count = 0;
    uint11_t seed_count = 0;
    uint11_t read_count = 0;
    uint4_t stage = 0;
    int stage_cnt = 0; // counter for substate

    int ch_wr_data[READ_PE_NUM];
    uint1_t ch_wr_en[READ_PE_NUM];
    #pragma unroll
    init_wren : for (int i=0; i<READ_PE_NUM; i++) {
        ch_wr_en[i] = 0;
    }

    //int cnt = 0;
        uint11_t chain_num;
        uint11_t seed_num;
        int64_t rmax_0;
        int64_t rmax_1;
        int tmpval_hi;
        int tmpval_lo;
        int seed_len; 
        int read_seq_length_div8;
        int rseq_length_div8;
        int read_seq_length;
        int rseq_length;        
        int seed_index;
        int seed_qbeg;
        int64_t seed_rbeg;
    data_parse_loop :
    for (int idx = 0; idx < total_size; idx ++) {

        // get an input from dram
        int data_input = totalinp[idx];

        if (idx < readEndIndex) {
            if (stage == 0) { // parse read param
                read_seq_length = data_input;
                read_seq_length_div8 = (read_seq_length + 7) >> 3;

                // write read seq length
                if (read_pe_idx == 0)
                    write_channel_intel(read_seqlen_ch[k_idx][0],read_seq_length);
                else
                    write_channel_intel(read_seqlen_ch[k_idx][1],read_seq_length);
                //cnt++;
                //printf("data parse %d %d %d %d\n", idx, total_size,read_pe_idx,cnt);
                read_count = 0;
                stage = 1;
            }
            else if (stage == 1) { // parse read seqs
                if (read_count < read_seq_length_div8) {
                    // write read seq data
                    if (read_pe_idx == 0)
                        write_channel_intel(read_seq_ch[k_idx][0],data_input);
                    else
                        write_channel_intel(read_seq_ch[k_idx][1],data_input);
                    read_count += 1 ;
                }
                else {
                    chain_num = data_input;

                    // write chain_num
                    if (read_pe_idx == 0)
                        write_channel_intel(read_cnum_ch[k_idx][0],chain_num);
                    else
                        write_channel_intel(read_cnum_ch[k_idx][1],chain_num);

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
                        rseq_length = tmpval_lo;
                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 1: // read in tmpval_hi, set rmax_0
                        tmpval_hi = data_input;
                        rmax_0 = ((int64_t)tmpval_hi << 32) | tmpval_lo;
                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 2: // read in tmpval_lo
                        tmpval_lo = data_input;

                        // only need to substract lower 32bit
                        rseq_length = tmpval_lo - rseq_length;

                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 3: // read in tmpval_hi, set rmax_1
                        tmpval_hi = data_input;
                        rmax_1 = ((int64_t)tmpval_hi << 32) | tmpval_lo;

                        // only need to substract lower 32bit
                        //rseq_length = rmax_1 - rmax_0 ;
                        rseq_length_div8 = (rseq_length + 7) >> 3 ; 

                        
                        if (read_pe_idx == 0)
                            write_channel_intel(read_cseqlen_ch[k_idx][0],rseq_length);
                        else
                            write_channel_intel(read_cseqlen_ch[k_idx][1],rseq_length);

                        chain_count += 1;
                        stage = 3;

                        increment(&stage_cnt, state_cnt_max);
                        break;
                    default: ;
                }
            }
            else if (stage == 3) { // parse the chain seq
                if (read_count < rseq_length_div8){
                    if (read_pe_idx == 0)
                        write_channel_intel(read_cseq_ch[k_idx][0],data_input);
                    else
                        write_channel_intel(read_cseq_ch[k_idx][1],data_input);
                    read_count += 1 ;
                }
                else {
                    seed_num = data_input;
                    if (read_pe_idx == 0)
                        write_channel_intel(read_snum_ch[k_idx][0],seed_num);
                    else
                        write_channel_intel(read_snum_ch[k_idx][1],seed_num);
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
                    else {
                        stage = 4;
                    }
                }
            }
            else if (stage == 4){ // parse seed param
                const int state_cnt_max = 5;
                int leftRlen;
                int rightRlen;
                switch (stage_cnt) {
                    case 0:
                        seed_index = data_input;
                        ch_wr_en[read_pe_idx] = 1;
                        ch_wr_data[read_pe_idx] = seed_index;
                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 1:
                        tmpval_lo = data_input;
                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 2:
                        tmpval_hi = data_input;
                        seed_rbeg = ((int64_t)tmpval_hi << 32) | tmpval_lo;

                        leftRlen = seed_rbeg - rmax_0;
                        ch_wr_en[read_pe_idx] = 1;
                        ch_wr_data[read_pe_idx] = leftRlen;

                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 3:
                        seed_qbeg = data_input;
                        ch_wr_en[read_pe_idx] = 1;
                        ch_wr_data[read_pe_idx] = seed_qbeg;
                        increment(&stage_cnt, state_cnt_max);
                        break;
                    case 4:
                        seed_len = data_input;
                        rightRlen = rmax_1 - seed_rbeg - seed_len;

                        if (seed_len > 65535 || rightRlen > 65535) {
                            printf("seed_len %d or rightRlen overflow %d\n",seed_len,rightRlen);
                        }

                        ch_wr_en[read_pe_idx] = 1;
                        ch_wr_data[read_pe_idx] = seed_len << 16 | rightRlen;

                        seed_count += 1;
                        if (seed_count >= seed_num) {
                            seed_count = 0;
                            if (chain_count < chain_num) {
                                stage = 2;
                            }
                            else {
                                stage = 0;
                                chain_count = 0;
                            }
                        }
                        increment(&stage_cnt, state_cnt_max);
                        break;
                    default: ;
                }
            }
            if (ch_wr_en[read_pe_idx]) {
                if (read_pe_idx == 0)
                    write_channel_intel(read_sparams_ch[k_idx][0],ch_wr_data[read_pe_idx]);
                else
                    write_channel_intel(read_sparams_ch[k_idx][1],ch_wr_data[read_pe_idx]);
                ch_wr_en[read_pe_idx] = 0;
            }
        }
        else {
            // update readEndIndx
            readEndIndex = data_input; 
            
            stage = 0;

            // switch to a new PE
            if (idx) {
                nextPE(&read_pe_idx, READ_PE_NUM);
            }
        }
    }
}

__attribute__((max_global_work_dim(0)))
__attribute__((autorun))
__attribute__((num_compute_units(K_NUM,READ_PE_NUM)))
__kernel void read_proc() {

    int k_idx = get_compute_id(0);
    int read_pe_idx = get_compute_id(1);

    //int cnt = 0;
    char seed_pe_idx = 0; 

    read_proc_outer: while (1) {
        uint16_t read_seq_len;
        uint3_t read_seq[512];

        read_seq_len = read_channel_intel(read_seqlen_ch[k_idx][read_pe_idx]);
        //cnt++;
    
        read_proc_read_seq: 
        for (uint16_t i = 0; i < (uint16_t)((read_seq_len + 7) >> 3); i++) {
            int data_input = read_channel_intel(read_seq_ch[k_idx][read_pe_idx]);
            #pragma unroll
            extract_rseq: 
            for (int k = 0; k < 8; k++) {
                read_seq[i*8 + k] = (data_input & 0xF0000000) >> 28 ; 
                data_input <<= 4;
            }
        }
    
        uint11_t chain_num;
        chain_num = read_channel_intel(read_cnum_ch[k_idx][read_pe_idx]);
        read_proc_chain_outer: 
        for (uint11_t chain = 0; chain < chain_num; chain++) {
            uint11_t chain_seq_len;
            chain_seq_len = read_channel_intel(read_cseqlen_ch[k_idx][read_pe_idx]);
      
            uint3_t chain_rseqs[2048];
            int data_input = 0;
            read_proc_chain_seq: 
            for (uint11_t i = 0; i < ((chain_seq_len + 7) >> 3); i++) {
                data_input = read_channel_intel(read_cseq_ch[k_idx][read_pe_idx]);
                #pragma unroll
                extract_chrseq: 
                for (int k = 0; k < 8; k++) {
                    chain_rseqs[i*8 + k] = (data_input & 0xF0000000) >> 28 ; 
                    data_input <<= 4;
                }
            }
     
            uint11_t seed_num;
            seed_num = read_channel_intel(read_snum_ch[k_idx][read_pe_idx]);
            read_proc_seed_outer: for (uint11_t seed = 0; seed < seed_num; seed ++) {
                int read_sparams[4];
                for (int m=0; m<4; m++) {
                    read_sparams[m] = read_channel_intel(read_sparams_ch[k_idx][read_pe_idx]);
                }
                uint16_t seed_index = read_sparams[0];
                uint16_t leftRlen = read_sparams[1];
                uint16_t seed_qbeg = read_sparams[2];
                int tmp_data = read_sparams[3];
                uint16_t rightRlen = tmp_data & 0x0000FFFF;
                uint16_t seed_len = (ap_uint<32>)tmp_data >> 16;

                uint16_t leftQlen = seed_qbeg;
                uint16_t rightQlen = read_seq_len - seed_qbeg - seed_len;

                uint16_t seed_params[7];
                seed_params[0] = leftQlen;
                seed_params[1] = leftRlen;
                seed_params[2] = rightQlen;
                seed_params[3] = rightRlen;
                seed_params[4] = seed_len;
                seed_params[5] = seed_qbeg;
                seed_params[6] = seed_index;

                // send the params
                for (int i=0; i<7; i++) {
                    if (seed_pe_idx == 0)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][0], seed_params[i]);
                    else if (seed_pe_idx == 1)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][1], seed_params[i]);
                    else if (seed_pe_idx == 2)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][2], seed_params[i]);
                    else if (seed_pe_idx == 3)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][3], seed_params[i]);
                    else if (seed_pe_idx == 4)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][4], seed_params[i]);
                    else if (seed_pe_idx == 5)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][5], seed_params[i]);
                    else if (seed_pe_idx == 6)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][6], seed_params[i]);
                    else if (seed_pe_idx == 7)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][7], seed_params[i]);
                    else if (seed_pe_idx == 8)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][8], seed_params[i]);
                    else if (seed_pe_idx == 9)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][9], seed_params[i]);
                    #if SEED_PE_NUM > 19
                    else if (seed_pe_idx == 10)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][10], seed_params[i]);
                    else if (seed_pe_idx == 11)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][11], seed_params[i]);
                    else if (seed_pe_idx == 12)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][12], seed_params[i]);
                    else if (seed_pe_idx == 13)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][13], seed_params[i]);
                    else if (seed_pe_idx == 14)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][14], seed_params[i]);
                    else if (seed_pe_idx == 15)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][15], seed_params[i]);
                    else if (seed_pe_idx == 16)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][16], seed_params[i]);
                    else if (seed_pe_idx == 17)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][17], seed_params[i]);
                    else if (seed_pe_idx == 18)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][18], seed_params[i]);
                    else if (seed_pe_idx == 19)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][19], seed_params[i]);
                    #endif
                    #if SEED_PE_NUM == 30
                    else if (seed_pe_idx == 20)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][20], seed_params[i]);
                    else if (seed_pe_idx == 21)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][21], seed_params[i]);
                    else if (seed_pe_idx == 22)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][22], seed_params[i]);
                    else if (seed_pe_idx == 23)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][23], seed_params[i]);
                    else if (seed_pe_idx == 24)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][24], seed_params[i]);
                    else if (seed_pe_idx == 25)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][25], seed_params[i]);
                    else if (seed_pe_idx == 26)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][26], seed_params[i]);
                    else if (seed_pe_idx == 27)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][27], seed_params[i]);
                    else if (seed_pe_idx == 28)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][28], seed_params[i]);
                    else if (seed_pe_idx == 29)
                        write_channel_intel(seed_param_ch[k_idx][read_pe_idx][29], seed_params[i]);
                    #endif
                }
                mem_fence(CLK_CHANNEL_MEM_FENCE);
                if (seed_pe_idx == 0)
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][0], 1);
                else if (seed_pe_idx == 1)
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][1], 1);
                else if (seed_pe_idx == 2)
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][2], 1);
                else if (seed_pe_idx == 3)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][3], 1);
                else if (seed_pe_idx == 4)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][4], 1);
                else if (seed_pe_idx == 5)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][5], 1);
                else if (seed_pe_idx == 6)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][6], 1);
                else if (seed_pe_idx == 7)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][7], 1);
                else if (seed_pe_idx == 8)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][8], 1);
                else if (seed_pe_idx == 9)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][9], 1);
                #if SEED_PE_NUM > 19
                else if (seed_pe_idx == 10)
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][10], 1);
                else if (seed_pe_idx == 11)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][11], 1);
                else if (seed_pe_idx == 12)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][12], 1);
                else if (seed_pe_idx == 13)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][13], 1);
                else if (seed_pe_idx == 14)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][14], 1);
                else if (seed_pe_idx == 15)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][15], 1);
                else if (seed_pe_idx == 16)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][16], 1);
                else if (seed_pe_idx == 17)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][17], 1);
                else if (seed_pe_idx == 18)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][18], 1);
                else if (seed_pe_idx == 19)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][19], 1);
                #endif                                                           
                #if SEED_PE_NUM == 30                                            
                else if (seed_pe_idx == 20)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][20], 1);
                else if (seed_pe_idx == 21)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][21], 1);
                else if (seed_pe_idx == 22)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][22], 1);
                else if (seed_pe_idx == 23)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][23], 1);
                else if (seed_pe_idx == 24)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][24], 1);
                else if (seed_pe_idx == 25)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][25], 1);
                else if (seed_pe_idx == 26)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][26], 1);
                else if (seed_pe_idx == 27)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][27], 1);
                else if (seed_pe_idx == 28)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][28], 1);
                else if (seed_pe_idx == 29)                                      
                    write_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][29], 1);
                #endif
                //cnt++;
                //if (k_idx==0 && cnt > 60) printf("read_proc %d %d %d\n",k_idx, read_pe_idx, cnt);

                uint8_t leftQlen_u8 = leftQlen;
                uint11_t leftRlen_u11 = leftRlen;
                uint8_t rightQlen_u8 = rightQlen;
                uint11_t rightRlen_u11 = rightRlen;
                // send the seqs
                uint12_t seed_seq_length = leftQlen_u8 + leftRlen_u11 + rightQlen_u8 + rightRlen_u11;

                uint8_t wr_data;
                read_proc_seed_seq: 
                for (uint12_t i = 0; i < seed_seq_length; i++) {
                    if ( i< leftQlen_u8){
                        wr_data = read_seq[leftQlen_u8-1-i];
                    }
                    else if (i < leftQlen_u8 + rightQlen_u8){
                        // seed_qbeg + seed_len +i - leftQlen
                        wr_data = read_seq[seed_len+i];
                    }
                    else if (i < leftQlen_u8 + leftRlen_u11 + rightQlen_u8){
                        // leftRlen -1 -(i- leftQlen - rightQlen) 
                        wr_data = chain_rseqs[leftRlen_u11 +leftQlen_u8 +rightQlen_u8-1-i];
                    }
                    else {
                        // seed_rbeg + seed_len - rmax0 + i - leftQlen - rightQlen - leftRlen
                        wr_data = chain_rseqs[2*seed_len - read_seq_len + i];
                    }
                    if (seed_pe_idx == 0)
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][0], wr_data);
                    else if (seed_pe_idx == 1)
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][1], wr_data);
                    else if (seed_pe_idx == 2)
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][2], wr_data);
                    else if (seed_pe_idx == 3)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][3], wr_data);
                    else if (seed_pe_idx == 4)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][4], wr_data);
                    else if (seed_pe_idx == 5)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][5], wr_data);
                    else if (seed_pe_idx == 6)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][6], wr_data);
                    else if (seed_pe_idx == 7)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][7], wr_data);
                    else if (seed_pe_idx == 8)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][8], wr_data);
                    else if (seed_pe_idx == 9)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][9], wr_data);
                    #if SEED_PE_NUM > 19
                    else if (seed_pe_idx == 10)
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][10], wr_data);
                    else if (seed_pe_idx == 11)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][11], wr_data);
                    else if (seed_pe_idx == 12)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][12], wr_data);
                    else if (seed_pe_idx == 13)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][13], wr_data);
                    else if (seed_pe_idx == 14)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][14], wr_data);
                    else if (seed_pe_idx == 15)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][15], wr_data);
                    else if (seed_pe_idx == 16)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][16], wr_data);
                    else if (seed_pe_idx == 17)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][17], wr_data);
                    else if (seed_pe_idx == 18)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][18], wr_data);
                    else if (seed_pe_idx == 19)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][19], wr_data);
                    #endif                                                             
                    #if SEED_PE_NUM == 30                                              
                    else if (seed_pe_idx == 20)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][20], wr_data);
                    else if (seed_pe_idx == 21)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][21], wr_data);
                    else if (seed_pe_idx == 22)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][22], wr_data);
                    else if (seed_pe_idx == 23)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][23], wr_data);
                    else if (seed_pe_idx == 24)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][24], wr_data);
                    else if (seed_pe_idx == 25)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][25], wr_data);
                    else if (seed_pe_idx == 26)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][26], wr_data);
                    else if (seed_pe_idx == 27)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][27], wr_data);
                    else if (seed_pe_idx == 28)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][28], wr_data);
                    else if (seed_pe_idx == 29)                               
                        write_channel_intel(seed_seq_ch[k_idx][read_pe_idx][29], wr_data);
                    #endif                                                           
                }
                nextPE(&seed_pe_idx, SEED_PE_NUM);
            }
        }
    }
}

__attribute__((max_global_work_dim(0)))
__attribute__((autorun))
__attribute__((num_compute_units(K_NUM,READ_PE_NUM,SEED_PE_NUM)))
__kernel void seed_proc()
{

    int k_idx = get_compute_id(0);
    int read_pe_idx = get_compute_id(1);
    int seed_pe_idx = get_compute_id(2);

    //int cnt = 0;

    seed_proc_outer:
    while (1) {
        uint8_t flag;
        flag = read_channel_intel(seed_param_ctrl_ch[k_idx][read_pe_idx][seed_pe_idx]);
        //cnt++;
        //if (k_idx==0 && cnt>5) printf("seed_proc %d %d %d %d\n",k_idx, read_pe_idx, seed_pe_idx, cnt);

        uint16_t seed_params[7];
        seed_proc_params: 
        for (int i=0; i<7; i++) {
            seed_params[i] = read_channel_intel(seed_param_ch[k_idx][read_pe_idx][seed_pe_idx]);
        }

        uint16_t leftQlen;
        uint16_t leftRlen;
        uint16_t rightQlen;
        uint16_t rightRlen;
        uint16_t seed_len;
        uint16_t seed_qbeg;
        uint16_t seed_index;
        uint8_t seq_mem[2048];
        
        leftQlen = seed_params[0];
        leftRlen = seed_params[1];
        rightQlen = seed_params[2];
        rightRlen = seed_params[3];
        seed_len = seed_params[4];
        seed_qbeg = seed_params[5];
        seed_index = seed_params[6];

        seed_proc_seq: 
        for (uint16_t k = 0; k<(leftQlen + leftRlen + rightQlen + rightRlen); k++) {
            seq_mem[k] = read_channel_intel(seed_seq_ch[k_idx][read_pe_idx][seed_pe_idx]);
        }

        uint8_t qlen[2];
        uint11_t tlen[2];
        int16_t h0;
        int16_t regScore;
        int16_t max_ins[2];
        int16_t max_del[2];
        int16_t qle;
        int16_t tle;
        int16_t gtle;
        int16_t gscore;
        int16_t maxoff;
        int16_t qBeg;
        int16_t rBeg;
        int16_t qEnd;
        int16_t rEnd;
        int16_t score;
        int16_t trueScore;
        int16_t aw[2];

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

        uint8_t qs_baddr;
        uint11_t ts_baddr;
        
        qs_baddr = 0;
        ts_baddr = qlen[0] + qlen[1];

        // ****** CAN THIS BE UNROLLED? OR IS THERE A DEP ON regScore?
        seed_proc_sw_lr: 
        for (uint2_t i=0; i<2; i++) {
            int16_t sc0;
            int16_t h0_arr[2];
            
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
            qs_baddr += qlen[i];
            ts_baddr += tlen[i];
        }
        
        int16_t width;
        if (aw[0] > aw[1]) width = aw[0];
        else width = aw[1];
        // get the results
        int wr_data[5];
        wr_data[0] = seed_index;
        wr_data[1] = (qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000);
        wr_data[2] = (rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000);
        wr_data[3] = (score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000);
        wr_data[4] = width & 0xFFFF;
        /*printf("%d %d %d %d %d\n",seed_index,
            (qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000),
            (rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000),
            (score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000),
            width & 0xFFFF);*/
        for (int i=0; i<5; i++) {
            write_channel_intel(match_ch[k_idx][read_pe_idx][seed_pe_idx], wr_data[i]);
        }
        // write to control ch to indicate data is avail in data ch
        mem_fence(CLK_CHANNEL_MEM_FENCE);
        write_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][seed_pe_idx], 1);
        //if (k_idx == 0 && cnt>5) printf("seed proc %d %d %d %d\n",k_idx,read_pe_idx,seed_pe_idx,cnt);
        //cnt++;
    }
}

__attribute__((max_global_work_dim(0)))
__attribute__((autorun))
__attribute__((num_compute_units(K_NUM,READ_PE_NUM,READ_BLOCK_NUM)))
__kernel void match_assemble() 
{
    int k_idx = get_compute_id(0);
    int read_pe_idx = get_compute_id(1);
    int read_blk_idx = get_compute_id(2);

    //int cnt=0;
    
    match_assemble_outer: 
    while (1) {
        match_assemble_inner: 
        for (int i=0; i<BLOCK_PE_NUM; i++) {
            bool valid0, valid1;
            uint8_t flag;
            //flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+i]);
            
            if (i == 0)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][read_blk_idx*BLOCK_PE_NUM]);
            else if (i == 1)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+1]);
            else if (i == 2)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+2]);
            else if (i == 3)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+3]);
            else if (i == 4)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+4]);
            else if (i == 5)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+5]);
            else if (i == 6)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+6]);
            else if (i == 7)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+7]);
            else if (i == 8)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+8]);
            else if (i == 9)
                flag = read_channel_intel(match_ctrl1_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+9]);
            
            for (int k = 0; k < 5; k++) {
                int one_match_buffer[5];
                //one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+i]);
                
                if (i == 0)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][read_blk_idx*BLOCK_PE_NUM]);
                else if (i == 1)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+1]);
                else if (i == 2)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+2]);
                else if (i == 3)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+3]);
                else if (i == 4)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+4]);
                else if (i == 5)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+5]);
                else if (i == 6)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+6]);
                else if (i == 7)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+7]);
                else if (i == 8)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+8]);
                else if (i == 9)
                    one_match_buffer[k] = read_channel_intel(match_ch[k_idx][read_pe_idx][(read_blk_idx*BLOCK_PE_NUM)+9]);
                
                write_channel_intel(block_ch[k_idx][(read_pe_idx*READ_BLOCK_NUM)+read_blk_idx],one_match_buffer[k]);
            }
            mem_fence(CLK_CHANNEL_MEM_FENCE);
            write_channel_intel(block_ctrl1_ch[k_idx][(read_pe_idx*READ_BLOCK_NUM)+read_blk_idx],1);
            //cnt++;
            //if (k_idx == 0 && ((read_pe_idx == 0 && cnt>61) || (read_pe_idx==1 && cnt>37))) 
                //printf("match asm %d %d %d %d \n",k_idx, read_pe_idx,i,cnt);
        }
    }
}

__attribute__((max_global_work_dim(0)))
__attribute__((autorun))
__attribute__((num_compute_units(K_NUM)))
__kernel void block_assemble() {
    int k_idx = get_compute_id(0);

    /*        
    int cnt = 0;
    int cnt0 = 0;
    int cnt1 = 0;
    */
    
    block_assemble_outer: 
    while (1) {
        block_assemble_inner: 
        for (int i=0; i<BLOCK_NUM; i++) {
            bool valid0;
            uint8_t flag;

            //flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][i],&valid0);
            if (i == 0)
                flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][0],&valid0);
            else if (i == 1)
                flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][1],&valid0);
        #if BLOCK_NUM > 2
            else if (i == 2)
                flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][2],&valid0);
            else if (i == 3)                                        
                flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][3],&valid0);
        #endif
        #if BLOCK_NUM >4 
            else if (i == 4)                                        
                flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][4],&valid0);
            else if (i == 5)                                        
                flag = read_channel_nb_intel(block_ctrl1_ch[k_idx][5],&valid0);
        #endif
            
            if (valid0) {
                for (int k = 0; k < 5; k++) {
                    int one_block_buffer[5];
                    //one_block_buffer[k] = read_channel_intel(block_ch[k_idx][i]);
                    if (i == 0)
                        one_block_buffer[k] = read_channel_intel(block_ch[k_idx][0]);
                    else if (i == 1)
                        one_block_buffer[k] = read_channel_intel(block_ch[k_idx][1]);
                #if BLOCK_NUM > 2
                    else if (i == 2)
                        one_block_buffer[k] = read_channel_intel(block_ch[k_idx][2]);
                    else if (i == 3)
                        one_block_buffer[k] = read_channel_intel(block_ch[k_idx][3]);
                #endif
                #if BLOCK_NUM >4 
                    else if (i == 4)
                        one_block_buffer[k] = read_channel_intel(block_ch[k_idx][4]);
                    else if (i == 5)
                        one_block_buffer[k] = read_channel_intel(block_ch[k_idx][5]);
                #endif

                    write_channel_intel(results_ch[k_idx],one_block_buffer[k]);
                }
                /*
                cnt++;
                if (i==0) cnt0++;
                if (i==1) cnt1++;
                if (k_idx == 0 && cnt>75) printf("block asm %d %d %d %d %d \n",k_idx, i, cnt0, cnt1, cnt);
                */
            }
        }
    }
}

void upload_results(__global int *output, int out_size, int k_idx) {
    for (int i=0; i<out_size; i++) {
        int tmp;
        tmp = read_channel_intel(results_ch[k_idx]);
        output[i] = tmp;
        //if (i<20) printf("%d \n",output[i]);
        //if (k_idx==0 && i>445 && i<510) printf("%d %d %d\n",k_idx,i,tmp);
    }
    //printf("upload results done %d\n",k_idx);
}


#define SPAWN_KERNEL(n) \
    __attribute__((max_global_work_dim(0))) \
    __kernel void data_parse ## n (__global int * restrict totalinp, int total_size) { \
        data_parse(totalinp, total_size, n); \
    } \
\
    __attribute__((max_global_work_dim(0))) \
    __kernel void upload_results ## n (__global int *output, int out_size) {\
        upload_results(output, out_size, n); \
    } 


    SPAWN_KERNEL(0)
#if K_NUM > 1
    SPAWN_KERNEL(1)
#endif
#if K_NUM > 2 
    SPAWN_KERNEL(2)
#endif
#if K_NUM > 3 
    SPAWN_KERNEL(3)
#endif

