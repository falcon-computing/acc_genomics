#ifdef CPP
#include <ap_int.h>
#include <ap_utils.h>
#include <hls_stream.h>
#endif

#define CAT_NX(A, B) A ## B
#define CAT(A, B) CAT_NX(A, B)

#define BWT_PRIMARY (long)2637708960
#define CNT_TABLE_CONTENT { 4, 259, 65539, 16777219, 259, 514, 65794, 16777474, 65539, 65794, \
131074, 16842754, 16777219, 16777474, 16842754, 33554434, 259, 514, 65794, 16777474, \
514, 769, 66049, 16777729, 65794, 66049, 131329, 16843009, 16777474, 16777729, \
16843009, 33554689, 65539, 65794, 131074, 16842754, 65794, 66049, 131329, 16843009, \
131074, 131329, 196609, 16908289, 16842754, 16843009, 16908289, 33619969, 16777219, 16777474, \
16842754, 33554434, 16777474, 16777729, 16843009, 33554689, 16842754, 16843009, 16908289, 33619969, \
33554434, 33554689, 33619969, 50331649, 259, 514, 65794, 16777474, 514, 769, \
66049, 16777729, 65794, 66049, 131329, 16843009, 16777474, 16777729, 16843009, 33554689, \
514, 769, 66049, 16777729, 769, 1024, 66304, 16777984, 66049, 66304, \
131584, 16843264, 16777729, 16777984, 16843264, 33554944, 65794, 66049, 131329, 16843009, \
66049, 66304, 131584, 16843264, 131329, 131584, 196864, 16908544, 16843009, 16843264, \
16908544, 33620224, 16777474, 16777729, 16843009, 33554689, 16777729, 16777984, 16843264, 33554944, \
16843009, 16843264, 16908544, 33620224, 33554689, 33554944, 33620224, 50331904, 65539, 65794, \
131074, 16842754, 65794, 66049, 131329, 16843009, 131074, 131329, 196609, 16908289, \
16842754, 16843009, 16908289, 33619969, 65794, 66049, 131329, 16843009, 66049, 66304, \
131584, 16843264, 131329, 131584, 196864, 16908544, 16843009, 16843264, 16908544, 33620224, \
131074, 131329, 196609, 16908289, 131329, 131584, 196864, 16908544, 196609, 196864, \
262144, 16973824, 16908289, 16908544, 16973824, 33685504, 16842754, 16843009, 16908289, 33619969, \
16843009, 16843264, 16908544, 33620224, 16908289, 16908544, 16973824, 33685504, 33619969, 33620224, \
33685504, 50397184, 16777219, 16777474, 16842754, 33554434, 16777474, 16777729, 16843009, 33554689, \
16842754, 16843009, 16908289, 33619969, 33554434, 33554689, 33619969, 50331649, 16777474, 16777729, \
16843009, 33554689, 16777729, 16777984, 16843264, 33554944, 16843009, 16843264, 16908544, 33620224, \
33554689, 33554944, 33620224, 50331904, 16842754, 16843009, 16908289, 33619969, 16843009, 16843264, \
16908544, 33620224, 16908289, 16908544, 16973824, 33685504, 33619969, 33620224, 33685504, 50397184, \
33554434, 33554689, 33619969, 50331649, 33554689, 33554944, 33620224, 50331904, 33619969, 33620224, \
33685504, 50397184, 50331649, 50331904, 50397184, 67108864 } 

#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#define MIN_SEED_LEN 19
#define SPLIT_WIDTH 10
#define SPLIT_LEN 28 //(int)(opt->min_seed_len * opt->split_factor + .499)
#define MAX_MEM_INTV 20
#define MAX_INTV_ALLOC 256
#define SEQ_LENGTH 256
#define MAX_BATCH_SIZE 1024
#define MAX_TILE_SIZE 16

#ifdef CL
#define BATCH_SIZE 16
#endif

#ifdef CL
typedef ulong  bwtint_t;
typedef ulong4 smem_t;
typedef uint   addr_bundle;
typedef uint16 data_bundle;
typedef uchar  byte_t;
typedef uint   uint_t;
#else 
typedef ap_uint<64>  bwtint_t;
typedef ap_uint<256> smem_t;
typedef ap_uint<32>  addr_bundle;
typedef ap_uint<512> data_bundle;
typedef unsigned char   byte_t;
typedef ap_uint<32>  uint_t;
#endif

#ifdef CL
#define addr_channel pipe addr_bundle
#define data_channel pipe data_bundle
#define ctrl_channel pipe byte_t 
#else
#define addr_channel hls::stream<addr_bundle> 
#define data_channel hls::stream<data_bundle>
#define ctrl_channel hls::stream<byte_t> 
#endif

typedef struct {
    bwtint_t x0;
    bwtint_t x1;
    bwtint_t x2;
    bwtint_t info;
} bwtintv_t_fpga;


#ifdef CL
addr_channel  fe_addr_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  fe_data_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel  be_addr_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  be_data_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lfe_addr_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel lfe_data_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lbe_addr_channel0 __attribute__((xcl_reqd_pipe_depth(1024)));
data_channel lbe_data_channel0 __attribute__((xcl_reqd_pipe_depth(1024)));
addr_channel afe_addr_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel afe_data_channel0 __attribute__((xcl_reqd_pipe_depth(256)));
ctrl_channel     ctrl_channel0 __attribute__((xcl_reqd_pipe_depth(16)));

addr_channel  fe_addr_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  fe_data_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel  be_addr_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  be_data_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lfe_addr_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel lfe_data_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lbe_addr_channel1 __attribute__((xcl_reqd_pipe_depth(1024)));
data_channel lbe_data_channel1 __attribute__((xcl_reqd_pipe_depth(1024)));
addr_channel afe_addr_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel afe_data_channel1 __attribute__((xcl_reqd_pipe_depth(256)));
ctrl_channel     ctrl_channel1 __attribute__((xcl_reqd_pipe_depth(16)));

addr_channel  fe_addr_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  fe_data_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel  be_addr_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  be_data_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lfe_addr_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel lfe_data_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lbe_addr_channel2 __attribute__((xcl_reqd_pipe_depth(1024)));
data_channel lbe_data_channel2 __attribute__((xcl_reqd_pipe_depth(1024)));
addr_channel afe_addr_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel afe_data_channel2 __attribute__((xcl_reqd_pipe_depth(256)));
ctrl_channel     ctrl_channel2 __attribute__((xcl_reqd_pipe_depth(16)));

addr_channel  fe_addr_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  fe_data_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel  be_addr_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel  be_data_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lfe_addr_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel lfe_data_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
addr_channel lbe_addr_channel3 __attribute__((xcl_reqd_pipe_depth(1024)));
data_channel lbe_data_channel3 __attribute__((xcl_reqd_pipe_depth(1024)));
addr_channel afe_addr_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
data_channel afe_data_channel3 __attribute__((xcl_reqd_pipe_depth(256)));
ctrl_channel     ctrl_channel3 __attribute__((xcl_reqd_pipe_depth(16)));
#endif

inline void write_data_ch_blk(data_channel* ch, data_bundle* val) {
#ifdef CPP
#pragma HLS inline
    ch->write(*val);
#else
    write_pipe_block(*ch, val);
#endif
}

inline void read_data_ch_blk(data_channel* ch, data_bundle* val) {
#ifdef CPP
#pragma HLS inline
    *val = ch->read();
#else
    read_pipe_block(*ch, val);
#endif
}

inline void write_addr_ch_blk(addr_channel* ch, addr_bundle* val) {
#ifdef CPP
#pragma HLS inline
    ch->write(*val);
#else
    write_pipe_block(*ch, val);
#endif
}

inline bool read_addr_ch_nblk(addr_channel* ch, addr_bundle* val) {
#ifdef CPP
#pragma HLS inline
    return ch->read_nb(*val);
#else
    int ret = read_pipe(*ch, val);
    return (ret == 0);
#endif
}

inline bool write_ctrl_ch_nblk(ctrl_channel* ch, byte_t* val) {
#ifdef CPP
#pragma HLS inline
    return ch->write_nb(*val);
#else
    int ret = write_pipe(*ch, val);
    return (ret == 0);
#endif
}

inline bool read_ctrl_ch_nblk(ctrl_channel* ch, byte_t* val) {
#ifdef CPP
#pragma HLS inline
    return ch->read_nb(*val);
#else
    int ret = read_pipe(*ch, val);
    return (ret == 0);
#endif
}

uint_t inline vector_uin16_sel(data_bundle v, byte_t idx){
#ifdef CPP
#pragma HLS inline
    return v(32 * idx + 31, 32 * idx);
#else
    uint_t ret;
    if(idx == 0)
        ret = v.s0;
    else if(idx == 1)
        ret = v.s1;
    else if(idx == 2)
        ret = v.s2;
    else if(idx == 3)
        ret = v.s3;
    else if(idx == 4)
        ret = v.s4;
    else if(idx == 5)
        ret = v.s5;
    else if(idx == 6)
        ret = v.s6;
    else if(idx == 7)
        ret = v.s7;
    else if(idx == 8)
        ret = v.s8;
    else if(idx == 9)
        ret = v.s9;
    else if(idx == 10)
        ret = v.sa;
    else if(idx == 11)
        ret = v.sb;
    else if(idx == 12)
        ret = v.sc;
    else if(idx == 13)
        ret = v.sd;
    else if(idx == 14)
        ret = v.se;
    else if(idx == 15)
        ret = v.sf;
    else 
        ret = 0;
    return ret;
#endif
}

bwtint_t cnt_table_lookup(uint_t tmp) {
#pragma HLS inline off
    uint_t cnt_table0[256] = CNT_TABLE_CONTENT;
    uint_t cnt_table1[256] = CNT_TABLE_CONTENT;
    return cnt_table0[tmp & 0xff] + cnt_table0[(tmp >> 8) & 0xff]\
        + cnt_table1[(tmp >> 16) & 0xff] + cnt_table1[tmp >> 24];
}

void inline intv_init(bwtintv_t_fpga* intv, bwtint_t L2[5], byte_t base) {
#pragma HLS inline
    intv->x0 = L2[base] + 1;
    intv->x2 = L2[base + 1] - L2[base];
    intv->x1 = L2[3 - base] + 1;
    intv->info = 0;
}

void inline bwt_occ4_fpga(data_bundle p_buf_long, bwtint_t k, bwtint_t cnt[4]) {
#pragma HLS inline
	bwtint_t x = 0;
    //byte_t end_addr = ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4)) & 0xf;
    byte_t end_addr = (k >> 4) & 0x7;
    bwtint_t lookup_results_s0[8];
#pragma HLS array_partition variable=lookup_results_s0 complete
    bwtint_t lookup_results_s1[4];
#pragma HLS array_partition variable=lookup_results_s1 complete
    bwtint_t lookup_results_s2[2];
#pragma HLS array_partition variable=lookup_results_s2 complete
    uint_t mask = (~((1U << ((~k & 15) << 1)) - 1));

    for (byte_t i = 0; i < 8; i++) {
#pragma HLS unroll
        uint_t tmp = vector_uin16_sel(p_buf_long, i + 8);
        if (i == end_addr) {
            tmp = tmp & mask;  
        }
        if(i <= end_addr){
            lookup_results_s0[i] = cnt_table_lookup(tmp);
        }
        else
            lookup_results_s0[i] = 0;
    }

    for (byte_t i = 0; i < 4; i++) {
#pragma HLS unroll
        lookup_results_s1[i] = lookup_results_s0[2 * i] + lookup_results_s0[2 * i + 1];
    }

    for (byte_t i = 0; i < 2; i++) {
#pragma HLS unroll
        lookup_results_s2[i] = lookup_results_s1[2 * i] + lookup_results_s1[2 * i + 1];
    }

    x = lookup_results_s2[0] + lookup_results_s2[1]; 
    x -= (~k & 15);

#ifdef CL
    cnt[0] = p_buf_long.s0 + ((bwtint_t)(p_buf_long.s1) << 32) + (x & 0xff);
    cnt[1] = p_buf_long.s2 + ((bwtint_t)(p_buf_long.s3) << 32) + ((x >> 8) & 0xff);
    cnt[2] = p_buf_long.s4 + ((bwtint_t)(p_buf_long.s5) << 32) + ((x >> 16) & 0xff);
    cnt[3] = p_buf_long.s6 + ((bwtint_t)(p_buf_long.s7) << 32) + (x >> 24);
#else
    cnt[0] = p_buf_long(63, 0) + (x & 0xff);
    cnt[1] = p_buf_long(127, 64) + ((x >> 8) & 0xff);
    cnt[2] = p_buf_long(191, 128) + ((x >> 16) & 0xff);
    cnt[3] = p_buf_long(255, 192) + (x >> 24);
#endif
}

inline void bwt_extend_fpga_addr(bwtintv_t_fpga *ik, bwtint_t bwt_primary, \
        int is_back, bwtint_t* _k, bwtint_t* _l) {
#pragma HLS inline
    bwtint_t k, l;
    if (!is_back) {
        k = ik->x1 - 1;
        l = ik->x1 - 1 + ik->x2;
    }
    else {
        k = ik->x0 - 1;
        l = ik->x0 - 1 + ik->x2;
    }
    *_k = k - (k >= bwt_primary);
    *_l = l - (l >= bwt_primary);
}

inline void bwt_extend_fpga_data(bwtint_t _k, bwtint_t _l, data_bundle k_data, \
        data_bundle l_data, bwtint_t bwt_primary, bwtint_t L2[5], bwtintv_t_fpga* ik, \
        bwtintv_t_fpga* final_ok, byte_t base, int is_back) {
#pragma HLS inline 
    bwtintv_t_fpga ok[4];
#pragma HLS array_partition variable=ok complete
	bwtint_t tk[4];
#pragma HLS array_partition variable=tk complete
    bwtint_t tl[4];
#pragma HLS array_partition variable=tl complete
    
    bwt_occ4_fpga(k_data, _k, tk);
    bwt_occ4_fpga(l_data, _l, tl);

    if (!is_back) {
        for (int i = 0; i < 4; ++i) {
#pragma HLS unroll
            ok[i].x1 = L2[i] + 1 + tk[i];
            ok[i].x2 = tl[i] - tk[i];
        }
        ok[3].x0 = ik->x0 + ((ik->x1 <= bwt_primary) \
                && (ik->x1 + ik->x2 - 1 >= bwt_primary));
        ok[2].x0 = ok[3].x0 + ok[3].x2;
        ok[1].x0 = ok[2].x0 + ok[2].x2;
        ok[0].x0 = ok[1].x0 + ok[1].x2;
    }
    else {
	    for (int i = 0; i < 4; ++i) {
#pragma HLS unroll
		    ok[i].x0 = L2[i] + 1 + tk[i];
		    ok[i].x2 = tl[i] - tk[i];
	    }
        ok[3].x1 = ik->x1 + ((ik->x0 <= bwt_primary) \
                && (ik->x0 + ik->x2 - 1 >= bwt_primary));
        ok[2].x1 = ok[3].x1 + ok[3].x2;
        ok[1].x1 = ok[2].x1 + ok[2].x2;
        ok[0].x1 = ok[1].x1 + ok[1].x2;
    }
    *final_ok = ok[base];
}


#ifdef CL
inline void addr_to_data(const global data_bundle *bwt, bwtint_t bwt_size, uint_t addr, data_bundle* data) {
#else
inline void addr_to_data(data_bundle *bwt, bwtint_t bwt_size, uint_t addr, data_bundle* data) {
#endif
#pragma HLS inline
    bwtint_t shift_addr = addr;
    if (shift_addr < bwt_size) 
        *data = bwt[shift_addr];
    else{
        *data = 0;
    }
}

inline void enqueue_addr(bwtint_t bwt_primary, bwtintv_t_fpga *ik, int is_back, addr_channel* addrs_pipe, \
        bwtint_t* k_addr, bwtint_t* l_addr) {
#pragma HLS inline 
    uint_t addr0, addr1;
    bwt_extend_fpga_addr(ik, bwt_primary, is_back, k_addr, l_addr);
    addr0 = (uint_t)((*k_addr) >> OCC_INTV_SHIFT);
    addr1 = (uint_t)((*l_addr) >> OCC_INTV_SHIFT);
    write_addr_ch_blk(addrs_pipe, &addr0);
#ifdef CPP
    ap_wait();
#endif
    if (addr1 != addr0)
        write_addr_ch_blk(addrs_pipe, &addr1);
}

inline void dequeue_data(bwtintv_t_fpga *ik, bwtint_t k_addr, bwtint_t l_addr, bwtint_t bwt_primary, bwtint_t L2[5], \
        bwtintv_t_fpga* final_ok, byte_t base, int is_back, data_channel* datas_pipe) {
#pragma HLS inline 
    uint_t addr0, addr1;
    data_bundle k_data, l_data;
    addr0 = (uint_t)(k_addr >> OCC_INTV_SHIFT);
    addr1 = (uint_t)(l_addr >> OCC_INTV_SHIFT);

    read_data_ch_blk(datas_pipe, &k_data);
#ifdef CPP
    ap_wait();
#endif
    if (addr1 != addr0) {
        read_data_ch_blk(datas_pipe, &l_data);
    }
    else
        l_data = k_data;
    bwt_extend_fpga_data(k_addr, l_addr, k_data, l_data, bwt_primary, L2, ik, final_ok, base, is_back);
}

inline void bwt_extend_fpga_opt(bwtint_t bwt_primary, bwtint_t L2[5], \
        bwtintv_t_fpga *ik, bwtintv_t_fpga* final_ok, byte_t base, \
        int is_back, addr_channel* addrs_pipe, data_channel* datas_pipe) {
#pragma HLS inline 
    bwtint_t k_addr, l_addr;
/*
    bwtint_t addrs;
    data_bundle k_data, l_data;
    bwt_extend_fpga_addr(ik, bwt_primary, is_back, &k_addr, &l_addr);
    addrs = (k_addr >> OCC_INTV_SHIFT) | ((l_addr >> OCC_INTV_SHIFT) << 32);
    write_addr_ch_blk(addrs_pipe, &addrs);*/
    enqueue_addr(bwt_primary, ik, is_back, addrs_pipe, &k_addr, &l_addr);
#ifdef CPP
    ap_wait();
#endif
    dequeue_data(ik, k_addr, l_addr, bwt_primary, L2, final_ok, base, is_back, datas_pipe);
    /*
    read_data_ch_blk(datas_pipe, &k_data);
#ifdef CPP
    ap_wait();
#endif
    read_data_ch_blk(datas_pipe, &l_data);
    bwt_extend_fpga_data(k_addr, l_addr, k_data, l_data, bwt_primary, L2, ik, final_ok, base, is_back);*/
}





#ifdef CL
void input_dup(const global byte_t *seq, const global byte_t* seq_len, int batch_id, \
        byte_t* seq_fe,  byte_t* seq_len_fe,\
        byte_t* seq_afe, byte_t* seq_len_afe, byte_t tile_size) {
#else
void input_dup(byte_t *seq, byte_t* seq_len, int batch_id, \
        hls::stream<byte_t>* seq_fe,  hls::stream<byte_t>* seq_len_fe,\
        hls::stream<byte_t>* seq_afe, hls::stream<byte_t>* seq_len_afe, byte_t tile_size) {
#endif
#pragma HLS inline
    int offset = batch_id * SEQ_LENGTH;
input_dup_seq_len:
    for (int i = 0; i < tile_size; i++) {
        byte_t seq_len_local = seq_len[batch_id + i];
#ifdef CL
        seq_len_fe[i]  = seq_len_local;
        seq_len_afe[i] = seq_len_local;
#else
        seq_len_fe->write(seq_len_local);
        seq_len_afe->write(seq_len_local);
#endif
    }

#ifdef CPP
    ap_wait();
#endif

input_dup_seq:
    for (int i = 0; i < tile_size * SEQ_LENGTH; i++) {
        byte_t cur_base = seq[i + offset];
#ifdef CL
        seq_fe[i] = cur_base;
        seq_afe[i] = cur_base;
#else
        seq_fe->write(cur_base);  
        seq_afe->write(cur_base);  
#endif
    }
}

#ifdef CL
void forward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        byte_t* seq, byte_t* seq_len, int* forward_intv_n, \
        bwtintv_t_fpga* forward, addr_channel* addrs_pipe, data_channel* datas_pipe, \
        byte_t* seq_be, byte_t* seq_len_be, byte_t tile_size) {
#else
void forward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        hls::stream<byte_t>* seq, hls::stream<byte_t>* seq_len, hls::stream<int>* forward_intv_n, \
        hls::stream<bwtintv_t_fpga>* forward, addr_channel* addrs_pipe, data_channel* datas_pipe, \
        hls::stream<byte_t>* seq_be, hls::stream<byte_t>* seq_len_be, byte_t tile_size) {
#endif
#pragma HLS inline
    int intv_n = 0;
    bwtintv_t_fpga ik, ok;
    int min_intv = 1;
    int i = 0;
    bwtintv_t_fpga forward_local[MAX_TILE_SIZE][MAX_INTV_ALLOC];
#pragma HLS resource variable=forward_local core=xpm_memory uram
    byte_t seq_local[MAX_TILE_SIZE][SEQ_LENGTH];
#pragma HLS resource variable=seq_local core=RAM_2P_BRAM
    byte_t seq_len_local[MAX_TILE_SIZE];
#pragma HLS resource variable=seq_len_local core=RAM_2P_LUTRAM
#ifdef CL
seq_len_fe_latch:
    for (int idx = 0; idx < tile_size; idx++) {
        seq_len_local[idx] = seq_len[idx];
    }
    int addr = 0;
seq_fe_latch:
    for (int idx = 0; idx < tile_size; idx++) {
        for (int jdx = 0; jdx < SEQ_LENGTH; jdx++) {
            seq_local[idx][jdx] = seq[addr++];
        }
    }
#else
    for (int idx = 0; idx < tile_size; idx++) {
        seq_len_local[idx] = seq_len->read();
    }
    ap_wait();
seq_fe_latch:
    for (int idx = 0; idx < tile_size; idx++) {
        for (int jdx = 0; jdx < SEQ_LENGTH; jdx++) {
            seq_local[idx][jdx] = seq->read();
        }
    }
#endif
    
    short intv_n_list[MAX_TILE_SIZE];
#pragma HLS resource variable=intv_n_list core=RAM_2P_LUTRAM
    bool raw_intv_list[MAX_TILE_SIZE];
#pragma HLS resource variable=raw_intv_list core=RAM_2P_LUTRAM
    bwtintv_t_fpga ik_list[MAX_TILE_SIZE];
#pragma HLS resource variable=ik_list core=RAM_2P_LUTRAM
    short start_pos_list[MAX_TILE_SIZE];
#pragma HLS resource variable=start_pos_list core=RAM_2P_LUTRAM
fe_init_raw_intv:
    for (i = 0; i < MAX_TILE_SIZE; i++) {
        raw_intv_list[i] = true;
        start_pos_list[i] = -1;
        intv_n_list[i] = 0;
    }
    int max_seq_len = 0;
fe_find_max_seq_len:
    for (int i = 0; i < tile_size; i++) {
        if(seq_len_local[i] > max_seq_len)
            max_seq_len = seq_len_local[i];
    }

forward_ext:
    for (i = 0; i < max_seq_len; ++i) {
        bwtint_t k_addr_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=k_addr_buf core=RAM_2P_LUTRAM
        bwtint_t l_addr_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=l_addr_buf core=RAM_2P_LUTRAM
fe_push_addr:
        for (int j = 0; j < tile_size; j++) {
            if (i < seq_len_local[j] && !raw_intv_list[j] && seq_local[j][i] <= 3) {
                enqueue_addr(bwt_primary, &ik_list[j], 0, addrs_pipe, &k_addr_buf[j], &l_addr_buf[j]); 
            }
        }
#ifdef CPP
        ap_wait();
#endif
fe_pop_data:
        for (int j = 0; j < tile_size; j++) {
            if (i < seq_len_local[j]) {
                if (seq_local[j][i] > 3) {
                    if (!raw_intv_list[j]) {
                        bwtintv_t_fpga ik_tmp = ik_list[j];
                        ik_tmp.info = i | ((bwtint_t)start_pos_list[j] << 32);
                        if (intv_n_list[j] < MAX_INTV_ALLOC)
                            forward_local[j][intv_n_list[j]] = ik_tmp;
                        intv_n_list[j]++;
                    }
                    raw_intv_list[j] = true;
                    start_pos_list[j] = i;
                }
                else if (raw_intv_list[j]) {
                    bwtintv_t_fpga ik_tmp;
                    intv_init(&ik_tmp, L2, seq_local[j][i]);
                    ik_list[j] = ik_tmp;
                    start_pos_list[j] = i - 1;
                    raw_intv_list[j] = false;
                }
                else {
                    bwtintv_t_fpga ok;
                    dequeue_data(&ik_list[j], k_addr_buf[j], l_addr_buf[j], bwt_primary, L2, &ok, 3 - seq_local[j][i], 0, datas_pipe);
                    bwtintv_t_fpga ik_tmp = ik_list[j];
                    if (ok.x2 != ik_tmp.x2) {
                        ik_tmp.info = i | ((bwtint_t)start_pos_list[j] << 32);
                        if (intv_n_list[j] < MAX_INTV_ALLOC)
                            forward_local[j][intv_n_list[j]] = ik_tmp;
                        intv_n_list[j]++;
                        if (ok.x2 < min_intv) {
                            start_pos_list[j] = i - 1;
                            intv_init(&ik_tmp, L2, seq_local[j][i]);
                            ik_list[j] = ik_tmp;
                        }
                        else {
                            ik_list[j] = ok;
                        }
                    }
                    else {
                        ik_list[j] = ok;
                    }
                }
            }
        }
    }
fe_intv_pad:
    for (int i = 0; i < tile_size; i++) {
        if (!raw_intv_list[i]) {
            bwtintv_t_fpga ik_tmp = ik_list[i];
            ik_tmp.info = seq_len_local[i] | ((bwtint_t)start_pos_list[i] << 32);
            if (intv_n_list[i] < MAX_INTV_ALLOC)
                forward_local[i][intv_n_list[i]] = ik_tmp;
            intv_n_list[i]++;
        }
    }


/*
    
forward_ext:
    for (i = 0; i < seq_len_local; ++i) {
        if (seq_local[i] > 3) {
            if (!raw_intv) {
                ik.info = i | ((bwtint_t)start_pos << 32);
                if (intv_n < MAX_INTV_ALLOC)
                    forward_local[intv_n] = ik;
                intv_n++;
            }
            raw_intv = true;
            start_pos = i;
        }
        else if (raw_intv) {
            intv_init(&ik, L2, seq_local[i]);
            start_pos = i - 1;
            raw_intv = false;
        }
        else {
            byte_t c = 3 - seq_local[i];
            bwt_extend_fpga_opt(bwt_primary, L2, &ik, &ok, c, 0, addrs_pipe, datas_pipe);
            if (ok.x2 != ik.x2) {
                ik.info = i | ((bwtint_t)start_pos << 32);
                if (intv_n < MAX_INTV_ALLOC)
                    forward_local[intv_n] = ik;
                intv_n++;
                if (ok.x2 < min_intv) {
                    start_pos = i - 1;
                    intv_init(&ik, L2, seq_local[i]);
                }
                else {
                    ik = ok;
                }
            }
            else {
                ik = ok;
            }
        }
    }
    if (!raw_intv) {
        ik.info = i | ((bwtint_t)start_pos << 32);
        if (intv_n < MAX_INTV_ALLOC)
            forward_local[intv_n] = ik;
        intv_n++;
    }*/
#ifdef CL
    int addr_seq = 0;
fe_output:
    for (int i = 0; i < tile_size; i++) {
        forward_intv_n[i] = intv_n_list[i];
        seq_len_be[i] = seq_len_local[i];
seq_be_relay:
        for (int idx = 0; idx < SEQ_LENGTH; idx++) {
            seq_be[addr_seq++] = seq_local[i][idx];
        }
        int bound = intv_n_list[i];
        if (bound > MAX_INTV_ALLOC)
            bound = MAX_INTV_ALLOC;
forward_intv_relay:
        for (int idx = 0; idx < bound; idx++) {
            forward[i * MAX_INTV_ALLOC + idx] = forward_local[i][idx];
        }
    }
#else
fe_output:
    for (int i = 0; i < tile_size; i++) {
        forward_intv_n->write(intv_n_list[i]);
        ap_wait();
        seq_len_be->write(seq_len_local[i]);
        ap_wait();
seq_be_relay:
        for (int idx = 0; idx < SEQ_LENGTH; idx++) {
            seq_be->write(seq_local[i][idx]);
        }   
        ap_wait();
        int bound = intv_n_list[i];
        if (bound > MAX_INTV_ALLOC)
            bound = MAX_INTV_ALLOC;
forward_intv_relay:
        for (int idx = 0; idx < bound; idx++) {
            forward->write(forward_local[i][idx]);
        }
    }
#endif
}

#ifdef CL
void backward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        byte_t seq[SEQ_LENGTH], byte_t seq_len, int forward_intv_n, \
        bwtintv_t_fpga forward[MAX_INTV_ALLOC], int* backward_intv_n, int* long_smem_n, \
        bwtintv_t_fpga backward[MAX_INTV_ALLOC], bwtintv_t_fpga long_smem[MAX_INTV_ALLOC], \
        addr_channel* addrs_pipe, data_channel* datas_pipe, \
        byte_t seq_lfe[SEQ_LENGTH], byte_t* seq_len_lfe, int core_id) {
#else
void backward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        hls::stream<byte_t>* seq, hls::stream<byte_t>& seq_len, hls::stream<int>& forward_intv_n, \
        hls::stream<bwtintv_t_fpga>* forward, hls::stream<int>* backward_intv_n, hls::stream<int>* long_smem_n, \
        hls::stream<bwtintv_t_fpga>* backward, hls::stream<bwtintv_t_fpga>* long_smem, \
        addr_channel* addrs_pipe, data_channel* datas_pipe, \
        hls::stream<byte_t>* seq_lfe, hls::stream<byte_t>* seq_len_lfe) {
#endif
#pragma HLS inline
    bwtintv_t_fpga forward_local[2][MAX_INTV_ALLOC];
#pragma HLS array_partition variable=forward_local complete dim=1
#pragma HLS resource variable=forward_local core=xpm_memory uram
    int forward_intv_n_local;
    bwtintv_t_fpga backward_local[MAX_INTV_ALLOC];
#pragma HLS resource variable=backward_local core=xpm_memory uram
    bwtintv_t_fpga long_smem_local[MAX_INTV_ALLOC];
#pragma HLS resource variable=long_smem_local core=xpm_memory uram
    byte_t seq_local[SEQ_LENGTH];
    byte_t seq_len_local;
    int last_intv_left_bound = 0;
    bwtintv_t_fpga p, ok;
    p.x0 = 0; p.x1 = 0; p.x2 = 0; p.info = 0;
    ok.x0 = 0; ok.x1 = 0; ok.x2 = 0; ok.info = 0;

    bool first_output = true;
    int intv_n = 0;
    int smem_n = 0;
    int min_intv = 1;


#ifdef CL
forward_intv_latch:
    forward_intv_n_local = forward_intv_n;
    seq_len_local = seq_len;
    int bound = forward_intv_n_local;
    if (bound > MAX_INTV_ALLOC)
        bound = MAX_INTV_ALLOC;

    for (int i = 0; i < bound; i++) {
        forward_local[0][i] = forward[i];
    }

seq_be_latch:
    for (int i = 0; i < SEQ_LENGTH; i++) {
        seq_local[i] = seq[i];
    }
    
#else
    forward_intv_n_local = forward_intv_n.read();
    ap_wait();
    seq_len_local = seq_len.read();
    ap_wait();
seq_be_latch:
    for (int i = 0; i < SEQ_LENGTH; i++) {
        seq_local[i] = seq->read();
    }
    ap_wait();
    int bound = forward_intv_n_local;
    if (bound > MAX_INTV_ALLOC)
        bound = MAX_INTV_ALLOC;

forward_intv_latch:
    for (int i = 0; i < bound; i++) {
        forward_local[0][i] = forward->read();
    }
    
#endif

    if (forward_intv_n_local > MAX_INTV_ALLOC || intv_n > MAX_INTV_ALLOC) {
        //exceeds capacity, kick it back to CPU
        forward_intv_n_local = 0;
        intv_n = MAX_INTV_ALLOC + 1;
        smem_n = MAX_INTV_ALLOC + 1;
    }
    
    bwtintv_t_fpga backward_intv_tmp[MAX_INTV_ALLOC];
#pragma HLS resource variable=backward_intv_tmp core=xpm_memory uram
    bool valid_intv_mask[MAX_INTV_ALLOC];
    int intv_id_buf[2][MAX_INTV_ALLOC];
#pragma HLS array_partition variable=intv_id_buf complete dim=1
be_init_intv_mask:
    for (int i = 0; i < MAX_INTV_ALLOC; i++) {
        valid_intv_mask[i] = false;
    }
    int intv_n_tmp = 0;
    int cur_intv_n[2];
#pragma HLS array_partition variable=cur_intv_n complete dim=1
    bool dir = true;
reverse_forward_local_intv:
    for (int i = 0; i < forward_intv_n_local; i++) {
        forward_local[1][forward_intv_n_local - 1 - i] = forward_local[0][i];
        intv_id_buf[1][i] = forward_intv_n_local - 1 - i;
    }
    cur_intv_n[1] = forward_intv_n_local;
    cur_intv_n[0] = 0;
    while (cur_intv_n[dir] > 0) {
        bwtint_t k_addr_buf[MAX_INTV_ALLOC];
        bwtint_t l_addr_buf[MAX_INTV_ALLOC];
        char base_buf[MAX_INTV_ALLOC];
be_push_addr:
        for (int i = 0; i < cur_intv_n[dir]; i++) {
            bwtintv_t_fpga p = forward_local[dir][i];
            int back_start_pos = (p.info >> 32);
            char c = back_start_pos < 0 ? -1 : seq_local[back_start_pos] < 4 ? seq_local[back_start_pos] : -1;
            base_buf[i] = c;
            if (c >= 0) {
                enqueue_addr(bwt_primary, &p, 1, addrs_pipe, &k_addr_buf[i], &l_addr_buf[i]); 
            }
        }
#ifdef CPP
        ap_wait();
#endif
        int cur_trip_count, nxt_trip_count;
        cur_trip_count = cur_intv_n[dir];
        nxt_trip_count = cur_intv_n[!dir];
be_pop_data:
        for (int i = 0; i < cur_trip_count; i++) {
#pragma HLS dependence variable=forward_local false
#pragma HLS dependence variable=intv_id_buf false
            char c = base_buf[i];
            int intv_id = intv_id_buf[dir][i];
            bwtintv_t_fpga ik = forward_local[dir][i];
            bwtintv_t_fpga ok;
            int left_end = (ik.info >> 32) + 1;
            if (c >= 0) {
                dequeue_data(&ik, k_addr_buf[i], l_addr_buf[i], bwt_primary, L2, &ok, c, 1, datas_pipe);
            }
#ifdef CPP
            ap_wait();
#endif
            if (c < 0 || ok.x2 < min_intv) {
                bwtintv_t_fpga tmp = ik;
                tmp.info &= 0xffffffff;
                tmp.info |= ((bwtint_t)left_end << 32);
                backward_intv_tmp[intv_id] = tmp;
                valid_intv_mask[intv_id] = true;
            }
            else {
                ok.info = (ik.info & 0xffffffff) | (((ik.info >> 32) - 1) << 32);
                forward_local[!dir][nxt_trip_count] = ok;
                intv_id_buf[!dir][nxt_trip_count] = intv_id; 
                nxt_trip_count++;
            }
            
        }
        cur_intv_n[!dir] = nxt_trip_count;
        cur_intv_n[dir] = 0;
        dir = !dir;
    }

be_filter_intv:
    for (int i = forward_intv_n_local - 1; i >= 0; i--) {
        if (valid_intv_mask[i]) {
            bwtintv_t_fpga ik = backward_intv_tmp[i];
            int left_end = (ik.info >> 32) + 1;
            int right_end = (uint_t)ik.info;
            if ((first_output || (left_end < last_intv_left_bound))) {
                last_intv_left_bound = left_end;
                if (first_output)
                    first_output = false;
                int seed_len = (uint_t)ik.info - (ik.info >> 32);
                if (seed_len >= MIN_SEED_LEN) {
                    if (intv_n < MAX_INTV_ALLOC) {
                        backward_local[intv_n] = ik;
                    }
                    intv_n++;
                    if (seed_len >= SPLIT_LEN && ik.x2 <= SPLIT_WIDTH) {
                        if (smem_n < MAX_INTV_ALLOC) {
                            long_smem_local[smem_n] = ik;
                        }
                        smem_n++;
                    }
                }
            }
        }
    }

/*
backward_ext:
    for (int j = forward_intv_n_local - 1; j >= 0; j--) {
        bwtintv_t_fpga p = forward_local[0][j];
        int back_start_pos = (p.info >> 32);
backward_ext_innner:
        for (int i = back_start_pos; i >= -1; i--) {
            char c = i < 0 ? -1 : seq_local[i] < 4 ? seq_local[i] : -1;
            if (c >= 0) {
                bwt_extend_fpga_opt(bwt_primary, L2, &p, &ok, c, 1, addrs_pipe, datas_pipe);
            }
            if (c < 0 || ok.x2 < min_intv) {
                if ((first_output || (i + 1 < last_intv_left_bound))) {
                    bwtintv_t_fpga tmp = p;
                    tmp.info &= 0xffffffff;
                    tmp.info |= ((bwtint_t)(i + 1) << 32);
                    last_intv_left_bound = i + 1;
                    if (first_output)
                        first_output = false;
                    int seed_len = (uint_t)tmp.info - (tmp.info >> 32);
                    if (seed_len >= MIN_SEED_LEN) {
                        if (intv_n < MAX_INTV_ALLOC) {
                            //printf("writing to backward list, intv_n = %d, intv info = %lx, intv.x2 = %lx\n", intv_n, tmp.info, tmp.x2);
                            backward_local[intv_n] = tmp;
                        }
                        intv_n++;
                        if (seed_len >= SPLIT_LEN && tmp.x2 <= SPLIT_WIDTH) {
                            if (smem_n < MAX_INTV_ALLOC) {
                                long_smem_local[smem_n] = tmp;
                            }
                            smem_n++;
                        }
                    }
                }
                break;
            }
            else {
                ok.info = p.info;
                p = ok;
            }
        }
    }*/
#ifdef CL
    *backward_intv_n = intv_n;
    *long_smem_n = smem_n;
    *seq_len_lfe = seq_len_local;
seq_be_relay:
    for (int i = 0; i < seq_len_local; i++) {
        seq_lfe[i] = seq_local[i];
    }
    int bound0;
    if (smem_n > MAX_INTV_ALLOC)
        bound0 = MAX_INTV_ALLOC;
    else
        bound0 = smem_n;
long_smem_relay:
    for (int i = 0; i < bound0; i++) {
        long_smem[i] = long_smem_local[i];
    }
    int bound1;
    if (intv_n > MAX_INTV_ALLOC)
        bound1 = MAX_INTV_ALLOC;
    else
        bound1 = intv_n;
backward_intv_relay:
    for (int i = 0; i < bound1; i++) {
        backward[i] = backward_local[i];
    }
#else
    backward_intv_n->write(intv_n);
    ap_wait();
    long_smem_n->write(smem_n);
    ap_wait();
    seq_len_lfe->write(seq_len_local);
    ap_wait();
seq_be_relay:
    for (int i = 0; i < seq_len_local; i++) {
        seq_lfe->write(seq_local[i]);
    }
    ap_wait();
    int bound0;
    if (smem_n > MAX_INTV_ALLOC)
        bound0 = MAX_INTV_ALLOC;
    else
        bound0 = smem_n;
long_smem_relay:
    for (int i = 0; i < bound0; i++) {
        long_smem->write(long_smem_local[i]);
    }
    ap_wait();
    int bound1;
    if (intv_n > MAX_INTV_ALLOC)
        bound1 = MAX_INTV_ALLOC;
    else
        bound1 = intv_n;
backward_intv_relay:
    for (int i = 0; i < bound1; i++) {
        backward->write(backward_local[i]);
    }
#endif
}

#ifdef CL
void long_smem_forward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        byte_t seq[SEQ_LENGTH], byte_t seq_len, int long_smem_n, \
        bwtintv_t_fpga long_smem[MAX_INTV_ALLOC], int* forward_split_intv_n, \
        bwtintv_t_fpga forward_split_intv[4 * MAX_INTV_ALLOC], int forward_split_min_intv[4 * MAX_INTV_ALLOC], \
        addr_channel* addrs_pipe, data_channel* datas_pipe, \
        byte_t seq_lbe[SEQ_LENGTH], \
        bwtintv_t_fpga be_intv_in[MAX_INTV_ALLOC], int be_intv_n_in, \
        bwtintv_t_fpga be_intv_out[MAX_INTV_ALLOC], int* be_intv_n_out, int core_id) {
#else
void long_smem_forward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        hls::stream<byte_t>* seq, hls::stream<byte_t>& seq_len, hls::stream<int>& long_smem_n, \
        hls::stream<bwtintv_t_fpga>* long_smem, hls::stream<int>* forward_split_intv_n, \
        hls::stream<bwtintv_t_fpga>* forward_split_intv, hls::stream<int>* forward_split_min_intv, \
        addr_channel* addrs_pipe, data_channel* datas_pipe, \
        hls::stream<byte_t>* seq_lbe, \
        hls::stream<bwtintv_t_fpga>* be_intv_in, hls::stream<int>& be_intv_n_in, \
        hls::stream<bwtintv_t_fpga>* be_intv_out, hls::stream<int>* be_intv_n_out) {
#endif
#pragma HLS inline
    bwtintv_t_fpga long_smem_local[MAX_INTV_ALLOC];
#pragma HLS resource variable=long_smem_local core=xpm_memory uram 
    int long_smem_n_local;
    bwtintv_t_fpga forward_split_intv_local[4 * MAX_INTV_ALLOC];
#pragma HLS resource variable=forward_split_intv_local core=xpm_memory uram 
    int forward_split_min_intv_local[4 * MAX_INTV_ALLOC];
#pragma HLS resource variable=forward_split_min_intv_local core=xpm_memory uram 
    byte_t seq_local[SEQ_LENGTH];
    byte_t seq_len_local;
    int intv_n = 0;

#ifdef CL
    long_smem_n_local = long_smem_n;
    seq_len_local = seq_len;

reseed_seq_latch:
    for (int i = 0; i < seq_len_local; i++) {
        seq_local[i] = seq[i];
    }
    int bound0;
    if (long_smem_n_local > MAX_INTV_ALLOC)
        bound0 = MAX_INTV_ALLOC;
    else
        bound0 = long_smem_n_local;
long_smem_latch:
    for (int i = 0; i < bound0; i++) {
        long_smem_local[i] = long_smem[i];
    }
#else
    long_smem_n_local = long_smem_n.read();
    ap_wait();
    seq_len_local = seq_len.read();
    ap_wait();
reseed_seq_latch:
    for (int i = 0; i < seq_len_local; i++) {
        seq_local[i] = seq->read();
    }
    ap_wait();
    int bound0;
    if (long_smem_n_local > MAX_INTV_ALLOC)
        bound0 = MAX_INTV_ALLOC;
    else
        bound0 = long_smem_n_local;
long_smem_latch:
    for (int i = 0; i < bound0; i++) {
        long_smem_local[i] = long_smem->read();
    }
#endif

    if (long_smem_n_local > MAX_INTV_ALLOC) {
        //exceeds capacity, kick it back to CPU
        long_smem_n_local = 0;
        intv_n = MAX_INTV_ALLOC + 1;
    }
   
reseed_forward_ext:
    for (int i = 0; i < long_smem_n_local; i += MAX_TILE_SIZE) {
        byte_t tile_size = i + MAX_TILE_SIZE < long_smem_n_local ? MAX_TILE_SIZE : long_smem_n_local - i;
        
        short min_intv_list[MAX_TILE_SIZE];
#pragma HLS resource variable=min_intv_list core=RAM_2P_LUTRAM
        bwtintv_t_fpga intv_list[2][MAX_TILE_SIZE];
#pragma HLS array_partition variable=intv_list complete dim=1
#pragma HLS resource variable=intv_list core=RAM_2P_LUTRAM
        short cur_intv_n[2];
#pragma HLS array_partition variable=cur_intv_n complete dim=1
        short intv_n_output[MAX_TILE_SIZE];
#pragma HLS resource variable=intv_n_output core=RAM_2P_LUTRAM
        bwtintv_t_fpga intv_output[MAX_TILE_SIZE][MAX_INTV_ALLOC];
#pragma HLS resource variable=intv_output core=xpm_memory uram
        short intv_id_buf[2][MAX_TILE_SIZE];
#pragma HLS array_partition variable=intv_id_buf complete dim=1
#pragma HLS resource variable=intv_id_buf core=RAM_2P_LUTRAM

        bool dir = true;
        cur_intv_n[1] = 0; 
        cur_intv_n[0] = 0; 
lfe_init_tiling:
        for (int j = 0; j < tile_size; j++) {
            int smem_start = (int)(long_smem_local[i + j].info >> 32);
            int smem_end = (int)(long_smem_local[i + j].info);
            int split = (smem_start + smem_end) >> 1;
            min_intv_list[j] = long_smem_local[i + j].x2 + 1;
            intv_n_output[j] = 0;
            byte_t c = seq_local[split];
            if (split < seq_len_local && c <= 3) {
                bwtintv_t_fpga tmp;
                intv_init(&tmp, L2, c);
                tmp.info = (bwtint_t)(split + 1) | ((bwtint_t)(split - 1) << 32);
                intv_list[1][cur_intv_n[1]] = tmp;
                intv_id_buf[1][cur_intv_n[1]] = j;
                cur_intv_n[1]++;
            }
        }

lfe_tiling:
        while (cur_intv_n[dir] > 0) {
            bwtint_t k_addr_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=k_addr_buf core=RAM_2P_LUTRAM
            bwtint_t l_addr_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=l_addr_buf core=RAM_2P_LUTRAM
            byte_t base_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=base_buf core=RAM_2P_LUTRAM

lfe_push_addr:
            for (int k = 0; k < cur_intv_n[dir]; k++) {
                bwtintv_t_fpga p = intv_list[dir][k];
                int forward_start_pos = (int)(p.info); 
                char c = forward_start_pos >= seq_len_local ? -1 : \
                         (seq_local[forward_start_pos] < 4 ? seq_local[forward_start_pos] : -1);  
                base_buf[k] = c;
                if (c >= 0) {
                    enqueue_addr(bwt_primary, &p, 0, addrs_pipe, &k_addr_buf[k], &l_addr_buf[k]); 
                }
            }
#ifdef CPP
            ap_wait();
#endif
            int cur_trip_count, nxt_trip_count;
            cur_trip_count = cur_intv_n[dir];
            nxt_trip_count = cur_intv_n[!dir];
lfe_pop_data:
            for (int k = 0; k < cur_trip_count; k++) {
#pragma HLS dependence variable=intv_list false
#pragma HLS dependence variable=intv_id_buf false

                char c = base_buf[k];
                int intv_id = intv_id_buf[dir][k];
                int min_intv = min_intv_list[intv_id];
                bwtintv_t_fpga ik = intv_list[dir][k];
                bwtintv_t_fpga ok;
                if (c >= 0) {
                    dequeue_data(&ik, k_addr_buf[k], l_addr_buf[k], bwt_primary, L2, &ok, 3 - c, 0, datas_pipe);
                }
#ifdef CPP
                ap_wait();
#endif
                if (c < 0 || ok.x2 != ik.x2) {
                    int output_intv_n = intv_n_output[intv_id];
                    intv_output[intv_id][output_intv_n] = ik;
                    intv_n_output[intv_id] = output_intv_n + 1;
                }
                if (c >= 0 && (!(ok.x2 != ik.x2 && (ok.x2 < min_intv)))) {   
                    ok.info = (ik.info & 0xffffffff00000000) | (bwtint_t)((int)(ik.info) + 1);
                    intv_list[!dir][nxt_trip_count] = ok;
                    intv_id_buf[!dir][nxt_trip_count] = intv_id;
                    nxt_trip_count++;
                }
            }
            cur_intv_n[!dir] = nxt_trip_count;
            cur_intv_n[dir] = 0;
            dir = !dir;
        }
        //mark the last forward intv generated from current long smem
lfe_concat_intv:
        for (int j = 0; j < tile_size; j++) {
            for (int k = 0; k < intv_n_output[j]; k++) {
                bwtintv_t_fpga tmp = intv_output[j][k]; 
                if (k == intv_n_output[j] - 1) {
                    tmp.x2 = tmp.x2 | (bwtint_t)(0x8000000000000000);
                }
                if (intv_n < 4 * MAX_INTV_ALLOC) {
                    forward_split_intv_local[intv_n] = tmp;
                    forward_split_min_intv_local[intv_n] = min_intv_list[j];
                }
                intv_n++;   
            }
        }
    }


/*
reseed_forward_ext:
    for (int i = 0; i < long_smem_n_local; i++) {
        int smem_start = (int)(long_smem_local[i].info >> 32);
        int smem_end = (int)(long_smem_local[i].info);
        int split = (smem_start + smem_end) >> 1;
        int start_pos;
        int min_intv = long_smem_local[i].x2 + 1;
        bool raw_intv = true;
        bwtintv_t_fpga ik, ok;
        int j = split;
reseed_forward_ext_inner:
        for (j = split; j < seq_len_local; j++) {
            if (seq_local[j] > 3) {
                if (!raw_intv) {
                    ik.info |= ((bwtint_t)(split - 1) << 32);
                    if (intv_n < 4 * MAX_INTV_ALLOC) {
                        forward_split_intv_local[intv_n] = ik;
                        forward_split_min_intv_local[intv_n] = min_intv;
                    }
                    intv_n++;
                }
                break;
            }
            else if (raw_intv) {
                intv_init(&ik, L2, seq_local[j]);
                ik.info = j + 1;
                raw_intv = false;
            }
            else {
                byte_t c = 3 - seq_local[j];
                bwt_extend_fpga_opt(bwt_primary, L2, &ik, &ok, c, 0, addrs_pipe, datas_pipe);
                if (ok.x2 != ik.x2) { // change of the interval size
                    ik.info |= ((bwtint_t)(split - 1) << 32);
                    if (intv_n < 4 * MAX_INTV_ALLOC) {
                        forward_split_intv_local[intv_n] = ik;
                        forward_split_min_intv_local[intv_n] = min_intv;
                    }
                    intv_n++;
                    if (ok.x2 < min_intv) {
                        break;
                    }
                }
                ik = ok;
                ik.info = j + 1;
            }
        }
        if (j == seq_len_local && !raw_intv) {
            ik.info |= ((bwtint_t)(split - 1) << 32);
            if (intv_n < 4 * MAX_INTV_ALLOC) {
                forward_split_intv_local[intv_n] = ik;
                forward_split_min_intv_local[intv_n] = min_intv;
            }
            intv_n++;
        }
        //mark the last forward intv generated from current long smem
        if (intv_n > 0 && intv_n < 4 * MAX_INTV_ALLOC) {
            bwtintv_t_fpga last = forward_split_intv_local[intv_n - 1]; 
            last.x2 = last.x2 | (bwtint_t)(0x8000000000000000);
            forward_split_intv_local[intv_n - 1] = last;
        }
    }*/
    
#ifdef CL
    *forward_split_intv_n = intv_n;
seq_lfe_relay:
    for (int i = 0; i < SEQ_LENGTH; i++) {
        seq_lbe[i] = seq_local[i];
    }
    int bound2;
    if (intv_n > 4 * MAX_INTV_ALLOC) {
        bound2 =  4 * MAX_INTV_ALLOC;
    }
    else
        bound2 = intv_n;
forward_split_intv_relay:
    for (int i = 0; i < bound2; ++i) {
        forward_split_intv[i] = forward_split_intv_local[i];
        forward_split_min_intv[i] = forward_split_min_intv_local[i];
    }
    *be_intv_n_out = be_intv_n_in;
    int bound1;
    if (be_intv_n_in > MAX_INTV_ALLOC)
        bound1 = MAX_INTV_ALLOC;
    else
        bound1 = be_intv_n_in;
be_intv_relay:
    for (int i = 0; i < bound1; i++) {
        be_intv_out[i] = be_intv_in[i];
    }
#else
    forward_split_intv_n->write(intv_n);
    ap_wait();
seq_lfe_relay:
    for (int i = 0; i < SEQ_LENGTH; i++) {
        seq_lbe->write(seq_local[i]);
    }
    ap_wait();
    int bound2;
    if (intv_n > 4 * MAX_INTV_ALLOC) {
        bound2 =  4 * MAX_INTV_ALLOC;
    }
    else
        bound2 = intv_n;
forward_split_intv_relay:
    for (int i = 0; i < bound2; ++i) {
        forward_split_intv->write(forward_split_intv_local[i]);
        forward_split_min_intv->write(forward_split_min_intv_local[i]);
    }
    ap_wait();
    int be_intv_n = be_intv_n_in.read();
    be_intv_n_out->write(be_intv_n);
    ap_wait();
    int bound1;
    if (be_intv_n > MAX_INTV_ALLOC)
        bound1 = MAX_INTV_ALLOC;
    else
        bound1 = be_intv_n;
be_intv_relay:
    for (int i = 0; i < bound1; i++) {
        be_intv_out->write(be_intv_in->read());
    }
#endif
}

#ifdef CL
void long_smem_backward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        byte_t seq[SEQ_LENGTH], int forward_intv_n, \
        bwtintv_t_fpga forward[4 * MAX_INTV_ALLOC], int forward_min_intv[4 * MAX_INTV_ALLOC], \
        int* backward_intv_n, bwtintv_t_fpga backward[MAX_INTV_ALLOC], \
        addr_channel* addrs_pipe, data_channel* datas_pipe, \
        bwtintv_t_fpga be_intv_in[MAX_INTV_ALLOC], int be_intv_n_in) {
#else
void long_smem_backward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        hls::stream<byte_t>* seq, hls::stream<int>& forward_intv_n, \
        hls::stream<bwtintv_t_fpga>* forward, hls::stream<int>* forward_min_intv, \
        hls::stream<int>* backward_intv_n, hls::stream<bwtintv_t_fpga>* backward, \
        addr_channel* addrs_pipe, data_channel* datas_pipe, \
        hls::stream<bwtintv_t_fpga>* be_intv_in, hls::stream<int>& be_intv_n_in) {
#endif
#pragma HLS inline
    bwtintv_t_fpga forward_local[2][4 * MAX_INTV_ALLOC];
#pragma HLS array_partition variable=forward_local complete dim=1
#pragma HLS resource variable=forward_local core=xpm_memory uram
    int forward_intv_n_local;
    int forward_min_intv_local[2][4 * MAX_INTV_ALLOC];
#pragma HLS array_partition variable=forward_min_intv_local complete dim=1
#pragma HLS resource variable=forward_min_intv_local core=xpm_memory uram
    bwtintv_t_fpga backward_local[MAX_INTV_ALLOC];
#pragma HLS resource variable=backward_local core=xpm_memory uram
    byte_t seq_local[SEQ_LENGTH];
#pragma HLS resource variable=seq_local core=RAM_2P_LUTRAM
    byte_t seq_len_local;
    int last_intv_left_bound = 0;
    bwtintv_t_fpga p, ok;
    p.x0 = 0; p.x1 = 0; p.x2 = 0; p.info = 0;
    ok.x0 = 0; ok.x1 = 0; ok.x2 = 0; ok.info = 0;
    //bool first_output = true;
    int intv_n = 0;

#ifdef CL
    forward_intv_n_local = forward_intv_n;
seq_lbe_latch:
    for (int i = 0; i < SEQ_LENGTH; i++) {
        seq_local[i] = seq[i];
    }
    int bound2;
    if (forward_intv_n_local > 4 * MAX_INTV_ALLOC) {
        bound2 =  4 * MAX_INTV_ALLOC;
    }
    else
        bound2 = forward_intv_n_local;
reseed_forward_intv_latch:
    for (int i = 0; i < bound2; i++) {
        forward_local[0][i] = forward[i];
        forward_min_intv_local[0][i] = forward_min_intv[i];
    }
    intv_n = be_intv_n_in;
    int bound1;
    if (intv_n > MAX_INTV_ALLOC) {
        bound1 = MAX_INTV_ALLOC;
    }
    else
        bound1 = intv_n;
reseed_backward_intv_prop:
    for (int i = 0; i < bound1; i++) {
        backward_local[i] = be_intv_in[i];
    }
#else
    forward_intv_n_local = forward_intv_n.read();
    ap_wait();
seq_lbe_latch:
    for (int i = 0; i < SEQ_LENGTH; i++) {
        seq_local[i] = seq->read();
    }
    ap_wait();
    int bound2;
    if (forward_intv_n_local > 4 * MAX_INTV_ALLOC) {
        bound2 =  4 * MAX_INTV_ALLOC;
    }
    else
        bound2 = forward_intv_n_local;
reseed_forward_intv_latch:
    for (int i = 0; i < bound2; i++) {
        forward_local[0][i] = forward->read();
        forward_min_intv_local[0][i] = forward_min_intv->read();
    }
    ap_wait();
    intv_n = be_intv_n_in.read();
    ap_wait();
    int bound1;
    if (intv_n > MAX_INTV_ALLOC) {
        bound1 = MAX_INTV_ALLOC;
    }
    else
        bound1 = intv_n;
reseed_backward_intv_prop:
    for (int i = 0; i < bound1; i++) {
        backward_local[i] = be_intv_in->read();
    }
#endif

    if (forward_intv_n_local > 4 * MAX_INTV_ALLOC || intv_n > MAX_INTV_ALLOC) {
        //exceeds capacity, kick it back to CPU
        forward_intv_n_local = 0;
        intv_n = MAX_INTV_ALLOC + 1;
    }

    bwtintv_t_fpga backward_intv_tmp[4 * MAX_INTV_ALLOC];
#pragma HLS resource variable=backward_intv_tmp core=xpm_memory uram
    bool valid_intv_mask[4 * MAX_INTV_ALLOC];
#pragma HLS resource variable=valid_intv_mask core=RAM_2P_LUTRAM
    bool first_intv_mask[4 * MAX_INTV_ALLOC];
#pragma HLS resource variable=first_intv_mask core=RAM_2P_LUTRAM
    int intv_id_buf[2][4 * MAX_INTV_ALLOC];
#pragma HLS array_partition variable=intv_id_buf complete dim=1
#pragma HLS resource variable=intv_id_buf core=xpm_memory uram
lbe_init_intv_mask:
    for (int i = 0; i < 4 * MAX_INTV_ALLOC; i++) {
        valid_intv_mask[i] = false;
        bwtintv_t_fpga p = forward_local[0][i];
        if((bool)(p.x2 >> 63) && i < forward_intv_n_local) {
            first_intv_mask[i] = true;
            p.x2 = p.x2 & (bwtint_t)(0x7fffffffffffffff);
            forward_local[0][i] = p;
        }
        else {
            first_intv_mask[i] = false;
        }
    }
    int intv_n_tmp = 0;
    int cur_intv_n[2];
#pragma HLS array_partition variable=cur_intv_n complete dim=1
    bool dir = true;
reverse_long_forward_local_intv:
    for (int i = 0; i < forward_intv_n_local; i++) {
        forward_local[1][forward_intv_n_local - 1 - i] = forward_local[0][i];
        forward_min_intv_local[1][forward_intv_n_local - 1 - i] = forward_min_intv_local[0][i]; 
        intv_id_buf[1][i] = forward_intv_n_local - 1 - i;
    }
    cur_intv_n[1] = forward_intv_n_local;
    cur_intv_n[0] = 0;

    while (cur_intv_n[dir] > 0) {
        bwtint_t k_addr_buf[4 * MAX_INTV_ALLOC];
        bwtint_t l_addr_buf[4 * MAX_INTV_ALLOC];
        char base_buf[4 * MAX_INTV_ALLOC];
lbe_push_addr:
        for (int i = 0; i < cur_intv_n[dir]; i++) {
            bwtintv_t_fpga p = forward_local[dir][i];
            int back_start_pos = (p.info >> 32);
            char c = back_start_pos < 0 ? -1 : seq_local[back_start_pos] < 4 ? seq_local[back_start_pos] : -1;
            base_buf[i] = c;
            if (c >= 0) {
                enqueue_addr(bwt_primary, &p, 1, addrs_pipe, &k_addr_buf[i], &l_addr_buf[i]); 
            }
        }
#ifdef CPP
        ap_wait();
#endif
        int cur_trip_count, nxt_trip_count;
        cur_trip_count = cur_intv_n[dir];
        nxt_trip_count = cur_intv_n[!dir];
lbe_pop_data:
        for (int i = 0; i < cur_trip_count; i++) {
#pragma HLS dependence variable=intv_id_buf false
#pragma HLS dependence variable=forward_min_intv_local false
#pragma HLS dependence variable=forward_local false
            char c = base_buf[i];
            int intv_id = intv_id_buf[dir][i];
            int min_intv = forward_min_intv_local[dir][i]; 
            bwtintv_t_fpga ik = forward_local[dir][i];
            bwtintv_t_fpga ok;
            int left_end = (ik.info >> 32) + 1;
            if (c >= 0) {
                dequeue_data(&ik, k_addr_buf[i], l_addr_buf[i], bwt_primary, L2, &ok, c, 1, datas_pipe);
            }
#ifdef CPP
            ap_wait();
#endif
            if (c < 0 || ok.x2 < min_intv) {
                bwtintv_t_fpga tmp = ik;
                tmp.info &= 0xffffffff;
                tmp.info |= ((bwtint_t)left_end << 32);
                backward_intv_tmp[intv_id] = tmp;
                valid_intv_mask[intv_id] = true;
            }
            else {
                ok.info = (ik.info & 0xffffffff) | (((ik.info >> 32) - 1) << 32);
                forward_local[!dir][nxt_trip_count] = ok;
                forward_min_intv_local[!dir][nxt_trip_count] = min_intv;
                intv_id_buf[!dir][nxt_trip_count] = intv_id; 
                nxt_trip_count++;
            }
        }
        cur_intv_n[!dir] = nxt_trip_count;
        cur_intv_n[dir] = 0;
        dir = !dir;
    }

lbe_filter_intv:
    for (int i = forward_intv_n_local - 1; i >= 0; i--) {
        if (valid_intv_mask[i]) {
            bwtintv_t_fpga ik = backward_intv_tmp[i];
            int left_end = (ik.info >> 32) + 1;
            int right_end = (uint_t)ik.info;
            bool first_output = first_intv_mask[i];
            if ((first_output || (left_end < last_intv_left_bound))) {
                last_intv_left_bound = left_end;
                int seed_len = (uint_t)ik.info - (ik.info >> 32);
                if (seed_len >= MIN_SEED_LEN) {
                    if (intv_n < MAX_INTV_ALLOC) {
                        backward_local[intv_n] = ik;
                    }
                    intv_n++;
                }
            }
        }
    }



/*
reseed_backward_ext:
    for (int j = forward_intv_n_local - 1; j >= 0; j--) {
        bwtintv_t_fpga p = forward_local[0][j];
        int min_intv = forward_min_intv_local[0][j];
        bool first_output = (bool)(p.x2 >> 63);
        p.x2 = p.x2 & (bwtint_t)(0x7fffffffffffffff);
        int back_start_pos = (p.info >> 32);
reseed_backward_inner:
        for (int i = back_start_pos; i >= -1; i--) {
            char c = i < 0 ? -1 : seq_local[i] < 4 ? seq_local[i] : -1;
            if (c >= 0) {
                bwt_extend_fpga_opt(bwt_primary, L2, &p, &ok, c, 1, addrs_pipe, datas_pipe);
            }
            if (c < 0 || ok.x2 < min_intv) {
                if ((first_output || (i + 1 < last_intv_left_bound))) {
                    bwtintv_t_fpga tmp = p;
                    tmp.info &= 0xffffffff;
                    tmp.info |= ((bwtint_t)(i + 1) << 32);
                    last_intv_left_bound = i + 1;
                    int seed_len = (uint_t)tmp.info - (tmp.info >> 32);
                    if (seed_len >= MIN_SEED_LEN) {
                        if (intv_n < MAX_INTV_ALLOC) {
                            backward_local[intv_n] = tmp;
                        }
                        intv_n++;
                    }
                }
                break;
            }
            else {
                ok.info = p.info;
                p = ok;
            }
        }
    }*/
#ifdef CL
    *backward_intv_n = intv_n;
    int bound0;
    if (intv_n > MAX_INTV_ALLOC) 
        bound0 = MAX_INTV_ALLOC;
    else
        bound0 = intv_n; 
final_intv_relay:
    for (int i = 0; i < bound0; i++) {
        backward[i] = backward_local[i];
    }
#else
    backward_intv_n->write(intv_n);
    int bound0;
    if (intv_n > MAX_INTV_ALLOC) 
        bound0 = MAX_INTV_ALLOC;
    else
        bound0 = intv_n; 
final_intv_relay:
    for (int i = 0; i < bound0; i++) {
        backward->write(backward_local[i]);
    }
#endif
}
#ifdef CL
void all_like_forward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        byte_t* seq, byte_t* seq_len, int* forward_intv_n, bwtintv_t_fpga* forward, \
        addr_channel* addrs_pipe, data_channel* datas_pipe, byte_t tile_size) {
#else
void all_like_forward_extension(bwtint_t bwt_primary, bwtint_t L2[5], \
        hls::stream<byte_t>* seq, hls::stream<byte_t>* seq_len, hls::stream<int>* forward_intv_n, \
        hls::stream<bwtintv_t_fpga>* forward, addr_channel* addrs_pipe, data_channel* datas_pipe, byte_t tile_size) {
#endif
#pragma HLS inline 
    bwtintv_t_fpga forward_local[MAX_TILE_SIZE][MAX_INTV_ALLOC];
#pragma HLS resource variable=forward_local core=xpm_memory uram
    byte_t seq_local[MAX_TILE_SIZE][SEQ_LENGTH];
#pragma HLS resource variable=seq_local core=RAM_2P_LUTRAM
    byte_t seq_len_local[MAX_TILE_SIZE];
#pragma HLS resource variable=seq_len_local core=RAM_2P_LUTRAM
#ifdef CL
seq_len_afe_latch:
    for (int idx = 0; idx < tile_size; idx++) {
        seq_len_local[idx] = seq_len[idx];
    }
    int addr = 0;
seq_afe_latch:
    for (int idx = 0; idx < tile_size; idx++) {
        for (int jdx = 0; jdx < SEQ_LENGTH; jdx++) {
            seq_local[idx][jdx] = seq[addr++];
        }
    }
#else
    for (int idx = 0; idx < tile_size; idx++) {
        seq_len_local[idx] = seq_len->read();
    }
    ap_wait();
seq_afe_latch:
    for (int idx = 0; idx < tile_size; idx++) {
        for (int jdx = 0; jdx < SEQ_LENGTH; jdx++) {
            seq_local[idx][jdx] = seq->read();
        }
    }
#endif
 
    short intv_n_list[MAX_TILE_SIZE];
#pragma HLS resource variable=intv_n_list core=RAM_2P_LUTRAM
    bool raw_intv_list[MAX_TILE_SIZE];
#pragma HLS resource variable=raw_intv_list core=RAM_2P_LUTRAM
    bwtintv_t_fpga ik_list[MAX_TILE_SIZE];
#pragma HLS resource variable=ik_list core=RAM_2P_LUTRAM
    short start_pos_list[MAX_TILE_SIZE];
#pragma HLS resource variable=start_pos_list core=RAM_2P_LUTRAM

afe_init_raw_intv:
    for (int i = 0; i < MAX_TILE_SIZE; i++) {
        raw_intv_list[i] = true;
        start_pos_list[i] = 0;
        intv_n_list[i] = 0;
    }
    int max_seq_len = 0;
afe_find_max_seq_len:
    for (int i = 0; i < tile_size; i++) {
        if(seq_len_local[i] > max_seq_len)
            max_seq_len = seq_len_local[i];
    }

all_like_forward_ext:
    for (int i = 0; i < max_seq_len; ++i) {
        bwtint_t k_addr_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=k_addr_buf core=RAM_2P_LUTRAM
        bwtint_t l_addr_buf[MAX_TILE_SIZE];
#pragma HLS resource variable=l_addr_buf core=RAM_2P_LUTRAM
afe_push_addr:
        for (int j = 0; j < tile_size; j++) {
            if (i < seq_len_local[j] && !raw_intv_list[j] && seq_local[j][i] <= 3) {
                enqueue_addr(bwt_primary, &ik_list[j], 0, addrs_pipe, &k_addr_buf[j], &l_addr_buf[j]); 
            }
        }
#ifdef CPP
        ap_wait();
#endif
afe_pop_data:
        for (int j = 0; j < tile_size; j++) {
            if (i < seq_len_local[j]) {
                if (seq_local[j][i] > 3) {
                    raw_intv_list[j] = true;
                }
                else if (raw_intv_list[j]) {
                    bwtintv_t_fpga ik_tmp;
                    intv_init(&ik_tmp, L2, seq_local[j][i]);
                    ik_list[j] = ik_tmp;
                    raw_intv_list[j] = false;
                    start_pos_list[j] = i;
                }
                else {
                    bwtintv_t_fpga ok;
                    dequeue_data(&ik_list[j], k_addr_buf[j], l_addr_buf[j], bwt_primary, L2, &ok, 3 - seq_local[j][i], 0, datas_pipe);
                    if (ok.x2 < MAX_MEM_INTV && i - start_pos_list[j] >= MIN_SEED_LEN) {
                        bwtintv_t_fpga tmp = ok;
                        tmp.info = (i + 1) | ((bwtint_t)start_pos_list[j] << 32);
                        if (tmp.x2 > 0) {
                            if (intv_n_list[j] < MAX_INTV_ALLOC)
                                forward_local[j][intv_n_list[j]] = tmp;
                            intv_n_list[j]++;
                        }
                        raw_intv_list[j] = true;
                    }
                    else {
                        ik_list[j] = ok;
                    }
                }
            }
        }
    }

/*
all_like_forward_ext:
    for (int i = 0; i < seq_len_local; ++i) {
        if (seq_local[i] > 3) {
            raw_intv = true;
        }
        else if (raw_intv) {
            intv_init(&ik, L2, seq_local[i]);
            raw_intv = false;
            start_pos = i;
        }
        else {
            byte_t c = 3 - seq_local[i];
            bwt_extend_fpga_opt(bwt_primary, L2, &ik, &ok, c, 0, addrs_pipe, datas_pipe);
            if (ok.x2 < MAX_MEM_INTV && i - start_pos >= MIN_SEED_LEN) {
                bwtintv_t_fpga tmp = ok;
                tmp.info = (i + 1) | ((bwtint_t)start_pos << 32);
                if (tmp.x2 > 0) {
                    if (intv_n < MAX_INTV_ALLOC) {
                        forward_local[intv_n] = tmp;
                    }
                    intv_n++;
                }
                raw_intv = true;
            }
            else {
                ik = ok;
            }
        }
    }*/

#ifdef CL
    int addr_mem = 0;
afe_output:
    for (int i = 0; i < tile_size; i++) {
        forward_intv_n[i] = intv_n_list[i];
        int bound = intv_n_list[i];
        if (bound > MAX_INTV_ALLOC)
            bound = MAX_INTV_ALLOC;
all_like_forward_intv_relay:
        for (int j = 0; j < bound; j++) {
            forward[i * MAX_INTV_ALLOC + j] = forward_local[i][j];
        }
    }
#else
afe_output:
    for (int i = 0; i < tile_size; i++) {
        forward_intv_n->write(intv_n_list[i]);
        ap_wait();
        int bound = intv_n_list[i];
        if (bound > MAX_INTV_ALLOC)
            bound = MAX_INTV_ALLOC;
all_like_forward_intv_relay:
        for (int j = 0; j < bound; j++) {
            forward->write(forward_local[i][j]);
        }
    }
#endif
}

#ifdef CL
void output_stream(bwtintv_t_fpga mem_afe[MAX_INTV_ALLOC], int mem_n_afe, bwtintv_t_fpga mem_lbe[MAX_INTV_ALLOC], int mem_n_lbe, global smem_t* mem, global uint_t* mem_num, int id, int batch_size, ctrl_channel* ch) {
#else
void output_stream(hls::stream<bwtintv_t_fpga>* mem_afe, hls::stream<int>& mem_n_afe, hls::stream<bwtintv_t_fpga>* mem_lbe, hls::stream<int>& mem_n_lbe, smem_t* mem, uint_t* mem_num, int id, int batch_size, ctrl_channel* ch) {
#endif
#pragma HLS inline
    int offset = id * MAX_INTV_ALLOC;
    int total_mem_n = 0;
    int bound0 = 0;
    int bound1 = 0;
#ifdef CL
    int mem_n_afe_local = mem_n_afe;
#else
    int mem_n_afe_local = mem_n_afe.read();
#endif
    if (mem_n_afe_local > MAX_INTV_ALLOC)
        bound0 = MAX_INTV_ALLOC;
    else 
        bound0 = mem_n_afe_local;
output_afe:
    for (int i = 0; i < bound0; i++) {
        smem_t tmp;
#ifdef CL
        tmp.s0 = mem_afe[i].x0;
        tmp.s1 = mem_afe[i].x1;
        tmp.s2 = mem_afe[i].x2;
        tmp.s3 = mem_afe[i].info;
#else
        bwtintv_t_fpga tmp_intv = mem_afe->read();
        tmp(63, 0)    = tmp_intv.x0;
        tmp(127, 64)  = tmp_intv.x1;
        tmp(191, 128) = tmp_intv.x2;
        tmp(255, 192) = tmp_intv.info;
#endif
        mem[offset + i] = tmp; 
    }
#ifdef CL
    int mem_n_lbe_local = mem_n_lbe;
#else
    int mem_n_lbe_local = mem_n_lbe.read();
#endif
    if (mem_n_lbe_local > MAX_INTV_ALLOC)
        bound1 = MAX_INTV_ALLOC;
    else 
        bound1 = mem_n_lbe_local;

output_lbe:
    for (int i = 0; i < bound1; i++) {
        smem_t tmp;
#ifdef CL
        tmp.s0 = mem_lbe[i].x0;
        tmp.s1 = mem_lbe[i].x1;
        tmp.s2 = mem_lbe[i].x2;
        tmp.s3 = mem_lbe[i].info;
#else
        bwtintv_t_fpga tmp_intv = mem_lbe->read();
        tmp(63, 0)    = tmp_intv.x0;
        tmp(127, 64)  = tmp_intv.x1;
        tmp(191, 128) = tmp_intv.x2;
        tmp(255, 192) = tmp_intv.info;
#endif
        if (i + mem_n_afe_local < MAX_INTV_ALLOC) {
            mem[offset + i + mem_n_afe_local] = tmp; 
        }
    }

    mem_num[id] = mem_n_afe_local + mem_n_lbe_local;
#ifdef CPP
    ap_wait();
#endif
    if (id == batch_size - 1) {
        byte_t end = 1;
        write_ctrl_ch_nblk(ch, &end);
    }
}

#ifdef CL
void smem_flow(bwtint_t bwt_primary, bwtint_t L2[5], const global byte_t *seq, const global byte_t* seq_len, \
        global smem_t *mem, global uint_t* mem_num, int batch_size, \
        addr_channel* fe_addr_ch, data_channel* fe_data_ch, addr_channel* be_addr_ch, data_channel* be_data_ch, \
        addr_channel* lfe_addr_ch, data_channel* lfe_data_ch, addr_channel* lbe_addr_ch, data_channel* lbe_data_ch, \
        addr_channel* afe_addr_ch, data_channel* afe_data_ch, ctrl_channel* ctrl_ch, int core_id) {
#else
void smem_flow(bwtint_t bwt_primary, bwtint_t L2[5], byte_t *seq, byte_t* seq_len, \
        smem_t *mem, uint_t* mem_num, int batch_size, \
        addr_channel* fe_addr_ch, data_channel* fe_data_ch, addr_channel* be_addr_ch, data_channel* be_data_ch, \
        addr_channel* lfe_addr_ch, data_channel* lfe_data_ch, addr_channel* lbe_addr_ch, data_channel* lbe_data_ch, \
        addr_channel* afe_addr_ch, data_channel* afe_data_ch, ctrl_channel* ctrl_ch) {
#endif

#pragma HLS inline off
 
#pragma HLS dataflow
    //there are five stages of mem extension
    //first stage.1 (abbr. fe):  forward extension to get all the forward intvs
    //first stage.2 (abbr, afe): all-like forward extension, the fe and afe stage are in parallel 
    //
    //second stage (abbr. be): backward extension all the foward intvs got 
    //                         in the fe stage to get SMEMs 
    //third stage (abbr, lfe): split all the long SMEMs found in the be stage,
    //                         and do the forward extension for them                         
    //forth stage (abbr, lbe): backward extension all the foward intvs got 
    //                         in the lfe stage
#ifdef CL
    byte_t seq_fe[BATCH_SIZE * SEQ_LENGTH];
    byte_t seq_be[BATCH_SIZE * SEQ_LENGTH];
    byte_t seq_lfe[BATCH_SIZE * SEQ_LENGTH];
    byte_t seq_lbe[BATCH_SIZE * SEQ_LENGTH];
    byte_t seq_afe[BATCH_SIZE * SEQ_LENGTH];
    byte_t seq_len_fe[BATCH_SIZE], seq_len_be[BATCH_SIZE], seq_len_lfe[BATCH_SIZE], seq_len_afe[BATCH_SIZE];
    bwtint_t bwt_primary_fe, bwt_primary_be, bwt_primary_lfe, bwt_primary_lbe, bwt_primary_afe;
    bwtint_t L2_fe[5], L2_be[5], L2_lfe[5], L2_lbe[5], L2_afe[5];

    bwtintv_t_fpga forward_intv[BATCH_SIZE * MAX_INTV_ALLOC]; //produced by fe stage, consumed by be stage
    int forward_intv_n[BATCH_SIZE];

    bwtintv_t_fpga backward_intv[BATCH_SIZE * MAX_INTV_ALLOC]; //produced by be stage, consumed by output stage        
    bwtintv_t_fpga long_smem[BATCH_SIZE * MAX_INTV_ALLOC]; //produced by be stage, consumed by lfe stage
    int backward_intv_n[BATCH_SIZE];
    int long_smem_n[BATCH_SIZE];

    bwtintv_t_fpga forward_split_intv[BATCH_SIZE * 4 * MAX_INTV_ALLOC]; //produced by lfe stage, consumed by bfe stage
    int forward_split_min_intv[BATCH_SIZE * 4 * MAX_INTV_ALLOC];
    int forward_split_intv_n[BATCH_SIZE];

    bwtintv_t_fpga backward_split_intv[BATCH_SIZE * MAX_INTV_ALLOC]; //produced by lbe stage, consumed by output stage
    int backward_split_intv_n[BATCH_SIZE];

    bwtintv_t_fpga forward_all_like_intv[BATCH_SIZE * MAX_INTV_ALLOC]; //produced by afe stage, consumed by output stage
    int forward_all_like_intv_n[BATCH_SIZE];
    
    smem_t mem_local[BATCH_SIZE * MAX_INTV_ALLOC];
    int mem_n_local[BATCH_SIZE];
    int id_fe, id_be, id_lfe, id_lbe, id_output;

    bwtintv_t_fpga be_intv_to_lbe[BATCH_SIZE * MAX_INTV_ALLOC];
    int be_intv_n_to_lbe[BATCH_SIZE];
#else
    hls::stream<byte_t> seq_fe;
#pragma HLS stream variable=seq_fe depth=512
    hls::stream<byte_t> seq_be;
#pragma HLS stream variable=seq_be depth=512
    hls::stream<byte_t> seq_lfe;
#pragma HLS stream variable=seq_lfe depth=512
    hls::stream<byte_t> seq_lbe;
#pragma HLS stream variable=seq_lbe depth=512
    hls::stream<byte_t> seq_afe;
#pragma HLS stream variable=seq_afe depth=512
    hls::stream<byte_t> seq_len_fe, seq_len_be, seq_len_lfe, seq_len_afe;
#pragma HLS stream variable=seq_len_fe depth=2
#pragma HLS stream variable=seq_len_be depth=2
#pragma HLS stream variable=seq_len_lfe depth=2
#pragma HLS stream variable=seq_len_afe depth=2

    hls::stream<bwtintv_t_fpga> forward_intv; //produced by fe stage, consumed by be stage
#pragma HLS stream variable=forward_intv depth=2048
#pragma HLS resource variable=forward_intv core=xpm_memory uram
    hls::stream<int> forward_intv_n;
#pragma HLS stream variable=forward_intv_n depth=32

    hls::stream<bwtintv_t_fpga> backward_intv; //produced by be stage, consumed by output stage        
#pragma HLS stream variable=backward_intv depth=512
    hls::stream<bwtintv_t_fpga> long_smem; //produced by be stage, consumed by lfe stage
#pragma HLS stream variable=long_smem depth=512
    hls::stream<int> backward_intv_n;
#pragma HLS stream variable=backward_intv_n depth=2
    hls::stream<int> long_smem_n;
#pragma HLS stream variable=long_smem_n depth=2

    hls::stream<bwtintv_t_fpga> forward_split_intv; //produced by lfe stage, consumed by bfe stage
#pragma HLS stream variable=forward_split_intv depth=2048
#pragma HLS resource variable=forward_split_intv core=xpm_memory uram
    hls::stream<int> forward_split_min_intv;
#pragma HLS stream variable=forward_split_min_intv depth=2048
    hls::stream<int> forward_split_intv_n;
#pragma HLS stream variable=long_smem_n depth=2

    hls::stream<bwtintv_t_fpga> backward_split_intv; //produced by lbe stage, consumed by output stage
#pragma HLS stream variable=backward_split_intv depth=512
    hls::stream<int> backward_split_intv_n;
#pragma HLS stream variable=backward_split_intv_n depth=2

    hls::stream<bwtintv_t_fpga> forward_all_like_intv; //produced by afe stage, consumed by output stage
#pragma HLS stream variable=forward_all_like_intv depth=2048
#pragma HLS resource variable=forward_all_like_intv core=xpm_memory uram
    hls::stream<int> forward_all_like_intv_n;
#pragma HLS stream variable=forward_all_like_intv_n depth=32
    hls::stream<smem_t> mem_local;
#pragma HLS stream variable=mem_local depth=512
    hls::stream<int> mem_n_local;
#pragma HLS stream variable=mem_n_local depth=2
    
    hls::stream<bwtintv_t_fpga> be_intv_to_lbe;
#pragma HLS stream variable=be_intv_to_lbe depth=512
    hls::stream<int> be_intv_n_to_lbe;
#pragma HLS stream variable=be_intv_n_to_lbe depth=2
#endif

        
input_dup: for (int i = 0; i < batch_size; i += MAX_TILE_SIZE) {
        byte_t tile_size = i + MAX_TILE_SIZE < batch_size ? MAX_TILE_SIZE : batch_size - i;
        //duplicate for fe and afe stages
#ifdef CL
        input_dup(seq, seq_len, i, &seq_fe[i * SEQ_LENGTH], &seq_len_fe[i],\
                &seq_afe[i * SEQ_LENGTH], &seq_len_afe[i], tile_size);
#else
        input_dup(seq, seq_len, i, &seq_fe, &seq_len_fe, &seq_afe, &seq_len_afe, tile_size);
#endif
    }

fe: for (int i = 0; i < batch_size; i += MAX_TILE_SIZE) {
        byte_t tile_size = i + MAX_TILE_SIZE < batch_size ? MAX_TILE_SIZE : batch_size - i;
#ifdef CL
        forward_extension(bwt_primary, L2, &seq_fe[i * SEQ_LENGTH], &seq_len_fe[i], \
                &forward_intv_n[i], &forward_intv[i * MAX_INTV_ALLOC], fe_addr_ch, fe_data_ch, \
                &seq_be[i * SEQ_LENGTH], &seq_len_be[i], tile_size);
#else 
        forward_extension(bwt_primary, L2, &seq_fe, &seq_len_fe, \
                &forward_intv_n, &forward_intv, fe_addr_ch, fe_data_ch, \
                &seq_be, &seq_len_be, tile_size);
#endif
    }

afe: for (int i = 0; i < batch_size; i += MAX_TILE_SIZE) {
        byte_t tile_size = i + MAX_TILE_SIZE < batch_size ? MAX_TILE_SIZE : batch_size - i;
#ifdef CL
        all_like_forward_extension(bwt_primary, L2, &seq_afe[i * SEQ_LENGTH], &seq_len_afe[i], \
                &forward_all_like_intv_n[i], &forward_all_like_intv[i * MAX_INTV_ALLOC], \
                afe_addr_ch, afe_data_ch, tile_size);
#else
        all_like_forward_extension(bwt_primary, L2, &seq_afe, &seq_len_afe, \
                &forward_all_like_intv_n, &forward_all_like_intv, \
                afe_addr_ch, afe_data_ch, tile_size);
#endif
    }

be: for (int i = 0; i < batch_size; i++) {
#ifdef CL
        backward_extension(bwt_primary, L2, &seq_be[i * SEQ_LENGTH], seq_len_be[i], \
                forward_intv_n[i], &forward_intv[i * MAX_INTV_ALLOC], &backward_intv_n[i], &long_smem_n[i], \
                &backward_intv[i * MAX_INTV_ALLOC], &long_smem[i * MAX_INTV_ALLOC], be_addr_ch, be_data_ch, \
                &seq_lfe[i * SEQ_LENGTH], &seq_len_lfe[i], core_id);
#else
        backward_extension(bwt_primary, L2, &seq_be, seq_len_be, \
                forward_intv_n, &forward_intv, &backward_intv_n, &long_smem_n, \
                &backward_intv, &long_smem, be_addr_ch, be_data_ch, \
                &seq_lfe, &seq_len_lfe);
#endif
    }

lfe: for (int i = 0; i < batch_size; i++) {
#ifdef CL
        long_smem_forward_extension(bwt_primary, L2, \
                &seq_lfe[i * SEQ_LENGTH], seq_len_lfe[i], long_smem_n[i], &long_smem[i * MAX_INTV_ALLOC], \
                &forward_split_intv_n[i], &forward_split_intv[i * 4 * MAX_INTV_ALLOC], &forward_split_min_intv[i * 4 * MAX_INTV_ALLOC], \
                lfe_addr_ch, lfe_data_ch, \
                &seq_lbe[i * SEQ_LENGTH], \
                &backward_intv[i * MAX_INTV_ALLOC], backward_intv_n[i], \
                &be_intv_to_lbe[i * MAX_INTV_ALLOC], &be_intv_n_to_lbe[i], core_id);
#else
        long_smem_forward_extension(bwt_primary, L2, \
                &seq_lfe, seq_len_lfe, long_smem_n, &long_smem, \
                &forward_split_intv_n, &forward_split_intv, &forward_split_min_intv, \
                lfe_addr_ch, lfe_data_ch, \
                &seq_lbe, \
                &backward_intv, backward_intv_n, \
                &be_intv_to_lbe, &be_intv_n_to_lbe);
#endif
    }
        
lbe: for (int i = 0; i < batch_size; i++) {
#ifdef CL
        long_smem_backward_extension(bwt_primary, L2, \
                &seq_lbe[i * SEQ_LENGTH], forward_split_intv_n[i], &forward_split_intv[i * 4 * MAX_INTV_ALLOC],\
                &forward_split_min_intv[i * 4 * MAX_INTV_ALLOC], &backward_split_intv_n[i], &backward_split_intv[i * MAX_INTV_ALLOC], \
                lbe_addr_ch, lbe_data_ch, &be_intv_to_lbe[i * MAX_INTV_ALLOC], be_intv_n_to_lbe[i]);
#else
        long_smem_backward_extension(bwt_primary, L2, \
                &seq_lbe, forward_split_intv_n, &forward_split_intv,\
                &forward_split_min_intv, &backward_split_intv_n, &backward_split_intv, \
                lbe_addr_ch, lbe_data_ch, &be_intv_to_lbe, be_intv_n_to_lbe);
#endif
    }

output: for (int i = 0; i < batch_size; i++) {
#ifdef CL
        output_stream(&forward_all_like_intv[i * MAX_INTV_ALLOC], forward_all_like_intv_n[i], \
                &backward_split_intv[i * MAX_INTV_ALLOC], backward_split_intv_n[i], mem, mem_num, i, batch_size, ctrl_ch);
#else
        output_stream(&forward_all_like_intv, forward_all_like_intv_n, &backward_split_intv, backward_split_intv_n, mem, mem_num, i, batch_size, ctrl_ch);
#endif
    }

}

#ifdef CPP

void bwt_manager(data_bundle *bwt, bwtint_t bwt_size, addr_channel* addr_ch0, data_channel* data_ch0, \
        addr_channel* addr_ch1, data_channel* data_ch1, \
        addr_channel* addr_ch2, data_channel* data_ch2, \
        addr_channel* addr_ch3, data_channel* data_ch3, \
        addr_channel* addr_ch4, data_channel* data_ch4, \
        ctrl_channel* ctrl_ch) {
#pragma HLS inline off
    
bwt_manage:
    while (1) {
        int idx;
        if (!addr_ch3->empty()) 
            idx = 3;
        else if (!addr_ch1->empty())
            idx = 1;
        else if (!addr_ch0->empty())
            idx = 0;
        else if (!addr_ch2->empty())
            idx = 2;
        else if (!addr_ch4->empty())
            idx = 4;
        else
            idx = -1;
        if (idx >= 0) {
            addr_bundle addrs;
            if(idx == 0) {
                read_addr_ch_nblk(addr_ch0, &addrs);
            }
            else if (idx == 1) {
                read_addr_ch_nblk(addr_ch1, &addrs);
            }
            else if (idx == 2) {
                read_addr_ch_nblk(addr_ch2, &addrs);
            }
            else if (idx == 3) {
                read_addr_ch_nblk(addr_ch3, &addrs);
            }
            else {
                read_addr_ch_nblk(addr_ch4, &addrs);
            }
            data_bundle data;
            addr_to_data(bwt, bwt_size, addrs, &data);
            if (idx == 0) {
                write_data_ch_blk(data_ch0, &data);
            }
            else if(idx == 1) {
                write_data_ch_blk(data_ch1, &data);
            }
            else if (idx == 2) {
                write_data_ch_blk(data_ch2, &data);
            }
            else if (idx == 3) {
                write_data_ch_blk(data_ch3, &data);
            }
            else {
                write_data_ch_blk(data_ch4, &data);
            }
        }
        ap_wait();
        byte_t isbreak;
        if (read_ctrl_ch_nblk(ctrl_ch, &isbreak)) {
            break;
        }
    }
}

void ddr_streaming(bwtint_t bwt_primary, bwtint_t L2[5], byte_t *seq, \
        byte_t* seq_len, smem_t *mem, uint_t* mem_num, int batch_size, data_bundle* bwt, bwtint_t bwt_size) {
#pragma HLS inline off
#pragma HLS dataflow
    static hls::stream<data_bundle> data_ch0, data_ch1, data_ch2, data_ch3, data_ch4;
#pragma HLS stream variable=data_ch0 depth=16
#pragma HLS resource variable=data_ch0 core=FIFO_SRL
#pragma HLS stream variable=data_ch1 depth=256
#pragma HLS stream variable=data_ch2 depth=16
#pragma HLS resource variable=data_ch2 core=FIFO_SRL
#pragma HLS stream variable=data_ch3 depth=1024
#pragma HLS stream variable=data_ch4 depth=16
#pragma HLS resource variable=data_ch4 core=FIFO_SRL
    static hls::stream<addr_bundle> addr_ch0, addr_ch1, addr_ch2, addr_ch3, addr_ch4;
#pragma HLS stream variable=addr_ch0 depth=16
#pragma HLS resource variable=addr_ch0 core=FIFO_SRL
#pragma HLS stream variable=addr_ch1 depth=256
#pragma HLS stream variable=addr_ch2 depth=16
#pragma HLS resource variable=addr_ch2 core=FIFO_SRL
#pragma HLS stream variable=addr_ch3 depth=1024
#pragma HLS stream variable=addr_ch4 depth=16
#pragma HLS resource variable=addr_ch4 core=FIFO_SRL
    static hls::stream<byte_t> ctrl_ch;
#pragma HLS stream variable=ctrl_ch depth=2
        smem_flow(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, \
            &addr_ch0, &data_ch0, &addr_ch1, &data_ch1, \
            &addr_ch2, &data_ch2, &addr_ch3, &data_ch3, \
            &addr_ch4, &data_ch4, &ctrl_ch);
        bwt_manager(bwt, bwt_size, &addr_ch0, &data_ch0, &addr_ch1, &data_ch1, \
            &addr_ch2, &data_ch2, &addr_ch3, &data_ch3, \
            &addr_ch4, &data_ch4, &ctrl_ch); 
}

extern "C" {
    void mem_collect_intv_core0(bwtint_t *bwt_para, byte_t *seq, \
            byte_t* seq_len, smem_t *mem, uint_t* mem_num, int batch_size, data_bundle* bwt) {
#pragma HLS interface m_axi port=bwt_para offset=slave bundle=bwt_para depth=7
#pragma HLS interface m_axi port=seq offset=slave bundle=seq depth=262144
#pragma HLS interface m_axi port=seq_len offset=slave bundle=seq_len depth=1024
#pragma HLS interface m_axi port=mem offset=slave bundle=mem depth=262144
#pragma HLS interface m_axi port=mem_num offset=slave bundle=mem_num depth=1024
#pragma HLS interface m_axi port=bwt offset=slave bundle=bwt depth=6000000 num_read_outstanding=32
#pragma HLS interface s_axilite port=bwt_para bundle=control
#pragma HLS interface s_axilite port=seq bundle=control
#pragma HLS interface s_axilite port=seq_len bundle=control
#pragma HLS interface s_axilite port=mem bundle=control
#pragma HLS interface s_axilite port=mem_num bundle=control
#pragma HLS interface s_axilite port=batch_size bundle=control
#pragma HLS interface s_axilite port=bwt bundle=control
#pragma HLS interface s_axilite port=return bundle=control
        bwtint_t bwt_primary = bwt_para[0];
        bwtint_t L2[5];
#pragma HLS array_partition variable=L2 complete
        for (int i = 0; i < 5; i++) {
#pragma HLS unroll
            L2[i] = bwt_para[i + 1];
        }
        bwtint_t bwt_size = bwt_para[6];
        
        ddr_streaming(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, bwt, bwt_size);
    }
}

extern "C" {
    void mem_collect_intv_core1(bwtint_t *bwt_para, byte_t *seq, \
            byte_t* seq_len, smem_t *mem, uint_t* mem_num, int batch_size, data_bundle* bwt) {
#pragma HLS interface m_axi port=bwt_para offset=slave bundle=bwt_para depth=7
#pragma HLS interface m_axi port=seq offset=slave bundle=seq depth=262144
#pragma HLS interface m_axi port=seq_len offset=slave bundle=seq_len depth=1024
#pragma HLS interface m_axi port=mem offset=slave bundle=mem depth=262144
#pragma HLS interface m_axi port=mem_num offset=slave bundle=mem_num depth=1024
#pragma HLS interface m_axi port=bwt offset=slave bundle=bwt depth=6000000 num_read_outstanding=32
#pragma HLS interface s_axilite port=bwt_para bundle=control
#pragma HLS interface s_axilite port=seq bundle=control
#pragma HLS interface s_axilite port=seq_len bundle=control
#pragma HLS interface s_axilite port=mem bundle=control
#pragma HLS interface s_axilite port=mem_num bundle=control
#pragma HLS interface s_axilite port=batch_size bundle=control
#pragma HLS interface s_axilite port=bwt bundle=control
#pragma HLS interface s_axilite port=return bundle=control
        bwtint_t bwt_primary = bwt_para[0];
        bwtint_t L2[5];
#pragma HLS array_partition variable=L2 complete
        for (int i = 0; i < 5; i++) {
#pragma HLS unroll
            L2[i] = bwt_para[i + 1];
        }
        bwtint_t bwt_size = bwt_para[6];
        
        ddr_streaming(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, bwt, bwt_size);
    }
}

extern "C" {
    void mem_collect_intv_core2(bwtint_t *bwt_para, byte_t *seq, \
            byte_t* seq_len, smem_t *mem, uint_t* mem_num, int batch_size, data_bundle* bwt) {
#pragma HLS interface m_axi port=bwt_para offset=slave bundle=bwt_para depth=7
#pragma HLS interface m_axi port=seq offset=slave bundle=seq depth=262144
#pragma HLS interface m_axi port=seq_len offset=slave bundle=seq_len depth=1024
#pragma HLS interface m_axi port=mem offset=slave bundle=mem depth=262144
#pragma HLS interface m_axi port=mem_num offset=slave bundle=mem_num depth=1024
#pragma HLS interface m_axi port=bwt offset=slave bundle=bwt depth=6000000 num_read_outstanding=32
#pragma HLS interface s_axilite port=bwt_para bundle=control
#pragma HLS interface s_axilite port=seq bundle=control
#pragma HLS interface s_axilite port=seq_len bundle=control
#pragma HLS interface s_axilite port=mem bundle=control
#pragma HLS interface s_axilite port=mem_num bundle=control
#pragma HLS interface s_axilite port=batch_size bundle=control
#pragma HLS interface s_axilite port=bwt bundle=control
#pragma HLS interface s_axilite port=return bundle=control
        bwtint_t bwt_primary = bwt_para[0];
        bwtint_t L2[5];
#pragma HLS array_partition variable=L2 complete
        for (int i = 0; i < 5; i++) {
#pragma HLS unroll
            L2[i] = bwt_para[i + 1];
        }
        bwtint_t bwt_size = bwt_para[6];
        
        ddr_streaming(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, bwt, bwt_size);
    }
}

extern "C" {
    void mem_collect_intv_core3(bwtint_t *bwt_para, byte_t *seq, \
            byte_t* seq_len, smem_t *mem, uint_t* mem_num, int batch_size, data_bundle* bwt) {
#pragma HLS interface m_axi port=bwt_para offset=slave bundle=bwt_para depth=7
#pragma HLS interface m_axi port=seq offset=slave bundle=seq depth=262144
#pragma HLS interface m_axi port=seq_len offset=slave bundle=seq_len depth=1024
#pragma HLS interface m_axi port=mem offset=slave bundle=mem depth=262144
#pragma HLS interface m_axi port=mem_num offset=slave bundle=mem_num depth=1024
#pragma HLS interface m_axi port=bwt offset=slave bundle=bwt depth=6000000 num_read_outstanding=32
#pragma HLS interface s_axilite port=bwt_para bundle=control
#pragma HLS interface s_axilite port=seq bundle=control
#pragma HLS interface s_axilite port=seq_len bundle=control
#pragma HLS interface s_axilite port=mem bundle=control
#pragma HLS interface s_axilite port=mem_num bundle=control
#pragma HLS interface s_axilite port=batch_size bundle=control
#pragma HLS interface s_axilite port=bwt bundle=control
#pragma HLS interface s_axilite port=return bundle=control
        bwtint_t bwt_primary = bwt_para[0];
        bwtint_t L2[5];
#pragma HLS array_partition variable=L2 complete
        for (int i = 0; i < 5; i++) {
#pragma HLS unroll
            L2[i] = bwt_para[i + 1];
        }
        bwtint_t bwt_size = bwt_para[6];
        
        ddr_streaming(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, bwt, bwt_size);
    }
}
#endif



#ifdef CL
kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void bwt_request_core0(const global data_bundle *bwt, bwtint_t bwt_size) {

    while (1) {
        addr_bundle addrs;
        int idx;
        if (get_pipe_num_packets(lbe_addr_channel0) > 0) {
            idx = 3;
        }
        else if (get_pipe_num_packets(be_addr_channel0) > 0) {
            idx = 1;
        }
        else if (get_pipe_num_packets(fe_addr_channel0) > 0) {
            idx = 0;
        }
        else if (get_pipe_num_packets(lfe_addr_channel0) > 0) {
            idx = 2;
        }
        else if (get_pipe_num_packets(afe_addr_channel0) > 0) {
            idx = 4;
        }
        else 
            idx = -1;
        if (idx >= 0) {
            if (idx == 0) {
                read_addr_ch_nblk(&fe_addr_channel0, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&fe_data_channel0, &data);   
            }
            else if (idx == 1) {
                read_addr_ch_nblk(&be_addr_channel0, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&be_data_channel0, &data);   
            }
            else if (idx == 2) {
                read_addr_ch_nblk(&lfe_addr_channel0, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lfe_data_channel0, &data);   
            }
            else if (idx == 3) {
                read_addr_ch_nblk(&lbe_addr_channel0, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lbe_data_channel0, &data);   
            }
            else if (idx == 4) {
                read_addr_ch_nblk(&afe_addr_channel0, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&afe_data_channel0, &data);
            }
        }
      
        byte_t isbreak;
        if (read_ctrl_ch_nblk(&ctrl_channel0, &isbreak)) {
            break;
        }
    }
}

kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void mem_collect_intv_core0(const global bwtint_t *bwt_para, const global byte_t *seq, const global byte_t* seq_len, global smem_t *mem, global uint* mem_num, int batch_size){
    bwtint_t bwt_primary = bwt_para[0];
    bwtint_t L2[5] __attribute__((xcl_array_partition(complete, 1)));
    for (int i = 0; i < 5; i++) {
        L2[i] = bwt_para[i + 1];
    }
    smem_flow(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, \
        &fe_addr_channel0, &fe_data_channel0, &be_addr_channel0, &be_data_channel0, \
        &lfe_addr_channel0, &lfe_data_channel0, &lbe_addr_channel0, &lbe_data_channel0, \
        &afe_addr_channel0, &afe_data_channel0, &ctrl_channel0, 0);
}

kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void bwt_request_core1(const global data_bundle *bwt, bwtint_t bwt_size) {

    while (1) {
        addr_bundle addrs;
        int idx;
        if (get_pipe_num_packets(lbe_addr_channel1) > 0) {
            idx = 3;
        }
        else if (get_pipe_num_packets(be_addr_channel1) > 0) {
            idx = 1;
        }
        else if (get_pipe_num_packets(fe_addr_channel1) > 0) {
            idx = 0;
        }
        else if (get_pipe_num_packets(lfe_addr_channel1) > 0) {
            idx = 2;
        }
        else if (get_pipe_num_packets(afe_addr_channel1) > 0) {
            idx = 4;
        }
        else 
            idx = -1;
        if (idx >= 0) {
            if (idx == 0) {
                read_addr_ch_nblk(&fe_addr_channel1, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&fe_data_channel1, &data);   
            }
            else if (idx == 1) {
                read_addr_ch_nblk(&be_addr_channel1, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&be_data_channel1, &data);   
            }
            else if (idx == 2) {
                read_addr_ch_nblk(&lfe_addr_channel1, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lfe_data_channel1, &data);   
            }
            else if (idx == 3) {
                read_addr_ch_nblk(&lbe_addr_channel1, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lbe_data_channel1, &data);   
            }
            else if (idx == 4) {
                read_addr_ch_nblk(&afe_addr_channel1, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&afe_data_channel1, &data);
            }
        }
      
        byte_t isbreak;
        if (read_ctrl_ch_nblk(&ctrl_channel1, &isbreak)) {
            break;
        }
    }
}


kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void mem_collect_intv_core1(const global bwtint_t *bwt_para, const global byte_t *seq, const global byte_t* seq_len, global smem_t *mem, global uint* mem_num, int batch_size){
    bwtint_t bwt_primary = bwt_para[0];
    bwtint_t L2[5] __attribute__((xcl_array_partition(complete, 1)));
    for (int i = 0; i < 5; i++) {
        L2[i] = bwt_para[i + 1];
    }
    smem_flow(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, \
        &fe_addr_channel1, &fe_data_channel1, &be_addr_channel1, &be_data_channel1, \
        &lfe_addr_channel1, &lfe_data_channel1, &lbe_addr_channel1, &lbe_data_channel1, \
        &afe_addr_channel1, &afe_data_channel1, &ctrl_channel1, 1);
}

kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void bwt_request_core2(const global data_bundle *bwt, bwtint_t bwt_size) {

    while (1) {
        addr_bundle addrs;
        int idx;
        if (get_pipe_num_packets(lbe_addr_channel2) > 0) {
            idx = 3;
        }
        else if (get_pipe_num_packets(be_addr_channel2) > 0) {
            idx = 1;
        }
        else if (get_pipe_num_packets(fe_addr_channel2) > 0) {
            idx = 0;
        }
        else if (get_pipe_num_packets(lfe_addr_channel2) > 0) {
            idx = 2;
        }
        else if (get_pipe_num_packets(afe_addr_channel2) > 0) {
            idx = 4;
        }
        else 
            idx = -1;
        if (idx >= 0) {
            if (idx == 0) {
                read_addr_ch_nblk(&fe_addr_channel2, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&fe_data_channel2, &data);   
            }
            else if (idx == 1) {
                read_addr_ch_nblk(&be_addr_channel2, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&be_data_channel2, &data);   
            }
            else if (idx == 2) {
                read_addr_ch_nblk(&lfe_addr_channel2, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lfe_data_channel2, &data);   
            }
            else if (idx == 3) {
                read_addr_ch_nblk(&lbe_addr_channel2, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lbe_data_channel2, &data);   
            }
            else if (idx == 4) {
                read_addr_ch_nblk(&afe_addr_channel2, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&afe_data_channel2, &data);
            }
        }
      
        byte_t isbreak;
        if (read_ctrl_ch_nblk(&ctrl_channel2, &isbreak)) {
            break;
        }
    }
}


kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void mem_collect_intv_core2(const global bwtint_t *bwt_para, const global byte_t *seq, const global byte_t* seq_len, global smem_t *mem, global uint* mem_num, int batch_size){
    bwtint_t bwt_primary = bwt_para[0];
    bwtint_t L2[5] __attribute__((xcl_array_partition(complete, 1)));
    for (int i = 0; i < 5; i++) {
        L2[i] = bwt_para[i + 1];
    }
    smem_flow(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, \
        &fe_addr_channel2, &fe_data_channel2, &be_addr_channel2, &be_data_channel2, \
        &lfe_addr_channel2, &lfe_data_channel2, &lbe_addr_channel2, &lbe_data_channel2, \
        &afe_addr_channel2, &afe_data_channel2, &ctrl_channel2, 2);
}


kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void bwt_request_core3(const global data_bundle *bwt, bwtint_t bwt_size) {

    while (1) {
        addr_bundle addrs;
        int idx;
        if (get_pipe_num_packets(lbe_addr_channel3) > 0) {
            idx = 3;
        }
        else if (get_pipe_num_packets(be_addr_channel3) > 0) {
            idx = 1;
        }
        else if (get_pipe_num_packets(fe_addr_channel3) > 0) {
            idx = 0;
        }
        else if (get_pipe_num_packets(lfe_addr_channel3) > 0) {
            idx = 2;
        }
        else if (get_pipe_num_packets(afe_addr_channel3) > 0) {
            idx = 4;
        }
        else 
            idx = -1;
        if (idx >= 0) {
            if (idx == 0) {
                read_addr_ch_nblk(&fe_addr_channel3, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&fe_data_channel3, &data);   
            }
            else if (idx == 1) {
                read_addr_ch_nblk(&be_addr_channel3, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&be_data_channel3, &data);   
            }
            else if (idx == 2) {
                read_addr_ch_nblk(&lfe_addr_channel3, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lfe_data_channel3, &data);   
            }
            else if (idx == 3) {
                read_addr_ch_nblk(&lbe_addr_channel3, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&lbe_data_channel3, &data);   
            }
            else if (idx == 4) {
                read_addr_ch_nblk(&afe_addr_channel3, &addrs);
                data_bundle data;
                addr_to_data(bwt, bwt_size, addrs, &data);
                write_data_ch_blk(&afe_data_channel3, &data);
            }
        }
      
        byte_t isbreak;
        if (read_ctrl_ch_nblk(&ctrl_channel3, &isbreak)) {
            break;
        }
    }
}


kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void mem_collect_intv_core3(const global bwtint_t *bwt_para, const global byte_t *seq, const global byte_t* seq_len, global smem_t *mem, global uint* mem_num, int batch_size){
    bwtint_t bwt_primary = bwt_para[0];
    bwtint_t L2[5] __attribute__((xcl_array_partition(complete, 1)));
    for (int i = 0; i < 5; i++) {
        L2[i] = bwt_para[i + 1];
    }
    smem_flow(bwt_primary, L2, seq, seq_len, mem, mem_num, batch_size, \
        &fe_addr_channel3, &fe_data_channel3, &be_addr_channel3, &be_data_channel3, \
        &lfe_addr_channel3, &lfe_data_channel3, &lbe_addr_channel3, &lbe_data_channel3, \
        &afe_addr_channel3, &afe_data_channel3, &ctrl_channel3, 3);
}
#endif
