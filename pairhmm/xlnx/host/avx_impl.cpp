#include "host/avx_impl.h"
#include "host/avx-pairhmm.h"

float (*compute_fp_avxs)(testcase*) = &compute_full_prob_avxs;
double (*compute_fp_avxd)(testcase*) = &compute_full_prob_avxd;

