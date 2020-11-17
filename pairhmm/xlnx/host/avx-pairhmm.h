#include <stdint.h>
//#include "host_type.h"
#include "host/Context.h"

#include "host/avx-types.h"

#undef SIMD_ENGINE
#define SIMD_ENGINE avx

#include "host/avx-functions-float.h"
#include "host/avx-vector-shift.h"
#include "host/avx-pairhmm-template.h"

#include "host/avx-functions-double.h"
#include "host/avx-vector-shift.h"
#include "host/avx-pairhmm-template.h"
