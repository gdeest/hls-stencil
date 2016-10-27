#ifndef COMMON_PARTIAL_SKEWING_H
#define COMMON_PARTIAL_SKEWING_H

#include <hls_stream.h>

// Numerical data type
typedef float data_t;

// Unrolling factor
// #define UF 2
/* typedef struct { */
/* 	data_t vals[UNROLL_FACTOR]; */
/* } pack_t; */

/* typedef struct { */
/*   int t; */
/*   int x0; */
/*   int x1; */
/* } tile_coord_t; */

#define COMPUTE_NODE(d0,d1,d2,d3,d4,d5,d6,d7,d8) compute_point(d0,d1,d2,d3,d4,d5,d6,d7,d8)


#endif
