#ifndef HLS_TYPES_H
#define HLS_TYPES_H

#include <hls_stream.h>
typedef float data_t;
#include "jacobi2d.h"

typedef struct {
  int t;
  int x0;
  int x1;
} coord_t;

typedef struct {
	data_t vals[UNROLL];
} pack_t;

typedef hls::stream<data_t> fdata_t;
typedef hls::stream<pack_t> fpack_t;
typedef hls::stream<coord_t> fcoord_t;

#endif
