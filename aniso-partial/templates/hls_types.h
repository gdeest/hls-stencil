#ifndef HLS_TYPES_H
#define HLS_TYPES_H

#include <hls_stream.h>
#include "jacobi2d.h"

typedef struct {
  int t;
  int x0;
  int x1;
} tile_coord_t;

typedef struct {
	data_t vals[{{unroll[1]}}];
} pack_t;

typedef hls::stream<pack_t> fifo_t;

#endif
