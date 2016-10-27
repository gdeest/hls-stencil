#ifndef INNER_TILE_H
#define INNER_TILE_H

#include "jacobi2d.h"
#include <cmath>
#include <hls_stream.h>

/* arr_in */
#pragma SDS data zero_copy(arr_in[0:((T+1)*N*N)])
#pragma SDS data sys_port(arr_in:AFI)
#pragma SDS data mem_attribute(arr_in:PHYSICAL_CONTIGUOUS|NON_CACHEABLE)
/* arr_out */
#pragma SDS data zero_copy(arr_out[0:((T+1)*N*N)])
#pragma SDS data sys_port(arr_out:AFI)
#pragma SDS data mem_attribute(arr_out:PHYSICAL_CONTIGUOUS|NON_CACHEABLE)
/* coords */
#pragma SDS data copy(coords[0:2*ntiles])
#pragma SDS data access_pattern(coords:SEQUENTIAL)
int inner_tiles(data_t *arr_in, data_t *arr_out, int T, int N, int ntiles, int *coords);

data_t compute_point(data_t NE, data_t E, data_t SE,
                     data_t N,  data_t X, data_t S,
                     data_t NW, data_t W, data_t SW);
#endif
