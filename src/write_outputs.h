#ifndef WRITE_OUTPUTS_H
#define WRITE_OUTPUTS_H

#include "hls_types.h"

void write_outputs(int ntiles, hls::stream<tile_coord_t> &coords, data_t *arr, int T, int N0, int N1, foutput_t &outputs);

#endif
