#ifndef READ_INPUTS_H
#define READ_INPUTS_H

#include "hls_types.h"

void read_inputs(int ntiles, hls::stream<tile_coord_t> &coords, data_t *arr, int T, int N0, int N1,
                 finput_t &inputs0,
                 finput_t &inputs1,
                 finput_t &inputs2,
                 finput_t &inputs3,
                 finput_t &inputs4,
                 finput_t &inputs5,
                 finput_t &inputs6,
                 finput_t &inputs7
                 );

#endif
