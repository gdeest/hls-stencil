#ifndef INNER_TILE_H
#define INNER_TILE_H

#include "jacobi2d.h"
#include "hls_types.h"
#include <hls_stream.h>

/* arr_in */
#pragma SDS data zero_copy(arr_in[0:104857600])
#pragma SDS data sys_port(arr_in:AFI)
#pragma SDS data mem_attribute(arr_in:PHYSICAL_CONTIGUOUS)

/* arr_out */
#pragma SDS data zero_copy(arr_out[0:104857600])
#pragma SDS data sys_port(arr_out:AFI)
#pragma SDS data mem_attribute(arr_out:PHYSICAL_CONTIGUOUS)

/* coords */
/* #pragma SDS data zero_copy(coords[0:104857600]) */
/* #pragma SDS data sys_port(coords:AFI) */
#pragma SDS data access_pattern(coords:SEQUENTIAL)
#pragma SDS data copy(coords[0:3*ntiles]) 
#pragma SDS data mem_attribute(coords:PHYSICAL_CONTIGUOUS|NON_CACHEABLE)
void inner_tiles(int *dummy, data_t *arr_in, data_t *arr_out,
                 int T, int N, int ntiles, int *coords);

void read_coords(int ntiles, int *coords,
                 fcoord_t &coords_load,
                 fcoord_t &coords_store);


void load(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in,
            fdata_t &inputs_0,
            fdata_t &inputs_1,
            fdata_t &inputs_2,
            fdata_t &inputs_aux);


void load_B0(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_0);
void load_B1(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_1);
void load_B2(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_2);
void load_Baux(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_aux);

void unshuffle_B1(int ntiles, fdata_t &in, fdata_t &out);

void pack_inputs(int ntiles,
                 fdata_t &inputs_0,
                 fdata_t &inputs_1,
                 fdata_t &inputs_2,
                 fdata_t &inputs_aux,
                 fpack_t &inputs);

data_t compute_point(data_t NE, data_t E, data_t SE,
                     data_t N,  data_t X, data_t S,
                     data_t NW, data_t W, data_t SW);

void compute(int ntiles, fpack_t &inputs, fpack_t &outputs);

void unpack_outputs(int ntiles,
                    fpack_t &outputs,
                    fdata_t &outputs_0,
                    fdata_t &inputs_1,
                    fdata_t &inputs_2,
                    fdata_t &inputs_aux);

/* void store(int ntiles, int T, int N, fcoord_t &coords, */
/*            fdata_t &outputs_0, */
/*            fdata_t &outputs_1, */
/*            fdata_t &outputs_2, */
/*            fdata_t &outputs_aux, */
/*            data_t *arr_out); */

void store(int ntiles, int T, int N, fcoord_t &coords,
           fdata_t &outputs_0,
           fdata_t &outputs_1,
           fdata_t &outputs_2,
           fdata_t &outputs_aux,
           data_t *arr_out);

void shuffle_B1(int ntiles, fdata_t &in, fdata_t &out);

void store_B0(int T, int N, int x0, int x1, int x2, fdata_t &outputs_0, data_t *arr_out);
void store_B1(int T, int N, int x0, int x1, int x2, fdata_t &outputs_1, data_t *arr_out);
void store_B2(int T, int N, int x0, int x1, int x2, fdata_t &outputs_2, data_t *arr_out);
void store_Baux(int T, int N, int x0, int x1, int x2, fdata_t &outputs_aux, data_t *arr_out);

#endif
