#ifndef INNER_TILE_H
#define INNER_TILE_H

#include "jacobi2d.h"
#include "hls_types.h"
#include <hls_stream.h>

/* /\* arr_in_0 *\/ */
/* #pragma SDS data zero_copy(arr_in_0[0:10485760]) */
/* #pragma SDS data sys_port(arr_in_0:ACP) */
/* #pragma SDS data mem_attribute(arr_in_0:PHYSICAL_CONTIGUOUS) */

/* /\* arr_in_1 *\/ */
/* #pragma SDS data zero_copy(arr_in_1[0:10485760]) */
/* #pragma SDS data sys_port(arr_in_1:ACP) */
/* #pragma SDS data mem_attribute(arr_in_1:PHYSICAL_CONTIGUOUS) */

/* /\* arr_in_2 *\/ */
/* #pragma SDS data zero_copy(arr_in_2[0:10485760]) */
/* #pragma SDS data sys_port(arr_in_2:ACP) */
/* #pragma SDS data mem_attribute(arr_in_2:PHYSICAL_CONTIGUOUS) */

/* /\* arr_in_aux *\/ */
/* #pragma SDS data zero_copy(arr_in_aux[0:10485760]) */
/* #pragma SDS data sys_port(arr_in_aux:ACP) */
/* #pragma SDS data mem_attribute(arr_in_aux:PHYSICAL_CONTIGUOUS) */

/* /\* arr_out_0 *\/ */
/* #pragma SDS data zero_copy(arr_out_0[0:10485760]) */
/* #pragma SDS data sys_port(arr_out_0:ACP) */
/* #pragma SDS data mem_attribute(arr_out_0:PHYSICAL_CONTIGUOUS) */

/* /\* arr_out_1 *\/ */
/* #pragma SDS data zero_copy(arr_out_1[0:10485760]) */
/* #pragma SDS data sys_port(arr_out_1:ACP) */
/* #pragma SDS data mem_attribute(arr_out_1:PHYSICAL_CONTIGUOUS) */

/* /\* arr_out_2 *\/ */
/* #pragma SDS data zero_copy(arr_out_2[0:10485760]) */
/* #pragma SDS data sys_port(arr_out_2:ACP) */
/* #pragma SDS data mem_attribute(arr_out_2:PHYSICAL_CONTIGUOUS) */

/* /\* arr_out_aux *\/ */
/* #pragma SDS data zero_copy(arr_out_aux[0:10485760]) */
/* #pragma SDS data sys_port(arr_out_aux:ACP) */
/* #pragma SDS data mem_attribute(arr_out_aux:PHYSICAL_CONTIGUOUS) */

/* /\* coords *\/ */
/* #pragma SDS data zero_copy(coords[0:3*ntiles]) */
/* #pragma SDS data sys_port(coords:ACP) */
/* #pragma SDS data mem_attribute(coords:PHYSICAL_CONTIGUOUS) */
int inner_tiles(data_t *arr_in_0,
                data_t *arr_in_1,
                data_t *arr_in_2,
                data_t *arr_in_aux,
                data_t *arr_out_0,
                data_t *arr_out_1,
                data_t *arr_out_2,
                data_t *arr_out_aux,
                int T, int N, int ntiles, int *coords);

void read_coords(int ntiles, int *coords,
                 fcoord_t &coords_load_0,
                 fcoord_t &coords_load_1,
                 fcoord_t &coords_load_2,
                 fcoord_t &coords_load_aux,
                 fcoord_t &coords_store_0,
                 fcoord_t &coords_store_1, 
                 fcoord_t &coords_store_2,
                 fcoord_t &coords_store_aux);

void load_B0(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in, fdata_t &inputs_0);
void load_B1(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in, fdata_t &inputs_1);
void load_B2(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in, fdata_t &inputs_2);
void load_Baux(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in, fdata_t &inputs_aux);

void pack_inputs(int ntiles,
                 fdata_t &inputs_0,
                 fdata_t &inputs_1,
                 fdata_t &inputs_2,
                 fdata_t &inputs_aux,
                 fpack_t &inputs);

void compute(int ntiles, fpack_t &inputs, fpack_t &outputs);

void unpack_outputs(int ntiles,
                    fpack_t &outputs,
                    fdata_t &outputs_0,
                    fdata_t &inputs_1,
                    fdata_t &inputs_2,
                    fdata_t &inputs_aux);

void store_B0(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_0, data_t *arr_out);
void store_B1(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_1, data_t *arr_out);
void store_B2(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_2, data_t *arr_out);
void store_Baux(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_aux, data_t *arr_out);

#endif
