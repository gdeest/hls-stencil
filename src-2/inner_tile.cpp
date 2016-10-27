#include "inner_tile.h"
#include "hls_types.h"
#include "access_macros.h"

void read_tile(data_t *arr, fifo_t &inputs, int T, int N0, int N1, int t, int x0, int x1) {
  // #pragma HLS INTERFACE m_axi port=arr depth=8000000

  int tt_offset  = N0*N1-N1-1;
  int xx0_offset = N1;

  data_t *curr_timestep = &AT_TILED(arr, T, N0, N1, t, x0, x1, -1, -2, -2);

  pack_t input;

#pragma HLS ARRAY_PARTITION variable=input.vals complete

  for (int tt=-1; tt<ST-1; tt++) {
    data_t *curr = curr_timestep;
    curr_timestep += tt_offset;

    for (int xx0=-2; xx0<S0; xx0++) {
      int len = (tt==-1 || xx0<0) ? S1+2 : 2;

      for (int i=0; i<len; i++) {
#pragma HLS PIPELINE II=1
        int j = (i+0)%2;
        input.vals[j] = curr[i];
        if (j==1) inputs.write(input);
      }

      curr += xx0_offset;
    }
  }
}

#define CEILING(x,y) (((x) + (y) - 1)/(y))

void write_tile(data_t *arr, fifo_t &outputs, int T, int N0, int N1, int t, int x0, int x1) {
  int tt_offset  = N0*N1-N1-1;
  int xx0_offset = N1;

  data_t *curr_timestep = &AT_TILED(arr, T, N0, N1, t, x0, x1, 0, 0, S1);
  pack_t output;

#pragma HLS ARRAY_PARTITION variable=output.vals complete

  for (int tt=0; tt<ST; tt++) {
    data_t *curr = curr_timestep;
    curr_timestep += tt_offset;

    for (int xx0=0; xx0<S0; xx0++) {
      bool long_burst = (tt==ST-1 || xx0>=S0-2);
      int len = long_burst ? S1 : 2;
      int pack_nums = long_burst ? S1/2 : CEILING(2, 2);
      int offset = long_burst ? 0 : 2 % 2;
      data_t *addr = curr - len;


      for (int i=0; i<len; i++) {
#pragma HLS PIPELINE II=1
        int j = (i + offset) % 2;
        if (j == offset) output = outputs.read();
        addr[i] = output.vals[j];
      }

      curr += xx0_offset;
    }
  }
}

#undef CEILING

pack_t compute_pack(pack_t win[3][2]) {
#pragma HLS INLINE
  pack_t result;

  result.vals[0] = COMPUTE(win[2][0].vals[1], win[1][1].vals[0], win[1][0].vals[1], win[1][0].vals[0], win[0][0].vals[1]);
  result.vals[1] = COMPUTE(win[2][1].vals[0], win[1][1].vals[1], win[1][1].vals[0], win[1][0].vals[1], win[0][1].vals[0]);
  return result;
}

#define BUFF(tt,xx0,xx1) buff[(tt)%2][xx0][xx1]

void compute_tile(fifo_t &inputs, fifo_t &outputs) {
#pragma HLS INLINE

  pack_t buff[2][S0/1][S1/2];
#pragma HLS RESOURCE variable=buff core=RAM_2P
#pragma HLS DATA_PACK variable=buff

  pack_t lines[2][(S1/2)+1];
#pragma HLS ARRAY_PARTITION variable=lines complete dim=1

  pack_t win[3][2];
#pragma HLS ARRAY_PARTITION variable=win complete dim=0


  for (int tt=0; tt<ST; tt++) {
    for (int xx0=-2; xx0<S0/1; xx0++) {
      for (int xx1=-1; xx1<S1/2; xx1++) {
#pragma HLS PIPELINE II=1
        bool outside = xx0<0 || xx1<0;
        bool read_fifo = outside || tt==0;
        bool output_result = (tt==ST-1 || xx0>=S0/1-2 || xx1==S1/2-1);

        pack_t input;

        {
#pragma HLS LATENCY max=16
          // Shift window down
          win[0][0] = win[0][1];
          win[1][0] = win[1][1];
          win[2][0] = win[2][1];
          
          
          /* Read top window row */
          // Retrieve top left and middle pack_t from line buffers
          win[0][1] = lines[0][xx1+1];
          win[1][1] = lines[1][xx1+1];
          

          // Retrieve top right pack_t from either previous timestep or input FIFO.
          if (read_fifo) {
            input = inputs.read();
          } else {
            input = BUFF(tt-1,xx0,xx1);
          }

          win[2][1] = input;

          // Shift lines right and store input.
          lines[0][xx1+1] = lines[1][xx1+1];
          lines[1][xx1+1] = input;
        }

        if (!outside) {
          pack_t result = compute_pack(win);

          BUFF(tt, xx0, xx1) = result;
          if (output_result)
            outputs.write(result);
        }
      }
    }
  }
}

#undef BUFF


void read_coords
( int ntiles,
  int *coords,
  hls::stream<tile_coord_t> &fi_ctrl2read,
  hls::stream<tile_coord_t> &fi_ctrl2write)
{
  tile_coord_t coord;
  for (int i=0; i<ntiles*3; i++) {
    #pragma HLS PIPELINE II=1
    int j = i%3;
    int v = coords[i];
    if (j==0) coord.t  = v;
    if (j==1) coord.x0 = v;
    if (j==2) {
      coord.x1 = v;
      fi_ctrl2read.write(coord);
      fi_ctrl2write.write(coord);
    }
  }
}

void read_actor(int ntiles, data_t *arr_in, fifo_t &inputs, int T, int N0, int N1, hls::stream<tile_coord_t> &coords) {
  for (int i=0; i<ntiles; i++) {
    tile_coord_t coord = coords.read();
    read_tile(arr_in, inputs, T, N0, N1, coord.t, coord.x0, coord.x1);
  }
}

void write_actor(int ntiles, data_t *arr_out, fifo_t &outputs, int T, int N0, int N1, hls::stream<tile_coord_t> &coords) {
  for (int i=0; i<ntiles; i++) {
    tile_coord_t coord = coords.read();
    write_tile(arr_out, outputs, T, N0, N1, coord.t, coord.x0, coord.x1);
  }
}

void compute_actor(int ntiles, fifo_t &inputs, fifo_t &outputs) {
  for (int i=0; i<ntiles; i++) {
    compute_tile(inputs, outputs);
  }
}

int inner_tiles(data_t *arr_in, data_t *arr_out, int T, int N0, int N1, int ntiles, int *coords) {
#pragma HLS DATAFLOW

  hls::stream<tile_coord_t> fi_ctrl2read;
  hls::stream<tile_coord_t> fi_ctrl2write;
#pragma HLS stream variable=fi_ctrl2read depth=10
#pragma HLS stream variable=fi_ctrl2write depth=10

  fifo_t inputs, outputs;
  #pragma HLS DATA_PACK variable=inputs
  #pragma HLS DATA_PACK variable=outputs
  #pragma HLS STREAM variable=inputs depth=256
  #pragma HLS STREAM variable=outputs depth=256

  read_coords(ntiles, coords, fi_ctrl2read, fi_ctrl2write);
  read_actor(ntiles, arr_in, inputs, T, N0, N1, fi_ctrl2read);
  compute_actor(ntiles, inputs, outputs);
  write_actor(ntiles, arr_out, outputs, T, N0, N1, fi_ctrl2write);

  return 0;
}
