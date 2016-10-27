#include "inner_tile.h"
#include "hls_types.h"
#include "access_macros_partial_skewing.h"

void read_line(data_t *addr, fifo_t &inputs) {
#pragma HLS INLINE off
  pack_t input;
#pragma HLS ARRAY_PARTITION variable=input.vals complete
  for (int xx1=0; xx1<S2; xx1++) {
#pragma HLS PIPELINE II=1
    int i = xx1%UNROLL;
    input.vals[i] = addr[i];
    if (i==UNROLL-1)
      inputs.write(input);
  }
}

void read_tile(data_t *arr, fifo_t &inputs, int T, int N, int t, int x0) {
  for (int tt=-1; tt<S0-1; tt++) {
    for (int xx0=-{{uhalo[0][0]}}; xx0<(tt==-1 ? S1 : 0); xx0++) {
      read_line(&AT_TILED(arr, T, N, t, x0, tt, xx0, 0), inputs);
    }
  }
}

void write_tile(data_t *arr, fifo_t &outputs, int T, int N, int t, int x0) {
  pack_t output;
#pragma HLS ARRAY_PARTITION variable=output.vals complete

  for (int tt=0; tt<S0; tt++) {
    for (int xx0=(tt==S0-1 ? 0 : S1-{{uhalo[0][0]}}); xx0<S1; xx0++) {
      for (int xx1=0; xx1<S2; xx1++) {
#pragma HLS PIPELINE II=1

        int i = xx1%UNROLL;
        if (i==0) {
          output = outputs.read();
        }
        AT_TILED(arr, T, N, t, x0, tt, xx0, xx1) = output.vals[i];
      }
    }
  }
}

pack_t compute_pack(pack_t win{% for i in winsize[0] -%}[{{i}}]{% endfor %}) {
#pragma HLS INLINE
  pack_t result;

  {% for num, args in updates -%}
  result.vals[{{num}}] = COMPUTE_NODE({% for idxs in args -%}
      win{% for j in idxs[0] -%}
        [{{j}}]
      {%- endfor %}.vals[{{idxs[1]}}]
      {%- if loop.index < loop.length -%}, {% endif %}
    {%- endfor %});
  {% endfor -%}


  return result;
}


#define BUFF(tt,xx0,xx1) buff[(tt)%2][xx0][xx1]

void compute_tile(fifo_t &inputs, fifo_t &outputs) {
#pragma HLS INLINE

  pack_t buff[2][S1][S2/UNROLL];
#pragma HLS RESOURCE variable=buff core=RAM_2P
#pragma HLS DATA_PACK variable=buff

  pack_t lines[2][(S2/UNROLL)+{{uhalo[1][0] + uhalo[1][1]}}];
#pragma HLS ARRAY_PARTITION variable=lines complete dim=1

  pack_t win[3][3];
#pragma HLS ARRAY_PARTITION variable=win complete dim=0

  pack_t zero;
  {% for i in range(unroll[1]) -%}
  zero.vals[{{i}}] = 0.0;
  {% endfor %}

  for (int tt=0; tt<S0; tt++) {
    for (int xx0=-2; xx0<S1; xx0++) {
      for (int xx1=-{{uhalo[1][0]}}; xx1<S2/UNROLL+{{uhalo[1][1]}}; xx1++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=buff inter=false
        bool boundary = xx1 < 0 || xx1 >= S2/UNROLL;
        bool read_fifo = tt == 0 || xx0<0;
        bool computation = xx0 >= 0 && xx1 >= {{uhalo[1][1]}};

        //bool output_result = (tt==ST-1 || xx0>=S0/1-2 || xx1==S1/-1);
        bool output_result= tt == S0-1 || xx0>=S1-{{uhalo[0][0]}};

        pack_t input;

        // Shift window down
        {% for j in range(winsize[0][1]-1) -%}
        {% for i in range(winsize[0][0]) -%}
        win[{{i}}][{{j}}] = win[{{i}}][{{j+1}}];
        {% endfor %}
        {%- endfor %}

        /* Read top window row */
        // Retrieve top left and middle pack_t from line buffers
        {% for i in range(winsize[0][0]-1) -%}
        win[{{i}}][{{winsize[0][1]-1}}] = lines[{{i}}][xx1+{{uhalo[1][0]}}];
        {% endfor %}

        // Retrieve top right pack_t from either previous timestep or input FIFO.
        if (boundary) {
          input = zero;
        } else if (read_fifo) {
          input = inputs.read();
        } else {
          input = BUFF(tt-1,xx0,xx1);
        }

        win[{{winsize[0][0]-1}}][{{winsize[0][1]-1}}] = input;

        // Shift lines right and store input.
        {% for i in range(winsize[0][0]-2) -%}
        lines[{{i}}][xx1+{{uhalo[1][0]}}] = lines[{{i+1}}][xx1+{{uhalo[1][0]}}];
        {% endfor -%}
        lines[{{winsize[0][0]-2}}][xx1+{{uhalo[1][0]}}] = input;

        // lines[0][xx1+1] = lines[1][xx1+1];
        // lines[1][xx1+1] = input;

        if (computation) {
          pack_t result = compute_pack(win);

          BUFF(tt, xx0, xx1-{{uhalo[1][1]}}) = result;
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
  for (int i=0; i<ntiles*2; i++) {
    #pragma HLS PIPELINE II=1
    int j = i%2;
    int v = coords[i];
    if (j==0) coord.t  = v;
    if (j==1) {
      coord.x0 = v;
      fi_ctrl2read.write(coord);
      fi_ctrl2write.write(coord);
    }
  }
}

void read_actor(int ntiles, data_t *arr_in, fifo_t &inputs, int T, int N, hls::stream<tile_coord_t> &coords) {
  for (int i=0; i<ntiles; i++) {
    tile_coord_t coord = coords.read();
    read_tile(arr_in, inputs, T, N, coord.t, coord.x0);
  }
}

void write_actor(int ntiles, data_t *arr_out, fifo_t &outputs, int T, int N, hls::stream<tile_coord_t> &coords) {
  for (int i=0; i<ntiles; i++) {
    tile_coord_t coord = coords.read();
    write_tile(arr_out, outputs, T, N, coord.t, coord.x0);
  }
}

void compute_actor(int ntiles, fifo_t &inputs, fifo_t &outputs) {
  for (int i=0; i<ntiles; i++) {
    compute_tile(inputs, outputs);
  }
}

int inner_tiles(data_t *arr_in, data_t *arr_out, int T, int N, int ntiles, int coords[2]) {
#pragma HLS INTERFACE m_axi port=arr_in offset=direct depth=9000
#pragma HLS INTERFACE m_axi port=arr_out offset=direct depth=9000//depth=2097152
#pragma HLS DATAFLOW

  hls::stream<tile_coord_t> fi_ctrl2read;
  hls::stream<tile_coord_t> fi_ctrl2write;
#pragma HLS stream variable=fi_ctrl2read depth=10
#pragma HLS stream variable=fi_ctrl2write depth=10

  fifo_t inputs, outputs;
#pragma HLS DATA_PACK variable=inputs
#pragma HLS DATA_PACK variable=outputs
#pragma HLS stream variable=inputs depth=1024
#pragma HLS stream variable=outputs depth=1024

  read_coords(ntiles, coords, fi_ctrl2read, fi_ctrl2write);
  read_actor(ntiles, arr_in, inputs, T, N, fi_ctrl2read);
  compute_actor(ntiles, inputs, outputs);
  write_actor(ntiles, arr_out, outputs, T, N, fi_ctrl2write);

#ifdef SW
  assert(fi_ctrl2write.empty());
  assert(inputs.empty());
  assert(outputs.empty());
#endif

  return 0;
}
