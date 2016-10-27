#include <assert.h>
#include "inner_tile.h"
#include "hls_types.h"
#include "access_macros.h"

#if not(defined(__SDSCC__) || defined(__SDSVHLS__))
#define DEBUG
#endif

#ifdef DEBUG
#define READ(fifo) (assert(!fifo.empty()), fifo.read())
#define ASSERT_EMPTY(fifo) assert(fifo.empty())
#else
#define READ(fifo) fifo.read()
#define ASSERT_EMPTY(fifo)
#endif

/*************************/
/******** TOPLEVEL********/
/*************************/

int inner_tiles(data_t *arr_in_0,
                data_t *arr_in_1,
                data_t *arr_in_2,
                data_t *arr_in_aux,
                data_t *arr_out_0,
                data_t *arr_out_1,
                data_t *arr_out_2,
                data_t *arr_out_aux,
                int T, int N, int ntiles, int *coords) {
#pragma HLS DATAFLOW

#pragma HLS INTERFACE m_axi port=arr_in_0 offset=slave bundle=m_axi_arr_in_0
#pragma HLS INTERFACE m_axi port=arr_in_1 offset=slave bundle=m_axi_arr_in_1
#pragma HLS INTERFACE m_axi port=arr_in_2 offset=slave bundle=m_axi_arr_in_2
#pragma HLS INTERFACE m_axi port=arr_in_aux offset=slave bundle=m_axi_arr_in_aux
#pragma HLS INTERFACE m_axi port=arr_out_0 offset=slave bundle=m_axi_arr_out_0
#pragma HLS INTERFACE m_axi port=arr_out_1 offset=slave bundle=m_axi_arr_out_1
#pragma HLS INTERFACE m_axi port=arr_out_2 offset=slave bundle=m_axi_arr_out_2
#pragma HLS INTERFACE m_axi port=arr_out_aux offset=slave bundle=m_axi_arr_out_aux
#pragma HLS INTERFACE m_axi port=coords offset=slave bundle=m_axi_coords

  fcoord_t coords_load_0;
  fcoord_t coords_load_1;
  fcoord_t coords_load_2;
  fcoord_t coords_load_aux;

  fcoord_t coords_store_0;
  fcoord_t coords_store_1;
  fcoord_t coords_store_2;
  fcoord_t coords_store_aux;

#pragma HLS STREAM variable=coords_load_0 depth=8
#pragma HLS STREAM variable=coords_load_1 depth=8
#pragma HLS STREAM variable=coords_load_2 depth=8
#pragma HLS STREAM variable=coords_load_aux depth=8

#pragma HLS STREAM variable=coords_store_0 depth=8
#pragma HLS STREAM variable=coords_store_1 depth=8
#pragma HLS STREAM variable=coords_store_2 depth=8
#pragma HLS STREAM variable=coords_store_aux depth=8

  fdata_t inputs_0;
  fdata_t inputs_1;
  fdata_t inputs_2;
  fdata_t inputs_aux;

#pragma HLS STREAM variable=inputs_0 depth=2176
#pragma HLS STREAM variable=inputs_1 depth=4352
#pragma HLS STREAM variable=inputs_2 depth=4352
#pragma HLS STREAM variable=inputs_aux depth=8

  fdata_t outputs_0;
  fdata_t outputs_1;
  fdata_t outputs_2;
  fdata_t outputs_aux;

#pragma HLS STREAM variable=outputs_0 depth=2048
#pragma HLS STREAM variable=outputs_1 depth=4096
#pragma HLS STREAM variable=outputs_2 depth=4096
#pragma HLS STREAM variable=outputs_aux depth=8

  fpack_t inputs;
  fpack_t outputs;

#pragma HLS STREAM variable=inputs depth=2048
#pragma HLS STREAM variable=outputs depth=2048

  read_coords(ntiles, coords,
              coords_load_0, coords_load_1, coords_load_2, coords_load_aux,
              coords_store_0, coords_store_1, coords_store_2, coords_store_aux);

  load_B0(ntiles, T, N, coords_load_0, arr_in_0, inputs_0);
  load_B1(ntiles, T, N, coords_load_1, arr_in_1, inputs_1);
  load_B2(ntiles, T, N, coords_load_2, arr_in_2, inputs_2);
  load_Baux(ntiles, T, N, coords_load_aux, arr_in_aux, inputs_aux);

  pack_inputs(ntiles, inputs_0, inputs_1, inputs_2, inputs_aux, inputs);
  compute(ntiles, inputs, outputs);
  unpack_outputs(ntiles, outputs, outputs_0, outputs_1, outputs_2, outputs_aux);

  store_B0(ntiles, T, N, coords_store_0, outputs_0, arr_out_0);
  store_B1(ntiles, T, N, coords_store_1, outputs_1, arr_out_1);
  store_B2(ntiles, T, N, coords_store_2, outputs_2, arr_out_2);
  store_Baux(ntiles, T, N, coords_store_aux, outputs_aux, arr_out_aux);

  return 0;
}

/*****************************/
/******** READ COORDS ********/
/*****************************/

void read_coords
( int ntiles,
  int *coords,
  fcoord_t &coords_load_0,
  fcoord_t &coords_load_1,
  fcoord_t &coords_load_2,
  fcoord_t &coords_load_aux,
  fcoord_t &coords_store_0,
  fcoord_t &coords_store_1,
  fcoord_t &coords_store_2,
  fcoord_t &coords_store_aux)
{
  coord_t coord;
  for (int i=0; i<ntiles*3; i++) {
#pragma HLS PIPELINE II=1
    int j = i%3;
    int v = coords[i];
    if (j==0) coord.t  = v;
    if (j==1) coord.x0 = v;
    if (j==2) {
      coord.x1 = v;
      coords_load_0.write(coord);
      coords_load_1.write(coord);
      coords_load_2.write(coord);
      coords_load_aux.write(coord);
      coords_store_0.write(coord);
      coords_store_1.write(coord);
      coords_store_2.write(coord);
      coords_store_aux.write(coord);
    }
  }
}

/*****************************/
/******** COMPUTE ************/
/*****************************/

pack_t compute_pack(pack_t win[3][2]) {
#pragma HLS INLINE
  pack_t result;

  result.vals[0] = COMPUTE(win[2][0].vals[1], win[1][1].vals[0], win[1][0].vals[1], win[1][0].vals[0], win[0][0].vals[1]);
  result.vals[1] = COMPUTE(win[2][1].vals[0], win[1][1].vals[1], win[1][1].vals[0], win[1][0].vals[1], win[0][1].vals[0]);

  return result;
}

#define BUFF(tt,xx0,xx1) buff[(tt)%2][xx0][xx1]

void compute_tile(fpack_t &inputs, fpack_t &outputs) {
#pragma HLS INLINE

  pack_t buff[2][S1/1][S2/2];
#pragma HLS RESOURCE variable=buff core=RAM_2P
#pragma HLS DATA_PACK variable=buff

  pack_t lines[2][(S2/2)+1];
#pragma HLS ARRAY_PARTITION variable=lines complete dim=1

  pack_t win[3][2];
#pragma HLS ARRAY_PARTITION variable=win complete dim=0


  for (int tt=0; tt<S0; tt++) {
    for (int xx0=-2; xx0<S1/1; xx0++) {
      for (int xx1=-1; xx1<S2/2; xx1++) {
#pragma HLS PIPELINE II=1
        bool outside = xx0<0 || xx1<0;
        bool read_fifo = outside || tt==0;
        bool output_result = (tt==S0-1 || xx0>=S1/1-2 || xx1==S2/2-1);

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
            input = READ(inputs);
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

void compute(int ntiles, fpack_t &inputs, fpack_t &outputs) {
  for (int i=0; i<ntiles; i++) {
    compute_tile(inputs, outputs);
  }

  ASSERT_EMPTY(inputs);
}

/********************************/
/************** LOAD ************/
/********************************/
void load_B0(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in_0, fdata_t &inputs_0) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_P0_acc(T, N, arr_in_0, coord.t, coord.x0, coord.x1, 0, -HALO_1, 0);
    int count = 0;

    // Assuming HALO_0=1, inputs are already in lexicographic order.
    for (int x=0; x<(S1+HALO_1)*S2; x++) {
#pragma HLS PIPELINE II=1
      inputs_0.write(base[x]);
    }
  }
}

void load_B1(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in_1, fdata_t &inputs_1) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_P1_acc(T, N, arr_in_1, coord.t, coord.x0, coord.x1, 0, 0, -HALO_2);
    data_t tmp[S0][HALO_1][S2+HALO_2];
    int count = 0;

    // Inputs are in order: j,t,i. We hence read the whole face before feeding
    // inputs in lexicographic order in a second swipe.
    //
    // We use a single increasing counter as offset and recover j,t,i using
    // integer modulo / division operations. We can't iterate over j,t,i and use
    // a counter here, as it prevents burst inference from Vivado HLS.
    for (int x=0; x<(S2+HALO_2)*S0*HALO_1; x++) {
#pragma HLS PIPELINE II=1
      int j = x / (S0*HALO_1);
      int r = x % (S0*HALO_1);
      int t = r / HALO_1;
      int i = r % HALO_1;
      tmp[t][i][j] = base[x];
    }

    for (int t=0; t<S0; t++) {
      for (int i=0; i<HALO_1; i++) {
        for (int j=0; j<S2+HALO_2; j++) {
#pragma HLS PIPELINE II=1
          inputs_1.write(tmp[t][i][j]);
        }
      }
    }
  }
}

void load_B2(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in_2, fdata_t &inputs_2) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_P2_acc(T, N, arr_in_2, coord.t, coord.x0, coord.x1, -HALO_0, 0, 0);
    int count = 0;

    // Inputs are already in lexicographic order.
    for (int x=0; x<(S0+HALO_0)*S1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
      inputs_2.write(base[x]);
    }
  }
}

void load_Baux(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in_aux, fdata_t &inputs_aux) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_AUX_acc(T, N, arr_in_aux, coord.t, coord.x0, coord.x1, 0);
    int count = 0;

    // Inputs are already in lexicographic order.
    for (int x=0; x<HALO_0*HALO_1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
      inputs_aux.write(base[x]);
    }
  }
}

/********************************/
/************* STORE ************/
/********************************/
void store_B0(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_0, data_t *arr_out_0) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_P0_acc(T, N, arr_out_0, coord.t+1, coord.x0, coord.x1, 0, 0, 0);

    // Assuming HALO_0=1, outputs are already in the correct order.
    for (int x=0; x<S1*S2; x++) {
#pragma HLS PIPELINE II=1
      base[x] = READ(outputs_0);
    }
  }

  ASSERT_EMPTY(outputs_0);
}

void store_B1(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_1, data_t *arr_out_1) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_P1_acc(T, N, arr_out_1, coord.t, coord.x0+1, coord.x1, 0, 0, 0);

    data_t tmp[S0][HALO_1][S2];

    // Outputs are in lexicographic order, while the required order is j,t,i We
    // store outputs in lexicographic order before storing them in the correct
    // order in a second swipe.
    for (int t=0; t<S0; t++) {
      for (int i=0; i<HALO_1; i++) {
        for (int j=0; j<S2; j++) {
#pragma HLS PIPELINE II=1
          tmp[t][i][j] = READ(outputs_1);
        }
      }
    }

    // See load_B1 for the rationale of using a single counter instead of a loop
    // nest.
    for (int x=0; x<S2*S0*HALO_1; x++) {
#pragma HLS PIPELINE II=1
      int j = x / (S0*HALO_1);
      int r = x % (S0*HALO_1);
      int t = r / HALO_1;
      int i = r % HALO_1;
      base[x] = tmp[t][i][j];
    }
  }

  ASSERT_EMPTY(outputs_1);
}

void store_B2(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_2, data_t *arr_out_2) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_P2_acc(T, N, arr_out_2, coord.t, coord.x0, coord.x1+1, 0, 0, 0);

    // Outputs are already in the correct order.
    for (int x=0; x<S0*S1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
      base[x] = READ(outputs_2);
    }
  }

  ASSERT_EMPTY(outputs_2);
}

void store_Baux(int ntiles, int T, int N, fcoord_t &coords, fdata_t &outputs_aux, data_t *arr_out_aux) {
  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    data_t *base = &B_AUX_acc(T, N, arr_out_aux, coord.t, coord.x0, coord.x1, 0);

    // Outputs are already in the correct order.
    for (int x=0; x<HALO_0*HALO_1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
      base[x] = READ(outputs_aux);
    }
  }

  ASSERT_EMPTY(outputs_aux);
}

/********************************/
/*********** PACK ***************/
/********************************/
void pack_inputs(int ntiles,
                 fdata_t &inputs_0,
                 fdata_t &inputs_1,
                 fdata_t &inputs_2,
                 fdata_t &inputs_aux,
                 fpack_t &inputs)
{
  pack_t input;
#pragma HLS ARRAY_PARTITION variable=input.vals complete

  int offset = UNROLL % HALO_2;
  for (int n=0; n<ntiles; n++) {
    for (int t=-HALO_0; t<S0; t++) {
      for (int i=-HALO_1; i<S1; i++) {
        for (int j=-HALO_2; j<(t<0 || i<0 ? S2 : 0); j++) {
#pragma HLS PIPELINE II=1
          int index = (offset + j + HALO_2) % UNROLL;
          data_t v = 0.0; // TODO: Replace by read from the correct FIFO.

          if (t<0) {
            if (j<0) {
              if (i<0) v = READ(inputs_aux);
              else v = READ(inputs_2);
            } else v = READ(inputs_0);
          } else {
            if (i<0) v = READ(inputs_1);
            else v = READ(inputs_2);
          }

          input.vals[index] = v;
          if (t < S0-1 && index == UNROLL-1) inputs.write(input);
        }
      }

    }
  }

  ASSERT_EMPTY(inputs_0);
  ASSERT_EMPTY(inputs_1);
  ASSERT_EMPTY(inputs_2);
  ASSERT_EMPTY(inputs_aux);
}

/********************************/
/*********** UNPACK *************/
/********************************/
void unpack_outputs(int ntiles,
                    fpack_t &outputs,
                    fdata_t &outputs_0,
                    fdata_t &outputs_1,
                    fdata_t &outputs_2,
                    fdata_t &outputs_aux)
{
  pack_t output;
#pragma HLS ARRAY_PARTITION variable=output.vals complete

  for (int n=0; n < ntiles; n++) {
    for (int t=0; t<S0; t++) {
      for (int i=0; i<S1; i++) {
        int startj = (t==S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;
        for (int j=startj; j<S2; j++) {
#pragma HLS PIPELINE II=1

          int index = j % UNROLL;
          if (j == startj || index == 0) output = READ(outputs);
          data_t v = output.vals[index];

          bool i1border = (i  >= S1-HALO_1);
          bool i2border = (j  >= S2-HALO_2);

          if (t >= S0-HALO_0) {
            outputs_0.write(v);

            if (i1border && i2border) {
              outputs_aux.write(v);
            }

            if (i1border) {
              outputs_1.write(v);
            }
            if (i2border) {
              outputs_2.write(v);
            }
          } else {
            if (i1border) {
              outputs_1.write(v);
            }

            if (i2border) {
              outputs_2.write(v);

            }
          }
        }
      }
    }
  }

  ASSERT_EMPTY(outputs);
}
