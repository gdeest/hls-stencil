#include <assert.h>
#include "inner_tile.h"
#include "hls_types.h"

#define K 4

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

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x) \



/*************************/
/******** TOPLEVEL********/
/*************************/

void inner_tiles(int *dummy, data_t *arr_in, data_t *arr_out,
                 int T, int N, int ntiles, int *coords) {
#pragma HLS DATAFLOW

  fcoord_t coords_load;
  fcoord_t coords_store;

#pragma HLS STREAM variable=coords_load depth=128
#pragma HLS STREAM variable=coords_store depth=128

  fdata_t inputs_0;
  fdata_t inputs_1_shuffled;
  fdata_t inputs_1;
  fdata_t inputs_2;
  fdata_t inputs_aux;

  DO_PRAGMA(HLS STREAM variable=inputs_0 depth=K*HALO_0*(S1+HALO_1)*S2);
  DO_PRAGMA(HLS STREAM variable=inputs_1 depth=K*S0*HALO_1*(S2+HALO_2));
  DO_PRAGMA(HLS STREAM variable=inputs_1_shuffled depth=K*S0*HALO_1*(S2+HALO_2));
  DO_PRAGMA(HLS STREAM variable=inputs_2 depth=K*(S0+HALO_0)*S1*HALO_2);
  DO_PRAGMA(HLS STREAM variable=inputs_aux depth=K*HALO_0*HALO_1*HALO_2);

  fdata_t outputs_0;
  fdata_t outputs_1;
  fdata_t outputs_1_shuffled;
  fdata_t outputs_2;
  fdata_t outputs_aux;

  DO_PRAGMA(HLS STREAM variable=outputs_0 depth=K*HALO_0*S1*S2);
  DO_PRAGMA(HLS STREAM variable=outputs_1 depth=K*S0*HALO_1*S2);
  DO_PRAGMA(HLS STREAM variable=outputs_1_shuffled depth=K*S0*HALO_1*S2);
  DO_PRAGMA(HLS STREAM variable=outputs_2 depth=K*S0*S1*HALO_2);
  DO_PRAGMA(HLS STREAM variable=outputs_aux depth=K*HALO_0*HALO_1*HALO_2);

  fpack_t inputs;
  fpack_t outputs;

  DO_PRAGMA(HLS STREAM variable=inputs depth=1024/UNROLL);
  DO_PRAGMA(HLS STREAM variable=outputs depth=1024/UNROLL);

  read_coords(ntiles, coords, coords_load, coords_store);

  load(ntiles, T, N, coords_load, arr_in, inputs_0, inputs_1_shuffled, inputs_2, inputs_aux);

  unshuffle_B1(ntiles, inputs_1_shuffled, inputs_1);

  pack_inputs(ntiles, inputs_0, inputs_1, inputs_2, inputs_aux, inputs);
  compute(ntiles, inputs, outputs);
  unpack_outputs(ntiles, outputs, outputs_0, outputs_1, outputs_2, outputs_aux);

  shuffle_B1(ntiles, outputs_1, outputs_1_shuffled);
  store(ntiles, T, N, coords_store, outputs_0, outputs_1_shuffled, outputs_2, outputs_aux, arr_out);

  *dummy = 0;
}

/*****************************/
/******** READ COORDS ********/
/*****************************/

void read_coords
( int ntiles,
  int *coords,
  fcoord_t &coords_load,
  fcoord_t &coords_store)
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
      coords_load.write(coord);
      coords_store.write(coord);
    }
  }
}

/*****************************/
/******** COMPUTE ************/
/*****************************/

pack_t compute_pack(pack_t win{% for i in winsize -%}[{{i}}]{% endfor %}) {
#pragma HLS INLINE
  pack_t result;

  {% for num, args in updates -%}
    result.vals[{{num}}] = COMPUTE({% for idxs in args -%}
      win{% for j in idxs[0] -%}
        [{{j}}]
      {%- endfor %}.vals[{{idxs[1]}}]
      {%- if loop.index < loop.length -%}, {% endif %}
    {%- endfor %});
  {% endfor -%}


  return result;
}

#define BUFF(tt,xx0,xx1) buff[(tt)%2][xx0][xx1]

void compute_tile(fpack_t &inputs, fpack_t &outputs) {
#pragma HLS INLINE

  pack_t buff[2][S1/{{unroll[0]}}][S2/{{unroll[1]}}];
#pragma HLS RESOURCE variable=buff core=RAM_2P
#pragma HLS DATA_PACK variable=buff

  pack_t lines[{{winsize[0]-1}}][(S2/{{unroll[1]}})+{{uhalo[1]}}];
#pragma HLS ARRAY_PARTITION variable=lines complete dim=1

  pack_t win[{{winsize[0]}}][{{winsize[1]}}];
#pragma HLS ARRAY_PARTITION variable=win complete dim=0


  for (int tt=0; tt<S0; tt++) {
    for (int xx0=-{{uhalo[0]}}; xx0<S1/{{unroll[0]}}; xx0++) {
      for (int xx1=-{{uhalo[1]}}; xx1<S2/{{unroll[1]}}; xx1++) {
#pragma HLS PIPELINE II=1
        bool outside = xx0<0 || xx1<0;
        bool read_fifo = outside || tt==0;
        bool output_result = (tt==S0-1 || xx0>=S1/{{unroll[0]}}-{{uhalo[0]}} || xx1>=S2/{{unroll[1]}}-{{uhalo[1]}});

        pack_t input;

        // Shift window down
        {% for j in range(winsize[1]-1) -%}
        {% for i in range(winsize[0]) -%}
        win[{{i}}][{{j}}] = win[{{i}}][{{j+1}}];
        {% endfor -%}
        {% endfor %}

        /* Read top window row */
        // Retrieve top left and middle pack_t from line buffers
        {% for i in range(winsize[0]-1) -%}
        win[{{i}}][{{winsize[1]-1}}] = lines[{{i}}][xx1+{{uhalo[1]}}];
        {% endfor %}

        // Retrieve top right pack_t from either previous timestep or input FIFO.
        if (read_fifo) {
          input = READ(inputs);
        } else {
          input = BUFF(tt-1,xx0,xx1);
        }

        win[{{winsize[0]-1}}][{{winsize[1]-1}}] = input;

        // Shift lines right and store input.
        {% for i in range(winsize[0]-2) -%}
        lines[{{i}}][xx1+{{uhalo[1]}}] = lines[{{i+1}}][xx1+{{uhalo[1]}}];
        {% endfor -%}
        lines[{{winsize[0]-2}}][xx1+{{uhalo[1]}}] = input;

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

void compute(int ntiles, fpack_t &inputs, fpack_t &outputs) {
  for (int i=0; i<ntiles; i++) {
    compute_tile(inputs, outputs);
  }

  ASSERT_EMPTY(inputs);
}

#undef BUFF

/********************************/
/************** LOAD ************/
/********************************/
// void load_0(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in,
//           fdata_t &inputs_0,
//           fdata_t &inputs_1)
// {
//   data_t *B0 = arr_in;
//   data_t *B1 = B0 + B_P0_SIZE(T,N);

//   for (int n=0; n<ntiles; n++) {
//     coord_t coord = coords.read();
//     load_B0(T, N, coord.t, coord.x0, coord.x1, B0, inputs_0);
//     load_B1(T, N, coord.t, coord.x0, coord.x1, B1, inputs_1);
//   }
// }

// void load_1(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in,
//           fdata_t &inputs_2,
//           fdata_t &inputs_aux)
// {
//   data_t *B0 = arr_in;
//   data_t *B1 = B0 + B_P0_SIZE(T,N);
//   data_t *B2 = B1 + B_P1_SIZE(T,N);
//   data_t *Baux = B2 + B_P2_SIZE(T,N);

//   for (int n=0; n<ntiles; n++) {
//     coord_t coord = coords.read();
//     load_Baux(T, N, coord.t, coord.x0, coord.x1, Baux, inputs_aux);
//     load_B2(T, N, coord.t, coord.x0, coord.x1, B2, inputs_2);
//   }
// }

void load(int ntiles, int T, int N, fcoord_t &coords, data_t *arr_in,
          fdata_t &inputs_0,
          fdata_t &inputs_1,
          fdata_t &inputs_2,
          fdata_t &inputs_aux)
{
  data_t *B0 = arr_in;
  data_t *B1 = B0 + B_P0_SIZE(T,N);
  data_t *B2 = B1 + B_P1_SIZE(T,N);
  data_t *Baux = B2 + B_P2_SIZE(T,N);

  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();
    load_Baux(T, N, coord.t, coord.x0, coord.x1, Baux, inputs_aux);
    load_B0(T, N, coord.t, coord.x0, coord.x1, B0, inputs_0);
    load_B2(T, N, coord.t, coord.x0, coord.x1, B2, inputs_2);
    load_B1(T, N, coord.t, coord.x0, coord.x1, B1, inputs_1);
  }
}

void load_B0(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_0) {
  data_t *base = &B_P0_acc(T, N, arr_in, x0, x1, x2, 0, -HALO_1, 0);
  int count = 0;

  // Assuming HALO_0=1, inputs are already in lexicographic order.
  for (int x=0; x<(S1+HALO_1)*S2; x++) {
#pragma HLS PIPELINE II=1
    inputs_0.write(base[x]);
  }
}

void load_B1(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_1) {
  data_t *base = &B_P1_acc(T, N, arr_in, x0, x1, x2, 0, 0, -HALO_2);
  // data_t tmp[S0][HALO_1][S2+HALO_2];
  // int count = 0;

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
    inputs_1.write(base[x]);
  }

//   for (int t=0; t<S0; t++) {
//     for (int i=0; i<HALO_1; i++) {
//       for (int j=0; j<S2+HALO_2; j++) {
// #pragma HLS PIPELINE II=1
//         inputs_1.write(tmp[t][i][j]);
//       }
//     }
//   }
}

void load_B2(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_2) {
  data_t *base = &B_P2_acc(T, N, arr_in, x0, x1, x2, -HALO_0, 0, 0);
  int count = 0;

  // Inputs are already in lexicographic order.
  for (int x=0; x<(S0+HALO_0)*S1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
    inputs_2.write(base[x]);
  }
}

void load_Baux(int T, int N, int x0, int x1, int x2, data_t *arr_in, fdata_t &inputs_aux) {
  data_t *base = &B_AUX_acc(T, N, arr_in, x0, x1, x2, 0);
  int count = 0;

  // Inputs are already in lexicographic order.
  for (int x=0; x<HALO_0*HALO_1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
    inputs_aux.write(base[x]);
  }
}

/********************************/
/************* STORE ************/
/********************************/
// void store_0(int ntiles, int T, int N, fcoord_t &coords,
//            fdata_t &outputs_0,
//            fdata_t &outputs_1,
//            data_t *arr_out) {

//   data_t *B0 = arr_out;
//   data_t *B1 = B0 + B_P0_SIZE(T,N);

//   for (int n=0; n<ntiles; n++) {
//     coord_t coord = coords.read();

//     store_B1(T, N, coord.t, coord.x0, coord.x1, outputs_1, B1);
//     store_B0(T, N, coord.t, coord.x0, coord.x1, outputs_0, B0);
//   }

//   ASSERT_EMPTY(outputs_0);
//   ASSERT_EMPTY(outputs_1);
// }

// void store_1(int ntiles, int T, int N, fcoord_t &coords,
//            fdata_t &outputs_2,
//            fdata_t &outputs_aux,
//            data_t *arr_out) {

//   data_t *B0 = arr_out;
//   data_t *B1 = B0 + B_P0_SIZE(T,N);
//   data_t *B2 = B1 + B_P1_SIZE(T,N);
//   data_t *Baux = B2 + B_P2_SIZE(T,N);

//   for (int n=0; n<ntiles; n++) {
//     coord_t coord = coords.read();

//     store_B2(T, N, coord.t, coord.x0, coord.x1, outputs_2, B2);
//     store_Baux(T, N, coord.t, coord.x0, coord.x1, outputs_aux, Baux);
//   }

//   ASSERT_EMPTY(outputs_2);
//   ASSERT_EMPTY(outputs_aux);
// }

void store(int ntiles, int T, int N, fcoord_t &coords,
           fdata_t &outputs_0,
           fdata_t &outputs_1,
           fdata_t &outputs_2,
           fdata_t &outputs_aux,
           data_t *arr_out) {

  data_t *B0 = arr_out;
  data_t *B1 = B0 + B_P0_SIZE(T,N);
  data_t *B2 = B1 + B_P1_SIZE(T,N);
  data_t *Baux = B2 + B_P2_SIZE(T,N);

  for (int n=0; n<ntiles; n++) {
    coord_t coord = coords.read();

    store_B2(T, N, coord.t, coord.x0, coord.x1, outputs_2, B2);
    store_B1(T, N, coord.t, coord.x0, coord.x1, outputs_1, B1);
    store_B0(T, N, coord.t, coord.x0, coord.x1, outputs_0, B0);
    store_Baux(T, N, coord.t, coord.x0, coord.x1, outputs_aux, Baux);
  }

  ASSERT_EMPTY(outputs_0);
  ASSERT_EMPTY(outputs_1);
  ASSERT_EMPTY(outputs_2);
  ASSERT_EMPTY(outputs_aux);
}

void store_B0(int T, int N, int x0, int x1, int x2, fdata_t &outputs_0, data_t *arr_out) {
  data_t *base = &B_P0_acc(T, N, arr_out, x0+1, x1, x2, 0, 0, 0);

  // Assuming HALO_0=1, outputs are already in the correct order.
  for (int x=0; x<S1*S2; x++) {
#pragma HLS PIPELINE II=1
    base[x] = READ(outputs_0);
  }
}

void store_B1(int T, int N, int x0, int x1, int x2, fdata_t &outputs_1, data_t *arr_out) {
  data_t *base = &B_P1_acc(T, N, arr_out, x0, x1+1, x2, 0, 0, 0);

//   data_t tmp[S0][HALO_1][S2];

//   // Outputs are in lexicographic order, while the required order is j,t,i We
//   // store outputs in lexicographic order before storing them in the correct
//   // order in a second swipe.
//   for (int t=0; t<S0; t++) {
//     for (int i=0; i<HALO_1; i++) {
//       for (int j=0; j<S2; j++) {
// #pragma HLS PIPELINE II=1
//         tmp[t][i][j] = READ(outputs_1);
//       }
//     }
//   }

  // See load_B1 for the rationale of using a single counter instead of a loop
  // nest.
  for (int x=0; x<S2*S0*HALO_1; x++) {
#pragma HLS PIPELINE II=1
    int j = x / (S0*HALO_1);
    int r = x % (S0*HALO_1);
    int t = r / HALO_1;
    int i = r % HALO_1;
    base[x] = READ(outputs_1);
  }
}

void store_B2(int T, int N, int x0, int x1, int x2, fdata_t &outputs_2, data_t *arr_out) {
  data_t *base = &B_P2_acc(T, N, arr_out, x0, x1, x2+1, 0, 0, 0);

  // Outputs are already in the correct order.
  for (int x=0; x<S0*S1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
    base[x] = READ(outputs_2);
  }
}

void store_Baux(int T, int N, int x0, int x1, int x2, fdata_t &outputs_aux, data_t *arr_out) {
  data_t *base = &B_AUX_acc(T, N, arr_out, x0, x1, x2, 0);

  // Outputs are already in the correct order.
  for (int x=0; x<HALO_0*HALO_1*HALO_2; x++) {
#pragma HLS PIPELINE II=1
    base[x] = READ(outputs_aux);
  }
}

/***********************************/
/*********** UNSHUFFLE ***************/
/***********************************/
void unshuffle_B1(int ntiles, fdata_t &in, fdata_t &out) {
  data_t tmp[2][S0][HALO_1][S2+HALO_2];
  // Using dual-port BRAM for simultaneous R/W
#pragma HLS RESOURCE variable=tmp core=RAM_2P

  int j_in, t_in, i_in;
  int t_out, i_out, j_out;

  for (int n=0; n<ntiles+1; n++) {
    for (int x=0; x<(S2+HALO_2)*S0*HALO_1; x++) {
#pragma HLS PIPELINE II=1
      // Asserting false dependency, as pipeline depth is much smaller than
      // buffer size.
#pragma HLS DEPENDENCE variable=tmp inter=false

      // Compute input coordinates
      bool t_in_max = t_in == S0-1;
      bool i_in_max = i_in == HALO_1-1;

      j_in = (x == 0) ? 0 : (t_in_max && i_in_max ? j_in+1 : j_in);
      t_in = (x == 0 || (t_in_max && i_in_max)) ? 0 : (i_in_max ? t_in+1 : t_in);
      i_in = (x == 0 || i_in_max) ? 0 : i_in+1;

#ifdef DEBUG
      // Verify
      int j_in2 = x / (S0*HALO_1);
      int r_in2 = x % (S0*HALO_1);
      int t_in2 = r_in2 / HALO_1;
      int i_in2 = r_in2 % HALO_1;

      assert(j_in == j_in2);
      assert(t_in == t_in2);
      assert(i_in == i_in2);
#endif

      // Compute output coordinates

      bool i_out_max = i_out == HALO_1-1;
      bool j_out_max = j_out == S2+HALO_2-1;

      t_out = (x == 0) ? 0 : (i_out_max && j_out_max ? t_out+1 : t_out);
      i_out = (x == 0 || (i_out_max && j_out_max)) ? 0 : (j_out_max ? i_out+1 : i_out);
      j_out = (x == 0 || j_out_max) ? 0 : j_out+1;

#ifdef DEBUG
      int t_out2 = x / (HALO_1*(S2+HALO_2));
      int r_out2 = x % (HALO_1*(S2+HALO_2));
      int i_out2 = r_out2 / (S2+HALO_2);
      int j_out2 = r_out2 % (S2+HALO_2);

      assert(t_out == t_out2);
      assert(i_out == i_out2);
      assert(j_out == j_out2);
#endif

      if (n<ntiles)
        tmp[(n+1)%2][t_in][i_in][j_in] = READ(in);
      if (n>0)
        out.write(tmp[n%2][t_out][i_out][j_out]);
    }
  }

  ASSERT_EMPTY(in);
}

/***********************************/
/*********** SHUFFLE *************/
/***********************************/
void shuffle_B1(int ntiles, fdata_t &in, fdata_t &out) {
  data_t tmp[2][S0][HALO_1][S2];
  // Using dual-port BRAM for simultaneous R/W
#pragma HLS RESOURCE variable=tmp core=RAM_2P

  for (int n=0; n<ntiles+1; n++) {
    for (int x=0; x<S2*S0*HALO_1; x++) {
#pragma HLS PIPELINE II=1
      // Asserting false dependency, as pipeline depth is much smaller than
      // buffer size.
#pragma HLS DEPENDENCE variable=tmp inter=false
      // Compute input coordinates
      int t_in = x / (HALO_1*S2);
      int r_in = x % (HALO_1*S2);
      int i_in = r_in / S2;
      int j_in = r_in % S2;

      // Compute output coordinates
      int j_out = x / (S0*HALO_1);
      int r_out = x % (S0*HALO_1);
      int t_out = r_out / HALO_1;
      int i_out = r_out % HALO_1;

      if (n<ntiles)
        tmp[(n+1)%2][t_in][i_in][j_in] = READ(in);
      if (n>0)
        out.write(tmp[n%2][t_out][i_out][j_out]);
    }
  }

  ASSERT_EMPTY(in);
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

  int offset = UNROLL - HALO_2 % UNROLL;

  // Total number of iterations in a tile.
  int per_tile = HALO_0*(S1+HALO_1)*(S2+HALO_2)+S0*(HALO_1*(S2+HALO_2)+S1*HALO_2);
  int t=-HALO_0, i=-HALO_1, j=-HALO_2;

  for (int x=0; x<ntiles*per_tile; x++) {
#pragma HLS PIPELINE II=1
    int index = (offset + j + HALO_2) % UNROLL;
    data_t v;

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

    // Compute new loop iterators.
    int jbound = (t < 0 || i < 0) ? S2 : 0;

    bool next_line = (j == jbound-1);
    bool next_plane = (i==S1-1 && j==jbound-1);
    bool next_tile = (t==S0-1 && i==S1-1 && j==jbound-1);

    if (next_tile) {
      t=-HALO_0; i=-HALO_1; j=-HALO_2;
    } else {
      if (next_plane) {
        t = t+1;
        i=-HALO_1;
        j=-HALO_2;
      } else {
        if (next_line) {
          i=i+1;
          j=-HALO_2;
        } else {
          j=j+1;
        }
      }
    }
  }


//   for (int n=0; n<ntiles; n++) {
//     for (int t=-HALO_0; t<S0; t++) {
//       for (int i=-HALO_1; i<S1; i++) {
//         for (int j=-HALO_2; j<(t<0 || i<0 ? S2 : 0); j++) {
// #pragma HLS PIPELINE II=1
//           int index = (offset + j + HALO_2) % UNROLL;
//           data_t v;

//           if (t<0) {
//             if (j<0) {
//               if (i<0) v = READ(inputs_aux);
//               else v = READ(inputs_2);
//             } else v = READ(inputs_0);
//           } else {
//             if (i<0) v = READ(inputs_1);
//             else v = READ(inputs_2);
//           }

//           input.vals[index] = v;
//           if (t < S0-1 && index == UNROLL-1) inputs.write(input);
//         }
//       }
//     }
//   }

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

  int per_tile=S2*(HALO_0*S1+(S0-HALO_0)*HALO_1)+HALO_2*(S0-HALO_0)*(S1-HALO_1);

  int t=0, i=0, j=(t>=S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;

  for (int x=0; x<ntiles*per_tile; x++) {
#pragma HLS PIPELINE II=1
    int index  = j % UNROLL;
    int startj = (t>=S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;

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

    // Compute new loop iterators
    bool next_line = (j==S2-1);
    bool next_plane = (j==S2-1 && i==S1-1);
    bool next_tile = (t== S0-1 && j==S2-1 && i==S1-1);

    if (next_tile) {
      t=0; i=0; j=(t>=S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;
    } else {
      if (next_plane) {
        t=t+1;
        i=0;
        j=(t>=S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;
      } else {
        if (next_line) {
          i=i+1;
          j=(t>=S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;
        } else {
          j = j+1;
        }
      }
    }
  }

//   for (int n=0; n < ntiles; n++) {
//     for (int t=0; t<S0; t++) {
//       for (int i=0; i<S1; i++) {
//         int startj = (t>=S0-HALO_0 || i>=S1-HALO_1) ? 0 : S2-HALO_2;
//         for (int j=startj; j<S2; j++) {
// #pragma HLS PIPELINE II=1

//           int index = j % UNROLL;
//           if (j == startj || index == 0) output = READ(outputs);
//           data_t v = output.vals[index];

//           bool i1border = (i  >= S1-HALO_1);
//           bool i2border = (j  >= S2-HALO_2);

//           if (t >= S0-HALO_0) {
//             outputs_0.write(v);

//             if (i1border && i2border) {
//               outputs_aux.write(v);
//             }

//             if (i1border) {
//               outputs_1.write(v);
//             }
//             if (i2border) {
//               outputs_2.write(v);
//             }
//           } else {
//             if (i1border) {
//               outputs_1.write(v);
//             }

//             if (i2border) {
//               outputs_2.write(v);

//             }
//           }
//         }
//       }
//     }
//   }

  ASSERT_EMPTY(outputs);
}
