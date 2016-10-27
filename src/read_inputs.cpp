#include "read_inputs.h"
#include "access_macros.h"
#include "params.h"
#include "access_macros.h"

#define FIFO_NAME(i) inputs ## i
#define WRITE(tag, i, x) FIFO_NAME(i).write(x)

void loop1(data_t *base_addr,
           finput_t &inputs6,
           finput_t &inputs7)
{
  #pragma HLS INTERFACE m_axi port=base_addr

  for (int i=0; i<S1; i+=1) {
#pragma HLS PIPELINE II=1
    data_t v = base_addr[i];
    if (i%2 == 0)  WRITE("0", 7, v);
    if (i%2 == 1)  WRITE("1", 6, v);
  }
}

void loop2(data_t *base_addr,
           finput_t &inputs2,
           finput_t &inputs3,
           finput_t &inputs4,
           finput_t &inputs5)
{
#pragma HLS INTERFACE m_axi port=base_addr

  for (int i=0; i<S1+2; i++) {
#pragma HLS PIPELINE II=1
    data_t v = base_addr[i];
    if (i == 0)           WRITE("2", 5, v);
    if (i == 1)           WRITE("3", 4, v);
    if (i>1 && i%2 == 0)  WRITE("4", 3, v);
    if (i>1 && i%2 == 1)  WRITE("5", 2, v);
  }
}

void loop3(data_t *base_addr,
           int xx0,
           finput_t &inputs0,
           finput_t &inputs1,
           finput_t &inputs2,
           finput_t &inputs5)
{
// #pragma HLS LATENCY min=50
#pragma HLS INTERFACE m_axi port=base_addr

  for (int i=0; i<S1+2; i++) {
#pragma HLS PIPELINE II=1
    data_t v = base_addr[i];
    if (i==0    &&  xx0<S0-1)  WRITE("6", 5, v);
    if (i%2==1  &&  i!=S1+1)   WRITE("7", 1, v);
    if (i>0     &&  i%2==0)    WRITE("8", 0, v);
    if (i==S1+1 &&  xx0<S0-1)  WRITE("9", 2, v);
  }
}

void loop4(data_t *base_addr,
           int xx0,
           finput_t &inputs1,
           finput_t &inputs5)
{
// #pragma HLS LATENCY min=50
#pragma HLS INTERFACE m_axi port=base_addr

  for (int i=0; i<2; i++) {
#pragma HLS PIPELINE II=1
    data_t v = base_addr[i];
    if (i==0    &&  xx0<S0-1)  WRITE("10", 5, v);
    if (i==1)                  WRITE("11", 1, v);
  }
}

void read_inputs(int ntiles, hls::stream<tile_coord_t> &coords, data_t *arr, int T, int N0, int N1,
                 finput_t &inputs0,
                 finput_t &inputs1,
                 finput_t &inputs2,
                 finput_t &inputs3,
                 finput_t &inputs4,
                 finput_t &inputs5,
                 finput_t &inputs6,
                 finput_t &inputs7
                 ) {
  data_t tmp[S1+2];

  int tt_offset  = N0*N1-N1-1;
  int xx0_offset = N1;
  int xx1_offset = 1;

  int offset_1   = -2*xx0_offset - xx1_offset;
  int offset_2   = -1*xx0_offset - 2*xx1_offset;
  int offset_3   = -2*xx1_offset;

  for (int n=0; n<ntiles; n++) {
    tile_coord_t coord = coords.read();
    int t  = coord.t;
    int x0 = coord.x0;
    int x1 = coord.x1;

    int tt=-1;

    data_t *timestep_addr = &AT_TILED(arr, T, N0, N1, t, x0, x1, -1, 0, 0);
    for (int tt=-1; tt<ST-1; tt++) {

      data_t *addr_m2, *addr_m1, *addr;
      {
#pragma HLS LATENCY max=1
        addr_m2 = timestep_addr + offset_1;
        addr_m1 = timestep_addr + offset_2;
        addr    = timestep_addr + offset_3;
      }

      // -2
      loop1(addr_m2, inputs6, inputs7);

      // -1
      loop2(addr_m1, inputs2, inputs3, inputs4, inputs5);

      // >= 0
      if (tt==-1) {
        for (int xx0=0; xx0<S0; xx0++) {
          loop3(addr, xx0, inputs0, inputs1, inputs2, inputs5);
          addr += N1;
        }
      } else {
        for (int xx0=0; xx0<S0; xx0++) {
          loop4(addr, xx0, inputs1, inputs5);
          addr += N1;
        }
      }
      timestep_addr += tt_offset;
    }
  }
}
#undef WRITE
#undef FIFO_NAME
