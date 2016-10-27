#include "assert.h"
#include "compute.h"


#define SHIFT_REGISTERS

#ifdef SHIFT_REGISTERS
#define REUSE(d0,d1) reuse[(d0)+2][(xx1_)+(d1)+2]
#else
#define REUSE(d0,d1) reuse[((xx0_)+(d0)+2)%3][(xx1_+(d1))+2]
#endif

// #define REUSE(xx0,xx1) prev[prev_tt][(xx0)+2][(xx1)+2]
#define PREV(xx0,xx1) prev[prev_tt][(xx0)+2][(xx1)+2]
// #define REUSE(xx0,xx1) prev[prev_tt][(xx0)+2][(xx1)+2]
#define FIFO_NAME(i) inputs ## i
#define READ_FIFO(i) deps[i] = FIFO_NAME(i).read()

void compute_timestep(int tt,
                      finput_t &inputs0,
                      finput_t &inputs1,
                      finput_t &inputs2,
                      finput_t &inputs3,
                      finput_t &inputs4,
                      finput_t &inputs5,
                      finput_t &inputs6,
                      finput_t &inputs7,
                      data_t prev[2][S0+2][S1+2],
                      foutput_t &output) {

  data_t reuse[3][S1+2];

#ifdef SHIFT_REGISTERS
#pragma HLS ARRAY_PARTITION variable=reuse dim=0 complete
#else
#pragma HLS RESOURCE variable=reuse core=RAM_T2P_BRAM
#pragma HLS ARRAY_PARTITION variable=reuse dim=1 complete
#pragma HLS ARRAY_PARTITION variable=reuse dim=2 cyclic   factor=4
#endif

  data_t deps[8];
#pragma HLS ARRAY_PARTITION variable=deps complete

  int curr_tt = tt%2;
  int prev_tt = (tt+1)%2;

  assert(curr_tt != prev_tt);

  for (int xx0=0; xx0<S0; xx0++) {
    for (int xx1=0; xx1<S1; xx1+=2) {
      #pragma HLS PIPELINE II=1

#pragma HLS DEPENDENCE variable=prev inter=false
#pragma HLS DEPENDENCE variable=reuse inter=false
#pragma HLS DEPENDENCE variable=deps inter=false
      // for (int i=0; i<8; i+=2) {
      // Each iteration of this loops performs 2 computations.
      // #pragma HLS UNROLL
      // int xx0_ = xx0, xx1_ = xx1+i;
      int xx0_ = xx0, xx1_ = xx1;

      data_t v0;
      data_t v1;

      #ifdef SHIFT_REGISTERS
      if (xx1_==0) {
        for (int i=0; i<2; i++) {
          #pragma HLS UNROLL
          for (int j=0; j<S1+1; j++) {
            #pragma HLS UNROLL
            reuse[i][j] = reuse[i+1][j];
          }
        }
      }
      #endif

      // Fetch dependencies
      // No need to store those either, they will never be reused.
      if (xx0_ == 0) {
        READ_FIFO(6);
        READ_FIFO(7);
      } else {
        deps[6] = REUSE(-2, 0);
        deps[7] = REUSE(-2, -1);
      }

      {
#pragma HLS LATENCY max=3
        if (xx1_ == 0) {                                                                  // C3
          if (xx0_==0) {
            READ_FIFO(4);
            REUSE(-1,-1) = deps[4];
          } else {
            deps[4] = REUSE(-1,-1);
          }
          // deps[5] = inputs5.read(); // No need to store it, it will never be reused
          READ_FIFO(5);
        } else {
          deps[5] = deps[3];
          deps[4] = deps[2];
        }
      }

      {
#pragma HLS LATENCY max=3
        if (xx0_==0) {
          READ_FIFO(3);                                                                   // C2
          REUSE(-1,0) = deps[3];
        } else {
          deps[3] = REUSE(-1,0);                                                   // C1
        }
      }

      {
#pragma HLS LATENCY max=3
        if (xx0_==0 || (tt==0 && xx1_==S1-2)) {
          READ_FIFO(2);
          REUSE(-1,+1) = deps[2];
        } else {
          data_t tmp;
          if (xx1_ == S1-2) {
            tmp = PREV(xx0_-1,xx1_+1);
            REUSE(-1,1) = tmp;
          } else {
            tmp = REUSE(-1,1);
          }
          deps[2] = tmp;
        }
      }

      {
#pragma HLS LATENCY max=3
        if (tt==0 || xx1_==0) {
          READ_FIFO(1);
          REUSE(0,-1) = deps[1];
        } else {
          data_t tmp = PREV(xx0_,xx1_-1);
          deps[1] = tmp;
          REUSE(0,-1) = tmp;
        }
      }

      {
#pragma HLS LATENCY max=3
        if (tt==0) {
          READ_FIFO(0);
          REUSE(0,0) = deps[0];
        } else {
          data_t tmp = PREV(xx0_, xx1_);
          deps[0] = tmp;
          REUSE(0,0) = tmp;
        }
      }

      // Compute two values
      v0 = COMPUTE(deps[1], deps[3], deps[4], deps[5], deps[7]);
      v1 = COMPUTE(deps[0], deps[2], deps[3], deps[4], deps[6]);

      output_t out;
      out.v[0] = v0;
      out.v[1] = v1;

      // deps[5] = deps[3];
      // deps[4] = deps[2];

      prev[curr_tt][xx0_+2][xx1_+2]   = v0;
      prev[curr_tt][xx0_+2][xx1_+1+2] = v1;

      if (tt==ST-1 || xx0_>=S0-2 || xx1_==S1-2)
        output.write(out);
    }
  }
}

#undef READ_FIFO
#undef FIFO_NAME
#undef PREV
#undef REUSE

void compute(int ntiles,
             finput_t &inputs0,
             finput_t &inputs1,
             finput_t &inputs2,
             finput_t &inputs3,
             finput_t &inputs4,
             finput_t &inputs5,
             finput_t &inputs6,
             finput_t &inputs7,
             foutput_t &output) {
  // Array to store previous results.
  data_t prev[2][S0+2][S1+2];
#pragma HLS RESOURCE variable=prev core=RAM_T2P_BRAM
#pragma HLS ARRAY_PARTITION variable=prev cyclic   dim=1 factor=2
#pragma HLS ARRAY_PARTITION variable=prev cyclic   dim=2 factor=2
#pragma HLS ARRAY_PARTITION variable=prev dim=3 cyclic factor=2
  // PRAGMA_ARRAY_PARTITION(prev, cyclic, 3+(UNROLL_FACTOR-1), 3)

  for (int n=0; n<ntiles; n++) {
    for (int tt=0; tt<ST; tt++) {
      compute_timestep(tt,
                       inputs0,
                       inputs1,
                       inputs2,
                       inputs3,
                       inputs4,
                       inputs5,
                       inputs6,
                       inputs7,
                       prev,
                       output);
    }
  }
}
