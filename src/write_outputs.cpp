#include "write_outputs.h"
#include "access_macros.h"
#include "params.h"
#include <string.h>

void write_outputs(int ntiles, hls::stream<tile_coord_t> &coords, data_t *arr, int T, int N0, int N1, foutput_t &outputs) {
  data_t tmp[S1+2]; // Size = longest burst.

  for (int i=0; i<ntiles; i++) {
    tile_coord_t coord = coords.read();
    int t  = coord.t;
    int x0 = coord.x0;
    int x1 = coord.x1;

    // First ST-1 timesteps
    for (int tt=0; tt<ST-1; tt++) {
      // First S0-2 lines
      for (int xx0=0; xx0<S0-2; xx0++) {
        // Discard S1/2-1 outputs.
//         for (int xx1=0; xx1<S1-2; xx1+=2)
// #pragma HLS PIPELINE II=1
//           output_t garbage = outputs.read();

        // Read end of the line and write in tmp.
        output_t out = outputs.read();
        tmp[0] = out.v[0];
        tmp[1] = out.v[1];

        // Copy to main memory.
        memcpy(&AT_TILED(arr, T, N0, N1, t, x0, x1, tt, xx0, S1-2), tmp, 2*sizeof(data_t));
      }

      // Last 2 lines
      for (int xx0=S0-2; xx0<S0; xx0++) {
        for (int xx1=0; xx1<S1; xx1+=2) {
#pragma HLS PIPELINE II=1
          output_t out = outputs.read();
          tmp[xx1  ] = out.v[0];
          tmp[xx1+1] = out.v[1];
        }

        // Copy to main memory.
        memcpy(&AT_TILED(arr, T, N0, N1, t, x0, x1, tt, xx0, 0), tmp, S1*sizeof(data_t));
      }
    }

    // Last timestep - copy everything.
    for (int xx0=0; xx0<S0; xx0++) {
      for (int xx1=0; xx1<S1; xx1+=2) {
#pragma HLS PIPELINE II=1
        output_t out = outputs.read();
        tmp[xx1  ] = out.v[0];
        tmp[xx1+1] = out.v[1];
      }

      memcpy(&AT_TILED(arr, T, N0, N1, t, x0, x1, ST-1, xx0, 0), tmp, S1*sizeof(data_t));
    }
  }
}
