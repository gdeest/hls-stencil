#include "jacobi2d.h"
#include "utils.h"


void outer_tile(data_t *arr, int T, int N0, int N1, int t, int x0, int x1) {

  int t_orig  = ST * t,
    x0_orig = S0 * x0,
    x1_orig = S1 * x1;

  for (int tt=0; tt < ST; tt++) {
    for (int xx0=0; xx0 < S0; xx0++) {
      for (int xx1=0; xx1 < S1; xx1++) {
        int _t, _x0, _x1;
        unskew_point(t_orig + tt, x0_orig + xx0, x1_orig + xx1, &_t, &_x0, &_x1);

        if ((_t >= 0) && in_array(T, N0, N1, _t, _x0, _x1)) {
          data_t d0, d1, d2, d3, d4;

#ifdef VALIDATE
          assert(!arr[AT(T, N0, N1, _t, _x0, _x1)]);
#endif

          // d0
          if (in_array(T, N0, N1, _t-1, _x0+1, _x1)) {
            d0 = arr[AT(T, N0, N1, _t-1, _x0+1, _x1)];
#ifdef VALIDATE
            assert(d0 == 1);
#endif
          } else d0 = 0;

          // d1
          if (in_array(T, N0, N1, _t-1, _x0, _x1+1)) {
            d1 = arr[AT(T, N0, N1, _t-1, _x0, _x1+1)];
#ifdef VALIDATE
            assert(d1 == 1);
#endif
          } else d1 = 0;

          // d2
          if (in_array(T, N0, N1, _t-1, _x0, _x1)) {
            d2 = arr[AT(T, N0, N1, _t-1, _x0, _x1)];
#ifdef VALIDATE
            assert(d2 == 1);
#endif
          } else d2 = 0;

          // d3
          if (in_array(T, N0, N1, _t-1, _x0, _x1-1)) {
            d3 = arr[AT(T, N0, N1, _t-1, _x0, _x1-1)];
#ifdef VALIDATE
            assert(d3 == 1);
#endif
          } else d3 = 0;

          // d4
          if (in_array(T, N0, N1, _t-1, _x0-1, _x1)) {
            d4 = arr[AT(T, N0, N1, _t-1, _x0-1, _x1)];
#ifdef VALIDATE
            // std::cout << "(tt,xx0,xx1) = (" << tt << "," << xx0 << "," << xx1 << ")" << std::endl;
            assert(d4 == 1);
#endif
          } else d4 = 0;

#ifdef VALIDATE
          arr[AT(T, N0, N1, _t, _x0, _x1)] = 1;
#else
          arr[AT(T, N0, N1, _t, _x0, _x1)] =
            COMPUTE(d0, d1, d2, d3, d4);
#endif

        }
      }

    }
  }
}
