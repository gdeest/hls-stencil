#include "jacobi2d.h"


int in_array(int T, int N0, int N1, int t, int x0, int x1) {
  return
    (t  >= -1) && (t < T)   &&
    (x0 >=  0) && (x0 < N0) &&
    (x1 >=  0) && (x1 < N1);
}

void unskew_tile(int t, int x0, int x1, int *_t, int *_x0, int *_x1) {
  *_t = t - x0 - x1;
  *_x0 = x0;
  *_x1 = x1;
}

void unskew_point(int t, int x0, int x1, int *_t, int *_x0, int *_x1) {
  *_t = t;
  *_x0 = x0 - t;
  *_x1 = x1 - t;
}
