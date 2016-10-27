#include "jacobi1d.h"

#undef max
#undef min

#include <vector>

#define INNER_TILE(t, x0) inner_tile_wrapper(arr, T, N0, (t), (x0), mutex_cout, cnt)
template<> void scan_hw_tiles<16,16>(int T,int N0, data_t *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt) {
if (T >= 1 && N0 >= 1)
  for (int c0 = 0; c0 <= (2 * T + N0 - 3) / 16; c0 += 1)
    if ((8 * c0 + 1 >= T && 2 * T + N0 >= ((T - 1) % 16) + 16 * c0 + 3) || (T >= 8 * c0 + 2 && T + N0 >= 8 * (c0 % 2) + 8 * c0 + 2 && N0 + 14 >= 16 * (c0 % 2)))
      {
        int d = c0;
        sem_wait(mutex_cout);
        cout << "HW blocked on barrier: " << (d&1) << endl;
sem_post(mutex_cout);
        pthread_barrier_wait(&barriers[d & 1]);
        sem_wait(mutex_cout);
        cout << "HW unblocked." << endl;
        sem_post(mutex_cout);

        int ntiles = 0;
        std::vector<int> coords;

        for (int c0 = max(floord(d, 2) + 1, d + floord(-T - 1, 16) + 2); c0 < min(d, floord(16 * d + N0 + 15, 32)); c0 += 1) {
          coords.push_back(d - c0);
          coords.push_back(c0);
          ntiles++;
          /* INNER_TILE(d - c0, c0); */
        }

        if (ntiles > 0) {
          int *coords_ = (int *) ALLOC(2 * ntiles * sizeof(int));
          for (int i=0; i<2*ntiles; i++)
            coords_[i] = coords[i];
          inner_tiles_wrapper(arr, T, N0, ntiles, coords_, mutex_cout, cnt);
          FREE(coords_);
        }
      }

}

#undef INNER_TILE
#define OUTER_TILE(t, x0) outer_tile_wrapper(arr, T, N0, (t), (x0), mutex_cout, cnt)
template<> void scan_sw_tiles<16,16>(int T,int N0, data_t *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt) {
if (T >= 1 && N0 >= 1)
  for (int c0 = 0; c0 <= (2 * T + N0 - 3) / 16; c0 += 1)
    if ((8 * c0 + 1 >= T && 2 * T + N0 >= ((T - 1) % 16) + 16 * c0 + 3) || (T >= 8 * c0 + 2 && T + N0 >= 8 * (c0 % 2) + 8 * c0 + 2 && N0 + 14 >= 16 * (c0 % 2)))
      {
        int d = c0;
        sem_wait(mutex_cout);
        cout << "SW blocked on barrier: " << (d&1) << endl;
sem_post(mutex_cout);
        pthread_barrier_wait(&barriers[d & 1]);
        sem_wait(mutex_cout);
        cout << "SW unblocked." << endl;
        sem_post(mutex_cout);
        {
          if (d >= 2 && T >= 8 * d + 16 && N0 >= 1 && d % 2 == 0) {
            OUTER_TILE(d / 2, d / 2);
          } else if (T >= 17 && N0 >= 1 && d - 2 * ((T + 15) / 16) >= -2 && 2 * T + N0 >= ((T - 1) % 16) + 16 * d + 3 && (T - 1) % 16 <= 14)
            OUTER_TILE((T + 16) / 16 - 1, d - (T + 16) / 16 + 1);
          if (d >= 0 && T >= 1 && T + N0 >= 16 * d + 2 && N0 >= 1 && N0 + 14 >= 16 * d) {
            OUTER_TILE(0, d);
          } else if (16 * d >= N0 + 16 && 2 * floord(16 * d + N0 + 15, 32) >= d + 1 && T + 16 * ((16 * d + N0 + 15) / 32) >= 16 * d + 16 && (16 * d - N0 - 16) % 32 <= 30)
            OUTER_TILE(d - (16 * d + N0 + 14) / 32, (16 * d + N0 + 14) / 32);
        }

      }

}

#undef OUTER_TILE
template<> void scan_wavefronts<16,16>(int T,int N0, pthread_barrier_t barriers[2], sem_t *mutex_cout) {
if (T >= 1 && N0 >= 1)
  for (int c0 = 0; c0 <= (2 * T + N0 - 3) / 16; c0 += 1)
    if ((8 * c0 + 1 >= T && 2 * T + N0 >= ((T - 1) % 16) + 16 * c0 + 3) || (T >= 8 * c0 + 2 && T + N0 >= 8 * (c0 % 2) + 8 * c0 + 2 && N0 + 14 >= 16 * (c0 % 2)))
      {
        int d = c0;
        pthread_barrier_wait(&barriers[d&1]);
        sem_wait(mutex_cout);
        cout << "Wavefront: " << d << endl;
        sem_post(mutex_cout);
        pthread_barrier_destroy(&barriers[d&1]);
        pthread_barrier_init(&barriers[d&1], 0, 3);
      }

}
