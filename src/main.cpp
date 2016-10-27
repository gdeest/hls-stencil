#include "jacobi2d.h"
#include "inner_tile.h"
#include "outer_tile.h"
#include "perf_counter.h"
#include "utils.h"

#include <math.h>

template<int _ST, int _S0, int _S1> void scan_hw_tiles
(int T, int N0, int N1, data_t *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt);

template<int _ST, int _S0, int _S1> void scan_sw_tiles
(int T, int N0, int N1, data_t *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt);

template<int _ST, int _S0, int _S1> void scan_wavefronts
(int T, int N0, int N1, pthread_barrier_t barriers[2], sem_t *mutex_cout);

#include "scan.h"

using namespace std;

inline int in_array(int T, int N0, int N1, int t, int x0, int x1) {
  return
    (t  >= -1) && (t < T)   &&
    (x0 >=  0) && (x0 < N0) &&
    (x1 >=  0) && (x1 < N1);
}

inline void unskew_tile(int t, int x0, int x1, int *_t, int *_x0, int *_x1) {
  *_t = t - x0 - x1;
  *_x0 = x0;
  *_x1 = x1;
}

inline void unskew_point(int t, int x0, int x1, int *_t, int *_x0, int *_x1) {
  *_t = t;
  *_x0 = x0 - t;
  *_x1 = x1 - t;
}

void inner_tiles_wrapper(data_t *arr, int T, int N0, int N1, int ntiles, int *coords, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "--- INNER TILES ---" << endl;
  for (int i=0; i<ntiles; i++) {
    cout << coords[3*i    ] << ", ";
    cout << coords[3*i + 1] << ", ";
    cout << coords[3*i + 2] << endl;
  }
  cout << "--- /INNER TILES ---" << endl;
  sem_post(mutex_cout);

#ifdef SW
  // Just call inner_tile_wrapper on each tile and let it handle the rest.
  for (int i=0; i<ntiles; i++) {
    inner_tile_wrapper(arr, T, N0, N1, coords[3*i], coords[3*i+1], coords[3*i+2], mutex_cout, cnt);
  }
#else
  // Actually call inner_tiles with multiple tile coordinates.

  // TODO: Performance counter does not really make sense when calling a varying number of tiles at once.
  cnt->start();
#ifdef __SDSCC__
  uint64_t before = sds_clock_counter();
#endif
  inner_tiles(arr, arr, T, N0, N1, ntiles, coords);
#ifdef __SDSCC__
  uint64_t after = sds_clock_counter();
  sem_wait(mutex_cout);
  cout << "MULTIPLE: " << ntiles << ", " << (after-before) << endl;
  sem_post(mutex_cout);
#endif
  cnt->stop();
#endif
}

void inner_tile_wrapper(data_t *arr, int T, int N0, int N1, int t, int x0, int x1, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "INNER TILE: " << t << ", " << x0 << ", " << x1 << endl;
  sem_post(mutex_cout);

#ifdef SW
  int
    t_orig  = ST * t,
    x0_orig = S0 * x0,
    x1_orig = S1 * x1;


  for (int tt=0; tt < ST; tt++) {
    for (int xx0=0; xx0 < S0; xx0++) {
      for (int xx1=0; xx1 < S1; xx1++) {
        int _t, _x0, _x1;
        unskew_point(t_orig + tt, x0_orig + xx0, x1_orig + xx1, &_t, &_x0, &_x1);
#ifdef VALIDATE
        assert((_t >= 0) && in_array(T, N0, N1, _t, _x0, _x1));
        assert(!arr[AT(T, N0, N1, _t, _x0, _x1)]);

        assert(in_array(T, N0, N1, _t-1, _x0-1, _x1));
        assert(in_array(T, N0, N1, _t-1, _x0+1, _x1));
        assert(in_array(T, N0, N1, _t-1, _x0, _x1-1));
        assert(in_array(T, N0, N1, _t-1, _x0, _x1+1));
        assert(in_array(T, N0, N1, _t-1, _x0, _x1));

        assert(arr[AT(T, N0, N1, _t-1, _x0-1, _x1)]);
        assert(arr[AT(T, N0, N1, _t-1, _x0+1, _x1)]);
        assert(arr[AT(T, N0, N1, _t-1, _x0, _x1-1)]);
        assert(arr[AT(T, N0, N1, _t-1, _x0, _x1+1)]);
        assert(arr[AT(T, N0, N1, _t-1, _x0, _x1)]);

        arr[AT(T, N0, N1, _t, _x0, _x1)] = 1;
#else
        arr[AT(T, N0, N1, _t, _x0, _x1)] =
          COMPUTE(arr[AT(T, N0, N1, _t-1,  _x0+1, _x1   )],
                  arr[AT(T, N0, N1, _t-1,  _x0,   _x1+1 )],
                  arr[AT(T, N0, N1, _t-1,  _x0,   _x1   )],
                  arr[AT(T, N0, N1, _t-1,  _x0,   _x1-1 )],
                  arr[AT(T, N0, N1, _t-1,  _x0-1, _x1   )]);
#endif
      }
    }
  }
#else

  int coords[3];
  coords[0] = t;
  coords[1] = x0;
  coords[2] = x1;

  cnt->start();
  inner_tiles(arr, arr, T, N0, N1, 1, coords);
  // inner_tile(arr, arr, T, N0, t, x0);
  cnt->stop();
#endif
}

void outer_tile_wrapper(data_t *arr, int T, int N0, int N1, int t, int x0, int x1, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "OUTER TILE: " << t << ", " << x0 << ", " << x1 << endl;
  sem_post(mutex_cout);

  cnt->start();

  outer_tile(arr, T, N0, N1, t, x0, x1);

  cnt->stop();
}


data_t *alloc_array(int T, int N0, int N1) {
  cout << "Allocating array of size: " << ARRAY_SIZE(T, N0, N1) << endl;
  data_t * arr = (data_t *) ALLOC(ARRAY_SIZE(T,N0,N1));

  if (!arr) {
    cerr << "Error allocating array." << endl;
    exit(-1);
  }

  return arr;
}

void init_array(data_t *arr, int T, int N0, int N1, unsigned int seed) {
  srand(seed);
  for (int x0=0; x0<N0; x0++) {
    for (int x1=0; x1<N1; x1++) {
#ifdef VALIDATE
      arr[AT(T, N0, N1, -1, x0, x1)] = 1;
#else
      // arr[AT(T, N0, N1, -1, x0, x1)] = 1;
      arr[AT(T, N0, N1, -1, x0, x1)] = (1.0*rand()/RAND_MAX)*300.0;
#endif
    }
  }

  // Initializing the array to 0 makes debugging (a lil bit) easier
  for (int t=0; t<T; t++) {
    for (int x0=0; x0<N0; x0++) {
      for (int x1=0; x1<N1; x1++) {
        arr[AT(T, N0, N1, t, x0, x1)] = 0;
      }
    }
  }
}

void jacobi2d_golden(data_t *arr, int T, int N0, int N1) {
  for (int t=0; t<T; t++) {
    for (int x0=0; x0<N0; x0++) {
      for (int x1=0; x1<N1; x1++) {
        // Update @(t,i,j)
#ifdef VALIDATE
        arr[AT(T,N0,N1,t,x0,x1)] = 1;
#else
        data_t d0, d1, d2, d3, d4;


          // d0
          if (in_array(T, N0, N1, t-1, x0+1, x1)) {
            d0 = arr[AT(T, N0, N1, t-1, x0+1, x1)];
          } else d0 = 0;

          // d1
          if (in_array(T, N0, N1, t-1, x0, x1+1)) {
            d1 = arr[AT(T, N0, N1, t-1, x0, x1+1)];
          } else d1 = 0;

          // d2
          if (in_array(T, N0, N1, t-1, x0, x1)) {
            d2 = arr[AT(T, N0, N1, t-1, x0, x1)];
          } else d2 = 0;

          // d3
          if (in_array(T, N0, N1, t-1, x0, x1-1)) {
            d3 = arr[AT(T, N0, N1, t-1, x0, x1-1)];
          } else d3 = 0;

          // d4
          if (in_array(T, N0, N1, t-1, x0-1, x1)) {
            d4 = arr[AT(T, N0, N1, t-1, x0-1, x1)];
          } else d4 = 0;

        arr[AT(T, N0, N1, t, x0, x1)] =
          COMPUTE(d0, d1, d2, d3, d4);
#endif
      }
    }
  }
}

void *sw_thread(void *_arg) {
  thread_arg *arg = (thread_arg *) _arg;
  scan_sw_tiles<ST,S0,S1>(arg->T, arg->N0, arg->N1, arg->arr, arg->barriers, arg->mutex_cout, arg->cnt);
  return 0;
}

void *hw_thread(void *_arg) {
  thread_arg *arg = (thread_arg *) _arg;
  scan_hw_tiles<ST,S0,S1>(arg->T, arg->N0, arg->N1, arg->arr, arg->barriers, arg->mutex_cout, arg->cnt);
  return 0;
}

void jacobi2d_tiled(data_t *arr, int T, int N0, int N1) {
  thread_arg arg, arg_hw, arg_sw;
  perf_counter *hw_counter, *sw_counter;
  hw_counter = new perf_counter;
  sw_counter = new perf_counter;

  pthread_barrier_t barriers[2];

  sw_counter->start();

  sem_t mutex_cout;

  sem_init(&mutex_cout, 0, 1);

  pthread_barrier_init(&barriers[0], 0, 3);
  pthread_barrier_init(&barriers[1], 0, 3);

  arg.T = T;
  arg.N0 = N0;
  arg.N1 = N1;
  arg.arr = arr;
  arg.barriers = barriers;
  arg.mutex_cout = &mutex_cout;

  arg_hw = arg;
  arg_hw.cnt = hw_counter;

  arg_sw = arg;
  arg_sw.cnt = sw_counter;

  pthread_t sw, hw;

  pthread_create(&sw, 0, sw_thread, &arg_sw);
  pthread_create(&hw, 0, hw_thread, &arg_hw);

  scan_wavefronts<ST,S0,S1>(T, N0, N1, arg.barriers, &mutex_cout);


  sem_wait(&mutex_cout);
  cout << "Waiting for threads to finish..." << endl;
  sem_post(&mutex_cout);

  pthread_join(hw, 0);
  pthread_join(sw, 0);

  cout << "Done !" << endl;
#ifdef __SDSCC__
  cout << "HW AVG: " << hw_counter->avg_cpu_cycles() << endl;
  cout << "SW AVG: " << sw_counter->avg_cpu_cycles() << endl;
#endif

  delete sw_counter;
  delete hw_counter;
}

int main(int argc, char **argv) {
  int T, N0, N1;

  data_t *A, *A_golden;

  if (argc < 3) {
    cout << "Usage: ./jacobi1d <T:uint> <N0:uint> <N1:uint>" << endl;
    exit(-1);
  } else {
    T = atoi(argv[1]);
    N0 = atoi(argv[2]);
    N1 = atoi(argv[3]);
  }

  A        = alloc_array(T, N0, N1);
  A_golden = alloc_array(T, N0, N1);

  init_array(A, T, N0, N1, SEED);
  init_array(A_golden, T, N0, N1, SEED);

  for (int x0=0; x0<N0; x0++) {
    for (int x1=0; x1<N1; x1++) {
      data_t init = A[AT(T,N0,N1,-1,x0,x1)];
      data_t init_golden = A_golden[AT(T,N0,N1,-1,x0,x1)];
      if (init != init_golden) {
        std::cout << "DIFF: " << init << "\t\t" << init_golden << "\t\t" << init_golden-init << std::endl;
      } else std::cout << "SAME: " << init << std::endl;
    }
  }

  jacobi2d_golden(A_golden, T, N0, N1);
  jacobi2d_tiled(A, T, N0, N1);

  int ok=1;

#ifndef SW
  // for (int t=0; t<T; t++) {
    for (int x0=0; x0<N0; x0++) {
      for (int x1=0; x1<N1; x1++) {
        data_t ret        = A[AT(T, N0, N1, T-1, x0, x1)];
        data_t ret_golden = A_golden[AT(T, N0, N1, T-1, x0, x1)];
        if (ret != ret_golden) {
          ok=0;
          cout << "Difference at: " << (T-1) << ", " << x0 << ", " << x1 << endl;
          cout << ret << "\t\t" << ret_golden << "\t\t" << (ret_golden-ret) << endl;
        }
        // assert(ret == ret_golden);
      }
    }
  // }
#else
    int t = ST-1;
    for (int x0=0; x0<N0; x0++) {
      for (int x1=0; x1<N1; x1++) {
        data_t ret        = A[AT(T, N0, N1, t, x0, x1)];
        data_t ret_golden = A_golden[AT(T, N0, N1, t, x0, x1)];
        if (ret != ret_golden) {
          ok=0;
          cout << "Difference at: " << t << ", " << x0 << ", " << x1 << endl;
          cout << ret << "\t\t" << ret_golden << endl;
        }
        // assert(ret == ret_golden);
      }
    }
#endif

  FREE(A);
  FREE(A_golden);

  return (1-ok);
}
