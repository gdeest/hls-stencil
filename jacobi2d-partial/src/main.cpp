#include <assert.h>
#include <stdio.h>
#include "jacobi2d.h"
#include "perf_counter.h"

#include "inner_tile.h"
#include "access_macros_partial_skewing.h"

#include <math.h>
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

template<int _S0, int _S1> void scan_hw_tiles
(int T, int N, data_t *da, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt);

template<int _S0, int _S1> void scan_sw_tiles
(int T, int N, data_t *da, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt);

template<int _S0, int _S1> void scan_wavefronts
(int T, int N, pthread_barrier_t barriers[2], sem_t *mutex_cout);

#include "scan.h"

using namespace std;

data_t safe_get(data_t *arr, int T, int N, int tt, int xx0, int xx1) {
  if (tt>=-1 && tt <T && xx0 >= 0 && xx0 < N && xx1 >=0 && xx1 < S2) {
    return AT_SIZE_UNSKEWED(arr, T, N, tt, xx0, xx1);
  } else {
    return 0;
  }
}

data_t safe_get_golden(int N, data_t *arr, int t, int xx0, int xx1) {
  if (xx0 >= 0 && xx0 < N && xx1 >= 0 && xx1<S2) {
    return arr[((t+1)%2)*N*S2 + xx0*S2 + xx1];
  } else {
    return 0;
  }
}

void inner_tiles_wrapper(data_t *arr, int T, int N, int ntiles, int *coords, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  for (int i=0; i<ntiles; i++) {
    cout << "INNER: ";
    cout << coords[2*i    ] << ", ";
    cout << coords[2*i + 1] << endl;
  }
  sem_post(mutex_cout);

#ifdef __SDSCC__
  uint64_t before = sds_clock_counter();
#endif

  inner_tiles(arr, arr, T, N, ntiles, coords);

#ifdef __SDSCC__
  uint64_t after = sds_clock_counter();
  sem_wait(mutex_cout);
  cout << "MULTIPLE: " << ntiles << ", " << (after-before) << endl;
  sem_post(mutex_cout);
#endif
}

void inner_tile_wrapper(data_t *arr, int T, int N, int t, int x0, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "INNER TILE: " << t << ", " << x0 << endl;
  sem_post(mutex_cout);

  compute_tile(arr, T, N, t, x0);

  // int coords[3];
  // coords[0] = t;
  // coords[1] = x0;
  // coords[2] = x1;

  // cnt->start();
  // inner_tiles(arr, arr, T, N0, N1, 1, coords);
  // // inner_tile(arr, arr, T, N0, t, x0);
  // cnt->stop();
}

void outer_tile_wrapper(data_t *arr, int T, int N, int t, int x0, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "OUTER: " << t << ", " << x0 << endl;
  sem_post(mutex_cout);

  cnt->start();

  compute_tile(arr, T, N, t, x0);
  // outer_tile(arr, T, N0, N1, t, x0, x1);

  cnt->stop();
}

data_t* alloc_array(int size) {
  // cout << "Allocating array of size: " <<  size*sizeof(data_t) << endl;
  data_t *arr = (data_t*) ALLOC(size*sizeof(data_t));

  if (!arr) {
    cout << "Error allocating array.\n" << endl;
    exit(-1);
  }

  return arr;
}


void compute_tile(data_t *arr, int T, int N, int t, int x0) {
  for (int tt=0; tt<S0; tt++) {
    for (int  xx0=0; xx0<S1; xx0++) {
      int t_ = t*S0 + tt, x0_ = x0*S1 + xx0 - t_;

      if (t_ >= 0 && t_ < T && x0_ >= 0 && x0_ < N) {
        for (int xx1=0; xx1<S2; xx1++) {
          // std::cout << t_ << ", " << x0_ << ", " << xx1 << std::endl;
          AT_SIZE_UNSKEWED(arr, T, N, t_, x0_, xx1) = 
            (safe_get(arr, T, N, t_-1, x0_+1, xx1) +
             safe_get(arr, T, N, t_-1, x0_, xx1+1) +
             safe_get(arr, T, N, t_-1, x0_, xx1) +
             safe_get(arr, T, N, t_-1, x0_, xx1-1) +
             safe_get(arr, T, N, t_-1, x0_-1, xx1)) * 0.4;
        }
      }
    }
  }
}

void *sw_thread(void *_arg) {
  thread_arg *arg = (thread_arg *) _arg;
  scan_sw_tiles<S0,S1>(arg->T, arg->N, arg->arr, arg->barriers, arg->mutex_cout, arg->cnt);
  return 0;
}

void *hw_thread(void *_arg) {
  thread_arg *arg = (thread_arg *) _arg;
  scan_hw_tiles<S0,S1>(arg->T, arg->N, arg->arr, arg->barriers, arg->mutex_cout, arg->cnt);
  return 0;
}

void jacobi2d_tiled(int T, int N, data_t *arr) {
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
  arg.N = N;
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

  scan_wavefronts<S0,S1>(T, N, arg.barriers, &mutex_cout);


  sem_wait(&mutex_cout);
  // cout << "Waiting for threads to finish..." << endl;
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
  int T, N;


  if (argc < 3) {
    cout << "Usage: ./jacobi2d <T:uint> <N:uint>" << endl;
    exit(-1);
  } else {
    T = atoi(argv[1]);
    N = atoi(argv[2]);
    cout << "Problem size: " << T << "x" << N << endl;
  }

  data_t *arr = alloc_array((T+1)*N*S2);
  data_t *arr2 = new data_t[2*N*S2];

  for (int i=0; i<T; i++) {
    for (int j=0; j<T; j++) {
      data_t v = (1.0*rand()/RAND_MAX)*300.0;
      AT_SIZE_UNSKEWED(arr, T, N, -1, i, j) = v;
      arr2[i*S2 + j] = v;
    }
  }

  for (int t=0; t<T; t++) {
    for (int x0=0; x0<N; x0++) {
      for (int x1=0; x1<S2; x1++) {
        arr2[((t+1)%2)*N*S2 + x0*S2 + x1] =
          (safe_get_golden(N, arr2, t-1, x0+1, x1) +
           safe_get_golden(N, arr2, t-1, x0, x1+1) +
           safe_get_golden(N, arr2, t-1, x0, x1) +
           safe_get_golden(N, arr2, t-1, x0, x1-1) +
           safe_get_golden(N, arr2, t-1, x0-1, x1)) * 0.4;
      }
    }
  }

  jacobi2d_tiled(T, N, arr);

  for (int x0=0; x0<N; x0++) {
    for (int x1=0; x1<S2; x1++) {
      data_t v1 = safe_get(arr, T, N, T-1, x0, x1);
      data_t v2 = safe_get_golden(N, arr2, T-1, x0, x1);
      assert(v1 == v2);
    }
  }

  int ok = 1;

  return (1-ok);
}
