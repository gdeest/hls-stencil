#ifndef JACOBI1D_H
#define JACOBI1D_H

#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <pthread.h>
#include <semaphore.h>
#include <math.h>

#include "perf_counter.h"

#include "params.h"

#if not(defined(ST) && defined(S0) && defined(S1))
#define ST 16
#define S0 16
#define S1 16
#endif

#ifndef UNROLLFAC
#define UNROLLFAC 2
#endif

#define AT(T,N0,N1,t,i,j) (((t)+1)*(N0*N1) + (i)*N1 + (j))
#define ARRAY_SIZE(T,N0,N1) (((T)+1)*(N0)*(N1)*(sizeof(data_t)))

#if  defined(__SDSCC__) || defined(__SDSVHLS__)
#define COMPUTE(d0,d1,d2,d3,d4) (((d0)+(d1)+(d2)+(d3)+(d4))/5.02)
#else
#define COMPUTE(d0,d1,d2,d3,d4) (cos((d0)-(d1)+(d2)-(d3)+(d4)))
#endif
/* #define COMPUTE(d0,d1,d2,d3,d4) d4 */

using namespace std;

// typedef unsigned long long data_t;
typedef float data_t;

#define SEED 4784

// Select alloc/free functions based on whether we're compiling with SDSoC or not.
#if defined(__SDSCC__) || defined(__SDSVHLS__)
#include "sds_lib.h"
#define ALLOC sds_alloc_non_cacheable
#define FREE sds_free
#else
#define ALLOC malloc
#define FREE  free
#endif


struct thread_arg {
  perf_counter *cnt;
  pthread_barrier_t *barriers;
  sem_t *mutex_cout;
  data_t *arr;
  int T, N0, N1;
};

void inner_tiles_wrapper(data_t *arr, int T, int N0, int N1, int ntiles, int *coords, sem_t *mutex_cout, perf_counter *cnt);
void inner_tile_wrapper(data_t *arr, int T, int N0, int N1, int t, int x0, int x1, sem_t *mutex_cout, perf_counter *cnt);
void outer_tile_wrapper(data_t *arr, int T, int N0, int N1, int t, int x0, int x1, sem_t *mutex_cout, perf_counter *cnt);

/* data_t f(data_t, data_t, data_t, data_t, data_t); */

/* #include "transforms.h" */

#endif
