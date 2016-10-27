#include <stdio.h>
#include "jacobi2d.h"
#include "perf_counter.h"
#include "inner_tile.h"

#include <math.h>

template<int _S0, int _S1, int _S2> void scan_hw_tiles
(int T, int N, data_array *da, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt);

template<int _S0, int _S1, int _S2> void scan_sw_tiles
(int T, int N, data_array *da, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt);

template<int _S0, int _S1, int _S2> void scan_wavefronts
(int T, int N, pthread_barrier_t barriers[2], sem_t *mutex_cout);

#include "scan.h"

using namespace std;

void inner_tiles_wrapper(data_array *arr, int ntiles, int *coords, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "--- INNER TILES ---" << endl;
  for (int i=0; i<ntiles; i++) {
    cout << coords[3*i    ] << ", ";
    cout << coords[3*i + 1] << ", ";
    cout << coords[3*i + 2] << endl;
  }
  cout << "--- /INNER TILES ---" << endl;
  sem_post(mutex_cout);

// #ifdef SW
  // Just call inner_tile_wrapper on each tile and let it handle the rest.

  inner_tiles(arr->B0,
              arr->B1,
              arr->B2,
              arr->Baux,
              arr->B0,
              arr->B1,
              arr->B2,
              arr->Baux,
              arr->T,
              arr->N,
              ntiles,
              coords);

  // for (int i=0; i<ntiles; i++) {
  //   inner_tile_wrapper(arr, coords[3*i], coords[3*i+1], coords[3*i+2], mutex_cout, cnt);
  // }
// #else
//   // Actually call inner_tiles with multiple tile coordinates.

//   // TODO: Performance counter does not really make sense when calling a varying number of tiles at once.
//   cnt->start();
// #ifdef __SDSCC__
//   uint64_t before = sds_clock_counter();
// #endif
//   inner_tiles(arr, arr, T, N0, N1, ntiles, coords);
// #ifdef __SDSCC__
//   uint64_t after = sds_clock_counter();
//   sem_wait(mutex_cout);
//   cout << "MULTIPLE: " << ntiles << ", " << (after-before) << endl;
//   sem_post(mutex_cout);
// #endif
//   cnt->stop();
// #endif
}

void inner_tile_wrapper(data_array *arr, int t, int x0, int x1, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "INNER TILE: " << t << ", " << x0 << ", " << x1 << endl;
  sem_post(mutex_cout);

  compute_tile(arr, t, x0, x1);

  // int coords[3];
  // coords[0] = t;
  // coords[1] = x0;
  // coords[2] = x1;

  // cnt->start();
  // inner_tiles(arr, arr, T, N0, N1, 1, coords);
  // // inner_tile(arr, arr, T, N0, t, x0);
  // cnt->stop();
}

void outer_tile_wrapper(data_array *arr, int t, int x0, int x1, sem_t *mutex_cout, perf_counter *cnt) {
  sem_wait(mutex_cout);
  cout << "OUTER TILE: " << t << ", " << x0 << ", " << x1 << endl;
  sem_post(mutex_cout);

  cnt->start();

  compute_tile(arr, t, x0, x1);
  // outer_tile(arr, T, N0, N1, t, x0, x1);

  cnt->stop();
}

data_t* alloc_array(int size) {
  cout << "Allocating array of size: " <<  size*sizeof(data_t) << endl;
  data_t *arr = (data_t*) ALLOC(size*sizeof(data_t));

  if (!arr) {
    cout << "Error allocating array.\n" << endl;
    exit(-1);
  }

  return arr;
}

void allocate_data_array(data_array *da, int T, int N) {
  da->T = T;
  da->N = N;
  da->B0 =  alloc_array(B_P0_SIZE(T,N));
  da->B1 =  alloc_array(B_P1_SIZE(T,N));
  da->B2 =  alloc_array(B_P2_SIZE(T,N));
  da->Baux =  alloc_array(B_AUX_SIZE(T,N));
}

void dump_input(data_array *da) {
	int t, i, j;
	t = 0;

	for (i=0; i<da->N; i++) {
		for (j=0; j<da->N; j++) {
			printf ("%.02f ", B_P0_accFull(da->T, da->N, da->B0, t, i, j));
		}
		printf("\n");
	}
}

void verify_output(data_array *da) {
	int t, i, j;
	t = da->T;

	printf("out(cmp)\n");

  int count=1;
	for (i=da->T; i<da->T+da->N; i++) {
		for (j=da->T; j<da->T+da->N; j++) {
      data_t v = B_P0_accFull(da->T, da->N, da->B0, t, i, j);
			// printf ("%.02f ", v);
      assert(v == count++);
		}
		// printf("\n");
	}
}

#define B_(t,i,j) B[(t)*(N*N)+(i)*N+(j)]
void golden(int T, int N, data_t *B, unsigned int seed) {
	srand(seed);
	int count = 1;

  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      data_t v;
      if (1<=i && i<N-1 && 1<=j && j<N-1)
        v = count;
      else
        v = count;
      B_(0,i,j) = v;//(1.0*rand()/RAND_MAX)*300.0;
      B_(1,i,j) = v;//(1.0*rand()/RAND_MAX)*300.0;
      count++;
    }
  }


	for (int t=0; t<T; t++) {
		for (int i=1; i<N-1; i++) {
			for (int j=1; j<N-1; j++) {
				B_((t+1)%2,i,j) =
          (B_(t%2,i,j) + B_(t%2,i+1,j) + B_(t%2,i-1,j) + B_(t%2,i,j+1) + B_(t%2,i,j-1))*0.2;
			}
		}
	}

	
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			fprintf(stdout, "%.2lf ", B_((T%2),i,j));
			//if (j%20==0) fprintf(stderr, "\n");
		}
		fprintf(stdout, "\n");
	}
}
#undef B_

data_t compute_point(data_t a, data_t b, data_t c, data_t d, data_t e) {
  // cout << "Compute: ";
  // cout << a << ", " << b << ", " << c << ", " << d << ", " << e << endl;
	data_t res =  (a + b + c + d + e)*0.2;
	//printf("%.02f = %.02f %.02f %.02f %.02f %.02f\n", res, a, b, c, d, e);
	return res;
}

void dump_reg(data_t reg[S0+HALO_0][S1+HALO_1][S2+HALO_2]) {

	for (int i0 = 0; i0<S0+HALO_0; i0++)  {
		for (int i1 = 0; i1<S1+HALO_1; i1++) {
			for (int i2 = 0; i2<S2+HALO_2; i2++) {
				printf("%.02f ", reg[i0][i1][i2]);
			}
			printf("\n");
		}
		printf("\n");
	}

}
void load_face(data_array *da, data_t reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
  int T = da->T, N = da->N;
	int p0, p1, p2;
	int x;
	data_t *offset;

  int st = t0*S0, si = t1*S1, sj = t2*S2;

	//load
	//P0
	offset = &B_P0_acc(T, N, da->B0, t0, t1, t2, 0, -HALO_1, 0);
	x = 0;
	for (p1=0; p1<S1+HALO_1; p1++)
		for (p2=HALO_2; p2<S2+HALO_2; p2++)
			for (p0=0; p0<HALO_0; p0++) {
        reg[p0][p1][p2] = offset[x];
				//printf("regP0[%d][%d][%d] = %.02f\n", p0, p1, p2, offset[x]);
				x++;
			}
	//P1
	offset = &B_P1_acc(T, N, da->B1, t0, t1, t2, 0, 0, -HALO_2);
//printf("offsetP1 = %d\n", B_P1_offset(t0, t1, t2, 0, 0, -HALO_2));
	x = 0;
	for (p2=0; p2<S2+HALO_2; p2++)
		for (p0=HALO_0; p0<S0+HALO_0; p0++)
			for (p1=0; p1<HALO_1; p1++) {
        reg[p0][p1][p2] = offset[x];
				//printf("regP1[%d][%d][%d] = %.02f\n", p0, p1, p2, offset[x]);
				x++;
			}

	//P2
	offset = &B_P2_acc(T, N, da->B2, t0, t1, t2, -HALO_0, 0, 0);
	x = 0;
	for (p0=0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2++) {
        reg[p0][p1][p2] = offset[x];
				//printf("regP2[%d][%d][%d] = %.02f\n", p0, p1, p2, offset[x]);
				x++;
			}

	//aux
	offset = &B_AUX_acc(T, N, da->Baux, t0, t1, t2, 0);
	x = 0;
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2++) {
        reg[p0][p1][p2] = offset[x];
				//printf("regAUX[%d][%d][%d] = %.02f @ %d\n", p0, p1, p2, offset[x], B_AUX_offset(t0, t1, t2, x));
				x++;
			}
}

void store_face(data_array *da, data_t reg[S0+HALO_0][S1+HALO_1][S2+HALO_2], int t0, int t1, int t2) {
  int T = da->T;
  int N = da->N;
	int p0, p1, p2, x;
	data_t *offset;

	//store
	//P0
	offset = &B_P0_acc(T, N, da->B0, t0+1, t1, t2, 0, 0, 0);
	x = 0;
	for (p1=HALO_1; p1<S1+HALO_1; p1++)
		for (p2=HALO_2; p2<S2+HALO_2; p2++)
			for (p0=S0; p0<S0+HALO_0; p0++) {
				offset[x] = reg[p0][p1][p2];
//printf("offsetP0[%d] = %.02f @ %d\n", x, offset[x], B_P0_offset(t0+1,t1,t2,0,0,0)+ x);
				x++;
			}
	//P1
	offset = &B_P1_acc(T, N, da->B1, t0, t1+1, t2, 0, 0, -HALO_2);
//printf("offsetP1 = %d\n", B_P1_offset(t0, t1+1, t2, 0, 0, -HALO_2));
	x = 0;
	for (p2=0; p2<S2+HALO_2; p2++)
		for (p0=HALO_0; p0<S0+HALO_0; p0++)
			for (p1=S1; p1<S1+HALO_1; p1++) {
				offset[x] = reg[p0][p1][p2];
//printf("offsetP1[%d] = %.02f @ %d\n", x, offset[x], B_P1_offset(t0, t1+1, t2, 0, 0, -HALO_2)+x);
				x++;
			}

	//P2
	offset = &B_P2_acc(T, N, da->B2, t0, t1, t2+1, -HALO_0, 0, 0);
	x = 0;
	for (p0=0; p0<S0+HALO_0; p0++)
		for (p1=HALO_1; p1<S1+HALO_1; p1++)
			for (p2=S2; p2<S2+HALO_2; p2++) {
				offset[x] = reg[p0][p1][p2];
//printf("offsetP2[%d] = %.02f\n", x, offset[x]);
				x++;
			}

	//aux
	offset = &B_AUX_acc(T, N, da->Baux, t0, t1, t2, 0);
	x = 0;
	for (p0=0; p0<HALO_0; p0++)
		for (p1=0; p1<HALO_1; p1++)
			for (p2=0; p2<HALO_2; p2++) {
				offset[x] = reg[p0+S0][p1+S1][p2+S2];
//printf("offsetAUX[%d] = %.02f @ %d\n", x, offset[x], B_AUX_offset(t0, t1, t2, x));
				x++;
			}


}

void compute_tile(data_array *da, int t0, int t1, int t2) {
  int T = da->T, N = da-> N;
	int p0, p1, p2;

	data_t reg[S0+HALO_0][S1+HALO_1][S2+HALO_2];

	//CMP for intra-tile; useful only for handling partial tiles
	//F out0[HALO_0][S1+HALO_1][S2];
	//F out1[S0][HALO_1][S2+HALO_2];
	//F out2[S0+HALO_0][S1][HALO_2];

	load_face(da, reg, t0, t1, t2);

#ifdef DUMP_REG_IN
	printf("tile (%d %d %d) input \n", t0, t1, t2);
	dump_reg(reg);
#endif

	//for (p0=HALO_0; p0<S0+HALO_0; p0++)
	//	for (p1=HALO_1; p1<S1+HALO_1; p1++)
	//		for (p2=HALO_2; p2<S2+HALO_2; p2++) {
	for (p0=1; p0<S0+HALO_0; p0++)
		for (p1=1; p1<S1+HALO_1; p1++)
			for (p2=1; p2<S2+HALO_2; p2++) {
				//guard for partial tiles
				int st = t0*S0+p0;
				int si = t1*S1+p1-HALO_1;
				int sj = t2*S2+p2-HALO_2;
				int t = st;
				int i = si - st;
				int j = sj - st;


				if (1<=t && t<T && 1<=i && i<N-1 && 1<=j && j<N-1 && HALO_0<p0 && HALO_1<p1 && HALO_2<p2) {
          // cout << t << "," << i << "," << j << endl;
//printf("%d %d %d -> %d %d %d -> %d %d %d\n", p0, p1, p2, st,si,sj,t,i,j);
					reg[p0][p1][p2] = compute_point(reg[p0-1][p1-1][p2], reg[p0-1][p1][p2-1], reg[p0-1][p1-1][p2-1], reg[p0-1][p1-1][p2-2], reg[p0-1][p1-2][p2-1]);
				} else {
          reg[p0][p1][p2] = reg[p0-1][p1-1][p2-1];
				}
			}

#ifdef DUMP_REG_OUT
	printf("tile (%d %d %d) output \n", t0, t1, t2);
	dump_reg(reg);
#endif

	store_face(da, reg, t0, t1, t2);
}

void init_data_array(data_array *da, unsigned int seed) {
  int T = da->T, N = da->N;
	srand(seed);
	data_t count = 1;

	//P0 
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			B_P0_accFull(T, N, da->B0, 0, i, j) = count;//(1.0*rand()/RAND_MAX)*300.0;
			count++;
		}
	}

	//P2 : might have some issue due to % of -1
	/*for (int i=0; i<N; i++)
		for (int j=0; j<N; j++) {
			B_P2_accFull(da->B2, -1, i, j+S2) = B_P0_accFull(da->B0, 0, i, j+S2-HALO_2);
		}

	for (int ti=0; ti<N/S1+1; ti++)
		for (int tj=1; tj<N/S2+1; tj++) {
			int x =0;
			for (int i=ti*S1; i<(ti+1)*S1; i++)
				for (int j=tj*S2+S2-HALO_2; j<(tj+1)*S2; j++) {
	B_P2_accFull(da->B2, 0, i, j) = B_P0_accFull(da->B0, 0, i, j);
				x++;
			}
		}
*/
	//aux
	for (int ti=1; ti<N/S1+1; ti++)
		for (int tj=1; tj<N/S2+1; tj++) {
			int x =0;
			for (int i=(ti-1)*S1+S1-HALO_1; i<(ti)*S1; i++)
				for (int j=(tj-1)*S2+S2-HALO_2; j<(tj)*S2; j++) {
          B_AUX_acc(T, N, da->Baux, 0, ti, tj, x) = B_P0_accFull(T, N, da->B0, 0, i, j);
					x++;
				}
	}

/*
	for (int i=0; i < B_P2_SIZE; i++) {
		printf("%.02f ", da->B2[i]);
	}
	printf("\n");
	for (int i=0; i < B_P0_SIZE; i++) {
		printf("%.02f ", da->B0[i]);
	}
	printf("\n");

	for (int i=0; i < B_P1_SIZE; i++) {
		printf("%.02f ", da->B1[i]);
	}
	printf("\n");

	for (int i=0; i < B_P2_SIZE; i++) {
		printf("%.02f ", da->B2[i]);
	}
	printf("\n");
*/

#ifdef DUMP_INPUT
	printf("input (cmp)");
	dump_input(da);
#endif
}

void *sw_thread(void *_arg) {
  thread_arg *arg = (thread_arg *) _arg;
  scan_sw_tiles<S0,S1,S2>(arg->da->T, arg->da->N, arg->da, arg->barriers, arg->mutex_cout, arg->cnt);
  return 0;
}

void *hw_thread(void *_arg) {
  thread_arg *arg = (thread_arg *) _arg;
  scan_hw_tiles<S0,S1,S2>(arg->da->T, arg->da->N, arg->da, arg->barriers, arg->mutex_cout, arg->cnt);
  return 0;
}

void jacobi2d_tiled(data_array *da) {
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

  arg.da = da;
  arg.barriers = barriers;
  arg.mutex_cout = &mutex_cout;

  arg_hw = arg;
  arg_hw.cnt = hw_counter;

  arg_sw = arg;
  arg_sw.cnt = sw_counter;

  pthread_t sw, hw;

  pthread_create(&sw, 0, sw_thread, &arg_sw);
  pthread_create(&hw, 0, hw_thread, &arg_hw);

  scan_wavefronts<S0,S1,S2>(da->T, da->N, arg.barriers, &mutex_cout);


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
  int T, N;

  data_array da;

  if (argc < 3) {
    cout << "Usage: ./jacobi1d <T:uint> <N:uint>" << endl;
    exit(-1);
  } else {
    T = atoi(argv[1]);
    N = atoi(argv[2]);
    cout << "Problem size: " << T << "x" << N << endl;
  }

  allocate_data_array(&da, T, N);
  init_data_array(&da, SEED);
  jacobi2d_tiled(&da);

  verify_output(&da);
  cout << "blihh" << endl;

  // data_t *arr = (data_t*)ALLOC(2*T*N*N*sizeof(data_t));

  // int ok=1;
  // golden(T, N, arr, SEED);

  // FREE(da.B0);
  // FREE(da.B1);
  // FREE(da.B2);
  // FREE(da.Baux);
  // FREE(A_golden);

  return 0;
}
