/* Domain:  [T, N] -> { D[i0, i1, i2] : 0 < i0 < T and 0 < i1 <= -2 + N and 0 < i2 <= -2 + N } */
#undef max
#undef min
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#include <vector>
#define INNER_TILE(t, x0, x1) \
  {\
    coords.push_back(t); \
    coords.push_back(x0); \
    coords.push_back(x1); \
    ntiles++; \
  }
template<> void scan_hw_tiles<32,32,32>(int T,int N, data_array *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt) {
if (T >= 2 && N >= 3)
  for (int c0 = 0; c0 <= (3 * T + 2 * N - 7) / 32; c0 += 1)
    if ((T + N >= 32 * c0 + 3 && N + 28 >= 32 * c0) || (16 * c0 >= T + 15 && 32 * c0 + 2 >= T + N && N + 28 >= 32 * c0) || (3 * T + N >= 32 * c0 + 65 && 32 * c0 >= N + 29 && 3 * ((N + 64 * c0 + 28) / 96) >= 2 * c0) || (32 * c0 >= N + 29 && 32 * c0 + 64 >= 3 * T + N && 3 * ((T + N + 32 * c0 + 61) / 64) >= 2 * c0 + 3 && T + 32 * ((T + N + 32 * c0 + 61) / 64) >= 32 * c0 + 33) || (32 * c0 + 4 >= 3 * T + N && 3 * T + 2 * N >= ((2 * T + N + 28) % 32) + 32 * c0 + 39 && (T - 1) % 32 >= ((2 * T + N + 28) % 32) + 1) || (32 * c0 >= N + 29 && 3 * T + N >= 32 * c0 + 5 && 32 * c0 + 64 >= 3 * T + N && T + N + 32 * ((T + N + 32 * c0 + 61) / 64) >= 32 * c0 + 67 && 32 * c0 + 32 >= T + 32 * ((T + N + 32 * c0 + 61) / 64)) || (32 * c0 >= N + 29 && 32 * c0 + 4 >= 3 * T + N && 3 * T + 2 * N >= ((2 * T + N + 28) % 32) + 32 * c0 + 7 && (2 * T + N + 28) % 32 >= (T - 1) % 32) || (N <= 35 && 3 * T + N >= 32 * c0 + 65 && (c0 - 2) % 3 == 0) || (T >= 18 && N >= 36 && T + N <= 66 && c0 == 2) || (3 * T + 3 * N >= 32 * c0 + 41 && 32 * c0 >= N + 29 && 3 * T + N >= 32 * c0 + 5 && 32 * c0 + 64 >= 3 * T + N && 32 * c0 + 134 >= 3 * T + 3 * N && (c0 - 2) % 3 == 0))
      {
      	int d = c0;
      	sem_wait(mutex_cout);
      	cout << "HW blocked on barrier: " << (d&1) << endl;
sem_post(mutex_cout);
      	pthread_barrier_wait(&barriers[d & 1]);
      	sem_wait(mutex_cout);
      	cout << "HW unblocked." << endl;
      	sem_post(mutex_cout);
      	// Before
      	int ntiles = 0;
      	std::vector<int> coords;
      	
      	for (int c0 = max(max(3, d + floord(-2 * T - N + 1, 32) + 4), floord(32 * d - N - 31, 96) + 3); c0 < min(min(d - 3, floord(T + N - 2, 32) - 1), floord(16 * d + N - 18, 48)); c0 += 1)
      	  for (int c1 = max(max(d - 2 * c0 + 2, d - c0 + floord(-T - 1, 32) + 2), -c0 + (d + c0 + 1) / 2 + 1); c1 < min(min(d - c0, d - 2 * c0 + (N - 2) / 32), -c0 + (32 * d + N + 32 * c0 + 30) / 64); c1 += 1)
      	    INNER_TILE(d - c0 - c1, c0, c1);
      	
      	// After WF
      	if (ntiles > 0) {
      	  int *coords_ = (int *) ALLOC(3 * ntiles * sizeof(int));
      	  for (int i=0; i<3*ntiles; i++)
      	    coords_[i] = coords[i];
      	  inner_tiles_wrapper(arr, ntiles, coords_, mutex_cout, cnt);
      	  FREE(coords_);
      	}
      	
      }

}

#undef INNER_TILE
#define OUTER_TILE(t, x0, x1) outer_tile_wrapper(arr, (t), (x0), (x1), mutex_cout, cnt)
template<> void scan_sw_tiles<32,32,32>(int T,int N, data_array *arr, pthread_barrier_t barriers[2], sem_t *mutex_cout, perf_counter *cnt) {
if (T >= 2 && N >= 3)
  for (int c0 = 0; c0 <= (3 * T + 2 * N - 7) / 32; c0 += 1)
    if ((T + N >= 32 * c0 + 3 && N + 28 >= 32 * c0) || (16 * c0 >= T + 15 && 32 * c0 + 2 >= T + N && N + 28 >= 32 * c0) || (3 * T + N >= 32 * c0 + 65 && 32 * c0 >= N + 29 && 3 * ((N + 64 * c0 + 28) / 96) >= 2 * c0) || (32 * c0 >= N + 29 && 32 * c0 + 64 >= 3 * T + N && 3 * ((T + N + 32 * c0 + 61) / 64) >= 2 * c0 + 3 && T + 32 * ((T + N + 32 * c0 + 61) / 64) >= 32 * c0 + 33) || (32 * c0 + 4 >= 3 * T + N && 3 * T + 2 * N >= ((2 * T + N + 28) % 32) + 32 * c0 + 39 && (T - 1) % 32 >= ((2 * T + N + 28) % 32) + 1) || (32 * c0 >= N + 29 && 3 * T + N >= 32 * c0 + 5 && 32 * c0 + 64 >= 3 * T + N && T + N + 32 * ((T + N + 32 * c0 + 61) / 64) >= 32 * c0 + 67 && 32 * c0 + 32 >= T + 32 * ((T + N + 32 * c0 + 61) / 64)) || (32 * c0 >= N + 29 && 32 * c0 + 4 >= 3 * T + N && 3 * T + 2 * N >= ((2 * T + N + 28) % 32) + 32 * c0 + 7 && (2 * T + N + 28) % 32 >= (T - 1) % 32) || (N <= 35 && 3 * T + N >= 32 * c0 + 65 && (c0 - 2) % 3 == 0) || (T >= 18 && N >= 36 && T + N <= 66 && c0 == 2) || (3 * T + 3 * N >= 32 * c0 + 41 && 32 * c0 >= N + 29 && 3 * T + N >= 32 * c0 + 5 && 32 * c0 + 64 >= 3 * T + N && 32 * c0 + 134 >= 3 * T + 3 * N && (c0 - 2) % 3 == 0))
      {
      	int d = c0;
      	sem_wait(mutex_cout);
      	cout << "SW blocked on barrier: " << (d&1) << endl;
sem_post(mutex_cout);
      	pthread_barrier_wait(&barriers[d & 1]);
      	sem_wait(mutex_cout);
      	cout << "SW unblocked." << endl;
      	sem_post(mutex_cout);
      	// Before
      	
      	if (T >= 33) {
      	  if (d >= 2 && T >= 16 * d + 16 && N + 28 >= 32 * d && (T - 1) % 32 <= 30) {
      	    OUTER_TILE(0, 0, d);
      	  } else if (16 * d + 15 >= T && N + 28 >= 32 * d && (T - 1) % 32 <= 30)
      	    OUTER_TILE(0, 0, d);
      	  if (N >= 3 && (T - 1) % 32 <= 30) {
      	    for (int c0 = floord(32 * d - N - 29, 96) + 1; c0 <= min(min(min(floord(T - 1, 32) - 1, d - (N + 30) / 32 - 1), floord(32 * d + N - 36, 96) + 1), floord(32 * d - N + 2, 96) + 1); c0 += 1)
      	      for (int c1 = max(d - 2 * c0, -c0 + (d + c0 + 1) / 2); c1 <= -c0 + (32 * d + N + 32 * c0 + 29) / 64; c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	    if (((T - 1) % 32) + 32 * d + 36 >= 3 * T + N && 2 * T + N + 32 * floord(T + N + 29, 32) >= ((T - 1) % 32) + 32 * d + 36 && ((T - 1) % 32) + ((T + N + 29) % 32) <= 30) {
      	      int c0 = 32 * d + 36 >= 3 * T + N ? d - (2 * T + N + 28) / 32 + 1 : T - 31 * T / 32 - 1;
      	      if (32 * c0 >= T || 1) {
      	        int c1 = 32 * c0 + 1 >= T ? d - c0 - floord(T + 31, 32) + 1 : d - 2 * c0;
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      }
      	    }
      	  } else if (N >= 3 && T % 32 == 0)
      	    for (int c0 = max(d + floord(-2 * T - N + 2, 32) + 2, floord(32 * d - N - 29, 96) + 1); c0 <= min(min(min(T / 32, d - (N + 30) / 32 - 1), floord(32 * d + N - 36, 96) + 1), floord(32 * d - N + 2, 96) + 1); c0 += 1)
      	      for (int c1 = max(max(d - 2 * c0, (-T / 32) + d - c0 + 1), -c0 + (d + c0 + 1) / 2); c1 <= -c0 + (32 * d + N + 32 * c0 + 29) / 64; c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	  if ((T - 1) % 32 <= 30)
      	    for (int c0 = max(floord(T - 1, 32), d + floord(-2 * T - N + 2, 32) + 2); c0 <= min(min(min(d + floord(-T - 1, 16) + 1, floord(T + N - 3, 32)), d + floord(-N + 1, 32) - 1), floord(32 * d - N + 2, 96) + 1); c0 += 1) {
      	      OUTER_TILE((T + 32) / 32 - 1, c0, d - c0 - (T + 32) / 32 + 1);
      	      if (T >= 32 * c0 + 1 && 32 * c0 + 31 >= T && N + 96 * c0 >= 32 * d + 35)
      	        OUTER_TILE(c0 - 1, c0, d - 2 * c0 + 1);
      	    }
      	  if (T >= 64 && N >= 98)
      	    for (int c0 = max(d - (2 * T + N + 29) / 32 + 2, floord(32 * d - N + 2, 96) + 2); c0 <= min(min(min(d - (T + 16) / 16 + 1, (T + N - 2) / 32 - 2), d - (N + 30) / 32 - 1), floord(32 * d + N - 35, 96)); c0 += 1) {
      	      if (32 * c0 + 32 >= T && (T - 1) % 32 <= 30)
      	        OUTER_TILE((T + 32) / 32 - 1, c0, d - c0 - (T + 32) / 32 + 1);
      	      for (int c1 = max(d - 2 * c0, d - c0 - (T + 32) / 32 + 2); c1 <= d - 2 * c0 + 1; c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (96 * c0 + 65 >= 32 * d + N) {
      	        OUTER_TILE(d - 2 * c0 - 1, c0, c0 + 1);
      	      } else if (T + 32 * c0 + 32 * ((32 * d + N - 32 * c0 + 30) / 64) >= 32 * d + 32 && (32 * d - N - 32 * c0 - 31) % 64 <= 62)
      	        OUTER_TILE(-c0 + (32 * d - N + 32 * c0 - 31) / 64 + 1, c0, d - (32 * d - N + 32 * c0 - 31) / 64 - 1);
      	    }
      	  if (N >= 98)
      	    for (int c0 = max(d + floord(-T - 1, 16) + 2, floord(32 * d - N + 2, 96) + 2); c0 <= min(min(floord(d - 1, 3) + 1, d - (N + 30) / 32 - 1), floord(32 * d + N - 35, 96)); c0 += 1) {
      	      if (3 * c0 >= d && 16 * d + 31 >= T + 16 * c0 && (d - c0) % 2 == 0)
      	        OUTER_TILE((d - c0) / 2, c0, (d - c0) / 2);
      	      for (int c1 = max(max(d - 2 * c0, d - c0 - (T + 32) / 32 + 2), -c0 + (d + c0 + 1) / 2); c1 <= d - 2 * c0 + 1; c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N <= 129 && 3 * c0 == d + 2) {
      	        OUTER_TILE((d - 7) / 3, (d + 2) / 3, (d + 5) / 3);
      	      } else if (32 * d + N >= 64 * c0 + 32 * ((32 * d + N - 32 * c0 + 30) / 64) + 34 && (32 * d - N - 32 * c0 - 31) % 64 <= 62)
      	        OUTER_TILE(-c0 + (32 * d - N + 32 * c0 - 31) / 64 + 1, c0, d - (32 * d - N + 32 * c0 - 31) / 64 - 1);
      	    }
      	  if ((T - 1) % 32 <= 30) {
      	    for (int c0 = max(1, d + floord(-N + 1, 32)); c0 <= min(min(floord(d - 1, 3) + 1, floord(T - 1, 32) - 1), d - 2); c0 += 1) {
      	      for (int c1 = max(d - 2 * c0, -c0 + (d + c0 + 1) / 2); c1 <= min(d - c0 - 1, d - 2 * c0 + 1); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N + 32 * c0 + 29 >= 32 * d)
      	        OUTER_TILE(0, c0, d - c0);
      	    }
      	    if (T >= 64 && d >= 3 * floord(d, 3) + 1 && T >= 32 * floord(d, 3) + 33 && N + 32 * floord(d, 3) + 62 >= 32 * d && 32 * floord(d, 3) + 63 >= T)
      	      if (3 * T >= 32 * d + 67 || 1) {
      	        OUTER_TILE((d + 3) / 3 - 1, (d + 3) / 3, d - 2 * ((d + 3) / 3) + 1);
      	        if (N + 32 * (d / 3) + 61 >= 32 * d)
      	          OUTER_TILE(0, (T + 32) / 32 - 1, d - (T + 32) / 32 + 1);
      	      }
      	    if (T >= 64 && 32 * d + 29 >= 3 * T && N <= 97 && 3 * T + 3 * N >= 32 * d + 230 && (d - 2) % 3 == 0) {
      	      OUTER_TILE((d - 2) / 3, (d + 1) / 3, (d + 1) / 3);
      	      OUTER_TILE((d - 5) / 3, (d + 1) / 3, (d + 4) / 3);
      	    } else if (T <= 63)
      	      for (int c0 = max(d + floord(-2 * T - N + 2, 32) + 2, floord(32 * d - N + 2, 96) + 2); c0 < min(floord(T + N - 2, 32) - 1, d + floord(-N + 1, 32)); c0 += 1)
      	        OUTER_TILE(1, c0, d - c0 - 1);
      	    for (int c0 = floord(32 * d + N - 35, 96) + 1; c0 <= min(min(d + floord(-T - 1, 16) + 1, floord(T + N - 2, 32) - 2), d + floord(-N + 1, 32) - 1); c0 += 1) {
      	      OUTER_TILE((T + 32) / 32 - 1, c0, d - c0 - (T + 32) / 32 + 1);
      	      for (int c1 = d - 2 * c0 + (N - 2) / 32; c1 <= min(d - 2 * c0 + (N - 3) / 32 + 1, -c0 + (32 * d + N + 32 * c0 + 29) / 64); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	    }
      	    for (int c0 = max(max(floord(T + N - 2, 32) - 1, d + floord(-2 * T - N + 2, 32) + 2), floord(32 * d - N + 2, 96) + 2); c0 <= min(min(d + floord(-T - 1, 16) + 1, floord(T + N - 3, 32)), d + floord(-N + 1, 32) - 1); c0 += 1) {
      	      OUTER_TILE((T + 32) / 32 - 1, c0, d - c0 - (T + 32) / 32 + 1);
      	      for (int c1 = d - c0 - (T + 32) / 32 + 2; c1 <= min(d - 2 * c0 + (N - 3) / 32 + 1, -c0 + (32 * d + N + 32 * c0 + 29) / 64); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	    }
      	    for (int c0 = max((T - 1) / 32, d + floord(-N + 1, 32)); c0 <= min(d - (T + 16) / 16 + 1, floord(T + N - 3, 32)); c0 += 1) {
      	      OUTER_TILE((T + 32) / 32 - 1, c0, d - c0 - (T + 32) / 32 + 1);
      	      if (c0 >= 2 && T >= 32 * c0 + 1 && 32 * c0 + 31 >= T)
      	        OUTER_TILE(c0 - 1, c0, d - 2 * c0 + 1);
      	      for (int c1 = max(d - c0 + floord(-T - 1, 32) + 2, d - 2 * c0 + floord(N - 2, 32)); c1 <= min(d - c0 - 1, d - 2 * c0 + (N + 29) / 32); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N + 32 * c0 + 29 >= 32 * d && N + 29 >= 32 * c0)
      	        OUTER_TILE(0, c0, d - c0);
      	    }
      	    if (d - 3 * ((T - 16) / 32) <= 1 && (T - 16) % 32 <= 15 && 2 * ((T - 16) % 32) + N + 62 >= 2 * T && 2 * ((T - 16) % 32) + 32 * d + 64 >= 3 * T) {
      	      OUTER_TILE(d / 3, d - 2 * (d / 3), d / 3);
      	      if (d >= 6 && d % 3 == 0)
      	        OUTER_TILE((d / 3) - 1, d / 3, (d / 3) + 1);
      	      if (2 * ((T - 16) % 32) + N + 61 >= 2 * T)
      	        OUTER_TILE(0, d - 2 * ((T + 32) / 32) + 2, 2 * ((T + 32) / 32) - 2);
      	    } else if (d >= 9 && 3 * T >= 32 * d + 48 && 32 * d + 93 >= 3 * T && N == 98 && d % 3 == 0) {
      	      OUTER_TILE(d / 3, d / 3, d / 3);
      	      OUTER_TILE((d / 3) - 1, d / 3, (d / 3) + 1);
      	    }
      	  } else {
      	    if (16 * d >= T && N + 28 >= 32 * d) {
      	      OUTER_TILE(0, 0, d);
      	    } else if (d >= 2 && T >= 16 * d + 16 && N + 28 >= 32 * d)
      	      OUTER_TILE(0, 0, d);
      	    for (int c0 = max(1, d + floord(-N + 1, 32)); c0 <= min(min(floord(d - 1, 3) + 1, d - 2), T / 32); c0 += 1) {
      	      for (int c1 = max(max(d - 2 * c0, (-T / 32) + d - c0 + 1), -c0 + (d + c0 + 1) / 2); c1 <= min(d - c0 - 1, d - 2 * c0 + 1); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N + 32 * c0 + 29 >= 32 * d)
      	        OUTER_TILE(0, c0, d - c0);
      	    }
      	    for (int c0 = max((T / 32) + 1, d + floord(-N + 2, 32)); c0 <= min((-T / 16) + d, floord(N - 2, 32)); c0 += 1)
      	      OUTER_TILE(0, c0, d - c0);
      	    for (int c0 = floord(32 * d + N - 35, 96) + 1; c0 <= min(min((-T / 16) + d, floord(T + N - 2, 32) - 2), d + floord(-N + 1, 32) - 1); c0 += 1)
      	      for (int c1 = d - 2 * c0 + (N - 2) / 32; c1 <= min(d - 2 * c0 + (N - 3) / 32 + 1, -c0 + (32 * d + N + 32 * c0 + 29) / 64); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	    for (int c0 = max(floord(N - 2, 32) + 1, d + floord(-N + 1, 32)); c0 <= min((-T / 16) + d, floord(T + N - 2, 32) - 2); c0 += 1) {
      	      for (int c1 = d - 2 * c0 + (N - 2) / 32; c1 <= min(d - c0 - 1, d - 2 * c0 + (N - 3) / 32 + 1); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N + 29 >= 32 * c0)
      	        OUTER_TILE(0, c0, d - c0);
      	    }
      	  }
      	  if (T >= 64 && N >= 98)
      	    for (int c0 = max(floord(d - 1, 3) + 2, d - (T + 16) / 16 + 2); c0 < min(min(d - 1, (T + N - 2) / 32 - 1), floord(16 * d + N + 14, 48)); c0 += 1) {
      	      if (16 * d + 31 >= T + 16 * c0 && (d - c0) % 2 == 0) {
      	        OUTER_TILE((d - c0) / 2, c0, (d - c0) / 2);
      	      } else if (16 * ((d - c0 - 1) % 2) + 16 * d + N >= 48 * c0 + 50 && T + 16 * c0 >= 16 * ((d - c0 - 1) % 2) + 16 * d + 16)
      	        OUTER_TILE(-c0 + (d + c0) / 2, c0, d - (d + c0) / 2);
      	      for (int c1 = d - 2 * c0 + (N - 2) / 32; c1 <= min(min(d - c0 - 1, d - 2 * c0 + (N - 3) / 32 + 1), -c0 + (32 * d + N + 32 * c0 + 29) / 64); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N + 32 * c0 + 29 >= 32 * d && N + 29 >= 32 * c0) {
      	        OUTER_TILE(0, c0, d - c0);
      	      } else if (32 * d >= N + 32 * c0 + 31 && 32 * d + N >= 64 * c0 + 32 * ((32 * d + N - 32 * c0 + 30) / 64) + 34 && (32 * d - N - 32 * c0 - 31) % 64 <= 62)
      	        OUTER_TILE(-c0 + (32 * d - N + 32 * c0 - 31) / 64 + 1, c0, d - (32 * d - N + 32 * c0 - 31) / 64 - 1);
      	    }
      	  if (N >= 3 && (T - 1) % 32 <= 30) {
      	    for (int c0 = max(max(0, d - 1), floord(32 * d - N - 29, 64) + 1); c0 <= min(min(min(floord(T - 1, 32) - 1, d), (N - 3) / 32 + 1), floord(32 * d + N + 28, 64)); c0 += 1)
      	      OUTER_TILE(0, c0, d - c0);
      	    if (d - 3 * ((T - 16) / 32) >= 2 && (T - 16) % 32 <= 15 && 3 * T + N >= 2 * ((T - 16) % 32) + 32 * d + 35 && 2 * ((T - 16) % 32) + 32 * d + 97 >= 3 * T + N) {
      	      OUTER_TILE((T + 32) / 32 - 1, d - 2 * ((T + 32) / 32) + 2, (T + 32) / 32 - 1);
      	      for (int c1 = -floord(-T - 1, 32); c1 < min(-2 * floord(-T - 1, 32) - 2, -d - 4 * floord(-T - 1, 32) + floord(N - 3, 32) - 2); c1 += 1)
      	        OUTER_TILE(-c1 + 2 * ((T + 32) / 32) - 2, d - 2 * ((T + 32) / 32) + 2, c1);
      	      if (2 * T + N >= 2 * ((T - 16) % 32) + 32 * d + 3)
      	        OUTER_TILE(0, d - 2 * ((T + 32) / 32) + 2, 2 * ((T + 32) / 32) - 2);
      	    } else if (d >= 5 && T >= 48 && T <= 63 && T + N >= 32 * d + 2)
      	      for (int c1 = 1; c1 <= 2; c1 += 1)
      	        OUTER_TILE(-c1 + 2, d - 2, c1);
      	  }
      	  if ((T - 1) % 32 <= 30) {
      	    for (int c0 = max((T - 1) / 32, d - 1); c0 <= min(min(d, floord(N - 3, 32) + 1), floord(32 * d + N + 28, 64)); c0 += 1)
      	      OUTER_TILE(0, c0, d - c0);
      	    if ((d <= 5 && T <= 95 && 32 * d >= N + 95 && N >= 19 && T + N + 29 >= 32 * d) || (3 * T + N >= 32 * d + 102 && 32 * d >= T + N + 30 && N >= 19 && d - 3 * ((T + 31) / 32) <= -4 && 3 * ((T + 31) % 32) + 32 * d + 101 >= 3 * T + N) || (32 * d >= T + N + 30 && N >= 3 && 32 * d + N + 63 >= 3 * T && 32 * d + 101 >= 3 * T + N && d - 3 * ((T + 31) / 32) <= -4) || (3 * T >= 32 * d + N + 64 && 32 * d >= T + N + 30 && N >= 3 && N <= 18 && 3 * T >= 3 * ((T - 1) % 32) + 32 * d + N + 16 && 3 * ((T + 31) % 32) + 32 * d + N + 63 >= 3 * T) || (d == 4 && T <= 95 && T + N >= 99 && N <= 18)) {
      	      OUTER_TILE((T + 32) / 32 - 2, (T + 32) / 32 - 1, d - 2 * ((T + 32) / 32) + 3);
      	    } else if (N >= 3 && N <= 97 && d + 3 * floord(-d + 2, 3) <= 1 && T + 32 * floord(-d + 2, 3) >= 16 && T + 32 * floord(-d + 2, 3) <= 31 && N - 64 * (d / 3) <= -31 && T + N >= 32 * d - 64 * (d / 3) + 3) {
      	      OUTER_TILE(d / 3, d - 2 * (d / 3), d / 3);
      	      if (N >= 35)
      	        OUTER_TILE(d / 3 - 1, d - 2 * (d / 3), d / 3 + 1);
      	    }
      	  }
      	  if (d >= 6 && N >= 82 && N <= 97 && 3 * T + 3 * N >= 32 * d + 294 && d % 3 == 0) {
      	    for (int c1 = d / 3; c1 <= min((2 * d / 3) - 2, (d / 3) + 1); c1 += 1)
      	      OUTER_TILE((2 * d / 3) - c1 - 1, (d / 3) + 1, c1);
      	    if (d == 6)
      	      OUTER_TILE(0, 3, 3);
      	  } else if ((d >= 11 && 3 * T == 32 * d + 32 && N >= 67 && N <= 97) || (d >= 7 && T + 32 >= 16 * d && N <= 97 && N + 189 >= 32 * d && T % 32 == 0) || (16 * d >= T + 48 && 6 * T + 2 >= 64 * d + N && N <= 97 && 3 * floord(64 * d + N + 93, 96) >= 2 * d + 4 && T % 32 == 0) || (T + 32 >= 16 * d && 32 * d >= N + 190 && N <= 97 && 3 * floord(64 * d + N + 93, 96) >= 2 * d + 4 && T % 32 == 0)) {
      	    int c0 = 32 * d >= N + 189 && 6 * T + 3 >= 64 * d + N && 3 * ((64 * d + N + 93) / 96) >= 2 * d + 4 ? d - (64 * d + N + 93) / 96 + 2 : 6 * T + 3 >= 64 * d + N && N + 188 >= 32 * d ? 3 : (d + 1) / 3;
      	    OUTER_TILE(c0 - 1, c0, d - 2 * c0 + 1);
      	    OUTER_TILE(c0 - 2, c0, d - 2 * c0 + 2);
      	  } else if (d >= 10 && N == 98 && d - 3 * ((T + 31) / 32) <= -5 && (T + 31) % 32 <= 30 && (d - 1) % 3 == 0) {
      	    OUTER_TILE((d - 1) / 3, (d + 2) / 3, (d - 1) / 3);
      	  } else if (d >= 10 && 3 * T >= 32 * d + 64 && N == 98 && T % 32 == 0 && (d - 1) % 3 == 0) {
      	    OUTER_TILE((d - 1) / 3, (d + 2) / 3, (d - 1) / 3);
      	  } else if ((d >= 7 && N <= 97 && N + 189 >= 32 * d && d - 2 * ((T + 31) / 32) <= -1 && (T + 31) % 32 <= 30) || (32 * d >= N + 190 && N <= 97 && d - 2 * ((T + 31) / 32) <= -1 && (T - 1) % 32 <= 30 && 3 * floord(64 * d + N + 93, 96) >= 2 * d + 4) || (N <= 97 && d - 2 * ((T + 31) / 32) >= 0 && (T - 1) % 32 <= 30 && 6 * T >= 6 * ((T - 1) % 32) + 64 * d + N + 4 && 3 * floord(64 * d + N + 93, 96) >= 2 * d + 4) || (d >= 8 && 3 * T >= 32 * d + 35 && 32 * d + 125 >= 3 * T && N >= 67 && N <= 97 && (d - 2) % 3 == 0)) {
      	    int c0 = N + 189 >= 32 * d && 2 * ((T + 31) / 32) >= d && (T - 1) % 32 <= 30 ? 3 : 32 * d >= N + 190 && 6 * T >= 64 * d + N + 3 && 3 * ((64 * d + N + 93) / 96) >= 2 * d + 4 && (T + 31) % 32 <= 30 && 6 * T + 92 >= 6 * ((T + 31) % 32) + 64 * d + N ? d - (64 * d + N + 93) / 96 + 2 : (d + 1) / 3;
      	    OUTER_TILE(c0 - 1, c0, d - 2 * c0 + 1);
      	    OUTER_TILE(c0 - 2, c0, d - 2 * c0 + 2);
      	  }
      	  if (N >= 35) {
      	    for (int c0 = max(floord(d - 1, 3) + 2, floord(16 * d + N + 14, 48)); c0 <= min(min(d - 2, floord(T + N - 2, 32) - 2), floord(16 * d + N - 20, 48) + 1); c0 += 1) {
      	      for (int c1 = -c0 + (d + c0 + 1) / 2; c1 <= min(d - c0 - 1, d - 2 * c0 + (N - 3) / 32 + 1); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (N + 29 >= 32 * c0)
      	        OUTER_TILE(0, c0, d - c0);
      	    }
      	    if (d >= 10 && 3 * T >= 32 * d + 16 && 32 * d + 61 >= 3 * T && N == 98 && (d - 1) % 3 == 0)
      	      OUTER_TILE((d - 1) / 3, (d + 2) / 3, (d - 1) / 3);
      	    if ((T - 1) % 32 <= 30)
      	      for (int c0 = max(max(floord(d - 1, 3) + 2, d + 2 * floord(-T, 32) + 3), floord(T + N - 2, 32) - 1); c0 <= min(min(d - 2, floord(T + N - 3, 32)), floord(16 * d + N - 20, 48) + 1); c0 += 1) {
      	        for (int c1 = -c0 + (d + c0 + 1) / 2; c1 <= min(d - c0 - 1, d - 2 * c0 + (N - 3) / 32 + 1); c1 += 1)
      	          OUTER_TILE(d - c0 - c1, c0, c1);
      	        if (N + 29 >= 32 * c0)
      	          OUTER_TILE(0, c0, d - c0);
      	      }
      	  }
      	  if (N >= 3 && T % 32 == 0)
      	    for (int c0 = max(max(0, d - 1), floord(32 * d - N - 29, 64) + 1); c0 <= min(min(d, (N - 3) / 32 + 1), floord(32 * d + N + 28, 64)); c0 += 1)
      	      OUTER_TILE(0, c0, d - c0);
      	  if (T % 32 == 0) {
      	    for (int c0 = max(max(floord(T + N - 2, 32) - 1, d + floord(-2 * T - N + 2, 32) + 2), floord(32 * d - N + 2, 96) + 2); c0 <= min(min((-T / 16) + d, floord(T + N - 3, 32)), d + floord(-N + 1, 32) - 1); c0 += 1)
      	      for (int c1 = (-T / 32) + d - c0 + 1; c1 <= min(d - 2 * c0 + (N - 3) / 32 + 1, -c0 + (32 * d + N + 32 * c0 + 29) / 64); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	    for (int c0 = max(floord(T + N - 2, 32) - 1, d + floord(-N + 1, 32)); c0 <= min((-T / 16) + d, floord(T + N - 3, 32)); c0 += 1) {
      	      for (int c1 = (-T / 32) + d - c0 + 1; c1 <= min(d - c0 - 1, d - 2 * c0 + floord(N - 3, 32) + 1); c1 += 1)
      	        OUTER_TILE(d - c0 - c1, c0, c1);
      	      if (T == 64 && N + 29 >= 32 * c0)
      	        OUTER_TILE(0, c0, d - c0);
      	    }
      	    if (N >= 35)
      	      for (int c0 = max(max((-T / 16) + d + 1, floord(d - 1, 3) + 2), floord(T + N - 2, 32) - 1); c0 <= min(min(d - 2, floord(T + N - 3, 32)), floord(16 * d + N - 20, 48) + 1); c0 += 1) {
      	        for (int c1 = -c0 + (d + c0 + 1) / 2; c1 <= min(d - c0 - 1, d - 2 * c0 + (N - 3) / 32 + 1); c1 += 1)
      	          OUTER_TILE(d - c0 - c1, c0, c1);
      	        if (T == 64 && N + 29 >= 32 * c0)
      	          OUTER_TILE(0, c0, d - c0);
      	      }
      	  }
      	} else if (T >= 2 && N >= 3) {
      	  for (int c0 = max(max(0, d + floord(-T - N + 2, 32) + 1), floord(32 * d - N - 29, 64) + 1); c0 <= min(d - 2, floord(T + N - 3, 32)); c0 += 1)
      	    OUTER_TILE(0, c0, d - c0);
      	  for (int c0 = max(max(max(0, d - 1), d - (T + N + 29) / 32 + 1), floord(32 * d - N - 29, 64) + 1); c0 <= min(min(d, (T + N - 3) / 32), floord(32 * d + N + 28, 64)); c0 += 1)
      	    OUTER_TILE(0, c0, d - c0);
      	}
      	
      	// After WF
      	
      }

}

#undef OUTER_TILE
template<> void scan_wavefronts<32,32,32>(int T,int N, pthread_barrier_t barriers[2], sem_t *mutex_cout) {
if (T >= 2 && N >= 3)
  for (int c0 = 0; c0 <= (3 * T + 2 * N - 7) / 32; c0 += 1)
    if ((T + N >= 32 * c0 + 3 && N + 28 >= 32 * c0) || (16 * c0 >= T + 15 && 32 * c0 + 2 >= T + N && N + 28 >= 32 * c0) || (3 * T + N >= 32 * c0 + 65 && 32 * c0 >= N + 29 && 3 * ((N + 64 * c0 + 28) / 96) >= 2 * c0) || (32 * c0 >= N + 29 && 32 * c0 + 64 >= 3 * T + N && 3 * ((T + N + 32 * c0 + 61) / 64) >= 2 * c0 + 3 && T + 32 * ((T + N + 32 * c0 + 61) / 64) >= 32 * c0 + 33) || (32 * c0 + 4 >= 3 * T + N && 3 * T + 2 * N >= ((2 * T + N + 28) % 32) + 32 * c0 + 39 && (T - 1) % 32 >= ((2 * T + N + 28) % 32) + 1) || (32 * c0 >= N + 29 && 3 * T + N >= 32 * c0 + 5 && 32 * c0 + 64 >= 3 * T + N && T + N + 32 * ((T + N + 32 * c0 + 61) / 64) >= 32 * c0 + 67 && 32 * c0 + 32 >= T + 32 * ((T + N + 32 * c0 + 61) / 64)) || (32 * c0 >= N + 29 && 32 * c0 + 4 >= 3 * T + N && 3 * T + 2 * N >= ((2 * T + N + 28) % 32) + 32 * c0 + 7 && (2 * T + N + 28) % 32 >= (T - 1) % 32) || (N <= 35 && 3 * T + N >= 32 * c0 + 65 && (c0 - 2) % 3 == 0) || (T >= 18 && N >= 36 && T + N <= 66 && c0 == 2) || (3 * T + 3 * N >= 32 * c0 + 41 && 32 * c0 >= N + 29 && 3 * T + N >= 32 * c0 + 5 && 32 * c0 + 64 >= 3 * T + N && 32 * c0 + 134 >= 3 * T + 3 * N && (c0 - 2) % 3 == 0))
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

