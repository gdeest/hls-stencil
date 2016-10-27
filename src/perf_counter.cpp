#include <iostream>
#include "perf_counter.h"

using namespace std;

#if defined(__SDSCC__) || defined(__SDSVHLS__)
#include "sds_lib.h"
#endif

perf_counter::perf_counter() {
  this->tot = 0;
  this->calls = 0;
  this->cnt = 0;
}

void perf_counter::start() {
  #ifdef __SDSCC__
  cnt = sds_clock_counter();
  #endif
  calls++;
}

void perf_counter::stop() {
#ifdef __SDSCC__
  tot += (sds_clock_counter() - cnt);
#endif 
}

uint64_t perf_counter::avg_cpu_cycles() {
  return ((tot+(calls>>1)) / calls);
}
