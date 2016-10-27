// Local Variables:
// mode: c++
// End:
#ifndef PERF_COUNTER_H
#define PERF_COUNTER_H

#include <stdint.h>

class perf_counter
{
private:
  uint64_t tot, cnt, calls;

public:
  perf_counter();
  void start();
  void stop();
  uint64_t avg_cpu_cycles();
};

#endif
