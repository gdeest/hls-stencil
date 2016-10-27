#ifndef COMPUTE_H
#define COMPUTE_H

#include "hls_types.h"

void compute(int ntiles,
             finput_t &inputs0,
             finput_t &inputs1,
             finput_t &inputs2,
             finput_t &inputs3,
             finput_t &inputs4,
             finput_t &inputs5,
             finput_t &inputs6,
             finput_t &inputs7,
             foutput_t &output);

#endif
