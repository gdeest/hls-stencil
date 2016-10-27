#ifndef ACCESS_MACROS_H
#define ACCESS_MACROS_H

#include "common_partial_skewing.h"

#define AT_SIZE_UNSKEWED(arr,T,N,t,x0,x1) arr[(t+1)*((N)*(S2))+(x0)*(S2)+(x1)]
#define AT_SIZE_SKEWED(arr,T,N,t,x0,x1)   AT_SIZE_UNSKEWED(arr,T,N,t,(x0)-(t),(x1))
#define AT_TILED(arr,T,N,t,x0,tt,xx0,xx1) AT_SIZE_SKEWED(arr,T,N,(t)*S0+(tt), (x0)*S1+(xx0), (xx1))

#endif
