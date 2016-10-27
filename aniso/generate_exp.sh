#!/usr/bin/bash -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

S0=$1
S1=$2
S2=$3
UNROLL=$4

DIRNAME="${S0}x${S1}x${S2}_${UNROLL}"

SRCS="               \
    inner_tile.h     \
    hls_types.h      \
    perf_counter.h   \
    perf_counter.cpp \
    jacobi2d.h       \
    main.cpp
    "



mkdir -p $DIRNAME/src

cp $DIR/Makefile $DIRNAME/

for f in $SRCS; do
    cp $DIR/src/$f $DIRNAME/src/$f
done

$DIR/generate_inner_tile.py $UNROLL > $DIRNAME/src/inner_tile.cpp

$DIR/codegen.py -d 2 -S "${S0},${S1},${S2}" \
                --dep "-1,0,0" \
                --dep "-1,-1,0" \
                --dep "-1,1,0" \
                --dep "-1,0,-1" \
                --dep "-1,0,1" \
                --dep "-1,-1,-1" \
                --dep "-1,1,-1" \
                --dep "-1,-1,-1" \
                --dep "-1,1,1" \
                --skew "1,0,0;1,1,0;1,0,1" \
                --tile-skew "1,1,1;0,1,0;0,0,1" > $DIRNAME/src/scan.h

cat >$DIRNAME/src/params.h <<EOF
#ifndef PARAMS_H
#define PARAMS_H

#define S0 ${S0}
#define S1 ${S1}
#define S2 ${S2}
#define UNROLL ${UNROLL}

#endif
EOF

