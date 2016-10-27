#!/bin/bash -e

source "/run/media/huginn/b85bf520-705c-4507-b6d8-262e11ce3f5b/soft/SDSoC/2016.2/settings64.sh"

# STs="32 64 128 256 512"
# S0s="32 64 128 256 512"
# STs="16 32 64"
# S0s="16 32 64"
# S1s="16 32 64"
# STs="128"
# S0s="128"
# S1s="128"
#STs="16 32 64"
#S0s="16 32 64"
#S1s="16 32 64"
# STs="32"
# S0s="32"
# S1s="128"
# STs="16 32 64"
# S0s="16 32 64"
# S1s="16 32 64"
# UNROLLFACs="2"


FIFO_TILE_NUM=1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR

INCDIR=/run/media/huginn/b85bf520-705c-4507-b6d8-262e11ce3f5b/soft/SDSoC/2016.2/Vivado_HLS/2016.2/include

cd src
CHECKSUM=$(cat read_coords.cpp read_inputs.cpp write_outputs.cpp compute.cpp | md5sum)
cat >checksum.h <<EOF
#ifndef CHECKSUM_H
#define CHECKSUM_H
#define CHECKSUM $CHECKSUM
#endif
EOF
cd ..

# for ST in $STs; do
#     for S0 in $S0s; do
#         for S1 in $S1s; do
#             for UNROLLFAC in $UNROLLFACs; do
regexp="([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)"
cat sizes.txt | while read line; do
    if [[ $line =~ $regexp ]]
    then
        ST="${BASH_REMATCH[1]}"
        S0="${BASH_REMATCH[2]}"
        S1="${BASH_REMATCH[3]}"
        UNROLLFAC="${BASH_REMATCH[4]}"
        echo $ST
        echo $S0
        echo $S1
        echo $UNROLLFAC

        NAME="${ST}x${S0}x${S1}_${UNROLLFAC}"
        echo "Generating $NAME..."
        if [ ! -d $NAME ]; then
            echo "Directory does not exist. Creating initial directory structure..."
            mkdir $NAME
            ln -s ../Makefile $NAME/Makefile
            cp -as $(pwd)/src $NAME/
            rm -f $NAME/src/scan.h
            rm -f $NAME/src/params.h
            rm -f $NAME/src/fifo_pragmas.h

            pushd $NAME
            echo "Code generation..."
            # cp ../scan.h src/scan.h
            cat >src/scan.h <<EOF
#include <math.h>
#define max(a,b) (((a)>=(b)) ? (a) : (b))
#define min(a,b) (((a)<=(b)) ? (a) : (b))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
EOF

            ../codegen.sh $ST $S0 $S1 >> src/scan.h

            cat >src/fifo_pragmas.h <<EOF
#pragma HLS STREAM variable=inputs depth=$(bc <<< "($FIFO_TILE_NUM)*(($S1+2+$S0)*$ST + $S0*$S1/2)")
#pragma HLS STREAM variable=outputs depth=$(bc <<< "($FIFO_TILE_NUM)*(($S0-2+$S1)*($ST-1) + $S0*$S1/2)")
EOF

            echo "#define ST $ST" >> src/params.h
            echo "#define S0 $S0" >> src/params.h
            echo "#define S1 $S1" >> src/params.h
            echo "#define UNROLLFAC $UNROLLFAC" >> src/params.h
        else
            pushd $NAME
        fi

        echo "Generating HW files"
        make
        popd
    else
        echo "ERREUR: $line"
    fi
done
#         done
#     done
# done
# done
