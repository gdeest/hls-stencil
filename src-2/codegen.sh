#!/bin/sh -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

$DIR/codegen.py -d 2 -S "$1,$2,$3" \
                --dep "-1,0,0" \
                --dep "-1,-1,0" \
                --dep "-1,1,0" \
                --dep "-1,0,-1" \
                --dep "-1,0,1" \
                --skew "1,0,0;1,1,0;1,0,1" \
                --tile-skew "1,1,1;0,1,0;0,0,1"
