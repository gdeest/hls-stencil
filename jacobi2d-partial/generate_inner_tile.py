#!/usr/bin/env python

from sys import argv
from math import ceil
from operator import mul
from functools import reduce
from jinja2 import Environment,PackageLoader

def dependence_depths(ndims, deps):
    """Computes positive and negative dependence depths for each dimension.

    Example:
    >>> dependence_depths(2, [(-1,-2), (-1,1)])
    [(1, 0), (2, 1)]

    """

    dep_pos = lambda dep: max(map(lambda v: max(0,  v[dep]), deps))
    dep_neg = lambda dep: max(map(lambda v: max(0, -v[dep]), deps))

    return [(dep_neg(dep), dep_pos(dep)) for dep in range(ndims)]


def window_size(ndims, deps, unroll_factors):
    """Assuming Jacobi-style dependencies, computes the (n-1)-dimensional window
    required to compute a set of points (pack), given unrolling factors for
    each spatial dimension.

    Returns a tuple containing:
    - The size of the window (in terms of unrolled blocks).
    - The "computation offset", from the last block in the window.

    This function does NOT assume that the tiling conditions are valid for each
    dimension, but the dimensions where this condition is invalid should not be tiled.

    Examples:
    >>> window_size(3, [(-1,0,-1), (-1,-1,0), (-1,-1,-1), (-1,-1,-2), (-1,-2,-1)], (1,1))
    ((3, 3), (0, 0))

    >>> window_size(3, [(-1,0,-1), (-1,-1,0), (-1,-1,-1), (-1,-1,-2), (-1,-2,-1)], (1,2))
    ((3, 2), (0, 0))

    >>> window_size(3, [(-1,0,0), (-1,-1,1), (-1,-1,0), (-1,-1,-1), (-1,-2,0)], (1,1))
    ((3, 3), (0, 1))

    >>> window_size(3, [(-1,0,0), (-1,-1,1), (-1,-1,0), (-1,-1,-1), (-1,-2,0)], (1,2))
    ((3, 3), (0, 1))
    """

    depths = dependence_depths(ndims, deps)
    size_neg = tuple(ceil(depths[i][0]/unroll_factors[i-1]) for i in range(1,ndims))
    size_pos = tuple(ceil(depths[i][1]/unroll_factors[i-1]) for i in range(1,ndims))
    size = tuple(1 + size_neg[i] + size_pos[i] for i in range(ndims-1))

    return (size, size_pos)

def lexscan_aux(i, ndims, unroll_factors):
    if i == ndims-1:
        return list(map(lambda x: (x,), range(unroll_factors[i-1])))

    inner_dims = lexscan_aux(i+1, ndims, unroll_factors)
    vecs = []

    for j in range(unroll_factors[i-1]):
        vecs = vecs + list(map(lambda v: (j,)+v, inner_dims))

    return vecs

def lexscan(ndims, unroll_factors):
    """Generates indices in a window in lexicographic order.

    Example:
    >>> lexscan(3, [2,2])
    [(0, 0), (0, 1), (1, 0), (1, 1)]
    """
    return list(map(tuple, lexscan_aux(1, ndims, unroll_factors)))


def compute_updates(ndims, deps, unroll_factors):
    """Compute updates in the pack.
    >>> deps = [(-1,0,-1), (-1,-1,0), (-1,-1,-1), (-1,-1,-2), (-1,-2,-1)]

    >>> compute_updates(3, deps, [1, 1])
    [(0, [([2, 1], 0), ([1, 2], 0), ([1, 1], 0), ([1, 0], 0), ([0, 1], 0)])]
    """
    # Compute window size based on unroll factors and dependencies
    size, compute_offset = window_size(ndims, deps, unroll_factors)

    # Compute function mapping point in DP-tile to its lexicographic index
    offsets = [reduce(mul, unroll_factors[(i+1):(ndims-1)], 1) for i in range(ndims-1)]
    compute_idx = lambda v: sum(map(lambda x: x[0]*x[1], zip(v,offsets)))

    # For each point in the pack, determine its accesses in the window and
    # generate the update arguments.
    updates = []
    for idx,v in enumerate(lexscan(ndims, unroll_factors)):
        # Add offset from the bottom-left of the window
        v2 = tuple(map(lambda x: x[1] + (size[x[0]]-compute_offset[x[0]]-1)*unroll_factors[x[0]], enumerate(v)))

        # Add v to each dependence
        accesses = map(lambda dep: tuple([sum(x) for x in zip(v2,dep[1:ndims])]), deps)

        args = []
        for access in accesses:
            divs = list(map(lambda x: x[1] // unroll_factors[x[0]], enumerate(access)))
            mods = list(map(lambda x: x[1]  % unroll_factors[x[0]], enumerate(access)))
            pos = compute_idx(mods)
            args = args + [(divs, pos)]

        updates = updates + [(idx, args)]

    return updates

def jacobi2d_partial_tiling(deps, innermost_unroll):
    for d in deps:
        assert(len(d) == 3)
        assert(d[0] == -1 and d[1] <=0)

    unroll= [1, innermost_unroll]

    winsize = window_size(3, deps, unroll)
    updates = compute_updates(3, deps, unroll)
    depths  = dependence_depths(3, deps)

    uhalo = [(ceil(depths[i][0] / unroll[i-1]), ceil(depths[i][1] / unroll[i-1])) for i in range(1,3)]

    params = dict(winsize=winsize,
                  updates=updates,
                  unroll=unroll,
                  depths=depths,
                  uhalo=uhalo)

    # Create template environment
    env = Environment(loader=PackageLoader(__name__, 'templates'))
    env.filters['ceil'] = ceil

    inner_tile = env.get_template("inner_tile.cpp")

    print(inner_tile.render(**params))

# Parameters to code generation
deps = [
    (-1,0,0),
    (-1,-1,1),
    (-1,-1,0),
    (-1,-1,-1),
    (-1,-2,0)
]

jacobi2d_partial_tiling(deps, int(argv[1]))

# ndims = 3
# winsize = window_size(ndims, deps, unroll)
# updates = compute_updates(ndims, deps, unroll)
# depths = dependence_depths(ndims, deps)
# uhalo = [ceil(depths[i] / unroll[i-1]) for i in range(1,ndims)]

# params = dict(winsize=winsize,
#               updates=updates,
#               unroll=unroll,
#               depths=depths,
#               uhalo=uhalo)

# # Create template environment
# env = Environment(loader=PackageLoader(__name__, 'templates'))
# env.filters['ceil'] = ceil

# # Load templates
# inner_tile = env.get_template("inner_tile.cpp")
# hls_types = env.get_template("hls_types.h")


# print(inner_tile.render(**params))
# #print(hls_types.render(**params))

import doctest
doctest.testmod()
