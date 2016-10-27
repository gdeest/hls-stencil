#!/usr/bin/env python

from sys import argv
from math import ceil
from operator import mul
from functools import reduce
from jinja2 import Environment,PackageLoader

def dependence_depths(ndims, deps):
    return {dep: max(map(lambda v: -v[dep], deps)) for dep in range(ndims)}

def window_size(ndims, deps, unroll_factors):
    depths = dependence_depths(ndims, deps)
    return tuple(ceil((depths[i] + unroll_factors[i-1])/ unroll_factors[i-1]) for i in range(1,ndims))

def lexscan_aux(i, ndims, unroll_factors):
    if i == ndims-1:
        return list(map(lambda x: [x], range(unroll_factors[i-1])))

    inner_dims = lexscan_aux(i+1, ndims, unroll_factors)
    vecs = []

    for j in range(unroll_factors[i-1]):
        vecs = vecs + list(map(lambda v: [j]+v, inner_dims))

    return vecs

def lexscan(ndims, unroll_factors):
    return list(map(tuple, lexscan_aux(1, ndims, unroll_factors)))


def compute_updates(ndims, deps, unroll_factors):
    # Compute window size based on unroll factors and dependencies
    sizes = window_size(ndims, deps, unroll_factors)

    # Compute function mapping point in DP-tile to its lexicographic index
    offsets = [reduce(mul, unroll_factors[(i+1):(ndims-1)], 1) for i in range(ndims-1)]
    compute_idx = lambda v: sum(map(lambda x: x[0]*x[1], zip(v,offsets)))

    # For each point in the top-right element of the window,
    # determine its accesses in the window and generate the update arguments.
    updates = []
    for idx,v in enumerate(lexscan(ndims, unroll_factors)):
        # Add offset from the bottom-left of the window
        v2 = tuple(map(lambda x: x[1] + (sizes[x[0]]-1)*unroll_factors[x[0]], enumerate(v)))

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

# Parameters to code generation
deps = [
    (-1,0,-1),
    (-1,-1,0),
    (-1,-1,-1),
    (-1,-1,-2),
    (-1,-2,-1)
]

unroll = (1,int(argv[1]))


ndims = 3
winsize = window_size(ndims, deps, unroll)
updates = compute_updates(ndims, deps, unroll)
depths = dependence_depths(ndims, deps)
uhalo = [ceil(depths[i] / unroll[i-1]) for i in range(1,ndims)]

params = dict(winsize=winsize,
              updates=updates,
              unroll=unroll,
              depths=depths,
              uhalo=uhalo)

# Create template environment
env = Environment(loader=PackageLoader(__name__, 'templates'))
env.filters['ceil'] = ceil

# Load templates
inner_tile = env.get_template("inner_tile.cpp")
hls_types = env.get_template("hls_types.h")


print(inner_tile.render(**params))
#print(hls_types.render(**params))

