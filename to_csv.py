#!/usr/bin/env python

import array
# from rpy2.robjects import IntVector, Formula
# from rpy2.robjects.packages import importr
# base = importr('base')
# stats = importr('stats')

import pickle
import csv

import numpy as np

def save_obj(obj, name):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

results = load_obj('results')

def linreg(pairs):
    ntiles = []
    cycles = []
    for ntile, cycle in pairs:
        ntiles.append(ntile)
        cycles.append(cycle)
    ntiles = np.array(ntiles)
    cycles = np.array(cycles)

    m, b = np.polyfit(ntiles, cycles, 1)
    return (b, m)

def compute_average(pairs):
    avgs = []
    for ntiles, cycles in pairs:
        avgs.append(cycles/ntiles)
    return np.mean(avgs)

with open('results.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    # writer.writerow(['ST', 'S0', 'UF', 'Offset', 'Per Tile'])
    writer.writerow(['ST', 'S0', 'UF', 'Per Tile'])

    for key, value in results.items():
        st, s0, s1, uf = key
        name = "%dx%dx%s_%d" % (st, s0, s1, uf)

        # m, b = linreg(sum(value, []))
        avg = compute_average(sum(value, []))

        writer.writerow([st, s0, s1, uf, avg])
