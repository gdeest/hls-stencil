#!/usr/bin/env python

from update_board import *
from os import path
from time import sleep

import pickle

# BURST_LENGTHS = [1, 2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 128, 256]
# # NBURSTS = [1, 100]
# NBURSTS = [100]
# #STRIDES = [0, 512, 1024, 2048]
# STRIDES = [2048]

# STs = [32, 48, 64, 80, 96]
# S0s = [32, 48, 64, 80, 96]
# STs = [32, 64, 128, 256]
# S0s = [32, 64, 128, 256]
UFs = [2]
STs = [16, 32, 64]
S0s = [16, 32, 64]
S1s = [16, 32, 64]

# HOSTNAME = "192.168.0.42"
HOSTNAME = "10.0.0.184"

def save_obj(obj, name):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

class DirectoryNotFound(Exception):
    pass

def measure_performance(st, s0, s1, uf, nexp=10):
    dirname = path.dirname(path.realpath(__file__))
    dirname = path.join(dirname, "%dx%dx%d_%d" % (st, s0, s1, uf))
    print(dirname)

    if path.isdir(dirname):
        copy_files(HOSTNAME, dirname)
        reboot(HOSTNAME)

        sleep(10)
        ssh = wait_for_reboot_and_ssh(HOSTNAME)

        ssh.exec_command("rm -f /mnt/timings.csv")
        ssh.exec_command("echo 'ntiles, cycles' > /mnt/timings.csv")

        timings = []
        for i in range(nexp):
            stdin, stdout, stderr = ssh.exec_command("/mnt/jacobi2d.elf 200 200 200 | grep MULTIPLE | sed 's/MULTIPLE: //'")
            output = stdout.readlines()
            ret = map(lambda line: line.split(','), output)
            ret = map(lambda lst: map(str.strip, lst), ret)
            ret = map(lambda lst: map(int, lst), ret)
            ret = map(tuple, ret)
            ret = list(ret)
            # print("Timings: %s" % ret)
            timings.append(ret)

        ssh.close()
    else:
       raise DirectoryNotFound()

    return(timings)

try:
    results = load_obj("results")
    print("Successfully loaded previous results from file.")
except FileNotFoundError:
    print("No previous results found.")
    results = {}

# for st in STs:
#     for s0 in S0s:
#         for s1 in S1s:
#             for uf in UFs:
with open("sizes.txt", 'r') as handle:
    for line in handle:
        st, s0, s1, uf = map(int, line.split(' '))
        key = (st, s0, s1, uf)

        if key in results:
            print("-- Skipping %s (already done) --" % str(key))
        else:
            try:
                print("----- EXPERIMENT %s -----" % str(key))
                timings = measure_performance(st, s0, s1, uf)
                print(timings)
                results[key] = timings
            except DirectoryNotFound:
                print("Directory not found - key: %s" % str(key))

                save_obj(results, "results")
