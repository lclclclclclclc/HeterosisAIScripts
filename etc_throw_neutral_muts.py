#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:51:27 2020

@author: egibson
"""

import pyslim
import tskit
import msprime

import numpy as np
import matplotlib.pyplot as plt

#%% read in tree seq from SLiM sim

DIR_tree = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/output/trees/"
tree_fn = "sgcb_m0_sexA"
# tree_fn = "ts_recording_test_full"

tree_file = DIR_tree + tree_fn + ".trees"
neg_ts = pyslim.load(tree_file)
# loaded with pyslim, so a SLiMTreeSeq object.
# This means I should be able to ts.recapitate via pyslim , which theoretically
# takes a recombination map (which msprime didn't recognize)

#%% messing around with original ts

print(f"There are {neg_ts.num_individuals} individuals in the SLiM tree.")
founders=neg_ts.first_generation_individuals()
initial_pops = [neg_ts.individual(f).population for f in founders]
print(f"Population(s) represented in the founders: {np.unique(initial_pops)}")
#numpy array of how long ago each individual was born
# SLiM automatically remembers the individuals that comprise the first generation of any new subpopulation created with addSubpop(), for easy recapitation and other analysis.
indiv_times = neg_ts.individual_times
plt.figure()
plt.hist(indiv_times)
recorded_times = np.unique(indiv_times)
print(f"When indivs. were added to the ts, in time ago: {recorded_times}")

root_num_spectrum = [t.num_roots for t in neg_ts.trees()]
plt.figure()
plt.hist(root_num_spectrum)
max_roots = max(root_num_spectrum)
print(f"Before recapitation, the max number of roots is {max_roots}")

#%% Recapitate

#TODO: sort out mutation rate variation / ratios
mut_rate = 10e-10

#TODO: sort out recombination rate variation / map
# https://pyslim.readthedocs.io/en/latest/tutorial.html#recapitation-with-a-nonuniform-recombination-map
# https://msprime.readthedocs.io/en/latest/api.html?highlight=map#recombination-map-limitations

# Let's use the sequence info doc with only recRate lines.  Feed to msprime.readHapMap
# https://msprime.readthedocs.io/en/stable/api.html?highlight=read_hapmap#msprime.RecombinationMap.read_hapmap

# Actually, let's not because it has to be subtlely changed.
# seq_info_parser.py is a protofunction returning a map.  Just need to integrate it.
L = int(neg_ts.sequence_length)
r1 = 10e-15
r2 = 10e-11
recomb_map = msprime.RecombinationMap([0, 100, L], [r1, r2, 0])
# I think I need to be doing msprime.simulation(from_ts), not msp.mutate...
# https://msprime.readthedocs.io/en/latest/api.html?highlight=recombination#msprime.simulate
# BUT  I now think I should be using the pyslim ts.recapitate
# https://pyslim.readthedocs.io/en/latest/python_api.html#pyslim.SlimTreeSequence.recapitate
# TODO: sort out reading in number of individuals.. necessary? why is default Ne=1?
n=10000
n_scale=10
n_p1 = n/n_scale
recap_ts = neg_ts.recapitate(recombination_map=recomb_map, Ne=1000)

# TODO: decide about simplifying here?  Essentially sampling here.
# I think this depends on how statistics are calculated

# throw mutations and output
ts = pyslim.SlimTreeSequence(msprime.mutate(recap_ts, rate=mut_rate, keep=True))
ts.dump(DIR_tree + tree_fn + "_with_neutrals" + ".trees")

# TODO: what format for existing code?  ms I think

#%% Mess around with full recap'd, mutated ts

recap_root_num_spectrum = [t.num_roots for t in ts.trees()]
recap_root = max(recap_root_num_spectrum)
print(f"After recapitation, the max number of roots is {recap_root}")

print("As for throwing mutations...")
print(f"Before there were {neg_ts.num_mutations}; after, {ts.num_mutations}. ")