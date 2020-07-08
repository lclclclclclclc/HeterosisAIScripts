#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:51:27 2020

@author: egibson
"""

import pyslim
import msprime

import numpy as np
import matplotlib.pyplot as plt

#%% set vars

dir_stem = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/"

## .trees output from SLiM sim
tree_fn = "sgcb_m0_sexA"
# tree_fn = "ts_recording_test_full"

## name of gene region simulated
region_name = "sgcb"


#%% read in tree seq from SLiM sim

DIR_tree = dir_stem + "output/trees/"
tree_file = DIR_tree + tree_fn + ".trees"
neg_ts = pyslim.load(tree_file)
# loaded with pyslim, so a SLiMTreeSeq object.
# This means I should be able to ts.recapitate via pyslim , which theoretically
# takes a recombination map (which msprime didn't recognize)


#%% messing around with original ts

print(f"There are {neg_ts.num_individuals} individuals in the SLiM tree.")
founders=neg_ts.first_generation_individuals()
initial_pops = [neg_ts.individual(f).population for f in founders]
print(f"Population(s) represented in the founders: {np.unique(initial_pops)}.")
#numpy array of how long ago each individual was born
# SLiM automatically remembers the individuals that comprise the first generation of any new subpopulation created with addSubpop(), for easy recapitation and other analysis.
indiv_times = neg_ts.individual_times
plt.figure()
plt.hist(indiv_times)
recorded_times = np.unique(indiv_times)
print(f"When indivs. were added to the ts, in time ago: {recorded_times}.")

root_num_spectrum = [t.num_roots for t in neg_ts.trees()]
# plt.figure()
# plt.hist(root_num_spectrum)
max_roots = max(root_num_spectrum)
print(f"Before recapitation, the max number of roots was {max_roots}.")


#%% Set recombination map

#TODO: sort out recombination rate variation / map
# https://pyslim.readthedocs.io/en/latest/tutorial.html#recapitation-with-a-nonuniform-recombination-map
# https://msprime.readthedocs.io/en/latest/api.html?highlight=map#recombination-map-limitations

# Let's use the sequence info doc with only recRate lines.  Feed to msprime.readHapMap
# https://msprime.readthedocs.io/en/stable/api.html?highlight=read_hapmap#msprime.RecombinationMap.read_hapmap

# Actually, let's not because it has to be subtlely changed.
# seq_info_parser.py is a protofunction returning a map.  Just need to integrate it.

# copied from run_slim_get_stats.main
DIR_region = dir_stem + "regions/"
region_file = DIR_region+"sim_seq_info_"+str(region_name)+".txt"

def make_region_recombination_map(region_filename):
    positions = []
    rates = []
    with open(region_filename, 'r') as file:
        for line in file:
            if "recRate" in line:
                components = line.split(" ")
                positions.append(int(components[1]))
                rates.append(float(components[2]))
    # adapted from https://pyslim.readthedocs.io/en/latest/tutorial.html#recapitation-with-a-nonuniform-recombination-map
    # step 1
    positions.insert(0, 0)
    # step 2
    rates.append(0.0)
    # step 3
    positions[-1] += 1

    recomb_map = msprime.RecombinationMap(positions, rates)
    return recomb_map


# L = int(neg_ts.sequence_length)
# r1 = 10e-15
# r2 = 10e-11
# recomb_map = msprime.RecombinationMap([0, 100, L], [r1, r2, 0])

recomb_map = make_region_recombination_map(region_file)

#%% Recapitate

# https://pyslim.readthedocs.io/en/latest/python_api.html#pyslim.SlimTreeSequence.recapitate
# TODO: sort out reading in number of individuals.. necessary? why is default Ne=1?
n=10000
n_scale=10
n_p1 = n/n_scale
recap_ts = neg_ts.recapitate(recombination_map=recomb_map, Ne=1000)

# TODO: decide about simplifying here?  Essentially sampling here.
# How does the simplifying default to SLiM tree writing figure in?
# I think this depends on how statistics are calculated


#%% Throw mutations and output
# TODO: Wait for mutationMaps to magically appear?
# Need mutation rate variation to ensure equal density of mutations along chromosome
# Mutation rate maps are coming to msprime... I think they're already implemented but not doc'd.
# Should I use a development version of msprime?  Mine does not have the .MutationMap class.
# see: https://github.com/tskit-dev/msprime/pull/920 , also 902, 711
#       https://github.com/tskit-dev/msprime/blob/master/msprime/mutations.py line 700ish

# mutation rate scaling notes from SLiM sim code:
#   // exon uses a mixture of syn and nonsyn at a 1:2.31 ratio (Huber et al.)
# 	// to achieve an overall mutation rate of 1.5e-8, need 2.31/(1+2.31) fraction of all muts to be nonsyn
# 	// i.e. ~0.6979 fraction of all muts should be the deleterious ones simulated here as "m1".
# 	// the remaining ~0.3021 fraction of all mut.s should be thrown later as neutrals.
# HOWEVER, this logic applies if we were only simulating exon regions.


# mut_rate = 1.5e-8 * n_scale  # fine for neutral model;  will over-mutate exonic regions in neg
mut_rate = 1.5e-8 * (1/(1+2.31)) * n_scale  # will under-mutate rest of genome

ts = pyslim.SlimTreeSequence(msprime.mutate(recap_ts, rate=mut_rate, keep=True))
ts.dump(DIR_tree + tree_fn + "_with_neutrals" + ".trees")


#%% Mess around with full recap'd, mutated ts

recap_root_num_spectrum = [t.num_roots for t in ts.trees()]
recap_root = max(recap_root_num_spectrum)
print(f"After recapitation, the max number of roots is {recap_root}.")

print("As for throwing mutations...")
print(f"before there were {neg_ts.num_mutations}; after, {ts.num_mutations}. ")