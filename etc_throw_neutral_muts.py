#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:51:27 2020

@author: egibson
"""

import pyslim
import tskit
import msprime

# read in tree seq from SLiM sim
DIR_tree = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/output/trees/"
tree_fn = "neg_m0"
tree_file = DIR_tree + tree_fn + ".trees"
neg_ts = pyslim.load(tree_file)

#TODO: sort out mutation rate variation / ratios
mut_rate_map = 10e-8

#TODO: sort out recombination rate variation / map
# https://msprime.readthedocs.io/en/latest/api.html?highlight=map#recombination-map-limitations
# L = int(from_ts.sequence_length)
# recomb_map = msprime.RecombinationMap.uniform_map(L, r, L)
# final_ts = mpsrime.simulate(from_ts=from_ts, recomb_map=recomb_map)
L = int(neg_ts.sequence_length)
r1 = 10e-7
r2 = 10e-11
recomb_map = msprime.RecombinationMap([0, 100, L], [r1, r2, 0])
# I think I need to be doing msprime.simulation(from_ts), not msp.mutate...
# https://msprime.readthedocs.io/en/latest/api.html?highlight=recombination#msprime.simulate

# throw mutations and output
ts = msprime.mutate(neg_ts, rate=mut_rate_map, keep=True, recombination_map=recomb_map)
ts.dump(DIR_tree + tree_fn + "_with_neutrals" + ".trees")

# TODO: what format for existing code?  ms I think