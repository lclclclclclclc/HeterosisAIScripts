#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:26:09 2020

@author: egibson
"""

import pyslim
import numpy as np

import tree_tools as tt

tree_file = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/output/trees/sgcb_m0_sexA.trees"
region_file = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/regions/sim_seq_info_sgcb.txt"


# ts = tt.throw_neutral_muts(tree_file, region_file)
ts = pyslim.load(tree_file)

# sanity check that ts is recapped
recap_root_num_spectrum = [t.num_roots for t in ts.trees()]
recap_root = max(recap_root_num_spectrum)
print(f"After recapitation, the max number of roots is {recap_root}.")
print(f"Num muts is {ts.num_mutations}.")

# My goal here is to immitate MS output of SLiM.
# https://tskit.readthedocs.io/en/latest/python-api.html?highlight=haplotypes#tskit.TreeSequence.simplify
# I need a list of variant positions for each population at end of sim
# and a LoL of 0/1 for pres/abs of each variant in each of 100 samples from each pop'n at end of sim

# first I should simplify to only the 100 from each pop'n at the end time.
# yes, we're taking 100 samples aka 50 people
# Make sure this is a random sample
# later, make sure this is a sample with the correct sex ratio.
# https://tskit.readthedocs.io/en/latest/python-api.html?highlight=haplotypes#tskit.TreeSequence.simplify

p1_id = 1
p2_id = 2
p3_id = 3

def sample_extant_haplotypes(population, n_haps=100):
    all_samps = np.asarray([ts.node(n).id for n in ts.samples(population=population)
                         if (ts.node(n).time == 0)])
    return np.random.choice(all_samps, size=n_haps, replace=False)


extant_samples = np.concatenate((sample_extant_haplotypes(p1_id),
                                sample_extant_haplotypes(p2_id),
                                sample_extant_haplotypes(p3_id)))

sts = ts.simplify(samples=extant_samples, reduce_to_site_topology=True)

# In the returned tree sequence, the node with ID 0 corresponds to samples[0], node 1 corresponds to samples[1], and so on. Besides the samples, node IDs in the returned tree sequence are then allocated sequentially in time order.
#   just visually, can see three distinct population groupings.
# columns are haplotypes, rows are variants
haplotype_matrix = sts.genotype_matrix()  # I don't understand how .haplotypes is different, but this works :)

# **note that these are 0,1,2 not 1,2,3
population_assignments = [sts.node(s).population for s in sts.samples()]

# flip to match other code: rows are individual haplotypes, columns are vars
hap_mat = np.transpose(haplotype_matrix)

# sanity checks
counts = np.sum(hap_mat, axis=0)
max(counts)  # this can be 300, meaning can have mono-morphic site of all 1s
min(counts)  # is not zero.  No all-ancestral sites

# split into populations
p1_haps = hap_mat[:100, :]
p2_haps = hap_mat[100:200, :]
p3_haps = hap_mat[200:, :]
# check that we've broken up the array without missing or duplicating any haps
assert np.all(np.vstack((p1_haps, p2_haps, p3_haps)) == hap_mat)

# Still need POSITIONS in chr coordinates
sts.sequence_length  # is still 5MB so that's good
# a Site object has Site.position that should be what I want.
# just need to figure out which sites are in my hapmat.  All, I believe
assert sts.num_sites == np.shape(hap_mat)[1]
positions = np.asarray([s.position for s in sts.sites()])
# great, that works but they're mostly floats.
# just check that the AI var is at an int position: 2325542
# sorted_positions = np.sort(positions)
insert_ai = 2325542
assert insert_ai in positions
# yes, it's there and there are some other "ints" as well.  Must be deleterious vars.



