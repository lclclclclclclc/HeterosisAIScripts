#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 13:50:21 2020

@author: egibson
"""

import pyslim
import numpy as np

def ancestry_p_varies(ts, source_popn, time_since_adm, model): #pop=source pop
    if model == 1:
        recip_popn = 4
    elif model == 0:
        recip_popn = 3

    time_just_before_adm = time_since_adm + 1  # in "generations ago"
    time_just_after_adm = time_since_adm       # during the adm generation, SLiM remembers individuals after the admixture event happens

    # collect node ids for the relevant samples
    source_samps_just_before = [ts.node(n).id for n in ts.samples(population=source_popn)
                                if (ts.node(n).time == time_just_before_adm)]
    # TODO: [perf] consider not looping through same ts.samples() three times
    recip_samps_just_before = [ts.node(n).id for n in ts.samples(population=recip_popn)
                               if (ts.node(n).time == time_just_before_adm)]
    recip_samps_just_after = [ts.node(n).id for n in ts.samples(population=recip_popn)
                              if (ts.node(n).time == time_just_after_adm)]
    recip_samps_today = [ts.node(n).id for n in ts.samples(population=recip_popn)
                         if (ts.node(n).time == 0)]
    # Find the source samples that introgressed into recip by looking at the
    # parents of the recipient population right after admixture.
        #  this is a very slow way to find just a few remaining parents (vs. just looking at first tree as below)
        # p3_first_parents = [ts.first().parent(i) for i in recip_samps_just_after]
    recip_parents_whole_chr = [[t.parent(i) for i in recip_samps_just_after] for t in ts.trees(sample_lists=True)]
    recip_parents = np.unique(recip_parents_whole_chr)
    # all of the recip parents should be samples, not untracked nodes
    assert np.all([ts.node(p).is_sample() for p in recip_parents])
    recip_parents_from_source = [p for p in recip_parents if (p in source_samps_just_before)]
    recip_parents_from_recip = [p for p in recip_parents if (p in recip_samps_just_before)]
    # all the parents must come from either p1 or p3
    assert np.all(recip_parents_from_source + recip_parents_from_recip == recip_parents)

    tree_p = [sum([t.num_tracked_samples(u) for u in recip_parents_from_source]) / len(recip_samps_today)
              for t in ts.trees(tracked_samples=recip_samps_today, sample_lists=True)]
    # TODO: check on ability to make u a list.  suggested in docs

    # etc: in my small sim at least this is all zero all the time
    # What is this trying to do?
    # # my original guess:
    # Return the fraction of all the "living"
    # (i.e. alive at end of simulation) recipient haplotypes that trace their
    # ancestry to the introgressing population.
    # This should be 10% in a "null" model, as the introgression event makes
    # the recipient population 90% recipient, 10% donor (in model 0, 90% pop3, 10% pop1).
    # # what I think now is largely the same
        # tree_p has length of ts.num_trees
        # for each tree (chr interval), give fract [as above].
        # later in ancestry_position_writeout we can see where there is more/less
    return tree_p


dir_stem = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/"
DIR_tree = dir_stem + "output/trees/"
tree_fn = "sgcb_m0_sexA"
tsrc = pyslim.load(DIR_tree + tree_fn + "_with_neutrals" + ".trees")
trp = ancestry_p_varies(tsrc, 1, 1000, 0)