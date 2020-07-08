#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:51:27 2020

@author: egibson
"""

import pyslim
import msprime


def throw_neutral_muts(tree_file, region_info_file, neu_or_neg=0,
                                  n_scale=10, initial_Ne=10000, verbose=False):
    """
    Recapitates .trees output from SLiM simulation and overlays neutral mutations.
    Writes out new .trees file with the neutral mutations.

    Parameters
    ----------
    tree_file : str
        Path to .trees file output from SLiM simulation.
    region_info_file : str
        Path to "sim_seq_info_[region].txt" file with recRate info.
    neu_or_neg : int, optional
        Null model of neutral variation or alt. model of rec. del. var.
        neu = 2; neg = 0;  per "dominance" parameter in original code.
        The default is 0.
    n_scale : int, optional
        Scaling factor. The default is 10.
    initial_Ne : int, optional
        Initial population size of founding population (p1). The default is 10000.
    verbose : bool, optional
        Print info about treeseqs before and after recapitation. The default is False.

    Output
    ------
    Writes out .trees file, recapitated, with neutral mutations overlayed.
    Currently OVERWRITES the original .trees file

    Returns
    ------
    ts : treeSeq

    """
    ## Set recombination map
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

    recomb_map = make_region_recombination_map(region_info_file)

    ## Load ts and recapitate
    slim_ts = pyslim.load(tree_file)
    # loaded with pyslim, so a SLiMTreeSeq object.
    # https://pyslim.readthedocs.io/en/latest/python_api.html#pyslim.SlimTreeSequence.recapitate
    # why is default Ne=1?
    n_p1 = initial_Ne / n_scale
    recap_ts = slim_ts.recapitate(recombination_map=recomb_map, Ne=n_p1)
    # TODO: decide about simplifying here?  Essentially sampling here.
    # How does the simplifying default to SLiM tree writing figure in?
    # I think this depends on how statistics are calculated

    ## Throw mutations and output
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

    # TODO: set base mutation rate by variable?
    mut_rate = 1.5e-8 * n_scale  # fine for neutral model;  will over-mutate exonic regions in neg
    if neu_or_neg == 0:  # meaning the "neg" model of recessive deleterious var.
        mut_rate *= 1 / (1 + 2.31)  # will under-mutate non-exons

    ts = pyslim.SlimTreeSequence(msprime.mutate(recap_ts, rate=mut_rate, keep=True))
    # TODO: Consider not overwritting original file
    # default is still to overwrite...
    out_name = tree_file
    # ... unless SLiM ts file has extension .orig
    if tree_file[-5:] == '.orig':
        out_name = tree_file[:-5]
    ts.dump(out_name)

    ## Messing around
    if verbose:
        import matplotlib.pyplot as plt
        import numpy as np

        # Mess around with original ts
        print(f"There are {slim_ts.num_individuals} individuals in the SLiM tree.")
        founders=slim_ts.first_generation_individuals()
        initial_pops = [slim_ts.individual(f).population for f in founders]
        print(f"Population(s) represented in the founders: {np.unique(initial_pops)}.")
        # numpy array of how long ago each individual was born
        # SLiM automatically remembers the individuals that comprise the first generation of any new subpopulation created with addSubpop(), for easy recapitation and other analysis.
        indiv_times = slim_ts.individual_times
        plt.figure()
        plt.hist(indiv_times)
        recorded_times = np.unique(indiv_times)
        print(f"When indivs. were added to the ts, in time ago: {recorded_times}.")
        root_num_spectrum = [t.num_roots for t in slim_ts.trees()]
        # plt.figure()
        # plt.hist(root_num_spectrum)
        max_roots = max(root_num_spectrum)
        print(f"Before recapitation, the max number of roots was {max_roots}.")

        # Mess around with full recap'd, mutated ts
        recap_root_num_spectrum = [t.num_roots for t in ts.trees()]
        recap_root = max(recap_root_num_spectrum)
        print(f"After recapitation, the max number of roots is {recap_root}.")

        print("As for throwing mutations...")
        print(f"before there were {slim_ts.num_mutations}; after, {ts.num_mutations}. ")

    return ts
