#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:51:27 2020

@author: egibson
"""

import pyslim
import msprime
import numpy as np

#%%
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


#%%
def translate_from_slim_pop(tree_seq, s_ids=None):
    popkey = [(p.id, p.metadata.slim_id) for p in tree_seq.populations()
              if p.metadata is not None]
    pops = [p for (p, s) in popkey]
    slims = [s for (p, s) in popkey]
    if s_ids is None:  # just get all populations with a slim_id
        s_ids = slims
    indices_of_slim_pops = np.where(np.isin(slims, s_ids))[0]
    return tuple(np.asarray(pops)[indices_of_slim_pops])


#%%
def sample_population_haplotypes(ts, popn_ids=(1,2,3), n_haps=100, check_loc=None):
    """
    Using a tree sequence, takes extant samples from three populations
    and generates both haplotype matrices (n_haps x [number of variants]) for each
    population and the chromosome position for each variant.

    Parameters
    ----------
    ts : TreeSeq
        Tree sequence from SLiM simulation.
    popn_ids : tuple or list, optional
        Three population ids from each of which to sample n_haps haplotypes.
        The default is (1,2,3).
    n_haps : int, optional
        Number of haplotypes to sample from each of the three populations.
        The default is 100.
    check_loc : int or float, optional
        If provided, will throw assertion error if this chromosome coordinate
        is not found in the variant positions. The default is None.

    Returns
    -------
    tuple
        3-tuple of haplotype matrices for each population, in order of popn_ids.
    np array
        Locations of each variant (matrix column) in SLiM chromosome coordinates.

    """
    ts_popn_ids = translate_from_slim_pop(ts, s_ids=popn_ids)

    def sample_extant_haplotypes(population, n_haps):
        # TODO: later, make sure this is a sample with the correct sex ratio.
        all_samps = np.asarray([ts.node(n).id for n in ts.samples(population=population)
                             if (ts.node(n).time == 0)])
        return np.random.choice(all_samps, size=n_haps, replace=False)

    ## Sample down to n_haps extant nodes per population
    # TODO: should I be sampling individuals and keeping both their nodes?
    #       perhaps not if I want to do same for sex chromosomes
    extant_samples = np.concatenate((sample_extant_haplotypes(ts_popn_ids[0], n_haps),
                                    sample_extant_haplotypes(ts_popn_ids[1], n_haps),
                                    sample_extant_haplotypes(ts_popn_ids[2], n_haps)))
    # reduce ts to only these samples so the matrix isn't too big
    sts = ts.simplify(samples=extant_samples, reduce_to_site_topology=True)
        # In the returned tree sequence, the node with ID 0 corresponds to samples[0], node 1 corresponds to samples[1], and so on. Besides the samples, node IDs in the returned tree sequence are then allocated sequentially in time order.

    ## Generate the haplotypes
    haplotype_matrix = sts.genotype_matrix()
    # flip to match other code: make rows individual haplotypes, columns -> vars
    hap_mat = np.transpose(haplotype_matrix)

    ## Split haplotype matrix into population-specific matrices
    p1_haps = hap_mat[:n_haps, :]
    p2_haps = hap_mat[n_haps:2*n_haps, :]
    p3_haps = hap_mat[2*n_haps:, :]
    # check that we've broken up the array without missing or duplicating any haps
    assert np.all(np.vstack((p1_haps, p2_haps, p3_haps)) == hap_mat)

    ## Get positions of mutations in chr coordinates
    # check that hap_mat contains all the sites
    assert sts.num_sites == np.shape(hap_mat)[1]
    positions = np.asarray([s.position for s in sts.sites()])
    if check_loc is not None:  # just check that the AI var is in the positions
        assert check_loc in positions

    return (p1_haps, p2_haps, p3_haps), positions


#%%
def calc_ancestry_frac_over_region(ts, source_popn, recip_popn, time_since_adm):
    """
    Calculates the fraction of ancestry from the source population seen in all
    extant samples from the recipient population.

    Parameters
    ----------
    ts : TreeSequence
        TreeSeq with samples from recip_popn taken at "generations ago"
        times 0, time_since_adm, and time_since_adm + 1,
        and samples from source_popn taken at time_since_adm + 1
    source_popn : int
        Population id number of source/introgressing population from SLiM sim.
    recip_popn : int
        Population id number of recipient/extant population from SLiM sim..
    time_since_adm : int
        Time of admixture in generations ago.  Assumes samples were remembered
        in late() stage of corresponding admixture generation within SLiM.

    Returns
    -------
    list of floats
        source ancestry fractions in each chromosome interval.
    list of tuples
        chromosome intervals (start, end) to which ancestry fractions correspond.

    """
    # check that the populations exist as expected in the original, recaped tree
    assert (source_popn, recip_popn) == translate_from_slim_pop(ts, (source_popn, recip_popn))

    # check that we got the recapitated treeseq
    assert max(t.num_roots for t in ts.trees()) == 1

    # TODO: recip_popn, source_popn, and time_since_adm could all be informed by model
    time_just_before_adm = time_since_adm + 1  # in "generations ago"
    time_just_after_adm = time_since_adm       # during the adm generation, SLiM remembers individuals after the admixture event happens

    def peri_admixture_sample_lists(tree_seq, source_popn, recip_popn, just_before, just_after):
        # translate slim population ids into the ts frame of reference
        source, recip = translate_from_slim_pop(tree_seq, (source_popn, recip_popn))
        # collect node ids for the relevant samples
        source_samps_just_before = [tree_seq.node(n).id for
                                    n in tree_seq.samples(population=source)
                                    if (tree_seq.node(n).time == just_before)]
        recip_samps_just_before = [tree_seq.node(n).id for
                                   n in tree_seq.samples(population=recip)
                                   if (tree_seq.node(n).time == just_before)]
        recip_samps_just_after = [tree_seq.node(n).id for
                                  n in tree_seq.samples(population=recip)
                                  if (tree_seq.node(n).time == just_after)]
        recip_samps_today = [tree_seq.node(n).id for
                             n in tree_seq.samples(population=recip)
                             if (tree_seq.node(n).time == 0)]
        return source_samps_just_before, recip_samps_just_before, recip_samps_just_after, recip_samps_today

    # Find the source samples that introgressed into recip by looking at the
    # parents of the recipient population right after admixture.
        #  this is a very slow way to find just a few remaining parents (vs. just looking at first tree as below)
        # p3_first_parents = [ts.first().parent(i) for i in recip_samps_just_after]
    # because it's slow, I'm going to try to simplify and see if that helps
    s_b, r_b, r_a, r_t = peri_admixture_sample_lists(ts, source_popn, recip_popn,
                                                     time_just_before_adm,
                                                     time_just_after_adm)
    samples_to_keep = np.concatenate((s_b, r_b, r_a, r_t))

    sts = ts.simplify(samples=samples_to_keep)
    # TODO: what precisely does reduce_to_site_topology do?
        # It certainly has an effect:  simplified tree goes from 15k trees to 6.7k trees.
    # After simplifying, now have new node ids for samples.  Need to split them out again.
        # Could (1) go by length of array chunks
        # (2) something with map_nodes=True in simplify call
        # (3) do same listcomps as above
        # (4) take the whole issue to SLiM and try to flag the intogressors specifically
    # For now, do option (3) via helper function.

    # these are being overwritten into new node_ids in the simplified ts world
    s_b, r_b, r_a, r_t = peri_admixture_sample_lists(sts, source_popn, recip_popn,
                                                     time_just_before_adm,
                                                     time_just_after_adm)
    assert len(r_t) > 0

    recip_parents_whole_chr = [[t.parent(i) for i in r_a] for t in sts.trees(sample_lists=True)]
    recip_parents = np.unique(recip_parents_whole_chr)
    # all of the recip parents should be samples, not untracked nodes
    assert np.all([sts.node(p).is_sample() for p in recip_parents])
    recip_parents_from_source = [p for p in recip_parents if (p in s_b)]
    recip_parents_from_recip = [p for p in recip_parents if (p in r_b)]
    # all the parents must come from either p1 or p3
    assert np.all(recip_parents_from_source + recip_parents_from_recip == recip_parents)

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
    # TODO: check on ability to make u a list.  suggested in docs
    if len(recip_parents_from_source) < 1:
        print("lost all source ancestry")
    elif len(recip_parents_from_recip) < 1:
        print("lost all recip ancestry")
    source_anc_sum_by_tree = [sum([t.num_tracked_samples(u) for u in recip_parents_from_source])
                              for t in sts.trees(tracked_samples=r_t, sample_lists=True)]
    source_anc_frac_by_tree = np.asarray(source_anc_sum_by_tree) / len(r_t)
    assert sts.num_trees == len(source_anc_frac_by_tree)
    # grab the windows too
    source_anc_tree_intervals = [t.interval for t in sts.trees()]
    return source_anc_frac_by_tree, source_anc_tree_intervals