#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 22:33:20 2020

@author: xinjunzhang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:44:17 2019

This script takes inputs from command line that specify the following parameters:

genomic segment, demographic model, population growth pattern in recipient population,
dominance coefficient for deleterious mutations, whether a human hs relationship is used,
selection coefficient for the adaptive mutation, scaling factor, and the number of simulation replicates

The script updates a slim script according to the above parameters, runs slim program, and extracts
adaptive introgression summary statistics in non-overlapping 50kb windows

@author: xinjunzhang
"""

import pyslim, os, random, argparse
import numpy as np
from multiprocessing import Manager, Pool

import tree_tools as tt  # etc


parser = argparse.ArgumentParser(description="A script for running slim and computing summary statistics in 50kb windows across a given chromosome")
parser.add_argument('-g', '--gene', action="store", dest="gene_id",
                        help="which simulation batch, default: 1; range 1-26",
                        default=1, type=int)
parser.add_argument('-het', '--dominance', action="store", dest="dominance_id",
                        help="dominance index, default: 1; range 0-100 (h=0-1); if running neutral model, value=200;",
                        default=1, type=int)
parser.add_argument('-m', '--model', action="store", dest="model_id",
                        help="model index, default: 0; 0=m0, 1=mh",
                        default=0, type=int)
parser.add_argument('-p', '--popsize', action="store", dest="growth_id",
                        help="growth index, default: 4; range:1-4",
                        default=4, type=int)
parser.add_argument('-d', '--hs', action="store", dest="hs_id",
                        help="growth index, default: 0; 0=use dominance input, 1=use hs relationship",
                        default=0, type=int)
parser.add_argument('-n', '--nscale', action="store", dest="nscale_id",
                        help="scaling factor index, default: 10",
                        default=10, type=int)
parser.add_argument('-s', '--selcoeff', action="store", dest="selcoeff_id",
                        help="adaptive mutation selection strength index, default: 0; range: 0-1 (s=0-0.1)",
                        default=0, type=int)
parser.add_argument('-r', '--rep', action="store", dest="numrep_id",
                        help="number of simulation replicates, default: 200",
                        default=200, type=int)
args = parser.parse_args()

whichgene = int(args.gene_id)
dominance = round(float(args.dominance_id)/100,2 ) #convert h-index to h value: 50 -> 0.5
model = int(args.model_id)
growth = int(args.growth_id)
hs = int(args.hs_id)
nscale = int(args.nscale_id)
m4s = float(args.selcoeff_id/100) #convert s-index to s: 1 -> 0.01
num_reps = int(args.numrep_id)

#sample command: python3 run_slim_get_stats.py -g 1 -h 0 -m 1 -p 4 -d 0 -n 10 -s 1 -r 200


def calc_p1ancestry (treepath, source_popn, t_sinceadm, model):
    ts = pyslim.load(treepath)
    any_ancestry = ancestry_p_varies(ts, source_popn, t_sinceadm, model)
    meanp1 = sum(any_ancestry)/len(any_ancestry)
    # This doesn't make much sense to me.  len(any_ancestry) equals ts.num_trees
    # So this is a per tree average.
    return meanp1

def ancestry_p_varies(ts, source_popn, time_since_adm, model):
    # TODO: surely source_popn and time_since_adm could all be informed by model
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


def ancestry_position_writeout(treepath, ancestry_filename, source_popn, t_sinceadm, model):
    ts = pyslim.load(treepath)
    starts=[]
    ends=[]

    for x in ts.trees():
        starts.append(x.interval[0])
        ends.append(x.interval[1])

    outfile = open(ancestry_filename, 'w')
    outfile.write('start,end,ancestry\n')

    # etc: I zero'd this b/c it was the second call to anc_p_var which takes a while
    # TODO: make only one call to anc_p_var in whole code
    # p1ancestry = ancestry_p_varies(ts, source_popn, t_sinceadm, model)
    p1ancestry = np.zeros(len(starts))

    for start, end, anc in zip(starts, ends, p1ancestry):
        outfile.write('{0},{1},{2}\n'.format(start, end, anc))

    outfile.close()

# No earthly idea why this is implemented like this.
# TODO: so much
def calc_ancestry_window (ancestry_file,len_genome):
    infile = open(ancestry_file,'r')
    end_pos = []
    ancestry = []
    line_counter=0
    for line_counter, line in enumerate(infile):
        fields = line.split(',')
        if fields[0] != "start":
            end_pos.append(float(fields[1]))
            ancestry.append(float(fields[2]))
    infile.close()
    allpos_bin = np.linspace(0,len_genome,int(len_genome/50000)) #windows of every 50kb

    endpos_digitized = np.digitize(end_pos, allpos_bin)
    end_pos = np.array(end_pos)
    ancestry = np.array(ancestry)

    anc_window = []
    anc_pos = []
    for w in range(1,100):
        these_pos = end_pos[endpos_digitized==w]
        these_anc = ancestry[endpos_digitized==w]

        if(len(these_pos))>0:
            anc_window.append(np.mean(these_anc))
            anc_pos.append(these_pos)
        else:
            anc_window.append(float('nan'))
            anc_pos.append(these_pos)

    return anc_window


# # NEVER CALLED
# import glob
# import matplotlib.pyplot
# # just copied from SLiM manual 17.5
# def ancestry_local (treepath):
#     starts, ends, subpops = [], [], []
#     ts = pyslim.load(treepath)
#     for tree in ts.trees(sample_counts=True):
#         subpop_sum, subpop_weights = 0, 0
#         for root in tree.roots:
#             leaves_count = tree.num_samples(root) - 1
#             subpop_sum += tree.population(root) * leaves_count
#             subpop_weights += leaves_count
#         starts.append(tree.interval[0])
#         ends.append(tree.interval[1])
#         subpops.append(subpop_sum / float(subpop_weights))
#     x = [x for pair in zip(starts, ends) for x in pair]
#     y = [x for x in subpops for _ in (0, 1)]
#     matplotlib.pyplot.plot(x, y)
#     matplotlib.pyplot.show()
#     return x,y #x=genome positions; y = ancestry
# # NEVER CALLED
# def write_ancestry (DIR_tree, output_anc_file):
#     tree_all = glob.glob(DIR_tree+'*.trees')
#     with open(output_anc_file, 'w') as outfile:
#         for file in tree_all:
#             x,y = ancestry_local (file)
#             for item in x:
#                 outfile.write("%s\t" % item)
#             outfile.write("\n")
#             for item in y:
#                 outfile.write("%s\t" % item)
#             outfile.write("\n")


def load_data_slim(file_path,len_genome,adm_gen,end_gen): # load slim's output
# TODO: etc: why does freqp4_b/a just get overridden a bunch?.. It's zero'd out later ln. 184ish
    pos_den, hapMat_den,freqp4_before,freqp4_after = get_pos_hap(file_path,'p1',len_genome,end_gen)
    pos_afr, hapMat_afr,freqp4_before,freqp4_after = get_pos_hap(file_path,'p2',len_genome,end_gen)
    pos_nonafr, hapMat_nonafr,freqp4_before,freqp4_after = get_pos_hap(file_path,'p3',len_genome,end_gen)
    # pos_preadm, hapMat_preadm,freqp4_before,freqp4_after = get_pos_hap(file_path,'p3',len_genome,adm_gen)
# TODO: either remove or utilize preadm info.
    pos_preadm = 'pointless'
    hapMat_preadm = 'pointless2'

    # TODO: etc: just really wrong....?
    # pos_den, hapMat_den,freqp4_before,freqp4_after = get_pos_hap(file_path,'p2',len_genome,end_gen)
    # pos_afr, hapMat_afr,freqp4_before,freqp4_after = get_pos_hap(file_path,'p1',len_genome,end_gen)
    # pos_nonafr, hapMat_nonafr,freqp4_before,freqp4_after = get_pos_hap(file_path,'p4',len_genome,end_gen)
    # pos_preadm, hapMat_preadm,freqp4_before,freqp4_after = get_pos_hap(file_path,'p2',len_genome,adm_gen)

    return pos_den, hapMat_den,pos_afr, hapMat_afr,pos_nonafr, hapMat_nonafr,pos_preadm, hapMat_preadm,freqp4_before,freqp4_after


def get_pos_hap(file_path,pop_id,len_genome,gen_time): #get pos and hapMat for a given pop from slim output
# etc: no! not in slim output for mod0 code
#   #OUT and SM header not in the output file when filepath is provided.  see manual 25.2.2
#   The first line is a header in the same format as for outputSample(), as described in the previous section. The output type code here, SM, represents “sample, MS format”. The outputMSSample() method allows output to be sent to a file, with the optional filePath argument. In this case, the #OUT: header line is not emitted, since it would not be conformant with the MS data format specification.
    infile = open(file_path,'r')
    end=0
    while end==0:
        line = infile.readline()
        if line[0:5]=='#OUT:': #output lines  # TODO: etc: lol no lines are marked #OUT:
            fields = line.split()
            out_type = fields[2]
            pop = fields[3]
            gen = fields[1]
            if out_type=='SM' and pop==pop_id and int(gen) == gen_time: #ms lines
                num_indiv = int(fields[4])
                infile.readline() # skip //
                infile.readline() # skip segsites
                pos = (np.array(infile.readline().split()[1:]).astype(float) * len_genome).astype(int)
                mult_mut_pos = find_mult_mut_pos(pos)+1
                pos = np.delete(pos,mult_mut_pos)
                hapMat = np.zeros((num_indiv,len(pos)),dtype=int)
                for indiv in range(0,num_indiv):
                    hap = np.array(list(infile.readline())[:-1]).astype(int)
                    hap = np.delete(hap,mult_mut_pos)
                    hapMat[indiv] = hap
                freqp4_before = 0  # why is this just zero'd?  also overwritten in load_data_slim.. but then output??
                freqp4_after = 0
                end=1
    infile.close()
    return pos, hapMat,freqp4_before,freqp4_after


def find_mult_mut_pos(pos): #find repeating mutations and remove them
    dist = np.array([pos[i+1]-pos[i] for i in range(0,len(pos)-1)])
    mult_mut_pos = np.where(dist==0)[0]
    return mult_mut_pos

def find_ai_site (segfile): #find an exon in the mid-range of the segment to insert AI mutation
    segs = open(segfile)
    starts = []
    ends = []
    total = 0
    for line_counter,line in enumerate (segs):
        if line[0:4]=="exon":
            fields = line.split()
            if (int(fields[1]) >= 2200000) & (int(fields[2]) <=2800000):
                starts.append(fields[1])
                ends.append(fields[2])
                total+=1
    any_exon = random.choice(range(0,total))
    window_start = int(starts[any_exon])
    window_end = int(ends[any_exon])
    segs.close()
    return window_start,window_end #return exon start and end position


def insert_anc_alleles(allpos, pos, haps):
    new_haps = np.zeros((haps.shape[0], allpos.size))
    insertidc = np.isin(allpos, pos)
    new_haps[:, insertidc] = haps
    return allpos, new_haps

def calc_derived_freq (pop_hap):
    popfreq = np.sum(pop_hap, axis=0)
    popfreq = popfreq/ float(pop_hap.shape[0])
    return popfreq

def vSumFunc(other_hap, currentArchi,p1_hapw):
    current_hap = np.array([p1_hapw[currentArchi,]])
    div = np.zeros(other_hap.shape)
    ones = np.ones((other_hap.shape[0],1))
    current_hap = current_hap
    current_hap_extended = np.dot(ones, current_hap)
    div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
    return np.add.reduce(div, 1)


def calc_stats (file_path,len_genome,adm_gen,end_gen):
    pos_den, hapMat_den,pos_afr, hapMat_afr,pos_nonafr, hapMat_nonafr,pos_preadm, hapMat_preadm,freqp4_before,freqp4_after = load_data_slim(file_path,len_genome,adm_gen,end_gen)

    p1_pos = pos_den
    p2_pos = pos_afr
    p3_pos = pos_nonafr
    p1_hap = hapMat_den
    p2_hap = hapMat_afr
    p3_hap = hapMat_nonafr

    all_pos = np.unique(np.concatenate((p1_pos,p2_pos,p3_pos)))

    p1_pos,p1_hap = insert_anc_alleles(all_pos,p1_pos,p1_hap)
    p2_pos,p2_hap = insert_anc_alleles(all_pos,p2_pos,p2_hap)
    p3_pos,p3_hap = insert_anc_alleles(all_pos,p3_pos,p3_hap)

    allpos_bin = np.linspace(0,len_genome,int(len_genome/50000)) #windows of every 50kb

    allpos_digitized = np.digitize(all_pos, allpos_bin)

    Dstat_list = []
    fD_list = []
    Het_list = []
    divratioavg_list = []
    Q_1_100_q95_list =[]
    Q_1_100_q90_list=[]
    Q_1_100_max_list=[]
    U_1_0_100_list =[]
    U_1_20_100_list =[]
    U_1_50_100_list =[]
    U_1_80_100_list =[]

    pos_start = []
    pos_end = []


    for w in range(1,100):  # etc: hardcoded here, but above was int(len_genome/50000)
        these_pos = all_pos[allpos_digitized==w]

        if len(these_pos)>1:
            pos_start.append(min(these_pos))
            pos_end.append(max(these_pos))

            these_pos_idx = np.nonzero(np.in1d(all_pos,these_pos))[0]
            p1_hapw = p1_hap[:,these_pos_idx]
            p2_hapw = p2_hap[:,these_pos_idx]
            p3_hapw = p3_hap[:,these_pos_idx]

            p1_freqw = calc_derived_freq (p1_hapw)
            p2_freqw = calc_derived_freq (p2_hapw)
            p3_freqw = calc_derived_freq (p3_hapw)

        # D-stat
            abbavecw = (1.0 - p2_freqw)*p3_freqw*p1_freqw
            babavecw = p2_freqw*(1.0 - p3_freqw)*p1_freqw
            abbacountsw = np.sum(abbavecw)
            babacountsw = np.sum(babavecw)
            if (abbacountsw + babacountsw > 0):
                Dstatw = (abbacountsw - babacountsw) / (abbacountsw + babacountsw)
            else:
                Dstatw = float('nan')
            Dstat_list.append(Dstatw)

        # fD
            checkfd1 = (p3_freqw > p1_freqw)
            abbafd1 = (1.0 - p2_freqw)*p3_freqw*p3_freqw
            babafd1 = p2_freqw*(1.0 - p3_freqw)*p3_freqw
            checkfd2 = (p3_freqw < p1_freqw)
            abbafd2 = (1.0 - p2_freqw)*p1_freqw*p1_freqw
            babafd2 = p2_freqw*(1.0 - p1_freqw)*p1_freqw
            abbafd = checkfd1 * abbafd1 + checkfd2 * abbafd2
            babafd = checkfd1 * babafd1 + checkfd2 * babafd2
            abbafdcounts = np.sum(abbafd)
            babafdcounts = np.sum(babafd)
            if (abbafdcounts + babafdcounts > 0):
                fD = (abbacountsw - babacountsw) / (abbafdcounts - babafdcounts)
            else:
                fD = float('nan')
            fD_list.append(fD)

        # Heterozygosity
            hetvec = 2 * p3_freqw * (1.0 - p3_freqw)
            Het = np.sum(hetvec) /50000  # etc: again, hardcoded window length
            Het_list.append(Het)


            # divratio = []
            # # for archi in range(p1_hapw.shape[0]): #iterate over 0-99 haps; 100 total)
            # for archi in range(0, 3): # TODO: etc: put this back.  just smaller for testing p1_hapw.shape[0]): #iterate over 0-99 haps; 100 total
            #     divarchintro = vSumFunc(p3_hapw, archi,p1_hapw)
            #     divarchintro = divarchintro.astype("float")
            #     divarchnonintro = vSumFunc(p2_hapw, archi,p1_hapw)

            #     divarchnonintro = divarchnonintro.astype("float")
            #     for comb in itertools.product(divarchintro,divarchnonintro):
            #         if comb[1] != 0:
            #             divratio.append(comb[0]/comb[1])
            # divratioavg = float(sum(divratio)) / float(len(divratio))
            # divratioavg_list.append(divratioavg)
            divratioavg_list.append(float('nan'))


            ArcHomoDer = (p1_freqw == 1)
            NonAdm_1 = (p2_freqw < 0.01)
            ArcHomoDerANDNonAdm_1 = (ArcHomoDer & NonAdm_1)
            DerFreqs_NonAdm_1 = p3_freqw[np.where(ArcHomoDerANDNonAdm_1 == True)]
            if DerFreqs_NonAdm_1.size > 0:
                Q_1_100_q95 = np.percentile(DerFreqs_NonAdm_1,95)
                Q_1_100_q90 = np.percentile(DerFreqs_NonAdm_1,90)
                Q_1_100_max = np.max(DerFreqs_NonAdm_1)
            else:
                Q_1_100_q95 = float('nan')
                Q_1_100_q90 = float('nan')
                Q_1_100_max = float('nan')

            Q_1_100_q95_list.append(Q_1_100_q95)
            Q_1_100_q90_list.append(Q_1_100_q90)
            Q_1_100_max_list.append(Q_1_100_max)

            U_1_0_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0) )
            U_1_20_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.2) )
            U_1_50_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.5) )
            U_1_80_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.8) )

            U_1_0_100 = np.sum(U_1_0_100)
            U_1_20_100 = np.sum(U_1_20_100)
            U_1_50_100 = np.sum(U_1_50_100)
            U_1_80_100 = np.sum(U_1_80_100)

            U_1_0_100_list.append(U_1_0_100)
            U_1_20_100_list.append(U_1_20_100)
            U_1_50_100_list.append(U_1_50_100)
            U_1_80_100_list.append(U_1_80_100)
        else:
            Dstat_list.append(float('nan'))
            fD_list.append(float('nan'))
            Het_list.append(float('nan'))
            divratioavg_list.append(float('nan'))
            Q_1_100_q95_list.append(float('nan'))
            Q_1_100_q90_list.append(float('nan'))
            Q_1_100_max_list.append(float('nan'))
            U_1_0_100_list.append(float('nan'))
            U_1_20_100_list.append(float('nan'))
            U_1_50_100_list.append(float('nan'))
            U_1_80_100_list.append(float('nan'))


    return pos_start, pos_end, freqp4_before, freqp4_after, Dstat_list, fD_list, Het_list, divratioavg_list, Q_1_100_q95_list, Q_1_100_q90_list, Q_1_100_max_list, U_1_0_100_list, U_1_20_100_list, U_1_50_100_list, U_1_80_100_list


def update_par_file(temp_par, new_par, model, growth, dominance,
                    nscale, m4s, hs, insert_ai, sex, trees_filename, region_filename):

    oldfile = open(temp_par)
    newfile = open(new_par,'w')
    line_counter=0
    for line_counter, line in enumerate(oldfile):
        fields = line.split()

        if model == 0: # etc: only implementing m0 rn
        # TODO: include hs?
        # TODO: calculate timepoints using adm_gen and end_gen??
            if line_counter == 1:
                fields[1] = str(dominance)  # irrelevant in neutral model
            elif line_counter == 2:
                fields[1] = str(nscale)
            elif line_counter == 3:
                fields[1] = str(m4s)
            if m4s == 2:  # neutral model
                if line_counter == 15:  # region info file
                    fields[2] = 'readFile("' + region_filename + '");'
                elif line_counter == 29:  # initializeSex
                    if sex is None:  # comment out the call
                        fields[0] = '// ' + fields[0]
                    else:  # initializeSex as autosome or Xchr ("A" or "X")
                        fields[1] = '"' + str(sex) + '"'
                elif line_counter == 40:  # p1/p2 split
                    fields[0] = str(int(2))
                elif line_counter == 43:  # remember p1/p2 split
                    fields[0] = str(int(2))
                elif line_counter == 47:  # AI variant emerges
                    fields[0] = str(int(100/nscale + 2))
                elif line_counter == 50:  # locus for AI variant
                    fields[1] = str(int(insert_ai))
                elif line_counter == 55:  # loop to check on AI variant
                    fields[0] = str(int(100/nscale + 2)) + ":"
                elif line_counter == 68:  # locus for AI variant in loop
                    fields[1] = str(int(insert_ai))
                elif line_counter == 73:  # p2/p3 split
                    fields[0] = str(int(10000/nscale))
                elif line_counter == 76:  # remember p2/p3 split
                    fields[0] = str(int(10000/nscale))
                elif line_counter == 80:  # admixture generation early()
                    fields[0] = str(int(20000/nscale))  # TODO: replace by adm_gen
                elif line_counter == 90:  # admixture generation late()
                    fields[0] = str(int(20000/nscale)) # TODO: replace by adm_gen
                elif line_counter == 96:  # final generation
                    fields[0] = str(int(30000/nscale))  # TODO: replace by end_gen
                elif line_counter == 100:  # write out .trees
                    fields[0] = 'sim.treeSeqOutput("' + trees_filename + '");'

            elif m4s != 2:   # recessive deleterious "negative" background model
                if line_counter == 23:  # region info file
                    fields[2] = 'readFile("' + region_filename + '");'
                elif line_counter == 45:  # initializeSex
                    if sex is None:  # comment out the call
                        fields[0] = '// ' + fields[0]
                    else:  # initializeSex as autosome or Xchr ("A" or "X")
                        fields[1] = '"' + str(sex) + '"'
                elif line_counter == 57:  # p1/p2 split
                    fields[0] = str(int(100000/nscale))
                elif line_counter == 60:  # remember p1/p2 split
                    fields[0] = str(int(100000/nscale))
                elif line_counter == 64:  # AI variant emerges
                    fields[0] = str(int(100/nscale) + int(100000/nscale))
                elif line_counter == 67:  # locus for AI variant
                    fields[1] = str(int(insert_ai))
                elif line_counter == 73:  # loop to check on AI variant
                    fields[0] = str(int(100/nscale) + int(100000/nscale))+":"
                elif line_counter == 86:  # locus for AI variant in loop
                    fields[1] = str(int(insert_ai))
                elif line_counter == 91:  # p2/p3 split
                    fields[0] = str(int(110000/nscale))
                elif line_counter == 94:  # remember p2/p3 split
                    fields[0] = str(int(110000/nscale))
                elif line_counter == 98:  # admixture generation early()
                    fields[0] = str(int(120000/nscale))  # TODO: replace by adm_gen
                elif line_counter == 107:  # admixture generation late()
                    fields[0] = str(int(120000/nscale))  # TODO: replace by adm_gen
                elif line_counter == 113:  # final generation
                    fields[0] = str(int(130000/nscale))  # TODO: replace by end_gen
                elif line_counter == 117:  # write out .trees
                    fields[0] = 'sim.treeSeqOutput("' + trees_filename + '");'

        elif model ==1:   #modelh
            if line_counter==1:
                fields[1] = str(int(growth))+");"
            elif line_counter==2:
                fields[1] = str(dominance)+");"
            elif line_counter==3:
                fields[1] = str(hs)+");"
            elif line_counter==4:
                fields[1] = str(nscale)+");"
            elif line_counter==5:
                fields[1] = str(m4s)+"*n);"
            elif line_counter==23:
                fields[2] = 'readFile("' + region_filename + '");'
            elif line_counter==77:
                fields[0] = "1:"+str(int(89000/nscale))
            elif line_counter==90:
                fields[0] = str(int(73000/nscale))
            elif line_counter==98:
                fields[0] = str(int(74000/nscale))
            elif line_counter==102:
                fields[1] = str(int(insert_ai))+");"
            elif line_counter==106:
                fields[0] = str(int(74000/nscale)) + ":"
            elif line_counter==118:
                fields[1] = str(int(insert_ai))+");"
            elif line_counter==124:
                fields[0] = str(int(83400/nscale))
            elif line_counter==128:
                fields[0] = str(int(86960/nscale))
            elif line_counter==133:
                fields[0] = str(int(87400/nscale - 1))
            elif line_counter==144:
                fields[0] = str(int(87400/nscale))
            elif line_counter==151:
                fields[0] = str(int(88080/nscale))
            elif line_counter==158:
                fields[0] = str(int(88080/nscale))+ ":"+str(int(89000/nscale))
            elif line_counter==183:
                fields[0] = str(int(89000/nscale))
            elif line_counter==187:
                fields[0] = 'sim.treeSeqOutput("' + trees_filename + '");'

        new_line=str()
        for item in fields:
            new_line = new_line+item+" "
        newfile.write(new_line+'\n')

    newfile.close()
    oldfile.close()



def run_slim_variable(n,q,r,dominance,nscale,m4s,model,growth,hs,insert_ai, sex):

    # set filenames
    region_name = region_all[r]
    region_info_filename = dir_stem + 'regions/sim_seq_info_' + str(region_name) + '.txt'
    trees_filename = dir_stem + 'output/trees/'+str(region_name)+'_m'+str(model)+'_sex'+str(sex)+'.trees'
    new_par = DIR_par +"par_"+region_name+str(dominance)+str(model)+ str(sex)+str(n)+".txt"
    ancestry_filename = DIR_anc + region_name+str(dominance)+ "_"+str(model)+ "_"+str(growth)+ "_"+str(m4s)+ "_"+str(hs) + "_"+str(n) + '.ancestry'

    segsize = 5000000  # nice that this is here, but hardcoded everywhere else

    if model ==1:  # etc: not handling this model yet
        if dominance != 2:
            temp_par = dir_stem + "slim/modelh_neg.txt"
        elif dominance == 2:
            temp_par = dir_stem + "slim/modelh_neu.txt"
        adm_gen = (87400-1)/nscale
        end_gen = 89000/nscale
        t_end = 1600/nscale -1
        popsize=41080/nscale # recipient population size at the end of simulation

        if growth ==1:
            popsize = 550/nscale
        elif growth ==2:
            popsize = 41080/nscale
        elif growth ==3:
            popsize = 7300/nscale
        elif growth ==4:
            popsize = 41080/nscale

    elif model == 0:
        # etc: why were these timepoints not used in writing par files????
        if dominance !=2:
            temp_par = dir_stem + "slim/ts_model0_neg.slim"
            adm_gen = 120000/nscale
            end_gen = 130000/nscale
        elif dominance == 2:
            temp_par = dir_stem + "slim/ts_model0_neu.slim"
            # recap'ing obviates need for 10k of SLiM burn-in in neutral model
            adm_gen = 20000/nscale
            end_gen = 30000/nscale

        # etc: why were these values not used in writing par files????
        t_end = 10000 / nscale  # etc: this means generations elapsed from admixture to end of simulation (present)
        popsize = 1000 / nscale  # size of p3 (as split off from p2 in mod0)

    update_par_file(temp_par, new_par, model, growth, dominance,
                    nscale, m4s, hs, insert_ai, sex, trees_filename+'.orig',
                    region_info_filename)

    # this is the slim output file in MS format (0x and 1s for haplotypes)
    # TODO: Move away from writing/reading these.
    # See below to extract from ts after throwing neutrals?
    # https://tskit.readthedocs.io/en/latest/python-api.html#tskit.TreeSequence.haplotypes !!!!
    # https://tskit.readthedocs.io/en/latest/python-api.html#tskit.TreeSequence.variants
    slim_output = DIR_out +'OUT_'+region_name+str(sex)+str(m4s)+ str(n)+".txt"

    # Run the SLiM simulation!
    os.system('slim %s > %s' %(new_par, slim_output))

    # Overlay neutral mutations onto TreeSequence
    # TODO: variable-ize initial_Ne (size of p1 at beginning of sim)
    tt.throw_neutral_muts(trees_filename+'.orig', region_info_filename,
                           neu_or_neg=dominance, n_scale=nscale)

    # Load ts, write ancestry file, read ancestry file
        # meanp1 works to my satisfaction!  However, the rest is crazy afaict
        # TODO: figure out / improve the windows and the file writing/reading
    if model==1:
        meanp1 = calc_p1ancestry(trees_filename, 2, t_end, model)
        ancestry_position_writeout(trees_filename, ancestry_filename, 2, t_end, model) #write out ancestry info
    elif model==0:
        meanp1 = calc_p1ancestry(trees_filename, 1, t_end, model)
        ancestry_position_writeout(trees_filename, ancestry_filename, 1, t_end, model)
    anc_window = calc_ancestry_window(ancestry_filename, segsize) #get mean ancestry per 50kb window

    # Calculate other statistics via loading the std-output from SLiM sim
    pos_start,pos_end,freqp4_before,freqp4_after,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list = calc_stats(slim_output,segsize,adm_gen,end_gen)

    q.put([n,insert_ai,growth,meanp1,pos_start,pos_end,freqp4_before,freqp4_after,anc_window, Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list])
    #other parameter info are stored in the output file name

    # os.system('rm '+slim_output)
    # os.system('rm '+treepath)
    # os.system('rm '+new_par)
    # os.system('rm '+ancestry_file)



def write_to_file(windowfile_name,q):
    windowfile = open(windowfile_name,'w')

    while 1:
        q_elem = q.get()

        if q_elem=='kill': # break if end of queue
            print ('END OF SIMULATIONS')
            break

        [n,insert_ai,growth,meanp1,pos_start,pos_end,freqp4_before,freqp4_after,anc_window,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list] = q_elem
        for i in range(len(Dstat_list)):
            windowfile.write("%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (n,insert_ai,growth,meanp1,pos_start[i],pos_end[i],freqp4_before,freqp4_after,anc_window[i],Dstat_list[i], fD_list[i], Het_list[i], divratioavg_list[i],Q_1_100_q95_list[i],Q_1_100_q90_list[i],Q_1_100_max_list[i],U_1_0_100_list[i],U_1_20_100_list[i],U_1_50_100_list[i],U_1_80_100_list[i]))

        windowfile.flush()

    windowfile.close()


#################################################################################
if __name__=='__main__':
    # etc: sex param takes None, 'A', or 'X'
    sex = 'A' #'A'
    #TODO: etc: these were originally commented out.  use to change defaults.
    whichgene = 14#+10  #1  15 was for project.  X is 25
    # model = 1 # 1=modelh; 0=model0 #define these two with parseargument
    #growth = 4
    #hs = 0 #0 = recessive or neutral; 1 = hs relationship
    dominance = 0 #if 0, run the deleterious recessive model #if 2, run the neutral model
    nscale = 10 #define scaling factor
    m4s = 0.01 #adaptive selection strength
    num_reps=1 #number of simulations per region
    region_all = ["chr11max","chr19region","chr3region","galnt18","hla","hyal2",
                  "krt71","nlrc5","oca2","pde6c","pou2f3","rnf34","sema6d","sgcb",
                  "sgcz","sipa1l2","slc16a11","slc19a3","slc5a10","stat2","tbx15",
                  "tlr1610","tnfa1p3","txn", 'X-0-5M']

    dir_stem = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/"

    DIR_region = dir_stem + "regions/"
    DIR_anc = dir_stem + "output/ancestry/"
    DIR_out = dir_stem + "output/out/"
    DIR_tree = dir_stem + "output/trees/"
    DIR_par = dir_stem + "slim/"

    # or loop over genes here I suppose
    r = int(whichgene-1)

    region_name = region_all[r]

    # Find an exon in the middle-ish of the region...
    window_start, window_end = find_ai_site(DIR_region+"sim_seq_info_"+str(region_name)+".txt")
    # ...and put the AI variant in the middle of that exon.
    insert_ai = int((int(window_end)+int(window_start))/2)

    attempt_num = np.random.randint(5000)
    windowfile_name = dir_stem + "output/stats/20200630/"+region_name+"-dominance"+str(dominance)+"-model"+str(model)+"-sex"+str(sex)+"-hs"+str(hs)+"-ai"+str(m4s)+'-attempt' + str(attempt_num) + '_human_windows.txt'
    num_proc = 10
    manager = Manager()
    pool = Pool(processes=num_proc)
    q = manager.Queue()
    watcher = pool.apply_async(write_to_file,(windowfile_name,q))
    reps = range(0,num_reps)
    args_iterable = list(zip(reps,[q]*num_reps))
    for i in args_iterable:
        n=i[0]
        print(str(n))
        run_slim_variable(i[0],i[1],r,dominance,nscale,m4s,model,growth,hs,insert_ai, sex)

    q.put('kill')
    pool.close()
    pool.join()

    # peek at first line of stats file to check in on things
    with open(windowfile_name, 'r') as f:
        print("First line of stats file for this run is below.")
        print(f.readline())

    print("END OF SIMULATION")
