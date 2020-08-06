import os
import random
import argparse
import numpy as np
from multiprocessing import Manager, Pool

import tree_tools as tt  # etc

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
dominance = round(float(args.dominance_id) / 100, 2)  # convert h-index to h value: 50 -> 0.5
model = int(args.model_id)
growth = int(args.growth_id)
hs = int(args.hs_id)
nscale = int(args.nscale_id)
m4s = float(args.selcoeff_id / 100)  # convert s-index to s: 1 -> 0.01
num_reps = int(args.numrep_id)

#sample command: python3 run_slim_get_stats.py -g 1 -h 0 -m 1 -p 4 -d 0 -n 10 -s 1 -r 200


# No earthly idea why this is implemented like this.
# TODO: scrap the file i/o but keep the digitization
def calc_ancestry_window(ancestry_fracs, intervals, num_windows=100, len_genome=None):

    start_positions = np.asarray([i[0] for i in intervals])
    end_positions = np.asarray([i[1] for i in intervals])
    interval_spans = end_positions - start_positions

    # TODO: actually fix the window size hardcoding
    if len_genome is None:  # just figure it out from intervals
        len_genome = max(end_positions)
    allpos_bin = np.linspace(0, len_genome, num_windows)  # windows of every 50kb
    endpos_digitized = np.digitize(end_positions, allpos_bin)

    ancestry_fracs = np.asarray(ancestry_fracs)
    anc_frac_by_window = []
    anc_windows = []
    for w in range(1, num_windows):
        window_mask = (endpos_digitized == w)
        if np.any(window_mask):
            window_start = min(start_positions[window_mask])
            window_end = max(end_positions[window_mask])
            anc_windows.append((window_start, window_end))

            these_spans = interval_spans[window_mask]
            these_anc = ancestry_fracs[window_mask]

            anc_frac_by_window.append(np.average(these_anc, weights=these_spans))
        else:
            anc_windows.append((float('nan'), float('nan')))
            anc_frac_by_window.append(float('nan'))

    return anc_frac_by_window, anc_windows


# Find an exon in the mid-range of the segment to insert AI mutation
def find_ai_site(segfile):
    segs = open(segfile)
    starts = []
    ends = []
    total = 0
    for line_counter, line in enumerate(segs):
        if line[0:4] == "exon":
            fields = line.split()
            if (int(fields[1]) >= 2200000) & (int(fields[2]) <= 2800000):
                starts.append(fields[1])
                ends.append(fields[2])
                total += 1
    any_exon = random.choice(range(0, total))
    window_start = int(starts[any_exon])
    window_end = int(ends[any_exon])
    segs.close()
    return window_start, window_end  # return exon start and end position


def insert_anc_alleles(allpos, pos, haps):
    new_haps = np.zeros((haps.shape[0], allpos.size))
    insertidc = np.isin(allpos, pos)
    new_haps[:, insertidc] = haps
    return allpos, new_haps


def calc_derived_freq(pop_hap):
    popfreq = np.sum(pop_hap, axis=0)
    popfreq = popfreq / float(pop_hap.shape[0])
    return popfreq


# Keep for divratio re-implementation
# def vSumFunc(other_hap, currentArchi,p1_hapw):
#     current_hap = np.array([p1_hapw[currentArchi,]])
#     div = np.zeros(other_hap.shape)
#     ones = np.ones((other_hap.shape[0],1))
#     current_hap = current_hap
#     current_hap_extended = np.dot(ones, current_hap)
#     div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
#     return np.add.reduce(div, 1)


def calc_stats(ts, sample_size, num_windows=100):
    (p1_hap, p2_hap, p3_hap), all_pos = tt.sample_population_haplotypes(ts, n_haps=int(sample_size))

    len_genome = ts.sequence_length
    allpos_bin = np.linspace(0, len_genome, num_windows)  # windows of every 50kb
    allpos_digitized = np.digitize(all_pos, allpos_bin)

    Dstat_list = []
    fD_list = []
    Het_list = []
    divratioavg_list = []
    Q_1_100_q95_list = []
    Q_1_100_q90_list = []
    Q_1_100_max_list = []
    U_1_0_100_list = []
    U_1_20_100_list = []
    U_1_50_100_list = []
    U_1_80_100_list = []

    pos_start = []
    pos_end = []

    for w in range(1, num_windows):  # etc: hardcoded here, but above was int(len_genome/50000)
        these_pos = all_pos[allpos_digitized == w]

        if len(these_pos) > 1:
            pos_start.append(min(these_pos))
            pos_end.append(max(these_pos))

            these_pos_idx = np.nonzero(np.in1d(all_pos, these_pos))[0]
            p1_hapw = p1_hap[:, these_pos_idx]
            p2_hapw = p2_hap[:, these_pos_idx]
            p3_hapw = p3_hap[:, these_pos_idx]

            p1_freqw = calc_derived_freq(p1_hapw)
            p2_freqw = calc_derived_freq(p2_hapw)
            p3_freqw = calc_derived_freq(p3_hapw)

        # D-stat
            abbavecw = (1.0 - p2_freqw) * p3_freqw * p1_freqw
            babavecw = p2_freqw * (1.0 - p3_freqw) * p1_freqw
            abbacountsw = np.sum(abbavecw)
            babacountsw = np.sum(babavecw)
            if (abbacountsw + babacountsw > 0):
                Dstatw = (abbacountsw - babacountsw) / (abbacountsw + babacountsw)
            else:
                Dstatw = float('nan')
            Dstat_list.append(Dstatw)

        # fD
            checkfd1 = (p3_freqw > p1_freqw)
            abbafd1 = (1.0 - p2_freqw) * p3_freqw * p3_freqw
            babafd1 = p2_freqw * (1.0 - p3_freqw) * p3_freqw
            checkfd2 = (p3_freqw < p1_freqw)
            abbafd2 = (1.0 - p2_freqw) * p1_freqw * p1_freqw
            babafd2 = p2_freqw * (1.0 - p1_freqw) * p1_freqw
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
            Het = np.sum(hetvec) / 50000  # etc: again, hardcoded window length
            Het_list.append(Het)

            # TODO: Re-implement this
            # divratio = []
            # for archi in range(p1_hapw.shape[0]): #iterate over 0-99 haps; 100 total)
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
            DerFreqs_NonAdm_1 = p3_freqw[ArcHomoDerANDNonAdm_1]
            if DerFreqs_NonAdm_1.size > 0:
                Q_1_100_q95 = np.percentile(DerFreqs_NonAdm_1, 95)
                Q_1_100_q90 = np.percentile(DerFreqs_NonAdm_1, 90)
                Q_1_100_max = np.max(DerFreqs_NonAdm_1)
            else:
                Q_1_100_q95 = float('nan')
                Q_1_100_q90 = float('nan')
                Q_1_100_max = float('nan')

            Q_1_100_q95_list.append(Q_1_100_q95)
            Q_1_100_q90_list.append(Q_1_100_q90)
            Q_1_100_max_list.append(Q_1_100_max)

            U_1_0_100 = (ArcHomoDerANDNonAdm_1 & (p3_freqw > 0))
            U_1_20_100 = (ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.2))
            U_1_50_100 = (ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.5))
            U_1_80_100 = (ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.8))

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

    return pos_start, pos_end, Dstat_list, fD_list, Het_list, divratioavg_list, Q_1_100_q95_list, Q_1_100_q90_list, Q_1_100_max_list, U_1_0_100_list, U_1_20_100_list, U_1_50_100_list, U_1_80_100_list


def update_par_file(temp_par, new_par, model, growth, dominance,
                    nscale, m4s, hs, insert_ai, sex, uniform_recombination,
                    trees_filename, region_filename):

    oldfile = open(temp_par)
    newfile = open(new_par, 'w')
    line_counter = 0
    for line_counter, line in enumerate(oldfile):
        fields = line.split()

        if model == 0:  # etc: only implementing m0 rn
            # Set line numbers to change
            sim_lines = [40, 43, 47, 53, 57, 73, 78, 81, 85, 90, 96, 97]
            if dominance == 2:  # neutral model
                reg_line, rec_line, sex_line = (15, 25, 28)
            elif sex == 'X':  # deleterious, Xchr
                sim_lines = [ln + 23 for ln in sim_lines]
                reg_line, rec_line, sex_line = (20, 31, 50)
            else:  # deleterious, A or None
                sim_lines = [ln + 17 for ln in sim_lines]
                reg_line, rec_line, sex_line = (23, 34, 45)
            # Set content for simulation section (not initialization)
            # TODO: calculate timepoints using adm_gen and end_gen??
            if (dominance == 2) and (sex != 'X'):  # skip burn-in entirely
                time_points = [2, 2,
                               100 / nscale + 2, 100 / nscale + 2,
                               10000 / nscale, 10000 / nscale,
                               20000 / nscale, 20000 / nscale,
                               30000 / nscale]
            else:  # need to burn in, set generations accordingly
                time_points = [100000 / nscale, 100000 / nscale,
                               100 / nscale + 100000 / nscale, 100 / nscale + 100000 / nscale,
                               110000 / nscale, 110000 / nscale,
                               120000 / nscale, 120000 / nscale,
                               130000 / nscale]
            time_points.insert(3, insert_ai)
            time_points.insert(5, insert_ai)
            sim_content = [str(int(i)) for i in time_points]
            sim_content[4] += ':'
            sim_content.extend(['sim.treeSeqOutput("' + trees_filename + '");'])
            assert len(sim_lines) == len(sim_content)
            assert len(sim_lines) == 12
            # Write changes for simulation part
            if line_counter in sim_lines:
                idx = sim_lines.index(line_counter)
                if idx in [3, 5]:  # this is an insert_ai for field 1
                    fields[1] = sim_content[idx]
                else:
                    fields[0] = sim_content[idx]
            # Write changes for initialization part
            elif line_counter == 1:
                fields[1] = str(dominance)  # irrelevant in neutral model
            elif line_counter == 2:
                fields[1] = str(nscale)
            elif line_counter == 3:
                fields[1] = str(m4s)
            elif line_counter == reg_line:  # region info file
                fields[2] = 'readFile("' + region_filename + '");'
            elif uniform_recombination and (line_counter == rec_line):
                fields[1] = '1e-09'
                fields[3] = '); //'
            elif line_counter == sex_line:  # initializeSex
                if sex is None:  # comment out the call
                    fields[0] = '// ' + fields[0]
                else:  # initializeSex as autosome or Xchr ("A" or "X")
                    fields[1] = '"' + str(sex) + '"'

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

        new_line = str()
        for item in fields:
            new_line = new_line + item + " "
        newfile.write(new_line + '\n')

    newfile.close()
    oldfile.close()


def run_slim_variable(n, q, r, dominance, nscale, m4s, model, growth, hs,
                      insert_ai, sex, uniform_recombination):

    # set filenames
    region_name = region_all[r]
    region_info_filename = dir_stem + 'regions/sim_seq_info_' + str(region_name) + '.txt'
    trees_filename = dir_stem + 'output/trees/' + str(region_name) + '_m' + str(model) + '_sex' + str(sex) + '.trees'
    new_par = DIR_par + "par_" + region_name + '_m' + str(model) + '_sex' + str(sex) + '_' + str(dominance) + ".txt"

# etc: why were these vales not used in writing par files????
# TODO: send these to update_par_file
    if model ==1:  # etc: not handling this model yet
        if dominance != 2:
            temp_par = dir_stem + "slim/modelh_neg.txt"
        elif dominance == 2:
            temp_par = dir_stem + "slim/modelh_neu.txt"
        adm_gen = (87400-1)/nscale
        end_gen = 89000/nscale
        # t_end = 1600/nscale -1
        popsize=41080/nscale # recipient population size at the end of simulation
        if growth ==1:
            popsize = 550/nscale
        elif growth ==2:
            popsize = 41080/nscale
        elif growth ==3:
            popsize = 7300/nscale
        elif growth ==4:
            popsize = 41080/nscale
        source_popn = 2  # etc: not sure
        recip_popn = 4  # etc: not sure

    elif model == 0:
        popsize = 1000 / nscale  # extant size of p3 (as split off from p2 in mod0)
        source_popn = 1
        recip_popn = 3
        adm_gen = 120000 / nscale
        end_gen = 130000 / nscale
        if dominance != 2:  # deleterious model
            if sex == 'X':
                temp_par = dir_stem + "slim/ts_Xchr_model0_neg.slim"
            else:
                temp_par = dir_stem + "slim/ts_model0_neg.slim"

        elif dominance == 2:  # neutral model
            temp_par = dir_stem + "slim/ts_model0_neu.slim"
            if sex == 'X':
                temp_par = dir_stem + "slim/ts_Xchr_model0_neu.slim"
            else:
                temp_par = dir_stem + "slim/ts_model0_neu.slim"
                # recap'ing obviates need for 10k of SLiM burn-in in neutral model
                adm_gen = 20000 / nscale
                end_gen = 30000 / nscale
    adm_gens_ago = end_gen - adm_gen
    # segsize = 5000000  # nice that this is here, but hardcoded everywhere else

    update_par_file(temp_par, new_par, model, growth, dominance,
                    nscale, m4s, hs, insert_ai, sex, uniform_recombination,
                    trees_filename + '.orig', region_info_filename)

    slim_stdout = DIR_out + 'OUT_' + region_name + str(sex) + str(m4s) + str(n) + ".txt"

    # Run the SLiM simulation!
    os.system('slim %s > %s' % (new_par, slim_stdout))

    # Recapitate and process TreeSequence from SLiM:
    # overlay neutral mutations if applicable
    # remove Y chr if applicable
    # TODO: variable-ize initial_Ne (size of p1 at beginning of sim)
    trees_from_slim = trees_filename + '.orig'
    ts = tt.process_treeseq(trees_from_slim, region_info_filename,
                            neu_or_neg=dominance, sex=sex, n_scale=nscale,
                            unif_recomb=uniform_recombination, mut_rate=base_mut_rate)

    # Calculate how much source ancestry is present in today's recipient popn
    mean_source_anc, source_anc_fracs, intervals = tt.calc_ancestry_frac(ts, source_popn, recip_popn, adm_gens_ago)
    # Mean ancestry per 50kb window
    anc_by_window, anc_windows = calc_ancestry_window(source_anc_fracs, intervals)

    # Calculate other statistics from genotype matrices
    pos_start, pos_end, Dstat_list, fD_list, Het_list, divratioavg_list, Q_1_100_q95_list, Q_1_100_q90_list, Q_1_100_max_list, U_1_0_100_list, U_1_20_100_list, U_1_50_100_list, U_1_80_100_list = calc_stats(ts, sample_size=popsize)

    q.put([n, insert_ai, growth, mean_source_anc, anc_windows, anc_by_window, pos_start, pos_end, Dstat_list, fD_list, Het_list, divratioavg_list, Q_1_100_q95_list, Q_1_100_q90_list, Q_1_100_max_list, U_1_0_100_list, U_1_20_100_list, U_1_50_100_list, U_1_80_100_list])
    #other parameter info are stored in the output file name

    # os.system('rm '+slim_stdout)
    # os.system('rm '+treepath)
    # os.system('rm '+new_par)


def write_to_file(windowfile_name, q):
    windowfile = open(windowfile_name, 'w')
    while 1:  # etc: terrifying
        q_elem = q.get()
        if q_elem == 'kill':  # break if end of queue
            print('END OF SIMULATIONS')
            break
        [n, insert_ai, growth, mean_source_anc, anc_windows, anc_by_window, pos_start, pos_end, Dstat_list, fD_list, Het_list, divratioavg_list, Q_1_100_q95_list, Q_1_100_q90_list, Q_1_100_max_list, U_1_0_100_list, U_1_20_100_list, U_1_50_100_list, U_1_80_100_list] = q_elem
        for i in range(len(Dstat_list)):
            format_string = "%d\t%d\t%d\t%f\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
            items_to_write = (n, insert_ai, growth, mean_source_anc, anc_windows[i][0], anc_windows[i][1], anc_by_window[i], pos_start[i], pos_end[i], Dstat_list[i], fD_list[i], Het_list[i], divratioavg_list[i], Q_1_100_q95_list[i], Q_1_100_q90_list[i], Q_1_100_max_list[i], U_1_0_100_list[i], U_1_20_100_list[i], U_1_50_100_list[i], U_1_80_100_list[i])
            if i == 0:  # first line for each replicate
                # Check that each item has a formatted location to go into
                assert format_string.count('%') == len(items_to_write)
            # formatting check for ancestry windows... might just want to
            # plot ancestry separately with the finer tree intervals and skip
            # weirdly forcing it to 50k windows different to stats' windows
            if np.isnan(anc_windows[i][0]):
                format_string = "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
            windowfile.write(format_string % items_to_write)
        windowfile.flush()
    windowfile.close()


#################################################################################
if __name__ == '__main__':
    # Set params.  (See parser defaults.)
    sex = 'X'  # etc: sex param takes None, 'A', or 'X'
    whichgene = 15 + 10  # 15 was for project.  X is 25

    dominance = 0  # if 0, run the deleterious recessive model #if 2, run the neutral model
    m4s = 0.01  # adaptive selection strength
    uniform_recombination = 1e-09
    base_mut_rate = 1.5e-8

    nscale = 100  # define scaling factor
    num_reps = 1  # number of simulations per region

    # Set directories
    dir_stem = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/"
    DIR_region = dir_stem + "regions/"
    DIR_out = dir_stem + "output/out/"
    DIR_tree = dir_stem + "output/trees/"
    DIR_par = dir_stem + "slim/"

    # Collect info about region
    r = int(whichgene - 1)
    region_all = ["chr11max", "chr19region", "chr3region", "galnt18", "hla", "hyal2",
                  "krt71", "nlrc5", "oca2", "pde6c", "pou2f3", "rnf34", "sema6d", "sgcb",
                  "sgcz", "sipa1l2", "slc16a11", "slc19a3", "slc5a10", "stat2", "tbx15",
                  "tlr1610", "tnfa1p3", "txn", 'X-0-5M']
    region_name = region_all[r]
    region_info_file = DIR_region + "sim_seq_info_" + str(region_name) + ".txt"
    # Find an exon in the middle-ish of the region...
    window_start, window_end = find_ai_site(region_info_file)
    # ...and put the AI variant in the middle of that exon.
    insert_ai = int((int(window_end) + int(window_start)) / 2)

    # Run simulations and calculate statistics
    attempt_num = np.random.randint(5000)
    print(attempt_num)
    windowfile_name = dir_stem + "output/stats/20200806/" + region_name + "-dominance" + str(dominance) + "-model" + str(model) + "-sex" + str(sex) + "-hs" + str(hs) + "-ai" + str(m4s) + '-attempt' + str(attempt_num) + '_human_windows.txt'
    num_proc = 10
    manager = Manager()
    pool = Pool(processes=num_proc)
    q = manager.Queue()
    watcher = pool.apply_async(write_to_file, (windowfile_name, q))
    reps = range(0, num_reps)
    args_iterable = list(zip(reps, [q] * num_reps))
    for i in args_iterable:
        n = i[0]
        print(str(n))
        run_slim_variable(i[0], i[1], r, dominance, nscale, m4s, model, growth,
                          hs, insert_ai, sex, uniform_recombination)

    q.put('kill')
    pool.close()
    pool.join()

    print("END OF SIMULATION")
