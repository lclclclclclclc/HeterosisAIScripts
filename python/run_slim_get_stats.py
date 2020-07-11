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

# No earthly idea why this is implemented like this.
# TODO: scrap the file i/o but keep the digitization
def calc_ancestry_window (ancestry_fracs, intervals, num_windows=100, len_genome=None):
    # TODO: actually fix the window size hardcoding
    end_positions = np.asarray([i[1] for i in intervals])
    if len_genome is None:  # just figure it out from intervals
        len_genome = max(end_positions)
    allpos_bin = np.linspace(0,len_genome, num_windows) #windows of every 50kb

    endpos_digitized = np.digitize(end_positions, allpos_bin)
    ancestry = np.asarray(ancestry_fracs)

    anc_window = []
    anc_pos = []
    for w in range(1, num_windows):
        these_pos = end_positions[endpos_digitized==w]
        these_anc = ancestry[endpos_digitized==w]

        if(len(these_pos))>0:
            anc_window.append(np.mean(these_anc))
            anc_pos.append(these_pos)
        else:
            anc_window.append(float('nan'))
            anc_pos.append(these_pos)

    return anc_window


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

## Keep for divratio re-implementation
# def vSumFunc(other_hap, currentArchi,p1_hapw):
#     current_hap = np.array([p1_hapw[currentArchi,]])
#     div = np.zeros(other_hap.shape)
#     ones = np.ones((other_hap.shape[0],1))
#     current_hap = current_hap
#     current_hap_extended = np.dot(ones, current_hap)
#     div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
#     return np.add.reduce(div, 1)


def calc_stats (trees_filename, sample_size, num_windows=100):
    ts = pyslim.load(trees_filename)
    (p1_hap, p2_hap, p3_hap), all_pos = tt.sample_population_haplotypes(ts, n_haps=int(sample_size), check_loc=insert_ai)

    len_genome = ts.sequence_length
    allpos_bin = np.linspace(0,len_genome,num_windows) #windows of every 50kb
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


    for w in range(1,num_windows):  # etc: hardcoded here, but above was int(len_genome/50000)
        these_pos = all_pos[allpos_digitized==w]

        if len(these_pos)>1:
            pos_start.append(min(these_pos))
            pos_end.append(max(these_pos))

            these_pos_idx = np.nonzero(np.in1d(all_pos,these_pos))[0]
            p1_hapw = p1_hap[:,these_pos_idx]
            p2_hapw = p2_hap[:,these_pos_idx]
            p3_hapw = p3_hap[:,these_pos_idx]

            p1_freqw = calc_derived_freq (p1_hapw)  # why is this called "derived"? It's just a w/in pop'n allele freq
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

            # TODO: Re-implement this
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


    return pos_start, pos_end, Dstat_list, fD_list, Het_list, divratioavg_list, Q_1_100_q95_list, Q_1_100_q90_list, Q_1_100_max_list, U_1_0_100_list, U_1_20_100_list, U_1_50_100_list, U_1_80_100_list


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
                elif line_counter == 85:  # admixture generation late()
                    fields[0] = str(int(20000/nscale)) # TODO: replace by adm_gen
                elif line_counter == 91:  # final generation
                    fields[0] = str(int(30000/nscale))  # TODO: replace by end_gen
                elif line_counter == 92:  # write out .trees
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
                elif line_counter == 103:  # admixture generation late()
                    fields[0] = str(int(120000/nscale))  # TODO: replace by adm_gen
                elif line_counter == 109:  # final generation
                    fields[0] = str(int(130000/nscale))  # TODO: replace by end_gen
                elif line_counter == 110:  # write out .trees
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
        if dominance !=2:
            temp_par = dir_stem + "slim/ts_model0_neg.slim"
            adm_gen = 120000/nscale
            end_gen = 130000/nscale
        elif dominance == 2:
            temp_par = dir_stem + "slim/ts_model0_neu.slim"
            # recap'ing obviates need for 10k of SLiM burn-in in neutral model
            adm_gen = 20000/nscale
            end_gen = 30000/nscale
        popsize = 1000 / nscale  # extant size of p3 (as split off from p2 in mod0)
        source_popn = 1
        recip_popn = 3

    # segsize = 5000000  # nice that this is here, but hardcoded everywhere else
    adm_gens_ago = end_gen - adm_gen

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
    ts = tt.throw_neutral_muts(trees_filename+'.orig', region_info_filename,
                           neu_or_neg=dominance, n_scale=nscale)

    # Calculate how much source ancestry is present in today's recipient popn
    source_anc_fracs, intervals = tt.calc_ancestry_frac_over_region(ts, source_popn, recip_popn, adm_gens_ago)
    mean_source_anc = np.mean(source_anc_fracs)
        # Mean ancestry per 50kb window
    source_anc_by_window = calc_ancestry_window(source_anc_fracs, intervals)

    # Calculate other statistics from genotype matrices
    pos_start,pos_end,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list = calc_stats(trees_filename, sample_size=popsize)

    q.put([n,insert_ai,growth,mean_source_anc,pos_start,pos_end,source_anc_by_window, Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list])
    #other parameter info are stored in the output file name

    # os.system('rm '+slim_output)
    # os.system('rm '+treepath)
    # os.system('rm '+new_par)


def write_to_file(windowfile_name, q):
    windowfile = open(windowfile_name, 'w')
    while 1:  # etc: terrifying
        q_elem = q.get()
        if q_elem=='kill': # break if end of queue
            print ('END OF SIMULATIONS')
            break
        [n,insert_ai,growth,meanp1,pos_start,pos_end,anc_window,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list] = q_elem
        for i in range(len(Dstat_list)):
            windowfile.write("%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (n,insert_ai,growth,meanp1,pos_start[i],pos_end[i],anc_window[i],Dstat_list[i], fD_list[i], Het_list[i], divratioavg_list[i],Q_1_100_q95_list[i],Q_1_100_q90_list[i],Q_1_100_max_list[i],U_1_0_100_list[i],U_1_20_100_list[i],U_1_50_100_list[i],U_1_80_100_list[i]))
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
    nscale = 100 #define scaling factor
    m4s = 0.01 #adaptive selection strength
    num_reps=5 #number of simulations per region
    region_all = ["chr11max","chr19region","chr3region","galnt18","hla","hyal2",
                  "krt71","nlrc5","oca2","pde6c","pou2f3","rnf34","sema6d","sgcb",
                  "sgcz","sipa1l2","slc16a11","slc19a3","slc5a10","stat2","tbx15",
                  "tlr1610","tnfa1p3","txn", 'X-0-5M']

    dir_stem = "/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/"

    DIR_region = dir_stem + "regions/"
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
    print(attempt_num)
    windowfile_name = dir_stem + "output/stats/20200709/"+region_name+"-dominance"+str(dominance)+"-model"+str(model)+"-sex"+str(sex)+"-hs"+str(hs)+"-ai"+str(m4s)+'-attempt' + str(attempt_num) + '_human_windows.txt'
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

    # # peek at first line of stats file to check in on things
    # with open(windowfile_name, 'r') as f:
    #     print("First line of stats file for this run is below.")
    #     print(f.readline())

    print("END OF SIMULATION")
