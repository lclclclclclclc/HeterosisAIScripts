import numpy as np
import matplotlib.pyplot as plt
import pyslim

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:06:40 2020

@author: egibson
"""

#%% TODO: import treetools
import warnings


def translate_from_slim_pop(tree_seq, s_ids=None):
    popkey = [(p.id, p.metadata.slim_id) for p in tree_seq.populations()
              if p.metadata is not None]
    pops = [p for (p, s) in popkey]
    slims = [s for (p, s) in popkey]
    if s_ids is None:  # just get all populations with a slim_id
        s_ids = slims
    indices_of_slim_pops = np.where(np.isin(slims, s_ids))[0]
    return tuple(np.asarray(pops)[indices_of_slim_pops])


# ts_popn_ids = translate_from_slim_pop(ts, s_ids=popn_ids)

# Sample down to n_haps extant nodes per population.
def sample_haplotypes(ts, population, n_haps, time_ago=0, sex_ratio_fm=False):
    all_samps = np.asarray([ts.node(n).id for n in ts.samples(population=population)
                            if (ts.node(n).time == time_ago)])
    if sex_ratio_fm:
        num_female_haps = round(n_haps * (sex_ratio_fm[0] / sum(sex_ratio_fm)))
        num_male_haps = n_haps - num_female_haps

        # Separate extant males and females.
        def get_sex_of_nodes(node_ids):
            individuals = [ts.node(n).individual for n in node_ids]
            female_bool = [ts.individual(i).metadata.sex == pyslim.INDIVIDUAL_TYPE_FEMALE
                           for i in individuals]
            return np.asarray(female_bool)

        sex_bool = get_sex_of_nodes(all_samps)
        females = all_samps[sex_bool]
        males = all_samps[~sex_bool]
        # Sample randomly from each sex.
        if num_female_haps > len(females):
            warnings.warn("Trying to sample more extant female nodes than exist.")
            print(f"Sampling {len(females)} female nodes, rather than requested {num_female_haps}.")
            female_samps = females
        else:
            female_samps = np.random.choice(females, size=num_female_haps, replace=False)
        if num_male_haps > len(males):
            warnings.warn("Trying to sample more extant male nodes than exist.")
            print(f"Sampling {len(males)} male nodes, rather than requested {num_male_haps}.")
            male_samps = males
        else:
            male_samps = np.random.choice(males, size=num_male_haps, replace=False)
        sample_of_samps = np.append(female_samps, male_samps)
    else:
        sample_of_samps = np.random.choice(all_samps, size=n_haps, replace=False)
    return sample_of_samps

#%%  Get the SFS from the treeseq
# X_neg_statfile = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance0-model0-sexX-hs0-ai0-attempt555_human_windows.txt'
# X_neu_statfile = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance2-model0-sexX-hs0-ai0-attempt2280_human_windows.txt'
# A_neg_statfile = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance0-model0-sexA-hs0-ai0-attempt1976_human_windows.txt'
# A_neu_statfile = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance2-model0-sexA-hs0-ai0-attempt2272_human_windows.txt'

# X_neg_treepref = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance0-model0-sexX-hs0-ai0-attempt555-rep'
# X_neu_treepref = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance2-model0-sexX-hs0-ai0-attempt2280-rep'

# X_neg_treefiles = [X_neg_treepref + str(i) + '.trees' for i in range(10)]
# X_neu_treefiles = [X_neu_treepref + str(i) + '.trees' for i in range(10)]


# for i in range(10):
#     ts = pyslim.load(X_neg_treefiles[i])
#     afs = ts.allele_frequency_spectrum()
#     plt.figure()
#     plt.plot(afs, 'r')

# for i in range(10):
#     ts = pyslim.load(X_neu_treefiles[i])
#     afs = ts.allele_frequency_spectrum()
#     plt.figure()
#     plt.plot(afs)


# big_idx = np.where(afs > 0.001)


#%%  Show no sex difference between Autosome SFSs, yes diff. btwn Xchr SFSs

A_AI_neu_treepref = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance2-model0-sexA-hs0-ai0.01-attempt3801-rep'

A_AI_neu_treefiles = [A_AI_neu_treepref + str(i) + '.trees' for i in range(10)]

tsA = pyslim.load(A_AI_neu_treefiles[0])
p1, p2 = translate_from_slim_pop(tsA, s_ids=(1, 2))

sample_times = np.unique([tsA.node(s).time for s in tsA.samples()])
# 0 is today
# 1000 is after introgression
# 1001 is before introgression
# 2000 is at p3 split-off from p2
#  -> AI emerges at 2988
# 2998 is p1/p2 split.  (immediate in this case)
# 3000 is start of SLiM
# We can use time 0 because looking at equilibrium of p1 and p2, not p3 craziness

# Autosome
extant_female_A_noAI = sample_haplotypes(tsA, p1, 500, time_ago=0,
                                                sex_ratio_fm=(1,0))
extant_male_A_noAI = sample_haplotypes(tsA, p1, 500, time_ago=0,
                                                sex_ratio_fm=(0,1))

extant_female_A_AI = sample_haplotypes(tsA, p2, 500, time_ago=0,
                                                sex_ratio_fm=(1,0))
extant_male_A_AI = sample_haplotypes(tsA, p2, 500, time_ago=0,
                                                sex_ratio_fm=(0,1))

afs_female_A_noAI = tsA.allele_frequency_spectrum(sample_sets=[extant_female_A_noAI])
afs_male_A_noAI = tsA.allele_frequency_spectrum(sample_sets=[extant_male_A_noAI])

afs_female_A_AI = tsA.allele_frequency_spectrum(sample_sets=[extant_female_A_AI])
afs_male_A_AI = tsA.allele_frequency_spectrum(sample_sets=[extant_male_A_AI])

past_A_noAI = sample_haplotypes(tsA, p1, 500, time_ago=2000)
past_A_AI = sample_haplotypes(tsA, p2, 500, time_ago=2000)

past_afs_Ano = tsA.allele_frequency_spectrum(sample_sets=[past_A_noAI])
past_afs_A = tsA.allele_frequency_spectrum(sample_sets=[past_A_AI])

# Xchr
X_AI_neu_treepref = '/Users/egibson/Documents/science/Grad/demog20/proj/HeterosisAIScripts/results/Xchr_check/X-0-5M-dominance2-model0-sexX-hs0-ai0.01-attempt2004-rep'

X_AI_neu_treefiles = [X_AI_neu_treepref + str(i) + '.trees' for i in range(10)]

tsX = pyslim.load(X_AI_neu_treefiles[0])
p1, p2 = translate_from_slim_pop(tsX, s_ids=(1,2))


extant_female_X_noAI = sample_haplotypes(tsX, p1, 500, time_ago=0,
                                                sex_ratio_fm=(1,0))
extant_male_X_noAI = sample_haplotypes(tsX, p1, 500, time_ago=0,
                                                sex_ratio_fm=(0,1))

extant_female_X_AI = sample_haplotypes(tsX, p2, 500, time_ago=0,
                                                sex_ratio_fm=(1,0))
extant_male_X_AI = sample_haplotypes(tsX, p2, 500, time_ago=0,
                                                sex_ratio_fm=(0,1))

afs_female_X_noAI = tsX.allele_frequency_spectrum(sample_sets=[extant_female_X_noAI])
afs_male_X_noAI = tsX.allele_frequency_spectrum(sample_sets=[extant_male_X_noAI])

afs_female_X_AI = tsX.allele_frequency_spectrum(sample_sets=[extant_female_X_AI])
afs_male_X_AI = tsX.allele_frequency_spectrum(sample_sets=[extant_male_X_AI])

# not extant: right after sweep

past_male_X_noAI = sample_haplotypes(tsX, p1, 500, time_ago=2000,
                                                sex_ratio_fm=(0,1))
past_female_X_noAI = sample_haplotypes(tsX, p1, 500, time_ago=2000,
                                                sex_ratio_fm=(1, 0))
past_female_X_AI = sample_haplotypes(tsX, p2, 500, time_ago=2000,
                                                sex_ratio_fm=(1,0))
past_male_X_AI = sample_haplotypes(tsX, p2, 500, time_ago=2000,
                                                sex_ratio_fm=(0,1))

past_afs_mXno = tsX.allele_frequency_spectrum(sample_sets=[past_male_X_noAI])
past_afs_mX = tsX.allele_frequency_spectrum(sample_sets=[past_male_X_AI])
past_afs_fXno = tsX.allele_frequency_spectrum(sample_sets=[past_female_X_noAI])
past_afs_fX = tsX.allele_frequency_spectrum(sample_sets=[past_female_X_AI])

#%% AFS comparisons

# I am expecting no difference between sexes (red and blue) in autosome, and
# the AI case (poorly named for one positive var.) shifted left.
# But seeing as the positive var fixes quickly... maybe no shift.
# looks like a bump instead
# Am expecting some sex difference in Xchr, but hard to eyeball

# afs
fig1, axes = plt.subplots(nrows=2, ncols=2, sharey=True)

axes[0, 0].set_title('Autosome, no AI var.')
axes[0,0].plot(afs_male_A_noAI[1:150], 'b.', alpha=0.5, label='male, s0')
axes[0,0].plot(afs_female_A_noAI[1:150], 'r.', alpha=0.5, label='female, s0')

axes[1, 0].set_title('Autosome, AI var.')
axes[1,0].plot(afs_male_A_AI[1:150], 'g.', alpha=0.5, label='male, s0.01')
axes[1,0].plot(afs_female_A_AI[1:150], 'y.', alpha=0.5, label='female, s0.01')

axes[0, 1].set_title('Xchr, no AI var.')
axes[0,1].plot(afs_male_X_noAI[1:150], 'b.', alpha=0.5)
axes[0,1].plot(afs_female_X_noAI[1:150], 'r.', alpha=0.5)

axes[1, 1].set_title('Xchr, AI var.')
axes[1,1].plot(afs_male_X_AI[1:150], 'g.', alpha=0.5)
axes[1,1].plot(afs_female_X_AI[1:150], 'y.', alpha=0.5)

fig1.suptitle('Allele frequency spectra, 500 extant individuals')
fig1.legend()


# is it skewed?
fig2, axes = plt.subplots(nrows=2, ncols=3, sharey=True)

axes[0, 0].set_title('extant autosomes')
axes[0, 0].plot(afs_male_A_noAI[1:60], 'k.', alpha=0.5)
axes[0, 0].plot(afs_female_A_noAI[1:60], 'k.', alpha=0.5)
axes[0, 0].plot(afs_male_A_AI[1:60], 'r.', alpha=0.5)
axes[0, 0].plot(afs_female_A_AI[1:60], 'r.', alpha=0.5)

axes[0, 1].set_title('extant female Xchrs')
axes[0, 1].plot(afs_female_X_noAI[1:60], 'k.', alpha=0.5)
axes[0, 1].plot(afs_female_X_AI[1:60], 'r.', alpha=0.5)

axes[0, 2].set_title('extant male Xchrs')
axes[0, 2].plot(afs_male_X_noAI[1:60], 'k.', alpha=0.5)
axes[0, 2].plot(afs_male_X_AI[1:60], 'r.', alpha=0.5)


axes[1, 0].set_title('nearer sweep autosomes')
axes[1, 0].plot(past_afs_Ano[1:60], 'k.', alpha=0.5)
axes[1, 0].plot(past_afs_A[1:60], 'r.', alpha=0.5)

axes[1, 1].set_title('nearer sweep female Xchrs')
axes[1, 1].plot(past_afs_fXno[1:60], 'k.', alpha=0.5)
axes[1, 1].plot(past_afs_fX[1:60], 'r.', alpha=0.5)

axes[1, 2].set_title('nearer sweep male Xchrs')
axes[1, 2].plot(past_afs_mXno[1:60], 'k.', alpha=0.5, label='s0')
axes[1, 2].plot(past_afs_mX[1:60], 'r.', alpha=0.5, label='s0.01')

fig2.suptitle('Is there a skew due to the selective sweep?')
fig2.legend()


# bumps showing fixed pos. var.
fig3, (ax31, ax32) = plt.subplots(nrows=1, ncols=2, sharey=True)

ax31.plot(afs_male_A_noAI[60:100], 'b.', alpha=0.5, label='male, s0')
ax31.plot(afs_female_A_noAI[60:100], 'r.', alpha=0.5, label='female, s0')
ax31.plot(afs_male_A_AI[60:100], 'g.', alpha=0.5, label='male, s0.01')
ax31.plot(afs_female_A_AI[60:100], 'y.', alpha=0.5, label='female, s0.01')


ax32.plot(afs_female_X_noAI[60:100], 'r.', alpha=0.5)
ax32.plot(afs_male_X_noAI[60:100], 'b.', alpha=0.5)
ax32.plot(afs_female_X_AI[60:100], 'g.', alpha=0.5)
ax32.plot(afs_male_X_AI[60:100], 'y.', alpha=0.5)

fig3.suptitle('This bump shows up in sweep model, not in neutral')
fig3.legend()