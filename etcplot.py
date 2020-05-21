#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 21:59:55 2020

@author: egibson
"""

import numpy as np
import matplotlib.pyplot as plt


sexX = np.loadtxt("output/stats/sexX_nscale100/sgcz-dominance0-model0-sexX-hs0-ai0.01_human_windows.txt")
sexXnull = np.loadtxt("output/stats/sexXnull_nscale100/sgcz-dominance0-model0-sexX-hs0-ai0attempt645forItAll_human_windows.txt")
# sexXnull = np.loadtxt("output/stats/sgcz-dominance0-model0-sexX-hs0-ai0_human_windows.txt")
sexAnull = np.loadtxt("output/stats/sexA_oscar/sgcz-dominance0-model0-sexA-hs0-ai0_human_windows.txt")
sexA = np.loadtxt("output/stats/sexA_oscar/sgcz-dominance0-model0-sexA-hs0-ai0.01_human_windows.txt")

colnames = ["n","insert_ai", "growth","meanpI","start","end","freq_before","freq_after","pI","Dstat","fD","Het","divratioavg","Q_1_100_q95","Q_1_100_q90","Q_1_100_max","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100"]
fDind = 10
u80ind = 19

sexX_pos = np.mean(sexX[:,[4,5]], axis=1)
sexXnull_pos = np.mean(sexXnull[:,[4,5]], axis=1)
sexAnull_pos = np.mean(sexAnull[:,[4,5]], axis=1)
sexA_pos = np.mean(sexAnull[:,[4,5]], axis=1)
 #%%
# fig1, ax1 = plt.subplots()
# ax1.scatter(sexX_pos, sexX[:, fDind], c='r', marker='.', alpha='.2')
# ax1.scatter(sexXnull_pos, sexXnull[:, fDind], c='b', marker='.', alpha='.2')

# fig2, ax2 = plt.subplots()
# ax2.hist(sexX[:,fDind], histtype='step', density=True, cumulative=True, label='AI on Xchr')
# ax2.hist(sexXnull[:,fDind], histtype='step', density=True, cumulative=True, label='null on Xchr')

#%%  FDR in windows

fD_thresh_X = np.percentile(sexXnull[:,fDind], 95)
u80_thresh_X = np.percentile(sexXnull[:,u80ind], 95)

fD_thresh_A = np.percentile(sexAnull[:,fDind], 95)
u80_thresh_A = np.percentile(sexAnull[:,u80ind], 95)


def true_pos_vs_neg_null(null_data, ai_data, stat_ind):
    thresh = np.percentile(null_data[:,stat_ind], 95)
    tp = ai_data[:, stat_ind] > thresh
    tpr = sum(tp) / ai_data.shape[0]
    return tp, tpr

tpX_fD, tprX_fD = true_pos_vs_neg_null(sexXnull, sexX, fDind)
tpX_u80, tprX_u80 = true_pos_vs_neg_null(sexXnull, sexX, u80ind)
tpA_fD, tprA_fD = true_pos_vs_neg_null(sexAnull, sexA, fDind)
tpA_u80, tprA_u80 = true_pos_vs_neg_null(sexAnull, sexA, u80ind)


#%%
stats2plot = [fDind, u80ind]

def twopack(null_data, null_pos, ai_data, ai_pos, stat_ind, chr_ax, cdf_ax,
            chr_line=None, cdf_line=None):
    chr_ax.scatter(ai_pos, ai_data[:, stat_ind], c='r', marker='.',
                   label='AI', alpha='.2')
    chr_ax.scatter(null_pos, null_data[:, stat_ind], c='b', marker='.',
                   label='null', alpha='.2')

    cdf_ax.hist(ai_data[:,stat_ind], histtype='step', cumulative=False, label='AI')
    cdf_ax.hist(null_data[:,stat_ind], label='null', cumulative=False, histtype='step')

    if chr_line is not None:
        chr_ax.axhline(chr_line, c='g', label='5% tail of null')
    if cdf_line is not None:
        cdf_ax.axvline(cdf_line, c='g', label='5% tail of null')

    pass

threshX = [fD_thresh_X, u80_thresh_X]
threshA = [fD_thresh_A, u80_thresh_A]

def plot_stats(stat_inds, stat_names, thresholds, null_data, null_pos,
               ai_data, ai_pos):
    fig, axes = plt.subplots(nrows=len(stat_inds), ncols=2, dpi=800)
    linaxes = axes.ravel()
    for i in range(len(stat_inds)):
        twopack(null_data, null_pos, ai_data, ai_pos, stats2plot[i],
                linaxes[i*2], linaxes[i*2+1], thresholds[i], thresholds[i])

    for ax, row in zip(axes[:,0], stat_names):
        ax.set_ylabel(row, rotation=0, size='large')
    fig.tight_layout(pad=2.0)
    return fig, axes

fig3, axes3 = plot_stats(stats2plot, ['fD', 'U80'], threshX, sexXnull, sexXnull_pos,
                         sexX, sexX_pos)
fig3.legend()
fig3.suptitle("X chromosome: AI summary stats")
fig4, axes4 = plot_stats(stats2plot, ['fD', 'U80'], threshA, sexAnull, sexAnull_pos,
                         sexA, sexA_pos)
fig4.suptitle("Autosome: AI summary stats")
fig4.legend()

#%%

window_names = np.arange(99)

def just_tpr(ai_data, stat_ind, stat_thresh):
    this_data = ai_data[:, stat_ind]
    return sum(this_data > stat_thresh) / this_data.shape[0]

def window_stat_tpr(ai_data, ai_pos, stat_ind, stat_thresh):
    windows = np.linspace(0,5000000, 100)
    window_assignments = np.digitize(ai_pos, windows)
    window_tprs = []
    for w in range(1,99):
        window_data = ai_data[window_assignments==w, :]
        window_tprs.append(just_tpr(window_data, stat_ind, stat_thresh))
    return window_tprs

fD_tprs = window_stat_tpr(sexX, sexX_pos, fDind, fD_thresh_X)
fD_fprs = window_stat_tpr(sexXnull, sexXnull_pos, fDind, fD_thresh_X)

def plot_stat_tprs(ax, ai_data, ai_pos, null_data, null_pos, stat_ind, stat_name,
                   stat_thresh, colors=['r', 'b']):
    tprs = window_stat_tpr(ai_data, ai_pos, stat_ind, stat_thresh)
    fprs = window_stat_tpr(null_data, null_pos, stat_ind, stat_thresh)

    ax.plot(tprs, label='TPR', c=colors[0])
    ax.plot(fprs, label='FPR', c=colors[1])
    ax.axhline(np.mean(tprs), c=colors[0], ls='--')
    ax.axhline(np.mean(fprs), c=colors[1], ls='--')

    ax.set_ylabel(stat_name)
    ax.set_xlabel("windowed chr axis")
    return ax

fig5, (ax51, ax52) = plt.subplots(2, dpi=800)
plot_stat_tprs(ax51, sexX, sexX_pos, sexXnull, sexXnull_pos, fDind, "fD", fD_thresh_X)
plot_stat_tprs(ax52, sexX, sexX_pos, sexXnull, sexXnull_pos, u80ind, "U80", u80_thresh_X)
fig5.legend()
fig5.suptitle("X chromosome: detecting AI")

fig6, (ax61, ax62) = plt.subplots(2, dpi=800)
plot_stat_tprs(ax61, sexA, sexA_pos, sexAnull, sexAnull_pos, fDind, "fD", fD_thresh_A)
plot_stat_tprs(ax62, sexA, sexA_pos, sexAnull, sexAnull_pos, u80ind, "U80", u80_thresh_A)
fig6.legend()
fig6.suptitle("Autosome: detecting AI")
