#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:17:10 2020

@author: egibson
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

final = pd.read_csv('./regions/sim_seq_info_sgcz.txt', ' ', header=None,
                    names=['type', 'start', 'end_or_rate'])


# # deconvolute recombination rate and end position
# final['rate'] = np.where(final['type'] == 'recRate', final['end_or_rate'], np.NaN)


segments = final[final['type'] != 'recRate']
segments['length'] = segments['end_or_rate'] - segments['start']

segments['length_cs']=segments['length'].cumsum()

segments['length_cs'].plot()

#%%

xchr = pd.read_csv('./regions/mart_export-1.txt', '\t', skiprows=1, header=None, names = ['ts_start', 'ts_end', 'g_start', 'g_end', 'ex_start', 'ex_end'])

xchr['gene'] = list(zip(xchr['g_start'], xchr['g_end']))
xchr['ts'] = list(zip(xchr['ts_start'], xchr['ts_end']))

midx = pd.MultiIndex.from_frame(xchr[['gene', 'ts']])

xchr2 = pd.DataFrame(xchr.iloc[:, :6], index=midx)
