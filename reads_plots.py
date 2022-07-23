# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:28:25 2022

@author: marsh
"""

import pandas as pd
import numpy as np
import seaborn as sns
import logomaker
import matplotlib.pyplot as plt

# =============================================================================
# preprocess data
# =============================================================================
proteins = ['M','F','X','W','2']
gates = [i+1 for i in range(12)]

data = pd.read_csv('merged.csv')
data = data.fillna(0)
#most strict, only keep peptides where they were measured in every case. num_pep=3394
# data0 = data[~(data==0).any(axis=1)]
#less strict, keep where each library needs n/12.
#n=6: num_pep = 25000
#n=3: num_pep = 57000
# n=3
# data0 = data[((~(data[['M'+str(g) for g in gates]]==0)).sum(axis=1) >= n) & ((~(data[['F'+str(g) for g in gates]]==0)).sum(axis=1) >= n) & 
#              ((~(data[['X'+str(g) for g in gates]]==0)).sum(axis=1) >= n) & ((~(data[['W'+str(g) for g in gates]]==0)).sum(axis=1) >= n) &
#              ((~(data[['2'+str(g) for g in gates]]==0)).sum(axis=1) >= n)]

#strictest, based on both number of gates and number of reads
n = 100 #num reads
g = 3 #gates
data0 = data.loc[(data[['M'+str(g) for g in gates]].sum(axis=1) >= n) & (data[['F'+str(g) for g in gates]].sum(axis=1) >= n) &
                        (data[['X'+str(g) for g in gates]].sum(axis=1) >= n) & (data[['W'+str(g) for g in gates]].sum(axis=1) >= n) &
                        (data[['2'+str(g) for g in gates]].sum(axis=1) >= n) & ((~(data[['M'+str(g) for g in gates]]==0)).sum(axis=1) >= g) & ((~(data[['F'+str(g) for g in gates]]==0)).sum(axis=1) >= g) & 
                  ((~(data[['X'+str(g) for g in gates]]==0)).sum(axis=1) >= g) & ((~(data[['W'+str(g) for g in gates]]==0)).sum(axis=1) >= g) &
                  ((~(data[['2'+str(g) for g in gates]]==0)).sum(axis=1) >= g)]

# =============================================================================
# #get strictness info based on number of non-zero gates
# =============================================================================
# num_peps = np.zeros(12)
# for n in range(1,13):
#     num_peps[n-1] = len(data[((~(data[['M'+str(g) for g in gates]]==0)).sum(axis=1) >= n) & ((~(data[['F'+str(g) for g in gates]]==0)).sum(axis=1) >= n) & 
#                  ((~(data[['X'+str(g) for g in gates]]==0)).sum(axis=1) >= n) & ((~(data[['W'+str(g) for g in gates]]==0)).sum(axis=1) >= n) &
#                  ((~(data[['2'+str(g) for g in gates]]==0)).sum(axis=1) >= n)])
    
# plt.plot(range(1,12),num_peps)

# =============================================================================
# get stricness info based on the number of reads in all gates
# =============================================================================
# num_peps = np.zeros(12)
# ns = np.logspace(np.log10(10),np.log10(5000),12)
# for i,n in zip(range(1,13),ns):
#     num_peps[i-1] = len(data.loc[(data[['M'+str(g) for g in gates]].sum(axis=1) >= n) & (data[['F'+str(g) for g in gates]].sum(axis=1) >= n) &
#                         (data[['X'+str(g) for g in gates]].sum(axis=1) >= n) & (data[['W'+str(g) for g in gates]].sum(axis=1) >= n) &
#                         (data[['2'+str(g) for g in gates]].sum(axis=1) >= n)])
# plt.scatter(ns,num_peps)
# plt.ylabel('number of peptides')
# plt.xlabel('minimum read counts per protein')

# =============================================================================
# get strictness info based on # reads and # gates
# =============================================================================
# n_g = [3,6,9,12]
# n_r = np.logspace(np.log10(10),np.log10(5000),12)


# for g in n_g:
#     for i,n in zip(range(1,13),n_r):
#         num_peps[i-1] = len(data.loc[(data[['M'+str(g) for g in gates]].sum(axis=1) >= n) & (data[['F'+str(g) for g in gates]].sum(axis=1) >= n) &
#                         (data[['X'+str(g) for g in gates]].sum(axis=1) >= n) & (data[['W'+str(g) for g in gates]].sum(axis=1) >= n) &
#                         (data[['2'+str(g) for g in gates]].sum(axis=1) >= n) & ((~(data[['M'+str(g) for g in gates]]==0)).sum(axis=1) >= g) & ((~(data[['F'+str(g) for g in gates]]==0)).sum(axis=1) >= g) & 
#                   ((~(data[['X'+str(g) for g in gates]]==0)).sum(axis=1) >= g) & ((~(data[['W'+str(g) for g in gates]]==0)).sum(axis=1) >= g) &
#                   ((~(data[['2'+str(g) for g in gates]]==0)).sum(axis=1) >= g)])
#     plt.scatter(n_r,num_peps,label='min # gates: ' + str(g))
    
# plt.ylabel('number of peptides')
# plt.xlabel('minimum read counts per protein')
# plt.legend()
# plt.xscale('log')
# plt.yscale('log')

# =============================================================================
# #gate score
# =============================================================================
# for i,p in enumerate(proteins):
#     columns = [p+str(g) for g in gates]
#     data0[str(p) +'_gate_score'] = data.loc[:,columns].dot(gates)/data.loc[:,columns].sum(axis=1)

# scores = data0[[str(p)+'_gate_score' for p in proteins]]
# scores.index = data0['peptide']

# =============================================================================
# plot gate score distribution
# =============================================================================
# for p in proteins:
#     #gate score
#     plt.scatter(np.arange(len(scores)),scores[p+'_gate_score'].sort_values(ascending=False),label=p)
    
# plt.legend()
# plt.xlabel('peptide number')
# plt.ylabel('gate score')

# =============================================================================
# logopots of top peptides by gate score
# =============================================================================
# for p in proteins:
#     plot_data = scores.loc[scores[str(p)+'_gate_score'] > 10]
    # logomatrix = logomaker.alignment_to_matrix(sequences=list(plot_data.index),to_type='information',
    #                                             counts=plot_data[p+'_gate_score'].values,characters_to_ignore='_')
    # logoplot = logomaker.Logo(logomatrix,color_scheme='NajafabadiEtAl2017')
#     logoplot.ax.set_ylabel('information (bits)')
#     plt.xticks(np.arange(0,23))
#     plt.title('number of peptides: ' + str(len(plot_data))+', protein = ' + str(p))
#     plt.savefig(str(p) + '.svg')
#     plt.close()


# #specificity score
# for i,p in enumerate(proteins):
#     off_targets=np.delete(proteins,i)
    
#     # data[str(p) +'_score'] = 4*data[p]-data[off_targets].sum(axis=1).sort_values(ascending=False)
#     data0[str(p) +'_score'] = 4*data0[p]-data0[off_targets].sum(axis=1).sort_values(ascending=False)
    

# proteins = np.hstack(np.array(data.columns[1:]))


# for p in proteins:
    # s = [p+str(g) for g in gates]
    
#     plt.scatter(np.arange(len(data)),data[s].sum(axis=1).sort_values(ascending=False)/data[s].sum(axis=1).sum(axis=0),label=p)
    
# plt.xscale('log')
# plt.yscale('log')
# plt.legend()
# plt.xlabel('peptide')
# plt.ylabel('frequency in library')