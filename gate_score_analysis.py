# -*- coding: utf-8 -*-
"""
Created on Mon May  9 10:42:05 2022

@author: marsh
"""

import pandas as pd
import numpy as np
import seaborn as sns
import logomaker
import matplotlib.pyplot as plt

data = pd.read_csv('gate_scores.csv')
data0 = data[~(data==0).any(axis=1)]

proteins = np.hstack(np.array(data.columns[1:]))

for i,p in enumerate(proteins):
    off_targets=np.delete(proteins,i)
    
    # data[str(p) +'_score'] = 4*data[p]-data[off_targets].sum(axis=1).sort_values(ascending=False)
    data0[str(p) +'_score'] = 4*data0[p]-data0[off_targets].sum(axis=1).sort_values(ascending=False)
# 

# sns.pairplot(data0[[str(p) +'_score' for p in proteins]])

# for p in proteins:
#     plot_data = data.loc[data[p+'_score'] > 15]
#     logomatrix = logomaker.alignment_to_matrix(sequences=list(plot_data['peptide']),to_type='information',
#                                                 counts=np.hstack(plot_data[p+'_score'].values),characters_to_ignore='_')
#     logoplot = logomaker.Logo(logomatrix,color_scheme='NajafabadiEtAl2017')
#     logoplot.ax.set_ylabel('information (bits)')
#     plt.xticks(np.arange(0,23))
#     plt.title(p)
#     plt.savefig(str(p) + '.svg')
#     plt.close()


for p in proteins:
    #gate score
    plot = data.loc[(data[p]!=0) & (data[proteins] % 1 != 0 )]
    plt.scatter(np.arange(len(plot)),plot[p].sort_values(ascending=False),label=p)
    
    #specificity score
    p = p+'_score'
    plot = data.loc[(data[p]!=0) & (data[proteins] % 1 != 0 ).any(axis=1)]
    plt.scatter(np.arange(len(plot)),plot[p].sort_values(ascending=False),label=p)
    
# plt.xscale('log')
# plt.yscale('log')
plt.legend()
plt.xlabel('peptide')
plt.ylabel('specificity score')