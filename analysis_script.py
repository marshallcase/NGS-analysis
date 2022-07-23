# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 09:24:55 2022

@author: Marshall
"""


import pandas as pd
import numpy as np
import os
from pandas import DataFrame
import sys
import matplotlib.pyplot as plt

print('load data to memory')
data = pd.read_csv('merged_patterns_combined.csv')
print('change index')
data = data.set_index(['peptide'])
print('normalize to read count')
fractions = data.fillna(0) / data.sum(axis=0)


# =============================================================================
# plot enrichment trajectories of single peptides
# =============================================================================

def plotEnrichmentTrajectories(sequence,fractions):
    plt.figure()
    x_axis_labels = ['unsorted','HA','MACS 1','MACS 2','FACS 1','FACS 2',
                      'FACS 3h','FACS 4h']
    for protein in ['M','F','X','W','2']:
        df_labels = ['M1','M2','F1','F2','F3h','F4h']
        df_labels = [protein+s for s in df_labels]
        df_labels.insert(0,'HA')
        df_labels.insert(0,'Unsorted')
        plt.plot(x_axis_labels,fractions.loc[sequence,df_labels],
                  label=protein)
    
    plt.xlabel('Sorting Round')
    plt.ylabel('Fraction')
    plt.yscale('log')
    plt.legend()
    plt.title(sequence)
    
# =============================================================================
# make logos
# =============================================================================
# import logomaker
    
# #individual    
# for protein in ['M','F','X','W','2']:
#     df_labels = ['HA'] #choose individual round (skip for loop)
#     test = data[df_labels].dropna()
#     logomatrix = logomaker.alignment_to_matrix(sequences=list(test.index),to_type='information',
#                                                 counts=np.hstack(test[df_labels].values),characters_to_ignore='_')
#     logoplot = logomaker.Logo(logomatrix,color_scheme='NajafabadiEtAl2017')
#     logoplot.ax.set_ylabel('information (bits)')
#     plt.xticks(np.arange(0,23))
#     plt.title(df_labels)
#     plt.savefig(str(df_labels[0]) + '_eCPX.svg')


# protein = ['M','F','X','W','2']
# df_labels = 'F4h' #choose round of FACS or MACS
# df_labels = [s+df_labels for s in protein]

# # df_labels = ['MF4h']
# # df_labels.insert(0,'HA')
# # df_labels.insert(0,'Unsorted')
# x_ticks = np.arange(0,23)
# x_ticks_labels = ['1e','1f','1g','2a','2b','2c','2d','2e','2f','2g','3a','3b','3c','3d','3e',
#            '3f','3g','4a','4b','4c','4d','4e','4f']
    
# fig,axes = plt.subplots(len(df_labels),1,sharey=True,sharex=True)
# for l,ax in zip(df_labels,axes):
#     test = data[l].dropna()
#     logomatrix = logomaker.alignment_to_matrix(sequences=list(test.index),to_type='information',
#                                                 counts=np.hstack(test.values),characters_to_ignore='_')
#     logoplot = logomaker.Logo(logomatrix,color_scheme='NajafabadiEtAl2017',ax=ax)
#     ax.set_xticks(x_ticks)
#     ax.set_xticklabels(x_ticks_labels)
#     ax.set_ylabel(l)
    
# fig.text(0.5, 0.04, 'peptide position', ha='center')
# fig.text(0.04, 0.5, 'information (bits)', va='center', rotation='vertical')
# # plt.savefig('unsorted vs FACS 4h.svg')
        
        
# =============================================================================
# plot frequency distributions of libraries
# =============================================================================
print('plotting')
# #individual
# for protein in ['M','F','X','W','2']:
#     print(protein)
#     plt.figure()
#     df_labels = ['M1','M2','F1','F2','F3h','F4h']
#     df_labels = [protein+s for s in df_labels]
#     df_labels.insert(0,'HA')
#     df_labels.insert(0,'Unsorted')
#     for l in df_labels:
#         test = fractions[l].dropna().sort_values(ascending=False)
#         plt.plot(np.arange(1,len(test)+1),test,label=l)
    
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.legend(loc='upper right')
#     plt.xlabel('peptides ranked by read count')
#     plt.ylabel('fraction of peptide in library')
#     plt.title(protein)
#     plt.tight_layout()
#     plt.savefig(protein + '_library_enrichment_trajectories_eCPX.png',dpi=1200)
#     plt.close()
    
#subplots

fig,axes = plt.subplots(3,2,sharey=True,sharex=True)
for protein,ax in zip(['M','F','X','W','2','-'],axes.ravel()):
    if protein == '-':
        1
    else:
        print(protein)
        df_labels = ['M1','M2','F1','F2','F3h','F4h']
        df_labels = [protein+s for s in df_labels]
        df_labels.insert(0,'HA')
        df_labels.insert(0,'Unsorted')
        for l in df_labels:
            test = fractions[l].dropna().sort_values(ascending=False)
            ax.plot(np.arange(1,len(test)+1),test,label=l)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc='upper right')
        ax.set_title(protein)

fig.text(0.5, 0.04, 'peptides ranked by read count', ha='center')
fig.text(0.04, 0.5, 'fraction of peptide in library', va='center', rotation='vertical')
plt.savefig('enrichment trajectories.png',dpi=1200)
# =============================================================================
# analyze staple locations
# =============================================================================
# print('plotting')
# patterns_array = ['p1','p2','p5','p6','p12','p13','p14','pn']
# patterns = ['MR.....M.ELRR..D.F..YYA','GM......MELRR..D.F..YYA','GR..M....ELMR..D.F..YYA','GR...M...ELRM..D.F..YYA',
#             'GR.......ELMR..D.FM.YYA','GR.......ELRM..D.F.MYYA','GR.......ELRRM.D.F..MYA','.........EL....D.F...YA']
# staple_location = pd.DataFrame(index=patterns_array,columns=['Unsorted','HA','2M1','2M2','2F1','2F2',
#                                          '2F3h','2F4h','2F3n','2F4p','FM1','FM2',
#                                          'FF1','FF2','FF3h','FF4h','FF3n','FF4p',
#                                          'MM1','MM2','MF1','MF2','MF3h','MF4h',
#                                          'MF3n','MF4p','WM1','WM2','WF1','WF2',
#                                          'WF3h','WF4h','WF3n','WF4p','XM1','XM2',
#                                          'XF1','XF2','XF3h','XF4h','XF3n','XF4p'])
# for protein in ['M','F','X','W','2']:
#     print(protein)
#     # df_labels = ['M1','M2','F1','F2','F3h','F4h']
#     df_labels=['F3n','F4p']
#     df_labels = [protein+s for s in df_labels]
#     # df_labels.insert(0,'HA')
#     # df_labels.insert(0,'Unsorted')
#     for l in df_labels:
#         test = fractions[l].dropna()
#         test = test.loc[test!=0]
#         test = test.sort_values(ascending=False)
        
#         for p,p_a in zip(patterns,patterns_array):
#             staple_location.loc[p_a,l]=test.loc[test.index.str.contains(p)].sum()


# x_axis_labels = ['unsorted','HA','MACS 1','MACS 2','FACS 1','FACS 2',
#                   'FACS 3h','FACS 4h']
# for protein in ['M','F','X','W','2']:
#     df_labels = ['M1','M2','F1','F2','F3h','F4h']
#     df_labels = [protein+s for s in df_labels]
#     df_labels.insert(0,'HA')
#     df_labels.insert(0,'Unsorted')
#     plt.figure()
#     for p_a in patterns_array:
#         plt.plot(x_axis_labels,staple_location.loc[p_a,df_labels],
#                   label=p_a)

#     plt.xlabel('Sorting Round')
#     plt.ylabel('Fraction')
#     plt.yscale('log')
#     plt.legend()
#     plt.title(protein)
#     plt.savefig(protein+'_staple_location_enrichment.png',dpi=1200)

# =============================================================================
# make logos split by staple location
# =============================================================================
# import logomaker
# patterns_array = ['p1','p2','p5','p6','p12','p13','p14','pn']
# patterns = ['MR.....M.ELRR..D.F..YYA','GM......MELRR..D.F..YYA','GR..M....ELMR..D.F..YYA','GR...M...ELRM..D.F..YYA',
#             'GR.......ELMR..D.FM.YYA','GR.......ELRM..D.F.MYYA','GR.......ELRRM.D.F..MYA','.........EL....D.F...YA']
# staple_location = pd.DataFrame(index=patterns_array,columns=['M','F','X','W','2'])
# for protein in ['M','F','X','W','2']:
#     # df_labels = ['M1','M2','F1','F2','F3h','F4h']
#     df_labels = ['F4h'] #choose round of FACS or MACS
#     df_labels = [protein+s for s in df_labels]
#     # df_labels.insert(0,'HA')
#     # df_labels.insert(0,'Unsorted')
#     # df_labels = ['HA'] #choose individual round (skip for loop)
#     for p,p_a in zip(patterns,patterns_array):
#         test = data[df_labels].dropna()
#         test = test.loc[test.index.str.contains(p)]
#         if len(test.values > 0):
#             # logomatrix = logomaker.alignment_to_matrix(sequences=list(test.index),to_type='information',
#             #                                             counts=np.hstack(test.values),characters_to_ignore='_')
#             # logoplot = logomaker.Logo(logomatrix,color_scheme='NajafabadiEtAl2017')
#             # logoplot.ax.set_ylabel('information (bits)')
#             # plt.xticks(np.arange(0,23))
#             # plt.title(df_labels)
#             # plt.savefig(str(df_labels[0]) + '_' + p_a + '.svg')
#             # plt.close()
#             staple_location.loc[p_a,protein]=test.sum().values[0]
