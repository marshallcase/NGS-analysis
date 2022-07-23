# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:02:13 2022

@author: Marshall
"""
import pandas as pd
import numpy as np
import os
from translate_dna import translate_dna
import glob
from pandas import DataFrame
from Bio import SeqIO
import sys

rootdir = os.getcwd()
index=0

#filepath1 = 'merged.fastq'
filepath1 = sys.argv[1]
listOfDataFrames = [filepath1.split('.')[0]]

#create dictionary of dataframes, one for each NGS pool
DFdict = {name: pd.DataFrame() for name in listOfDataFrames}

if not os.path.exists(filepath1):
    print('%s doesn\'t exist' %filepath1)
    sys.exit()
    
#patterns based on DNA library
patterns_array = ['p1','p2','p5','p6','p12','p13','p14','pn']
patterns = ['MR.....M.ELRR..D.F..YYA','GM......MELRR..D.F..YYA','GR..M....ELMR..D.F..YYA','GR...M...ELRM..D.F..YYA',
            'GR.......ELMR..D.FM.YYA','GR.......ELRM..D.F.MYYA','GR.......ELRRM.D.F..MYA','.........EL....D.F...YA']
patterns_re = '|'.join(patterns)
#patterns based on eCPX scaffold
patterns_array_eCPX = ['linker_only','HA_linker','eCPX_3prime','complete_read3']
patterns_eCPX = [            'GGSGGSGGQS.......................','YPYDVPDYAAGGSGGSGGQS.......................',
            '.......................GGQSGQSGDYNKNQYYGITAGPAYRINDWASIYG',
            'GTSVAGQSGQYPYDVPDYAAGGSGGSGGQS.......................GGQSGQSGDYNKNQYYGITAGPAYRINDWASIYG']

patterns_re_eCPX = '|'.join(patterns_eCPX)
patterns_peptide_keys = [p.find('.') for p in patterns_eCPX]
patterns_dict = dict(zip(patterns_eCPX,patterns_peptide_keys))

#get patterns based on eCPX scaffold proteins
for subdir, dirs, files in os.walk(rootdir):
    if subdir.find('__pycache__') == -1:
        print(os.path.join(subdir, filepath1))
        print('forward DNA')
        #use .decode("utf-8") when on Great Lakes cluster and remove it when on local computer
        sequences = [seq_record.seq._data.decode("utf-8") for seq_record in SeqIO.parse(filepath1,"fastq")] #.decode("utf-8")
        print('make dataframe')
        DFdict[listOfDataFrames[index]] = pd.DataFrame(list(zip(sequences)),columns=['DNA'])
        print('translate fwd read 1')
        DFdict[listOfDataFrames[index]]['protein_read1']=[translate_dna(sample[0:]) for sample in DFdict[listOfDataFrames[index]]['DNA']]
        print('translate fwd read 2')
        DFdict[listOfDataFrames[index]]['protein_read2']=[translate_dna(sample[1:]) for sample in DFdict[listOfDataFrames[index]]['DNA']]
        print('translate fwd read 3')
        DFdict[listOfDataFrames[index]]['protein_read3']=[translate_dna(sample[2:]) for sample in DFdict[listOfDataFrames[index]]['DNA']]
        
        #eCPX patterns
        for r in range(1,4):
            DFdict[listOfDataFrames[index]]['sequences_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.findall(patterns_re_eCPX).str[0]
            # DFdict[listOfDataFrames[index]]['sequences_num_matches'+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.findall(patterns_re).str.len()
            print('take substrings if recognition is based on surrounding sequence, not library conservation patterns')
            for p,p_a in zip(patterns_eCPX,patterns_array_eCPX):
                DFdict[listOfDataFrames[index]].loc[:,'peptide_'+p_a+'_read_'+str(r)]=np.nan
                n = patterns_dict[p]
                test = DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['sequences_read_' + str(r)].astype('str').str.match(p).fillna(False)]
                if(len(test) > 0):
                    DFdict[listOfDataFrames[index]].loc[:,'peptide_'+p_a+'_read_'+str(r)]=\
                        test['sequences_read_'+str(r)].str[n:n+23]
         #peptide patterns       
        for r in range(1,4):
            DFdict[listOfDataFrames[index]]['sequences'+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.findall(patterns_re).str[0]
            DFdict[listOfDataFrames[index]]['sequences_num_matches'+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.findall(patterns_re).str.len()
            for p,p_a in zip(patterns,patterns_array):
                print('count the number of matching sequences (should be 1 for all three read frames, total)')
                DFdict[listOfDataFrames[index]][p_a+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.contains(p)
        print('catch rare instances where alignment or read frame issues create multiple peptide-like sequences and drop them for removal of ambiguity')
        DFdict[listOfDataFrames[index]] = DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[0]][['pn_read_1','pn_read_2','pn_read_3']].sum(axis=1) <= 1]
        print('combine peptides across all 3 read frames')
        DFdict[listOfDataFrames[index]]['peptide_patterns']=pd.concat([DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['pn_read_1']]['sequences_read_1'],
                                                              DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['pn_read_2']]['sequences_read_2'],
                                                              DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['pn_read_3']]['sequences_read_3']])
        print('combine peptides across all 3 read frames')
        columns_to_combine = np.append(np.hstack([['peptide_'+p_a+'_read_'+str(r) for p_a in patterns_array_eCPX] for r in range(1,4)]),'peptide_patterns')
        DFdict[listOfDataFrames[index]]['peptide']=DFdict[listOfDataFrames[index]][columns_to_combine].ffill().iloc[:,-1]
    
# #get patterns based on library conservation patterns
# for subdir, dirs, files in os.walk(rootdir):
#     if subdir.find('__pycache__') == -1:
#         print(os.path.join(subdir, filepath1))
#         print('forward DNA')
#         sequences = [seq_record.seq._data for seq_record in SeqIO.parse(filepath1,"fastq")] #.decode("utf-8")
#         print('make dataframe')
#         DFdict[listOfDataFrames[index]] = pd.DataFrame(list(zip(sequences)),columns=['DNA'])
#         print('translate fwd read 1')
#         DFdict[listOfDataFrames[index]]['protein_read1']=[translate_dna(sample[0:]) for sample in DFdict[listOfDataFrames[index]]['DNA']]
#         print('translate fwd read 2')
#         DFdict[listOfDataFrames[index]]['protein_read2']=[translate_dna(sample[1:]) for sample in DFdict[listOfDataFrames[index]]['DNA']]
#         print('translate fwd read 3')
#         DFdict[listOfDataFrames[index]]['protein_read3']=[translate_dna(sample[2:]) for sample in DFdict[listOfDataFrames[index]]['DNA']]
        
#         for r in range(1,4):
#             DFdict[listOfDataFrames[index]]['sequences'+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.findall(patterns_re).str[0]
#             DFdict[listOfDataFrames[index]]['sequences_num_matches'+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.findall(patterns_re).str.len()
#             for p,p_a in zip(patterns,patterns_array):
#                 print('count the number of matching sequences (should be 1 for all three read frames, total)')
#                 DFdict[listOfDataFrames[index]][p_a+'_read_'+str(r)]=DFdict[listOfDataFrames[index]]['protein_read'+str(r)].str.contains(p)
#         print('catch rare instances where alignment or read frame issues create multiple peptide-like sequences and drop them for removal of ambiguity')
#         DFdict[listOfDataFrames[index]] = DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[0]][['pn_read_1','pn_read_2','pn_read_3']].sum(axis=1) <= 1]
#         print('combine peptides across all 3 read frames')
#         DFdict[listOfDataFrames[index]]['peptide']=pd.concat([DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['pn_read_1']]['sequences_read_1'],
#                                                               DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['pn_read_2']]['sequences_read_2'],
#                                                               DFdict[listOfDataFrames[index]].loc[DFdict[listOfDataFrames[index]]['pn_read_3']]['sequences_read_3']])
print('count number of peptides and consolidate entries')
for df_number,df in enumerate(DFdict):
    if not DFdict[listOfDataFrames[df_number]].empty:
        df2 = DFdict[listOfDataFrames[df_number]].groupby('peptide').first()
        df2['count']=DFdict[listOfDataFrames[df_number]]['peptide'].value_counts()
        DFdict[listOfDataFrames[df_number]] = df2
        del df2

print('output data in csv format')
for df_number,df in enumerate(DFdict):
    if not DFdict[listOfDataFrames[df_number]].empty:
        DFdict[listOfDataFrames[df_number]].to_csv(str(listOfDataFrames[df_number]) + '_patterns_combined.csv')
