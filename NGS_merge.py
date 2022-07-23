# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 14:49:58 2022

@author: Marshall
"""
import pandas as pd
import numpy as np
import os
from pandas import DataFrame
import sys


#xlsx_files = ['2F2','2F3h','2F4h','2F3n','2F4p','FF1','FF2','FF3h','FF4h',
#              'FF3n','FF4p','MF2','MF3h','MF4h','MF3n','MF4p','MM1','MM2',
#              'WF1','WF2','WF3h','WF4h','WF3n','WF4p','XF1','XF2','XF3h',
#              'XF4h','XF3n','XF4p']

csv_files = ['21','22','23','24','25','26','27','28','29','210','211','212',
	'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12',
	'M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12',
	'W1','W2','W3','W4','W5','W6','W7','W8','W9','W10','W11','W12',
	'X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12']
file_titles = ['5829-MC-1_AGGAGTCC-GCGATCTA_patterns_combined.csv','5829-MC-1_AGGAGTCC-ATAGAGAG_patterns_combined.csv',
	'5829-MC-1_AGGAGTCC-AGAGGATA_patterns_combined.csv','5829-MC-1_AGGAGTCC-TCTACTCT_patterns_combined.csv',
	'5829-MC-1_AGGAGTCC-CTCCTTAC_patterns_combined.csv','5829-MC-1_AGGAGTCC-TATGCAGT_patterns_combined.csv',
	'5829-MC-1_AGGAGTCC-TACTCCTT_patterns_combined.csv','5829-MC-1_AGGAGTCC-AGGCTTAG_patterns_combined.csv',
	'5829-MC-1_AGGAGTCC-ATTAGACG_patterns_combined.csv','5829-MC-1_AGGAGTCC-CGGAGAGA_patterns_combined.csv',
	'5829-MC-1_AGGAGTCC-CTAGTCGA_patterns_combined.csv','5829-MC-1_AGGAGTCC-AGCTAGAA_patterns_combined.csv',
	'5829-MC-1_CTAGTACG-GCGATCTA_patterns_combined.csv','5829-MC-1_CTAGTACG-ATAGAGAG_patterns_combined.csv',
	'5829-MC-1_CTAGTACG-AGAGGATA_patterns_combined.csv','5829-MC-1_CTAGTACG-TCTACTCT_patterns_combined.csv',
	'5829-MC-1_CTAGTACG-CTCCTTAC_patterns_combined.csv','5829-MC-1_CTAGTACG-TATGCAGT_patterns_combined.csv',
	'5829-MC-1_CTAGTACG-TACTCCTT_patterns_combined.csv','5829-MC-1_CTAGTACG-AGGCTTAG_patterns_combined.csv',
	'5829-MC-1_CTAGTACG-ATTAGACG_patterns_combined.csv','5829-MC-1_CTAGTACG-CGGAGAGA_patterns_combined.csv',
	'5829-MC-1_CTAGTACG-CTAGTCGA_patterns_combined.csv','5829-MC-1_CTAGTACG-AGCTAGAA_patterns_combined.csv',
	'5829-MC-1_TCGCCTTA-GCGATCTA_patterns_combined.csv','5829-MC-1_TCGCCTTA-ATAGAGAG_patterns_combined.csv',
	'5829-MC-1_TCGCCTTA-AGAGGATA_patterns_combined.csv','5829-MC-1_TCGCCTTA-TCTACTCT_patterns_combined.csv',
	'5829-MC-1_TCGCCTTA-CTCCTTAC_patterns_combined.csv','5829-MC-1_TCGCCTTA-TATGCAGT_patterns_combined.csv',
	'5829-MC-1_TCGCCTTA-TACTCCTT_patterns_combined.csv','5829-MC-1_TCGCCTTA-AGGCTTAG_patterns_combined.csv',
	'5829-MC-1_TCGCCTTA-ATTAGACG_patterns_combined.csv','5829-MC-1_TCGCCTTA-CGGAGAGA_patterns_combined.csv',
	'5829-MC-1_TCGCCTTA-CTAGTCGA_patterns_combined.csv','5829-MC-1_TCGCCTTA-AGCTAGAA_patterns_combined.csv',
	'5829-MC-1_GCTCAGGA-GCGATCTA_patterns_combined.csv','5829-MC-1_GCTCAGGA-ATAGAGAG_patterns_combined.csv',
	'5829-MC-1_GCTCAGGA-AGAGGATA_patterns_combined.csv','5829-MC-1_GCTCAGGA-TCTACTCT_patterns_combined.csv',
	'5829-MC-1_GCTCAGGA-CTCCTTAC_patterns_combined.csv','5829-MC-1_GCTCAGGA-TATGCAGT_patterns_combined.csv',
	'5829-MC-1_GCTCAGGA-TACTCCTT_patterns_combined.csv','5829-MC-1_GCTCAGGA-AGGCTTAG_patterns_combined.csv',
	'5829-MC-1_GCTCAGGA-ATTAGACG_patterns_combined.csv','5829-MC-1_GCTCAGGA-CGGAGAGA_patterns_combined.csv',
	'5829-MC-1_GCTCAGGA-CTAGTCGA_patterns_combined.csv','5829-MC-1_GCTCAGGA-AGCTAGAA_patterns_combined.csv',
	'5829-MC-1_TTCTGCCT-GCGATCTA_patterns_combined.csv','5829-MC-1_TTCTGCCT-ATAGAGAG_patterns_combined.csv',
	'5829-MC-1_TTCTGCCT-AGAGGATA_patterns_combined.csv','5829-MC-1_TTCTGCCT-TCTACTCT_patterns_combined.csv',
	'5829-MC-1_TTCTGCCT-CTCCTTAC_patterns_combined.csv','5829-MC-1_TTCTGCCT-TATGCAGT_patterns_combined.csv',
	'5829-MC-1_TTCTGCCT-TACTCCTT_patterns_combined.csv','5829-MC-1_TTCTGCCT-AGGCTTAG_patterns_combined.csv',
	'5829-MC-1_TTCTGCCT-ATTAGACG_patterns_combined.csv','5829-MC-1_TTCTGCCT-CGGAGAGA_patterns_combined.csv',
	'5829-MC-1_TTCTGCCT-CTAGTCGA_patterns_combined.csv','5829-MC-1_TTCTGCCT-AGCTAGAA_patterns_combined.csv']

DFdict = {name: pd.DataFrame() for name in csv_files}

#print('read xlsx')
#for col in xlsx_files:
#    print(col)
#    df = pd.read_excel(col+'.xlsx')
#    df = df.set_index('peptide')
#    DFdict[col] = pd.DataFrame(index=df.index)
#    DFdict[col][col] = df['count']
#    del df

print('read csv')
for col,file_title in zip(csv_files,file_titles):
    print(col)
    df = pd.read_csv(file_title)
    df = df.set_index('peptide')
    DFdict[col] = pd.DataFrame(index=df.index)
    DFdict[col][col] = df['count']
    del df

merged = pd.concat(DFdict.values(),axis=1,join='outer')
merged.to_csv('merged.csv')
