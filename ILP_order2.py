#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 16:33:59 2022

@author: marshallcase
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:49:19 2022

@author: marsh
"""

#integer linear programming optimization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.multioutput import MultiOutputRegressor
from pulp import *
from plot_ML import plotKDE
from sklearn.preprocessing import PolynomialFeatures
from sklearn.svm import LinearSVR
import time


position_vector = ['1e','1f','1g','2a','2b','2c','2d','2e','2f','2g','3a','3b',
                   '3c','3d','3e','3f','3g','4a','4b','4c','4d','4e','4f']

#import data
scores_vis = pd.read_excel('scores.xlsx')
scores_vis.set_index('peptide',inplace=True)
scores_vis = scores_vis.drop(columns=scores_vis.columns[5:])
#create feature matrix
features = pd.DataFrame(index=scores_vis.index)
for i in range(23):
    features[i]=features.index.str[i]
    
features.columns = position_vector
#one hot encode feature matrix
enc = OneHotEncoder(sparse=False)
enc.fit(features)
X = enc.transform(features)

#define label vector
y = scores_vis.iloc[:,:] #scores_vis.iloc[:,[5,7,9,11,13]]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

poly = PolynomialFeatures(2,interaction_only=True)
X_2 = poly.fit_transform(pd.DataFrame(X,columns=enc.get_feature_names(input_features=position_vector)))
X_train_2, X_test_2, y_train, y_test = train_test_split(X_2, y, test_size=0.3, random_state=42)

features_titles = poly.get_feature_names(input_features=enc.get_feature_names(input_features=position_vector))
features_titles = [features_title.replace(' ','_') for features_title in features_titles]
# #debug: drop zero columns
# analysis = pd.DataFrame(X_2,columns=features_titles)
# analysis = analysis.loc[:,~(analysis.sum(axis=0)==0)]

# # =============================================================================
# # train Mcl-1 score predictor  linear regression
# # =============================================================================
svr_mcl1 = MultiOutputRegressor(estimator=LinearSVR(C=0.005001177728215577,
                                          epsilon=0.009344323060800772,
                                          tol=0.001))
svr_mcl1.fit(X_train_2,y_train)
# plotKDE(linreg_multi.predict(X_test),y_test)
coefs = [model.coef_ for model in svr_mcl1.estimators_]



# =============================================================================
# #define ILP constraints
# =============================================================================
position_constraints = enc.categories_
allVariables = pd.DataFrame(index=features_titles,columns=['Variables'])

#define all first order input variables (limited by sequence space of identified peptides)
for v in features_titles:
    allVariables.at[v,'Variables'] = LpVariable(v,0,1,cat='Binary')
    

#define model, set to minimize
peptide_model = LpProblem("Peptide minimize", LpMinimize)

#add constraints (each position needs one and exactly one amino acid)
for p_v in position_vector:
    # peptide_model += (lpSum(allVariables.at[p_v,'Variables'])==1,p_v)
    peptide_model += (lpSum(np.hstack(allVariables.loc[(allVariables.index.str[:2]==p_v) & (allVariables.index.str.len() < 5)].values))==1,p_v)
    
#add constraints that each combo of positions can only add to 1
# i.e. 1e_A,1e_C and 1f_A,1f_C can only have a single 1 value
for i,p_v in enumerate(position_vector):
    for j,p_v2 in enumerate(position_vector[i+1:]):
        peptide_model += (lpSum(np.hstack(allVariables.loc[(allVariables.index.str[:2]==p_v) & (allVariables.index.str[5:7]==p_v2)].values))==1,p_v+'_'+p_v2)

#add constraint that the peptide must be stapled
peptide_model += (lpSum(np.hstack(allVariables.loc[allVariables.index.str.contains('M') & (allVariables.index.str.len() < 5)].values))==2,'stapled')
#equivalent
#np.hstack(allVariables.loc[allVariables.index.str.count('M')==2].values)

#add constraint that staples have to be 7 residues apart (hard coded)
peptide_model += (lpSum(np.delete(np.hstack(allVariables.loc[allVariables.index.str.count('M')==2].values),[4,16,27,37,65,70,-4]))==0,'stapled correctly')

    
for c in allVariables.loc[allVariables.index.str.len() > 5].index:
    c1 = c[:4]
    c2 = c[5:]
    peptide_model += (lpSum(allVariables.loc[c] - allVariables.loc[c1] - allVariables.loc[c2]) >= -1, c + '_identity1')
    peptide_model += (lpSum(allVariables.loc[c] - allVariables.loc[c1]) <= 0, c + '_identity2')
    peptide_model += (lpSum(allVariables.loc[c] - allVariables.loc[c2]) <= 0, c + '_identity3')
    
#TODO: add looped constraint to find multiple solutions
# =============================================================================
# #add minimization criteria - single protein
# =============================================================================

#explicit formula for optimum specificity (Mcl-1=0,Bfl-1=1,Bcl-xL=2,Bcl-w=3,Bcl-2=4)
# t = 0
# alpha = 0.25
# off_target_objective = np.sum([-np.dot(coefs[i],np.hstack(allVariables.values.tolist())) for i in range(5) if i != t])
# on_target_objective = np.dot(coefs[t],np.hstack(allVariables.values.tolist()))
# objective = on_target_objective+alpha*off_target_objective
# peptide_model += (objective,'objective')

# =============================================================================
# # #solve the model
# =============================================================================
# peptide_model.solve()

# =============================================================================
# solve the model with adding multiple constraints - single protein
# =============================================================================
# num_solns = 5
# solutions = pd.DataFrame(index=allVariables[allVariables.index.str.len() < 5].index,columns=range(num_solns))
# for n_s in range(num_solns):
#     print('solution' + str(n_s) + 'start')
#     start = time.time()
#     peptide_model.solve()
#     end = time.time()
#     print('solution' + str(n_s) + ' end')
#     print('execution time: ' + str(end-start))
#     solution = pd.DataFrame(index=[str(v.name) for v in peptide_model.variables()],data=[str(v.varValue) for v in peptide_model.variables()])
#     pep_seq = solution.loc[solution[0].astype('float64')==1]
#     pep_seq = pep_seq.loc[pep_seq.index.str.len() < 5]
#     solutions.loc[pep_seq.index,n_s]=1
#     peptide_model += (lpSum(np.hstack(allVariables.loc[pep_seq.index].values))<=22,'solution #: ' + str(n_s+1))
    
# =============================================================================
# # multiple protein solutions
# =============================================================================
num_solns = 5
solution_dict = [pd.DataFrame(index=allVariables[allVariables.index.str.len() < 5].index,columns=range(num_solns)) for i in range(5)]
for t in range(5): #(Mcl-1=0,Bfl-1=1,Bcl-xL=2,Bcl-w=3,Bcl-2=4)
    #define objective - overwrite if t!=0
    print('protein: ' +str(t))
    alpha = 0.25
    off_target_objective = np.sum([-np.dot(coefs[i],np.hstack(allVariables.values.tolist())) for i in range(5) if i != t])
    on_target_objective = np.dot(coefs[t],np.hstack(allVariables.values.tolist()))
    objective = on_target_objective+alpha*off_target_objective
    peptide_model += (objective,'objective')
    num_solns = 5
    solutions = pd.DataFrame(index=allVariables[allVariables.index.str.len() < 5].index,columns=range(num_solns))
    for n_s in range(num_solns):
        print('solution' + str(n_s) + 'start')
        start = time.time()
        peptide_model.solve()
        end = time.time()
        print('solution' + str(n_s) + ' end')
        print('execution time: ' + str(end-start))
        solution = pd.DataFrame(index=[str(v.name) for v in peptide_model.variables()],data=[str(v.varValue) for v in peptide_model.variables()])
        pep_seq = solution.loc[solution[0].astype('float64')==1]
        pep_seq = pep_seq.loc[pep_seq.index.str.len() < 5]
        solution_dict[t].loc[pep_seq.index,n_s]=1
        peptide_model += (lpSum(np.hstack(allVariables.loc[pep_seq.index].values))<=22,'solution #: ' + str(n_s+1) + ' protein #: ' + str(t))

solutions = solutions.drop(labels=['1'])
# =============================================================================
# # look at its solution
# =============================================================================
# t=0
# solution = solutions[t].loc[solutions[t] == 1]
# inv = pd.DataFrame(index=allVariables.loc[allVariables.index.str.len()<5].iloc[1:,:].index)
# inv.loc[solution.index,'0']=1
# inv = inv.fillna(0)
# encoded = np.hstack(inv.values).reshape(1,-1)
# enc.inverse_transform(np.hstack(inv.values).reshape(1,-1))
# svr_mcl1.predict(poly.transform(encoded))

# pep_seq = solution.loc[solution[0].astype('float64')==1]
# pep_seq = pep_seq.loc[pep_seq.index.str.len() < 5]


# for v in peptide_model.variables():
#     test = [v.varValue for v in peptide_model.variables()]
# enc.inverse_transform(np.array(test).reshape(1,-1))

# =============================================================================
# look at dict solutions
# =============================================================================
for t in range(5):
    solution_dict[t].to_excel('mcl1_ilp_solution_'+t+'.xlsx')
    for s in range(5):
        solution = solution_dict[t][s]=solution_dict[t][s].loc[solution_dict[t][s]==1]
        if solution.index[0] == '1':
            solution.drop(index=['1'],inplace=True)
        inv = pd.DataFrame(index=allVariables.loc[allVariables.index.str.len()<5].iloc[1:,:].index)
        inv.loc[solution.index,'0']=1
        inv = inv.fillna(0)
        encoded = np.hstack(inv.values).reshape(1,-1)
        pep_seq = enc.inverse_transform(np.hstack(inv.values).reshape(1,-1))
        encoded_poly = poly.transform(encoded)
        predicted_scores = svr_mcl1.predict(encoded_poly)
        print('protein: ' + str(t) + 'solution #: ' + str(s))
        print('peptide sequence: ' + str(pep_seq))
        print('predicted scores: ' + str(predicted_scores))