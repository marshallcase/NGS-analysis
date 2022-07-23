# -*- coding: utf-8 -*-
"""
Created on Sat May 14 13:42:21 2022

@author: Marshall
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.cluster import KMeans
import logomaker
from plot_ML import *
from skopt import BayesSearchCV
from skopt.plots import plot_objective
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.svm import LinearSVR
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler


#import data
scores_vis = pd.read_excel('scores.xlsx')
scores_vis.set_index('peptide',inplace=True)
scores_vis = scores_vis.loc[~scores_vis.index.str.contains('_')]
scores_vis = scores_vis.drop(columns=scores_vis.columns[5:])
#create feature matrix
features = pd.DataFrame(index=scores_vis.index)
for i in range(23):
    features[i]=features.index.str[i]
    
#one hot encode feature matrix
enc = OneHotEncoder(sparse=False)
enc.fit(features)
X = enc.transform(features)

#define label vector
y = scores_vis.iloc[:,:] #scores_vis.iloc[:,[5,7,9,11,13]]
scaler = StandardScaler()
y_scaled = scaler.fit_transform(y)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
_, _, y_train_scaled, y_test_scaled = train_test_split(X, y_scaled, test_size=0.3, random_state=42)

# =============================================================================
# pick protein
# =============================================================================
i=0

# # =============================================================================
# # train Mcl-1 score predictor  linear regression
# # =============================================================================
# linreg_mcl1 = LinearRegression()
# linreg_mcl1.fit(X_train,y_train.iloc[:,i])
# plt.scatter(linreg_mcl1.predict(X_train),y_train.iloc[:,i])
# plt.scatter(linreg_mcl1.predict(X_test),y_test.iloc[:,i])

# # =============================================================================
# # train Mcl-1 score predictor Support vector regressor
# # =============================================================================
# svr_mcl1 = SVR(C=10,epsilon=0.1,kernel='linear')
# svr_mcl1.fit(X_train,y_train.iloc[:,i])
# plt.scatter(svr_mcl1.predict(X_train),y_train.iloc[:,i])
# plt.scatter(svr_mcl1.predict(X_test),y_test.iloc[:,i])

# =============================================================================
# train SVR with pair-wise features
# =============================================================================
poly = PolynomialFeatures(2,interaction_only=True)
X_2 = poly.fit_transform(pd.DataFrame(X,columns=enc.get_feature_names()))
X_train_2, X_test_2, y_train, y_test = train_test_split(X_2, y, test_size=0.3, random_state=42)

# features_titles = poly.get_feature_names(input_features=enc.get_feature_names())
# analysis = pd.DataFrame(X_2,columns=features_titles)
# #drop zero columns
# analysis = analysis.loc[:,~(analysis.sum(axis=0)==0)]


# # svr_mcl1 = SVR(C=10,epsilon=0.1,kernel='linear')
# # svr_mcl1.fit(X_train_2,y_train.iloc[:,i])

# svr_mcl1 = LinearSVR(epsilon=0.1,C=10)
# svr_mcl1.fit(X_train_2,y_train.iloc[:,i])

# plt.scatter(svr_mcl1.predict(X_train_2),y_train.iloc[:,i])
# plt.scatter(svr_mcl1.predict(X_test_2),y_test.iloc[:,i])

# =============================================================================
# evaluate model
# =============================================================================
# mean_squared_error(y_train.iloc[:,0],svr_mcl1.predict(X_train))
# mean_squared_error(y_test.iloc[:,0],svr_mcl1.predict(X_test))
# mean_squared_error(y_test.iloc[:,0],svr_mcl1.predict(X_test),squared=False)
# mean_squared_error(y_train.iloc[:,0],svr_mcl1.predict(X_train),squared=False)

# =============================================================================
# gridsearchCV - SVR
# =============================================================================
# param_grid = [
#   # {'C': [1, 10, 100, 1000], 'epsilon': [0.05,0.1,0.25],'kernel': ['linear']},
#   {'C': [1, 10], 'epsilon': [0.05,0.1,0.25],'gamma': [0.001, 0.0001], 'kernel': ['poly'],'degree':[2,3]}]

# clf = GridSearchCV(estimator=SVR(),param_grid=param_grid,n_jobs=-1,verbose=4)
# clf.fit(X_train,y_train.iloc[:,i])

# =============================================================================
# multioutput SVR regressor
# =============================================================================
# regr = MultiOutputRegressor(SVR(C=10,kernel='poly',degree=2),n_jobs=-1).fit(X_train,y_train)

# =============================================================================
# random forest regressor
# =============================================================================
# rfr_multi = MultiOutputRegressor(RandomForestRegressor(n_estimators=100)).fit(X_train,y_train)
# rfr = RandomForestRegressor(n_estimators=100)
# rfr.fit(X_train,y_train.iloc[:,0])

# =============================================================================
# grid search CV - RGR
# =============================================================================
# param_grid = [c
#   {'estimator__max_depth': [5, 10, 20], 'estimator__max_features': [2,4,8],'estimator__max_leaf_nodes':[100,1000,5000]}]
# clf = GridSearchCV(estimator=MultiOutputRegressor(RandomForestRegressor(n_estimators=100)),param_grid=param_grid,n_jobs=-1,
#                    verbose=4)
# clf.fit(X_train,y_train)

# =============================================================================
# look at the UMAP distribution of peptides
# =============================================================================
# reducer = umap.UMAP()
# reducer.fit(X)
# embeddings = reducer.transform(X)

# n_clusters = 5
# kmc = KMeans(n_clusters=n_clusters).fit(X)
# labels = kmc.labels_
# data = pd.DataFrame(index=features.index,data=labels)

# fig,ax=plt.subplots()
# scatter_x = embeddings[:,0]
# scatter_y = embeddings[:,1]
# group = labels
# for g in np.unique(group):
#     i = np.where(group == g)
#     ax.scatter(scatter_x[i], scatter_y[i], label=g)
# ax.legend()
# plt.show()

# fig,axs = plt.subplots(nrows=n_clusters,ncols=1,figsize=(15,10),sharex=True)
# for cluster,ax in enumerate(axs.ravel()):
#     plot_data = data.loc[data[0]==cluster]
#     logomatrix = logomaker.alignment_to_matrix(sequences=list(plot_data.index),to_type='information',
#                                                 characters_to_ignore='_')
#     logoplot = logomaker.Logo(logomatrix,ax=ax,color_scheme='NajafabadiEtAl2017')
#     logoplot.ax = ax
#     logoplot.ax.set_ylabel('information (bits)')
#     logoplot.ax.set_xticks(np.arange(0,23))
#     logoplot.ax.set_title('number of peptides: ' + str(len(plot_data))+', cluster = ' + str(cluster))
#     # plt.savefig(str(p) + '.svg')
#     # plt.close()

# =============================================================================
# bayesian optimization of hyperparameters
# =============================================================================
# opt = BayesSearchCV(
#     MultiOutputRegressor(LinearSVR()),
#     {
#         'estimator__C': (1e-6, 1e+6, 'log-uniform'),
#         'estimator__epsilon': (1e-6, 1e+1, 'log-uniform'),
#         'estimator__tol': (1e-6,1e-3,'log-uniform')
#     },
#     n_iter=32,
#     cv=5,
#     n_jobs=-1,
#     verbose=10
# )
# opt = BayesSearchCV(
#     RandomForestRegressor(),
#     {
#         'max_depth': (2,23),
#         'min_samples_split': (2,6),
#         'n_estimators': (10, 500, 'log-uniform'),
#         'max_features': (1, 5),  # integer valued parameter
#     },s
#     n_iter=32,
#     cv=5,
#     n_jobs=-1,
#     return_train_score=True,
#     verbose=10
# )

# print('fit')
# opt.fit(X_train_2, y_train_scaled)
# print("val. score: %s" % opt.best_score_)
# print("test score: %s" % opt.score(X_test, y_test.iloc[:,0]))


#bayessearchcv multi:
svr_mcl1 = MultiOutputRegressor(estimator=LinearSVR(C=0.005001177728215577,
                                          epsilon=0.009344323060800772,
                                          tol=0.001))
svr_mcl1.fit(X_train_2,y_train)
# =============================================================================
# plot bayesian optimization
# =============================================================================

# plot_objective(opt.optimizer_results_[0],
#                     dimensions=["C", "epsilon", "tol"])
# plt.tight_layout()

# =============================================================================
# multioutput from bayesian optimization
# =============================================================================
# rfr_multi = MultiOutputRegressor(opt.best_estimator_)
# rfr_multi.fit(X_train,y_train)

# =============================================================================
# scipy optimize input peptide - does not work well
# =============================================================================
# from scipy.optimize import minimize
# bnds = [(0,1)]*220
# x0 = enc.transform(np.array(features.iloc[1000]).reshape(1,-1))

# rfr = RandomForestRegressor(n_estimators=500,max_depth=23)
# rfr.fit(X_train,y_train.iloc[:,0])

# optimize = minimize(lambda x: rfr.predict(np.array([x])), x0, method='COBYLA',bounds=bnds,options = {'maxiter':200})
# output = enc.inverse_transform(np.array(optimize.x).reshape(1,-1))
