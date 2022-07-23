# -*- coding: utf-8 -*-
"""
Created on Wed May 18 13:51:12 2022

@author: marsh

"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats import pearsonr
# =============================================================================
# plot KDE of predictions
# =============================================================================
def plotKDE(pred,y_train,gate_score=True):
    fig,axs = plt.subplots(1,5,sharey=True,figsize=(15,5))
    if gate_score:
        xmin=1
        xmax=12
    else:
        xmin = min(pred.min(axis=0).min(),y_train.min(axis=0).min())
        xmax = max(pred.max(axis=0).max(),y_train.max(axis=0).max())
    X, Y = np.mgrid[xmin:xmax:100j, xmin:xmax:100j]
    proteins=['M','F','X','W','2']
    prot_dict = dict(zip(range(5),proteins))
    positions = np.vstack([X.ravel(), Y.ravel()])
    for i,ax in enumerate(axs.ravel()):
        m1=pred[:,i]
        m2=y_train.iloc[:,i]
        corr = pearsonr(m1,m2)[0]
        kernel = gaussian_kde((m1,m2))
        Z = np.reshape(kernel(positions).T, X.shape)
        ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
          extent=[xmin, xmax, xmin, xmax])
        ax.plot(m1, m2, 'k.', markersize=2)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([xmin, xmax])
        ax.set_title('Protein: ' + prot_dict[i] + ', pearson corr: '+ str(np.round(corr,3)))
    # fig.text(0.5, 0.04, '# of peptides', ha='center')
    # fig.text(0.04, 0.5, 'gate score (lower better)', va='center', rotation='vertical')
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('predicted score')
    plt.ylabel('actual score')
    return fig

def plot_loss(history):
  plt.plot(history.history['loss'], label='loss')
  plt.plot(history.history['val_loss'], label='val_loss')
  plt.xlabel('Epoch')
  plt.ylabel('Error')
  plt.legend()
  plt.grid(True)