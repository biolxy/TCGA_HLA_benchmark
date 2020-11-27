#coding=utf-8 
from __future__ import division
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as scips
import numpy as np
from matplotlib import rcParams
import os
from lifelines import CoxPHFitter
from numpy.random import exponential
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
import re
import statsmodels.stats.weightstats as ssw
sns.set(font_scale=1.8,style='white')
#import rpy2.robjects as robjects  
rcParams.update({'font.size': 16, 'font.family': 'serif'})
rcParams['figure.figsize'] = 10, 8
df = pd.read_csv(filepath_or_buffer=r"v2.GSE135222_infor.txt", header='infer', sep='\t', index_col=None,low_memory=False)
#sns.violinplot(x="Cancer Stage", y="Immune cytolytic activity (CYT)",hue = 'B44', hue_order = ['present','absent'],data=df,fliersize=0)

ax=sns.violinplot(x="Clinical benefit", y="Immune cytolytic activity (CYT)",data=df,fliersize=0)
#sns.boxplot(x="supertype:Condition", y="Immune cytolytic activity (CYT)",hue = 'supertype status', hue_order = ['present','absent'],data=df,fliersize=0)
#sns.swarmplot(x="Clinical benefit", y="Immune cytolytic activity (CYT)", size=3.5,data=df,dodge=True, color=".25")
#conditions = ['early SKCM','advanced READ','early BRCA','advanced OV']
dcb_df = df[df['Clinical benefit']=='DCB']
ndb_df = df[df['Clinical benefit']=='NDB']
textall = []
text0 = 'Jeong-Yeon Kim NSCLC cohort under anti-PD1 therapy'
textall.append(text0)
textall.append('DCB vs NDB:\n')
pos_record = dcb_df['Immune cytolytic activity (CYT)']
neg_record = ndb_df['Immune cytolytic activity (CYT)']
print('DCB vs NDB:')
textall.append('Samplesize: ' + str(len(pos_record)) + ' vs ' + str(len(neg_record)) + '\n')
print('Samplesize: ' + str(len(pos_record)) + ' vs ' + str(len(neg_record)))


mean_dcb = round(np.mean(dcb_df['Immune cytolytic activity (CYT)']),2)
mean_ndb = round(np.mean(ndb_df['Immune cytolytic activity (CYT)']),2)
std_dcb = round(np.std(dcb_df['Immune cytolytic activity (CYT)']),2)
std_ndb = round(np.std(ndb_df['Immune cytolytic activity (CYT)']),2)
fc1 = round(mean_dcb/mean_ndb,2)
mannwhitneyu_p1 = round(scips.mannwhitneyu(pos_record,neg_record, alternative='greater')[1],2)
r_call1 = round(ssw.ttest_ind(np.array(pos_record),np.array(neg_record),alternative='larger')[1],2)

print('mean CYT in DCB: ' + str(round(mean_dcb,2)))
print('mean CYT in NDB: ' + str(round(mean_ndb,2)))
print('foldchange: ' + str(round(fc1,2)))
print('mannwhitneyu pvalue: ' + str(round(mannwhitneyu_p1,3)))
print('t test pvalue: ' + str(round(r_call1,3)))
textall.append('mean CYT: ' + str(mean_dcb) + ' vs ' + str(mean_ndb) + '\n')
textall.append('Foldchange=' + str(fc1) + '; pvalue=' + str(mannwhitneyu_p1) + '\n')

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.263, 0.98, '\n'.join(textall), transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)

#plt.text(0.3, 0.98, textall, transform=plt.gca().transAxes, fontsize=13,verticalalignment='top', bbox=props)	
ax.set_ylim(-0.9,35)
plt.xticks(rotation=45,horizontalalignment="right")
plt.savefig("clinical_benefit.pdf", dpi=330)
plt.savefig("clinical_benefit.png", dpi=330)

