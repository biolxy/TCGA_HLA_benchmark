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
rcParams.update({'font.size': 16, 'font.family': 'serif'})
rcParams['figure.figsize'] = 10, 7

df = pd.read_csv(filepath_or_buffer="Allen_survival.txt",header='infer',sep='\t',index_col=None,low_memory=False)
pair1 = ['HighTMB','LowTMB']
pair2 = ['B44.present','B44.absent']
pair3 = ['B44.present&HighTMB','B44.absent&LowTMB']
tag1 = 'TMB.group'
tag2 = 'B44_status'
tag3 = 'group'


def cyt_plot(df,tag,groups):
	pos_data = df[df[tag] == groups[0]]
	neg_data = df[df[tag] == groups[1]]
	pos_cyt = pos_data['Immune cytolytic activity (CYT)']
	neg_cyt = neg_data['Immune cytolytic activity (CYT)']
	pos_mean = np.mean(pos_cyt)
	neg_mean = np.mean(neg_cyt)
	fc1 = pos_mean/neg_mean
	if fc1 >1:
		#r_cal1 = robjects.r(r_command_greater1)
		r_call1 = ssw.ttest_ind(np.array(pos_cyt),np.array(neg_cyt),alternative='larger')[1]		
		mannwhitneyu_p1 = scips.mannwhitneyu(pos_cyt,neg_cyt, alternative='greater')[1]	
	else:
		#r_cal1 = robjects.r(r_command_less1)
		r_call1 = ssw.ttest_ind(np.array(pos_cyt),np.array(neg_cyt),alternative='smaller')[1]			
		mannwhitneyu_p1 = scips.mannwhitneyu(pos_cyt,neg_cyt, alternative='less')[1]		
	ax = sns.violinplot(x=tag, y="Immune cytolytic activity (CYT)", order = groups, data=df)
	text0 = ' vs '.join(groups) +'\n'
	text1 = 'mean CYT: ' + '%.2f' % pos_mean + ' vs ' + '%.2f' % neg_mean
	text2 = 'FC = %.2f' % fc1 +'; ' + 'p value = %.3f' % r_call1
	textstr = text0 + '\n' +text1 + '\n' + text2
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.568, 0.98, textstr, transform=ax.transAxes, fontsize=13,horizontalalignment='center', verticalalignment='top', bbox=props)
	ax.set_ylim(0,32)
	#ax.set_xlim(0,54)
	#ax.set_ylabel('Probability of OS(Months)  ', fontsize = 16)
	ax.set_ylabel("Immune cytolytic activity (CYT)", fontsize = 18)
	ax.set_xlabel('Group',fontsize=18)
	plt.savefig(tag +'.CYT.pdf', dpi=330)
	plt.savefig(tag +'.CYT.png', dpi=330)
	plt.close()

cyt_plot(df,tag1,pair1)
cyt_plot(df,tag2,pair2)
cyt_plot(df,tag3,pair3)
	
	