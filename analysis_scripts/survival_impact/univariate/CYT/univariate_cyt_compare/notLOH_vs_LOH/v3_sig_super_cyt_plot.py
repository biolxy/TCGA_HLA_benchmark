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
#import rpy2.robjects as robjects  
rcParams.update({'font.size': 12, 'font.family': 'serif'})


def get_stats(df1, allele, condition, parameter):
	df1[parameter] = df1[parameter].astype('float64')
	df = df1[df1['Supertype:condition'] == condition]
	pos_record = df[df[allele] == 'notLOH'][parameter]
	neg_record = df[df[allele] == 'absent/LOH'][parameter]
	pos_len = len(pos_record)
	neg_len = len(neg_record)
	text0 = condition + ' notLOH ' + str(pos_len) + ' absent/LOH: ' + str(neg_len)
	print text0
	mean_pos = np.mean(pos_record)
	mean_neg = np.mean(neg_record)
	std_pos = np.std(pos_record)
	std_neg = np.std(neg_record)
	fc1 = mean_pos/mean_neg
	#r_command_greater1 = r't.test(c(' + ','.join([str(temp) for temp in pos_record]) + '), y=c(' + ','.join([str(i) for i in neg_record]) + '), alternative="greater")$p.value'
	#r_command_less1 = r't.test(c(' + ','.join([str(temp) for temp in pos_record]) + '), y=c(' + ','.join([str(i) for i in neg_record]) + '), alternative="less")$p.value'
	if fc1 >1:
		#r_cal1 = robjects.r(r_command_greater1)
		r_call1 = ssw.ttest_ind(np.array(pos_record),np.array(neg_record),alternative='larger')		
		mannwhitneyu_p1 = scips.mannwhitneyu(pos_record,neg_record, alternative='greater')[1]	
	else:
		#r_cal1 = robjects.r(r_command_less1)
		r_call1 = ssw.ttest_ind(np.array(pos_record),np.array(neg_record),alternative='smaller')			
		mannwhitneyu_p1 = scips.mannwhitneyu(pos_record,neg_record, alternative='less')[1]	
	#r_cal1_p = list(r_cal1)[0]	
	r_cal1_p = list(r_call1)[1]	
	infors = condition.split(r' ')
	superallele = infors[0].split(':')[0]
	infors_stage = infors[0].split(':')[1]
	infors_cancer = infors[1].strip()
	text1 = superallele + ': ' + infors_stage + ' stage ' + infors_cancer
	#text2 = allele + ' notLOH group: size = %.0f, mean CYT = %0.2f, std = %0.2f' % (pos_len, mean_pos,std_pos)
	#text3 = allele + ' absent/LOH group: size = %.0f, mean CYT = %0.2f, std = %0.2f' % (neg_len, mean_neg,std_neg)
	#text2 = superallele + ' notLOH group: size = %.0f, mean CYT = %0.2f' % (pos_len, mean_pos)
	#text3 = superallele + ' absent/LOH group: size = %.0f, mean CYT = %0.2f' % (neg_len, mean_neg)
	text2 = superallele + ' notLOH group: mean CYT = %0.2f' % mean_pos
	text3 = superallele + ' absent/LOH group: mean CYT = %0.2f' % mean_neg
	text4 = superallele + ' notLOH vs absent/LOH: FC = %.2f, p value = %.3f' % (fc1, r_cal1_p)
	text5 = ''
	return '\t'.join([text0, text1,text2, text3, text4,text5])
	
df = pd.read_csv(filepath_or_buffer=r"select_super.records.xls", header='infer', sep='\t', index_col=None,low_memory=False)
#sns.violinplot(x="Cancer Stage", y="Immune cytolytic activity (CYT)",hue = 'B44', hue_order = ['notLOH','absent/LOH'],data=df,fliersize=0)
ax = plt.subplot(111)
sns.violinplot(x="Supertype:condition", y="Immune cytolytic activity (CYT)",hue = 'Supertype status', hue_order = ['notLOH','absent/LOH'],data=df,fliersize=0)
#sns.boxplot(x="Supertype:condition", y="Immune cytolytic activity (CYT)",hue = 'Supertype status', hue_order = ['notLOH','absent/LOH'],data=df,fliersize=0)
#sns.swarmplot(x="Supertype:condition", y="Immune cytolytic activity (CYT)", hue = 'Supertype status',hue_order = ['notLOH','absent/LOH'],size=3.5,data=df,dodge=True, color=".25")
#conditions = ['early SKCM','advanced READ','early BRCA','advanced OV']
conditions = ['A03:advanced READ', 'A03:advanced KIRP', 'A24:early STAD', 'A02:advanced UCEC', 'A02:early LIHC', 'A03:early LIHC']
conditions2 = ['B44:early SKCM', 'B44:advanced OV', 'A24:advanced OV', 'A02:early SKCM', 'B44:early BRCA', 'B08:early UCEC']
conditions3 = ['A03:early CESC', 'A01:early HNSC', 'B08:early CESC', 'A02:advanced MESO', 'B07:advanced CESC', 'B08:early HNSC']

textall = []
textall2 = []
textall3 = []
for temp in conditions:
	textall.append(get_stats(df, 'Supertype status', temp, 'Immune cytolytic activity (CYT)'))

for temp2 in conditions2:
	textall2.append(get_stats(df, 'Supertype status', temp2, 'Immune cytolytic activity (CYT)'))
for temp3 in conditions3:
	textall3.append(get_stats(df, 'Supertype status', temp3, 'Immune cytolytic activity (CYT)'))	
	
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.25)
#ax.text(0.02, 0.99, '\n'.join(textall), transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
#ax.text(0.33, 0.99, '\n'.join(textall2), transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)	
#ax.text(0.65, 0.99, '\n'.join(textall3), transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)	

new_out = open(r'fig_stats.txt','w')
print >> new_out, '\n'.join(textall)
print >> new_out, '\n'.join(textall2)
print >> new_out, '\n'.join(textall3)
new_out.close()
	#plt.text(0.3, 0.98, textall, transform=plt.gca().transAxes, fontsize=13,verticalalignment='top', bbox=props)	
ax.set_ylim(-0.9,58)
plt.xticks(rotation=45,horizontalalignment="right")
plt.show()
