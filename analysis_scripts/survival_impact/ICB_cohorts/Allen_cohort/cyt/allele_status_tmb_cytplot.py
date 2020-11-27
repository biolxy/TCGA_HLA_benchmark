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
rcParams.update({'font.size': 13, 'font.family': 'serif'})

#def cya_plot(df, v1_len, v2_len, mean_v1, mean_v2, std_v1, std_v2, fc,pvalue, outdir, allele, out_prefix):
def cya_plot(df, outdir, out_prefix, x_name):
	#fig,ax = plt.subplots()
	#fig = plt.figure()
	#sns.set()
	#class_mapping = {'0':'absent','1':'present'}
	#df[allele] = df[allele].map(class_mapping)
	#df.apply(pd.to_numeric,errors='ignore')
	df['Immune cytolytic activity (CYT)'] = df['Immune cytolytic activity (CYT)'].astype('float64')
	ax = plt.subplot(111)
	sns.violinplot(x=allele, y="Immune cytolytic activity (CYT)", data=df)
	#sns.boxplot(x=x_name, y="Immune cytolytic activity (CYT)",data=df,fliersize=0)
	#sns.swarmplot(x=x_name, y="Immune cytolytic activity (CYT)",size=4, data=df,dodge=True, color=".25")
	#textstr0 = allele + ' present group: size = %.0f, mean CYT = %.2f, std = %.2f' % (v1_len, mean_v1,std_v1)
	#textstr1 = allele + ' absent group: size = %.0f, mean CYT = %.2f, std = %.2f' % (v2_len, mean_v2,std_v2)
	#textstr0 = allele + ' present group: size = %.0f, mean CYT = %.2f' % (v1_len, mean_v1)
	#textstr1 = allele + ' absent group: size = %.0f, mean CYT = %.2f' % (v2_len, mean_v2)
	#print textstr0
	#print textstr1
	#textstr2 = allele + ' present vs absent: FC =%.2f, p value = %.3f' % (fc,pvalue)
	#textall = '\n'.join([textstr0,textstr1,textstr2])
	#props = dict(boxstyle='round', facecolor='wheat', alpha=0.25)
	#ax.text(0.03, 0.98, textall, transform=ax.transAxes, fontsize=12,verticalalignment='top', bbox=props)	
	#plt.text(0.03, 0.98, textall, transform=plt.gca().transAxes, fontsize=12,verticalalignment='top', bbox=props)	
	ax.set_ylim(-1,55)
	#plt.ylim(0, 150)
	#plt.show()
	plt.savefig(os.path.join(outdir,out_prefix + '_cyt.pdf'))
	#plt.cla()
	#plt.clf()
	plt.close('all')

def get_data_plot(dataframe, cancer, stage, allele,outdir1, outdir2, out_prefix1, out_prefix2):
	project = cancer
	cancer_df = dataframe[dataframe['project_id'] == project] #select cancer record
	stage_df = cancer_df[cancer_df['Study'] == stage] # select_stage record
	os_df = stage_df[stage_df['OS_tag'] == 'valid'] # slect valid os record
	cyt_df = os_df[os_df['cya_tag'] == 'valid'] # select valid cyt record
	allele_pos_df = cyt_df[cyt_df[allele] == 'present']
	allele_neg_df = cyt_df[cyt_df[allele] == 'absent']
	allele_pos_htmb = allele_pos_df[allele_pos_df['TMB'] == 'High']
	allele_neg_htmb = allele_neg_df[allele_neg_df['TMB'] == 'High']
	allele_pos_ltmb = allele_pos_df[allele_pos_df['TMB'] == 'Low']
	allele_neg_ltmb = allele_neg_df[allele_neg_df['TMB'] == 'Low']
	print 'processing ' + allele
	allele_pos_htmb[allele + r'&HighTMB'] = allele + r'.present&HighTMB'
	allele_neg_ltmb[allele + r'&HighTMB'] = allele + r'.absent&LowTMB'		
	#comparing1 allele.present.HighTMB vs allele.absent.LowTMB
	frames1 = [allele_pos_htmb,allele_neg_ltmb]
	compare1 = pd.concat(frames1)
	cya_plot(compare1, outdir1, out_prefix1, allele + r'&HighTMB')
	allele_pos_ltmb[allele + r'&LowTMB'] = allele + r'.present&LowTMB'
	allele_neg_htmb[allele + r'&LowTMB'] = allele + r'.absent&HighTMB'
	frames2 = [allele_pos_ltmb,allele_neg_htmb]
	compare2 = pd.concat(frames2)
	cya_plot(compare2, outdir2, out_prefix2, allele + r'$LowTMB')
		


df = pd.read_csv(filepath_or_buffer="v4_allen_new_cyt_os.txt",header='infer',sep='\t',index_col=None,low_memory=False)
out_dir1 = os.path.join(os.getcwd(),'allele_hightmb')
out_dir2 = os.path.join(os.getcwd(),'allele_lowtmb')
if not os.path.exists(out_dir1):
	os.mkdir(out_dir1)
if not os.path.exists(out_dir2):
	os.mkdir(out_dir2)

indata =open(r'superhrresults.HR.xls','r')	
indata.readline()

for line in indata:
	entry = line.split('\t')
	allele = entry[0].strip()
	cancer = entry[1].strip()
	stage = entry[2].strip()
	out_prefix1 = allele + '_highTMB'
	out_prefix2 = allele + '_lowTMB'
	get_data_plot(df,cancer,stage, allele,out_dir1,out_dir2,out_prefix1,out_prefix2)
indata.close()


	
	
	