import seaborn as sns
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib 
import pandas as pd
import sys
import os
from numpy.random import exponential
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from matplotlib import pyplot as plt
from lifelines.statistics import logrank_test
from matplotlib import rcParams
rcParams.update({'font.size': 16, 'font.family': 'serif'})
rcParams['figure.figsize'] = 10, 7

def draw_KM_curve(virus,taxid,pos_os_event,neg_os_event,outpath):
	ax = plt.subplot(111)
	label1 = r'(+) ' + virus
	label2 = r'(-) ' + virus
	pos_os = []
	pos_event = []
	neg_os = []
	neg_event = []
	for item in pos_os_event.strip().split(';'):
		temp = item.strip().split('_')
		pos_os.append(float(temp[0]))
		pos_event.append(int(temp[1]))
		
	for item in neg_os_event.strip().split(';'):
		temp = item.strip().split('_')
		neg_os.append(float(temp[0]))
		neg_event.append(int(temp[1]))
	kmf_control = KaplanMeierFitter()
	ax = kmf_control.fit(np.array(neg_os), event_observed= np.array(neg_event), label=label2).plot(ax=ax,ci_show=True)
	kmf_exp = KaplanMeierFitter()
	ax = kmf_exp.fit(np.array(pos_os), event_observed= np.array(pos_event),label=label1).plot(ax=ax,ci_show=True)
	add_at_risk_counts(kmf_exp, kmf_control, ax=ax,labels=['Positive Group','Negative Group'])
	results = logrank_test(pos_os,neg_os,event_observed_A=pos_event,event_observed_B=neg_event)
	# logrank pvalue
	p = results.p_value
	textstr = r'$logRank\hspace{0.7}p\hspace{0.4}value=%.4f$'% p
	#textstr = r'$logRank\hspace{0.7}p\hspace{0.4}value=3.9E-21$'
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.618, 0.8, textstr, transform=ax.transAxes, horizontalalignment='center', verticalalignment='top', bbox=props)
	ax.set_ylim(0,1.18)
	ax.set_xlim(0,1900)
	#ax.set_ylabel('Percent progression-free')
	ax.set_ylabel('Percent Overall Survival')
	ax.set_xlabel('Day')
	plt.savefig(os.path.join(outpath,taxid+ '.png'))
	plt.close('all') 

data = open(r'tumor_plot.txt','r')
data.readline()
out_path = os.getcwd()
#for line in data:
#	entry = line.split('\t')
#	virus = entry[1].strip()
#	taxid = entry[4].strip()
#	pos_os_event = entry[11].strip()
#	neg_os_event = entry[12].strip()
#	draw_KM_curve(virus,taxid,pos_os_event,neg_os_event,out_path)
#data.close()
draw_KM_curve('CYT',taxid,pos_os_event,neg_os_event,out_path)


	
	