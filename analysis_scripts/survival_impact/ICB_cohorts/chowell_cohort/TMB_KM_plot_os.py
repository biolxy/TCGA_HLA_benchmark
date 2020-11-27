from numpy.random import exponential
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
import numpy as np 
from matplotlib import pyplot as plt
from lifelines.statistics import logrank_test
from matplotlib import rcParams
from lifelines.utils import median_survival_times
rcParams.update({'font.size': 16, 'font.family': 'serif'})
rcParams['figure.figsize'] = 10, 7

def coldict(strings):
	entry = strings.strip().split('\t')
	dicts = {}
	for i in range(0,len(entry),1):
		dicts[entry[i]] = i
	return dicts
	

ax = plt.subplot(111)
data = open(r'Chowell_melanoma_CTLA4.txt','r')
first = data.readline().strip()
colindex = coldict(first)
#tag1 = 'High TMB'
#tag2 = 'Low TMB'
os_index = colindex['OS(Months)']
event_index = colindex['OS_Event']
cyt_index = colindex['TMB.group']

tag1 = 'HighTMB'
tag2 = 'LowTMB'

#define PFS, OS
T_low = []
T_high = []
#define event such 0 ,1
E_low = []
E_high = []
pos_size = 0
pos_death = 0
neg_size = 0
neg_death = 0
for line in data:
	entry = line.strip().split('\t')
	if entry[cyt_index] == tag2:
		T_low.append(float(entry[os_index]))
		E_low.append(float(entry[event_index]))
		if int(entry[event_index]) == 1:
			neg_death +=1
		neg_size +=1
	else:
		T_high.append(float(entry[os_index]))
		E_high.append(float(entry[event_index]))
		if int(entry[event_index]) == 1:
			pos_death +=1
		pos_size +=1		
data.close()


kmf_exp = KaplanMeierFitter()
ax = kmf_exp.fit(np.array(T_high), label='HighTMB').plot(ax=ax,ci_show=True,linewidth=3)
#ax = kmf_exp.fit(np.array(T_high), event_observed= np.array(E_high),label=tag1).plot(ax=ax,ci_show=False)
exp_median = kmf_exp.median_survival_time_
#test_median = kmf_exp.median_
#print str(test_median)
#exp_median = kmf_exp.median_
#exp_median_confidence_interval_ = median_survival_times(kmf_exp.confidence_interval_)

print('Positive group median OS: ' + str(float(exp_median)))

kmf_control = KaplanMeierFitter()
ax = kmf_control.fit(np.array(T_low), label='LowTMB').plot(ax=ax,ci_show=True,linewidth=3)
control_median = kmf_control.median_survival_time_
#control_median = kmf_control.median_
#control_median_confidence_interval_ = median_survival_times(kmf_control.confidence_interval_)

print('negative group median OS: ' + str(float(control_median)))
#ax = kmf_control.fit(np.array(T_low), event_observed= np.array(E_low), label=tag2).plot(ax=ax,ci_show=False)
#add_at_risk_counts(kmf_exp, kmf_control, ax=ax)
results = logrank_test(T_high,T_low,event_observed_A=E_high,event_observed_B=E_low)
# logrank pvalue
p = results.p_value
text0 = 'Chowell Melanoma cohorts under anti-CTLA4 therapy' + '\n'
text0_1 = 'HighTMB vs LowTMB'
text0_2 = 'deaths/sample_size: ' + str(pos_death) + '/' + str(pos_size) + ' vs ' + str(neg_death) + '/' + str(neg_size)
text1 = r'Median OS: %0.1f'%exp_median + 'm vs ' + '%0.1f'%control_median + 'm'
#text2 = r'LowTMB Median OS: %0.1f'%control_median
text3 = r'logRank p value = %.4f'%p
#text3 = r'$logRank\hspace{0.7}p\hspace{0.4}value=%.4f$'%p
textstr = text0 + '\n' +text0_1 + '\n' + text0_2 + '\n'+ text1 + '\n' + text3
#textstr = r'$logRank\hspace{0.7}p\hspace{0.4}value=3.9E-21$'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.548, 0.8, textstr, transform=ax.transAxes, fontsize=13,horizontalalignment='center', verticalalignment='top', bbox=props)
ax.set_ylim(0,1.0)
ax.set_xlim(0,95)
#ax.set_ylabel('Probability of OS(Months)  ', fontsize = 16)
ax.set_ylabel('Probability of overall survival  ', fontsize = 18)
ax.set_xlabel('Month',fontsize=18)
plt.savefig("Chowell_melanoma_CTLA4_TMB.png", dpi=330)
#plt.show()