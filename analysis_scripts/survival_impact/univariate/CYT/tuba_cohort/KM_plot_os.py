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
	for i in xrange(len(entry)):
		dicts[entry[i]] = i
	return dicts
	

ax = plt.subplot(111)
data = open(r'tuba_cyt_os.txt','r')
first = data.readline().strip()
colindex = coldict(first)
#tag1 = 'High TMB'
#tag2 = 'Low TMB'
os_index = colindex['OS']
event_index = colindex['Event']
cyt_index = colindex['CYT_group']

tag1 = 'CYT-High'
tag2 = 'CYT-Low'

#define PFS, OS
T_low = []
T_high = []
#define event such 0 ,1
E_low = []
E_high = []
for line in data:
	entry = line.strip().split('\t')
	if entry[cyt_index] == tag2:
		T_low.append(float(entry[os_index]))
		E_low.append(float(entry[event_index]))
	else:
		T_high.append(float(entry[os_index]))
		E_high.append(float(entry[event_index]))
		
data.close()


kmf_exp = KaplanMeierFitter()
ax = kmf_exp.fit(np.array(T_high), label='CYT-High group').plot(ax=ax,ci_show=True,linewidth=3)
#ax = kmf_exp.fit(np.array(T_high), event_observed= np.array(E_high),label=tag1).plot(ax=ax,ci_show=False)
#exp_median_ = kmf_exp.median_survival_time_
#test_median = kmf_exp.median_
#print str(test_median)
exp_median = kmf_exp.median_
#exp_median_confidence_interval_ = median_survival_times(kmf_exp.confidence_interval_)

print 'Positive group median OS: ' + str(float(exp_median))

kmf_control = KaplanMeierFitter()
ax = kmf_control.fit(np.array(T_low), label='CYT-Low group').plot(ax=ax,ci_show=True,linewidth=3)

#control_median_ = kmf_control.median_survival_time_
control_median = kmf_control.median_
#control_median_confidence_interval_ = median_survival_times(kmf_control.confidence_interval_)

print 'negative group median OS: ' + str(float(control_median))
#ax = kmf_control.fit(np.array(T_low), event_observed= np.array(E_low), label=tag2).plot(ax=ax,ci_show=False)
add_at_risk_counts(kmf_exp, kmf_control, ax=ax)
results = logrank_test(T_high,T_low,event_observed_A=E_high,event_observed_B=E_low)
# logrank pvalue
p = results.p_value
text0 = 'Tuba Nur Gide melanoma cohort' + '\n'
text1 = r'CYT-High group Median OS: %0.1f'%exp_median
text2 = r'CYT-Low group Median OS: %0.1f'%control_median
text3 = r'logRank p value %.4f'%p
#text3 = r'$logRank\hspace{0.7}p\hspace{0.4}value=%.4f$'%p
textstr = text0 + '\n' + text1 + '\n' + text2 + '\n' + text3
#textstr = r'$logRank\hspace{0.7}p\hspace{0.4}value=3.9E-21$'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.718, 0.8, textstr, transform=ax.transAxes, fontsize=13,horizontalalignment='center', verticalalignment='top', bbox=props)
ax.set_ylim(0,1.0)
ax.set_xlim(0,1400)
#ax.set_ylabel('Percent progression-free')
ax.set_ylabel('Percent Overall Survival', fontsize = 16)
ax.set_xlabel('Day',fontsize=16)
plt.savefig("Tuba.pdf", dpi=300)
#plt.show()