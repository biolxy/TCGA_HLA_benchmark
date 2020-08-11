import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as scips
import numpy as np
from matplotlib import rcParams
import os
from lifelines import CoxPHFitter
rcParams.update({'font.size': 12, 'font.family': 'serif'})
sys.setrecursionlimit(6000)

def colindex(strings):
    dicts = {}
    entry = strings.rstrip().split('\t')
    for i in range(0,len(entry),1):
        dicts[entry[i]] = i
    return dicts
def multivariableCOXph(dataframe, variable, stage, cancer,outdir):
	cph = CoxPHFitter(penalizer=0.00001)
	cph.fit(dataframe, duration_col= 'OS', event_col= 'OS_event', show_progress=True)
	p = cph.summary
	csv_out = os.path.join(outdir,variable + '.' + stage + '.' + cancer + r'.txt')
	p.to_csv(csv_out, sep='\t', index=True)
	# print cph.hazards_
	#cph.plot(standardized=True)
	#plt.savefig(os.path.join(outdir, variable + '.' + stage + '.' + cancer + ".pdf"))
		#modify coef to HR
	indata = open(csv_out, 'r')
	out = open(os.path.join(outdir, variable + '.' + stage + '.' + cancer + ".HR.txt"), 'w')
	out.write('\t'.join(['multivariable', 'Cancer','stage', 'Pts_pos','Pts_neg','Pts_pos_death','Pts_neg_death','HR', '95% CI Low', '95% CI High', 'P.Value']) + '\n')
	
	title_index = colindex(indata.readline())
	hr_index = title_index['exp(coef)']
	low_index = title_index['exp(coef) lower 95%']
	high_index = title_index['exp(coef) upper 95%']
	p_index = title_index['p']
	for line in indata:
		entry = line.split('\t')
		if len(entry) >= 8:
		#hr = str(entry[hr_index])[0:5]
		#low = str(np.exp(float(entry[low_index])))[0:5]
		#high = str(np.exp(float(entry[high_index])))[0:5]
		#p = str(entry[p_index])[0:5]
			hr = str(entry[hr_index])
			#if entry[low_index].strip() == 'NaN' or entry[low_index].strip() == '':
				#low = 'NaN'
			#else:
				#low = str(np.exp(float(entry[low_index])))
			#if entry[high_index].strip() == 'NaN' or entry[high_index].strip() == '':	
			#	high = 'NaN'
			#else:
			#	high = str(np.exp(float(entry[high_index])))
			low = str(entry[low_index])
			high = str(entry[high_index])
			p = str(entry[p_index])	
			pos_all = dataframe[dataframe[entry[0].strip()] == 1]
			neg_all = dataframe[dataframe[entry[0].strip()] == 0]
			pos_all_death = pos_all[pos_all['OS_event'] == 1]
			neg_all_death = neg_all[neg_all['OS_event'] == 1]
			out.write(entry[0].strip() + '\t'+ '\t'.join([cancer,stage]) + '\t' + str(len(pos_all['OS'])) + '\t' + str(len(neg_all['OS'])) + '\t' + str(len(pos_all_death['OS'])) + '\t' + str(len(neg_all_death['OS'])) + '\t' + '\t'.join([hr, low, high, p]) + '\n')
		else:
			print('Error: short length in ' + cancer  + ' ' + stage + ' ' + variable)
	
	
	
df = pd.read_csv(filepath_or_buffer="notLOH_homo_uncertain_failcall_vs_LOH_absent.txt",
                 header='infer', sep='\t', index_col=None)
outpath = os.path.join(os.getcwd(),'py3_LOH_multivariable')
#print outpath

if not os.path.exists(outpath):
	os.mkdir(outpath)
	
#get parameter list from file
data = open(r'test_notLOH_vs_LOH_absent.HR.plot.txt','r')
data.readline()
for line in data:
	entry = line.split('\t')
	variable = entry[0].strip()
	#cancer = entry[1].strip()
	#stage = entry[2].strip()
	stage, cancer = entry[1].strip().split(' ')	
	print('processing ' + ' '.join([variable,cancer,stage]))
	df1 = df[df['project_id'] == 'TCGA-' + cancer]
	df2 = df1[df1['tumor_stage_type'] == stage]
	df_valid = df2[df2['OS_tag'].str.contains('valid')]
	df_valid2 = df_valid[df_valid['age_tag'].isin(['1','0'])]
	df_valid3 = df_valid2[df_valid2['v2_TMB_binary'].isin(['1','0'])]
	#df_valid4 = df_valid3[df_valid3['gender'].isin(['1','0'])]
	if variable.find(r'_') == -1:
		df_valid5 = df_valid3[df_valid3[variable].isin(['1','0'])]
		df3 = df_valid5[['OS','OS_event', variable, 'v2_TMB_binary','age_tag']].apply(pd.to_numeric)
		multivariableCOXph(df3, variable, stage, cancer,outpath)
	else:
		variable0, variable1 = variable.split(r'_')
		df_valid5 = df_valid3[df_valid3[variable0].isin(['1','0'])]
		df3 = df_valid5[['OS','OS_event', variable0, variable1, 'v2_TMB_binary','age_tag']].apply(pd.to_numeric)
		multivariableCOXph(df3, variable, stage, cancer,outpath)
data.close()
