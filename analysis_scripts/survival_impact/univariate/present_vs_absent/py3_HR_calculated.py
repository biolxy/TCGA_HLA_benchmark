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
    for i in range(1,len(entry),1):
        dicts[entry[i]] = i
    return dicts

def univariableCOXph(dataframe,event,time,variable,stage,cancer,outdir):
	print('processing ' + cancer + ' ' + stage + ' ' + variable)
	#print outdir
	#print os.path.join(outdir,variable)
	#df1_1 = dataframe[dataframe['OS_tag'].str.contains('valid')]
	#early stage
	#df_early = df1_1[df1_1['tumor_stage_type'].str.contains('early')]
	#df_advanced = df1_1[df1_1['tumor_stage_type'].str.contains('advanced')]
	df3 = dataframe[[time, event, variable]].apply(pd.to_numeric)
	tag1_df = df3[df3[variable]== 1]
	tag0_df = df3[df3[variable]== 0]
	tag1_num = str(len(tag1_df[variable]))
	tag0_num = str(len(tag0_df[variable]))
	#tag0 death event
	tag0_death = tag0_df[tag0_df[event] == 1]
	tag1_death = tag1_df[tag1_df[event] == 1]
	tag0_death_num = len(tag0_death[variable])
	tag1_death_num = len(tag1_death[variable])
	tag0_survival_num = int(tag0_num) - tag0_death_num
	tag1_survival_num = int(tag1_num) - tag1_death_num
	tag1_list = [tag1_death_num, tag1_survival_num]
	tag0_list = [tag0_death_num, tag0_survival_num]

	#chisquare_p = scips.chi2_contingency(np.array([tag1_list,tag0_list]))[1]
	if tag0_death_num >=3 and tag1_death_num >=3:
		fisher_p = scips.fisher_exact([tag1_list,tag0_list],alternative='two-sided')[1]
		#chisquare_p = scips.chi2_contingency(np.array([tag1_list,tag0_list]), lambda_="log-likelihood")[1]
		chisquare_p = scips.chi2_contingency(np.array([tag1_list,tag0_list]))[1]		
		#cph = CoxPHFitter(baseline_estimation_method = 'spline', n_baseline_knots = 1)
		cph = CoxPHFitter(penalizer=0.00001)
		cph.fit(df3, duration_col= time, event_col= event,show_progress=True)
		p = cph.summary
		csv_out = os.path.join(outdir,variable + '.' + stage + '.' + cancer + r'.txt')
		p.to_csv(csv_out, sep='\t', index=True)
		# print cph.hazards_
		#cph.plot(standardized=True)
		#plt.savefig(os.path.join(outdir, variable + '.' + stage + '.' + cancer + ".pdf"))
		#modify coef to HR
		indata = open(csv_out, 'r')
		out = open(os.path.join(outdir, variable + '.' + stage + '.' + cancer + ".HR.txt"), 'w')
		out.write('\t'.join(['variable', 'Cancer','stage',variable+r'_neg',variable + '_pos',variable+r'_neg_death',variable+r'_pos_death','HR', 'CI_Low', 'CI_High', 'pvalue','fisher_p','chisquare_p']) + '\n')
		title_index = colindex(indata.readline())
		hr_index = title_index['exp(coef)']
		low_index = title_index['exp(coef) lower 95%']
		high_index = title_index['exp(coef) upper 95%']
		p_index = title_index['p']
		for line in indata:
			entry = line.split('\t')
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
			out.write(entry[0].strip() + '\t'+ '\t'.join([cancer,stage,tag0_num,tag1_num,str(tag0_death_num),str(tag1_death_num)]) + '\t' + '\t'.join([hr, low, high, p]) + '\t' + str(fisher_p) + '\t' + str(chisquare_p) + '\n')
		out.close()
		indata.close()	
	
	
df = pd.read_csv(filepath_or_buffer="present_vs_absent.txt",
                 header='infer', sep='\t', index_col=None)
outpath = os.path.join(os.getcwd(),'absent_vs_present_univariablehrresults2')
#print outpath

if not os.path.exists(outpath):
	os.mkdir(outpath)
	
#get parameter list from file
parameters_list = []
data = open(r'superhla1.txt','r')
for line in data:
	entry = line.strip().split('\t')
	parameters_list.append(entry)
data.close()

cancer_data = open(r'cancer_type.txt','r')
cancer_data.readline()
for line in cancer_data:
	entry = line.split('\t')
	cancer = entry[0].strip()
	#select valid OS
	df_valid = df[df['OS_tag']=='valid']
	#print str(len(df_valid['OS']))	
	#select cancer
	df_cancer = df_valid[df_valid['project_id'].str.contains('TCGA-' + cancer)]
	#select stage
	df_early = df_cancer[df_cancer['tumor_stage_type'].str.contains('early')]
	df_advanced = df_cancer[df_cancer['tumor_stage_type'].str.contains('advanced')]	
	for temp in parameters_list:
		#print temp[0]
		df_cancer2 = df_cancer[df_cancer[temp[0].strip()].isin(['0','1'])]
		df_early2 = df_early[df_early[temp[0].strip()].isin(['0','1'])]
		df_advanced2 = df_advanced[df_advanced[temp[0].strip()].isin(['0','1'])]
		#print str(len(df_cancer2[temp[0].strip()]))
		death_all = df_cancer2[df_cancer2['OS_event'].str.contains('1')]
		death_early = df_early2[df_early2['OS_event'].str.contains('1')]
		death_advanced = df_advanced2[df_advanced2['OS_event'].str.contains('1')]
		#live_all = df_cancer2[df_cancer2['OS_event'].str.contains('0')]
		#live_early = df_early2[df_early2['OS_event'].str.contains('0')]
		#live_advanced = df_advanced2[df_advanced2['OS_event'].str.contains('0')]		
		#death number in Tag0 group
		#death_all_tag0 = death_all[death_all[temp[0].strip()].isin(['0'])]
		#death_all_tag1 = death_all[death_all[temp[0].strip()].isin(['1'])]
		#death_early_tag0 = death_early[death_early[temp[0].strip()].isin(['0'])]
		#death_early_tag1 = death_early[death_early[temp[0].strip()].isin(['1'])]		
		#death_advanced_tag0 = death_advanced[death_advanced[temp[0].strip()].isin(['0'])]
		#death_advanced_tag1 = death_advanced[death_advanced[temp[0].strip()].isin(['1'])]		
		
		if len(death_all['OS']) >=3:
			if len(df_cancer2[temp[0].strip()]) >= 1:
				univariableCOXph(df_cancer2,'OS_event','OS',temp[0].strip(),'all',cancer,outpath)
		if len(death_early['OS']) >=3:
			print(cancer + ' ' + temp[0].strip() + ' death num. ' + str(len(death_early['OS'])))
			if len(df_early2[temp[0].strip()]) >= 1:
				univariableCOXph(df_early2,'OS_event','OS',temp[0].strip(),'early',cancer,outpath)		
		if len(death_advanced['OS']) >=3:
			if len(df_advanced2[temp[0].strip()]) >= 1:		
				univariableCOXph(df_advanced2,'OS_event','OS',temp[0].strip(),'advanced',cancer,outpath)
		
cancer_data.close()
