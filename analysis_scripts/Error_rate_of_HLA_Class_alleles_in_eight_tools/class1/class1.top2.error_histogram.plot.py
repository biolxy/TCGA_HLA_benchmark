#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : plot.py .
Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-05-25 19:59:25
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
import sys
# rcParams.update({'font.family': 'serif'})
dict1 = { 
 'simulation': "combining results",
 'polysolver': "POLYSOLVER",
 'OptiType': "OptiType" ,
 'xHLA': 'xHLA',
 'hlagenotyper': "hla-genotyper",
 'hlahd': "HLA-HD",
 'kourami': "Kourami",
 'soaphla': "SOAP-HLA",
 'HLA-VBSeq': "HLA-VBSeq"
}

def setlist(x, indict):
    for i in indict:
        if x.startswith(i):
            x = x.replace(i,indict[i])
            break
    return x

# 设置 图片尺寸等
params = {'legend.fontsize': 10,
          'figure.figsize': (10, 8),
         'axes.labelsize': 10,
         'axes.titlesize':10,
         'xtick.labelsize':10,
         'ytick.labelsize':10}
plt.rcParams.update(params)


file1 = sys.argv[1] # "/Users/sheng/Desktop/yeb/allele_error_num_and_top_error/class2.top3.error2.txt"
#设置子图的大小
#导入数据集car crash dataset
data1 = pd.read_csv(file1, sep = "\t")
df = pd.DataFrame(data1 , columns=pd.Index(['soft_allele', 'top1_ratio', 'top2_ratio', 'other_ratio']))
df = df.rename(columns={'top1_ratio': 'top1', 'top2_ratio': 'top2', 'other_ratio':'other'})
df['soft_allele'] = df.apply(lambda x: setlist(x.soft_allele, dict1), axis=1)
print(df)
list1 = df.soft_allele.tolist()
print(list1)


df.index = list1
df2 = data1[['allerrorNum' ,'str3', 'str4']]
df2 = df2.rename(columns={'allerrorNum': 'Error fre.', 'str3': 'top1', 'str4': 'top2'})
print(df2.values)
df2.index = list1
fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(3,3,5)
df.plot(ax =ax, kind='barh', stacked=True, x = 'soft_allele', alpha=0.5, figsize = (15, 15))
ax.set(ylabel = "Tool True allel", xlabel="specific error rate")
ax.legend(['top1', 'top2', 'other'])
list2 = np.array(df2.values)
list2 = list2[::-1]
print(list2)
the_table = ax.table(cellText=list2, colLabels= df2.columns.tolist() , cellLoc ='center', loc='right', edges ='open') # ,  rowLabels=df2.index.tolist()
the_table.auto_set_font_size(False)
the_table.set_fontsize(10)
ax.autoscale()
# ax = fig.add_subplot(2,2,1+1)
# df.plot(ax =ax, kind='bar', x = 'soft_allele', alpha=0.5, figsize = (10, 10))


fig = ax.get_figure()
fig.savefig( str(file1) + '.pdf', dpi = 300)
fig.savefig( str(file1) + '.png', dpi = 300)
plt.close()