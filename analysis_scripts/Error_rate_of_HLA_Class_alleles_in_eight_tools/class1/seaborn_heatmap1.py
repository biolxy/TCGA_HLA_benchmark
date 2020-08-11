# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib import rcParams
import sys
rcParams.update({'font.size': 8})
# rcParams.update({'font.size': 9, 'font.family': 'serif'})


#load data matrix through pandas
df = pd.read_csv(filepath_or_buffer="class1.table3.txt",header='infer',sep='\t',index_col=0)
flights2 = df[[u'POLYSOLVER', u'OptiType', u'xHLA', u'hla-genotyper', u'HLA-HD', u'Kourami', u'SOAP-HLA', u'HLA-VBSeq']]
f, ax = plt.subplots(figsize=(5, 11))
sns.heatmap(flights2, cmap="YlGnBu", linewidths=.3, ax=ax)
ax.text(0.3, 1.02, 'General Error rate', transform=ax.transAxes, fontsize=9,verticalalignment='top')
plt.xticks(rotation=45,horizontalalignment="right")
plt.ylabel('HLA Class I Allele (Frequency)', fontsize = 12) # y-axis label with  fontsize = 12
f.savefig("class1.error.ratio.heatmap.pdf", bbox_inches='tight')
f.savefig("class1.error.ratio.heatmap.png", bbox_inches='tight')