import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib import rcParams
import sys
rcParams.update({'font.size': 8})
# rcParams.update({'font.size': 9, 'font.family': 'serif'})
# get_new_df

#load data matrix through pandas
df = pd.read_csv(filepath_or_buffer="class2.table3.txt",header='infer',sep='\t',index_col=0)
flights = df[['HLA-HD', 'SOAP-HLA', 'xHLA', 'hla-genotyper', 'Kourami', 'HLA-VBSeq']]
f, ax = plt.subplots(figsize=(5, 11))
sns.heatmap(flights, cmap="YlGnBu",linewidths=.3, ax=ax)
ax.text(0.3, 1.02, 'General Error rate', transform=ax.transAxes, fontsize=9,verticalalignment='top')
plt.xticks(rotation=45,horizontalalignment="right")
plt.ylabel('HLA Class II Allele (Frequency)', fontsize = 12) # y-axis label with  fontsize = 12
plt.tight_layout()
f.savefig("class2.error.ratio.heatmap.pdf", bbox_inches='tight')
f.savefig("class2.error.ratio.heatmap.png", bbox_inches='tight')
