import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sklearn.metrics as skms
from matplotlib import rcParams
#rcParams.update({'font.size': 28, 'font.family': 'serif'})
sns.set(font_scale=4)
sns.set(style="white")
colors = ["g","b","y","purple"]
# B-,C, C-, D
color_palette = sns.color_palette(colors)
#iris = sns.load_dataset("iris")
#iris.to_csv(path_or_buf="D:\project\cooperation\zhoujian\huangbo\iris_demo.txt",sep='\t',header=True,index=False)
demodata = pd.read_csv(filepath_or_buffer="log_and_fc.txt",header='infer',sep='\t',index_col=None)

#ax = sns.jointplot("FFPE_gDNA_BAF", "NGS_lib_BAF", data=demodata,kind="reg", truncate=False,xlim=(0, 1), ylim=(0, 1),color="m", height=7)
#plt.figure(figsize=(14, 12), dpi=300)	
#ax = plt.subplots()
#ax = sns.kdeplot(demodata['FFPE_gDNA_BAF'], demodata['NGS_lib_BAF'], cmap="Reds", shade=True, shade_lowest=False)
ax = sns.kdeplot(demodata['FC (mean IC50_B44/IC50_nonB44)'], demodata['-log10(p)'], cmap="Reds", shade=True, shade_lowest=False)
					  
#ax = sns.jointplot("FFPE_gDNA_BAF", "NGS_lib_BAF", data=demodata, kind="reg", truncate=False,xlim=(0, 1), ylim=(0, 1),color="b", height=7)
#ax = sns.jointplot(x="Expected_AF", y="P81_AF", data=demodata, kind="reg", truncate=False,xlim=(0, 0.5), ylim=(0, 0.5),color="b", height=7)
#ax.scatter(demodata['FFPE_gDNA_BAF'], demodata['NGS_lib_BAF'],s=4,c='red',marker='o',edgecolors='face')
	
#ax = sns.scatterplot(x="mean IC50_B44/IC50_nonB44", y="-log10(p)", hue = 'DNA_Quality', size_norm = (100,100),palette = color_palette, data=demodata)
ax = sns.scatterplot(x="FC (mean IC50_B44/IC50_nonB44)", y="-log10(p)",data=demodata)

#ax.scatter(demodata['FFPE_gDNA_BAF'], demodata['NGS_lib_BAF'],s=4,c='red',marker='o',edgecolors='face')
	
#ax = sns.kdeplot(cluster1.Pt1_Primary_cancer_cell_fraction, cluster1.Pt1_Relapse_cancer_cell_fraction, cmap="Reds",

#r2 = skms.r2_score(demodata['Expected_AF'], demodata['P81_AF'])
#print str(r2)
#textstr = r'R square = %.2f' % (r2)
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.681)
#ax.text(0.05, 0.75, textstr, transform=ax.transAxes, fontsize=18,verticalalignment='top', bbox=props)
ax.set_xlim(0.72,1.25)
ax.set_ylim(0,17)	
#plt.plot([0.72,1.3],[1.25,1.3],linewidth='2',color='black')
plt.hlines(1.3, 0.72, 1.25)
plt.vlines(1, 0, 17)
#plt.plot([0,0.45],[0,0.45],linewidth='4',color='purple')
#plt.show()
plt.savefig('mean_IC50_B44_SKCM.png')

