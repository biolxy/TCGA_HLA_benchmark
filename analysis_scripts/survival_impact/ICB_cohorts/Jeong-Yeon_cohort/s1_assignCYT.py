data=open(r'CYT_gene.txt','r')
data.readline()
data.readline()
dicts = {}
for line in data:
	entry = line.split('\t')
	pt = entry[0].strip()
	gzma = entry[2].strip()
	prf1 = entry[-2].strip()
	cyt = entry[-1].strip()
	dicts[pt] = [gzma, prf1, cyt]
data.close()

indata = open(r'GSE135222_infor.txt','r')
out = open(r'v2.GSE135222_infor.txt','w')
out.write(indata.readline().strip() + '\t' +'GZMA' +'\t' + 'PRF1' + '\t' + 'CYT' +'\n')
for line in indata:
	entry = line.split('\t')
	if entry[0].strip() in dicts.keys():
		out.write(line.strip() + '\t' + '\t'.join(dicts[entry[0].strip()]) +'\n')
	else:
		out.write(line.strip() + '\t' +'\n')
out.close()
indata.close()
