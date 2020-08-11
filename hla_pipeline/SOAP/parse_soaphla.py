import os
import re
import sys
file_path = sys.argv[1]
dicts = {}
out = open(r'soaphla.result.xls', 'w')
hla_list = ['HLA-A', 'HLA-B', 'HLA-C']
print >> out, 'Sample' + '\t' + '\t'.join(hla_list)
for item in os.listdir(file_path):
    data_path = os.path.join(file_path, item)
    dicts[item] = [[], [], []]
    data = open(os.path.join(data_path, item + r'.type'), 'r')
    # data.readline()
    for line in data:
        entry = line.split('\t')
        a = re.match(r'A', entry[0].strip())
        if a:
            if entry[0] not in dicts[item][0]:
                dicts[item][0].append(entry[0])
        b = re.match(r'B', entry[0].strip())
        if b:
            if entry[0] not in dicts[item][1]:
                dicts[item][1].append(entry[0])
        c = re.match(r'C', entry[0].strip())
        if c:
            if entry[0] not in dicts[item][2]:
                dicts[item][2].append(entry[0])
    data.close()
for item in dicts.keys():
    print >> out, item + '\t' + \
        '\t'.join([', '.join(dicts[item][0]), ', '.join(
            dicts[item][1]), ', '.join(dicts[item][2])])
out.close()
