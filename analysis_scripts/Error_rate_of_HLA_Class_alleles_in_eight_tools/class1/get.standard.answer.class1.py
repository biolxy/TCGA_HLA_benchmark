#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : 2.get.standard.answer.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2018-11-23 14:27:16
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import sys
import os
import json
from collections import Counter


class MagicDict(dict):
    """Implementation of perl's autovivification feature."""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value



# list_key = ['DRB1-1', 'DRB1-2', 'DPB1-1', 'DPB1-2', 'DQA1-1', 'DQA1-2', 'DQB1-1', 'DQB1-2']
list_key = ['A-1', 'A-2', 'B-1', 'B-2', 'C-1', 'C-2']

# def get_highest_frequency_item(list_i):
#     # return list
#     # [('DRB1*09:01', 3), ('DRB1*09:19', 1), ('None', 0)]
#     hla_count = Counter(list_i)
#     top_list = hla_count.most_common(len(list_i))
#     # remove none
#     if top_list[0][0] == 'None' and top_list[0][1] != len(list_i):
#         top_list.remove(top_list[0])
#     #
#     list_hla = []
#     list_hla.append(top_list[0][0])
#     for i in range(1, len(top_list)):
#         if top_list[0][1] == top_list[i][1] and top_list[i][0] != 'None':
#             list_hla.append(top_list[i][0])
#     str_return = "|".join(list_hla) + "," + str(top_list[0][1])
#     return str_return


def get_highest_frequency_item(list_i):
    # return list
    hla_count = Counter(list_i)
    top_list = hla_count.most_common()
    # top_list = [('DRB1*09:01', 3), ('DRB1*09:19', 1), ('None', 0)]

    # remove none
    if top_list[0][0] == 'None' and top_list[0][1] != len(list_i):
        top_list.remove(top_list[0])

    list_hla = []
    list_hla.append(top_list[0][0])

    if int(top_list[0][1]) < 2:
        str_return = "None" + "," + str(top_list[0][1])
    else:
        for i in range(1, len(top_list)):
            if  top_list[i][1] == top_list[0][1] and top_list[i][0] != 'None':
                list_hla.append(top_list[i][0])
        str_return = "|".join(list_hla) + "," + str(top_list[0][1])
    return str_return


print(sys.argv)
alldict = MagicDict()
for i in sys.argv[1:]:
    # print i
    with open(i, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("A1\t"):
                continue
            linelist = line.split("\t")
            sample = linelist[0]
            for i in range(1, 7):
                typehla = list_key[i-1]
                hla = linelist[i]
                # print("{} {} {}".format(sample, typehla, hla))
                if not alldict[sample][typehla]:
                    alldict[sample][typehla] = []
                alldict[sample][typehla].append(hla)
                 # print alldict[sample][typehla], hla

# a = json.dumps(alldict, indent=4, separators=(',', ':'))
with open("result.json", 'w') as ff:
    json.dump(alldict, ff, indent=4, separators=(',', ':'))


out = open("result.t", 'w')
for key in alldict:
    linelist = []
    linelist.append(key)
    for item in list_key:
        hlatype = get_highest_frequency_item(alldict[key][item])
        linelist.append(hlatype)
    line = "\t".join(linelist)
    out.write("{}\n".format(line))
    print("{}".format(line))
out.close()
