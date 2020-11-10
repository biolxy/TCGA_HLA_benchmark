#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : get_simulation_result_by_three_software.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-06-09 10:55:33
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""

import sys, os
import numpy as np
from collections import Counter

class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def get_highest_frequency_item(indict):
    # return list
    software = ['polysolver', 'OptiType', 'xHLA']
    listtmp = []
    for item in software:
        if not indict[item]:
            indict[item] = ['None', 'None', 'None', 'None', 'None', 'None']
        listtmp.append(indict[item])
    a = np.array(listtmp, dtype=str)
    return_list = []
    for i in range(0,6):
        listtmp2 = a[:,i] # 第0列
        hla_count = Counter(listtmp2)
        top_list = hla_count.most_common()
        # print(top_list) [('C17:01', 2), ('C17:03', 1)]
        if int(top_list[0][1]) < 2:
            str_hla = a[:,0][0] # 如果没有没有cover的hla, 则取polysolver的结果作为答案
        else:
            str_hla = top_list[0][0] # C17:01
        return_list.append(str_hla)
    return return_list


def addsoftwareResult(dictResult, infile, softwareName):
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip()
            linelist = line.split("\t")
            sample = linelist[0]
            if not dictResult[sample][softwareName]:
                dictResult[sample][softwareName] = []
                dictResult[sample][softwareName].extend(linelist[1:])
            else:
                print(sample, softwareName)
    return dictResult


outfile = sys.argv[-1]

dictResult = MagicDict()
list1 = ['A-1', 'A-2', 'B-1', 'B-2', 'C-1', 'C-2']
software = ['polysolver', 'OptiType', 'xHLA']

addsoftwareResult(dictResult, '6.polysolver.o.standard.class1', 'polysolver')
addsoftwareResult(dictResult, '5.OptiType.o.standard.class1', 'OptiType')
addsoftwareResult(dictResult, '8.xHLA.o.standard.class1', 'xHLA')


out = open(outfile, 'w')
for sampleid in dictResult:
    linelist = []
    linelist.append(sampleid)
    try:
        hlatypelist = get_highest_frequency_item(dictResult[sampleid])
    except:
        print(dictResult[sampleid], sampleid)
        break
    linelist.extend(hlatypelist)
    line = "\t".join(linelist)
    out.write("{}\n".format(line))
    print("{}".format(line))
out.close()