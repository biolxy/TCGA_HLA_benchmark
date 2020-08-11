#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : top3.error.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-05-23 19:27:22
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""

import sys
import json
from collections import Counter


class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def sortedDictValues2(adict):
    list1 = sorted(adict.items(), key=lambda adict:adict[1])
    return list1


def get_allele_list(infile):
    list_1 = []
    with open(infile, 'r') as ff:
        for line in ff:
            line = line.strip()
            tag = line.split("\t")
            if tag[0] != "":
                allele = tag[0].split("(")[0]
                list_1.append(allele)
    return list_1


if __name__ == "__main__":
    tableS3 = 'class2.table3.txt'
    allele_list = get_allele_list(tableS3)
    print(allele_list)
    for infile in sys.argv[1:]:
        with open(infile, 'r') as ff:
            for line in ff:
                line = line.rstrip()
                if line.startswith("\t"):
                    headline = line
                    continue
                linelist = line.split("\t")
                error_num_total = int(linelist[-3])
                allele_all_num = int(linelist[-2])
                error_ratio = float(linelist[-1])
                right_allele = linelist[0]
                # if error_num_total >= 100 and error_ratio >= 0.05:
                if right_allele in allele_list:
                    headlinelist = headline.split("\t")[1:]
                    linelistNumD = {}
                    for x,item in enumerate(linelist[1:-1]):
                        allele = headlinelist[x]
                        linelistNumD[allele] = int(item)
                    list1 = sortedDictValues2(linelistNumD)
                    listNumtemp = []
                    listAlleletemp = []
                    for i in list1[-4:]:
                        number = i[-1]
                        allele = i[0]
                        listNumtemp.append(str(number))
                        listAlleletemp.append(allele)
                    listNumtemp.reverse()
                    other = int(listNumtemp[1]) - int(listNumtemp[2]) - int(listNumtemp[3])
                    top1_ratio = int(listNumtemp[2]) * 1.0 / int(listNumtemp[1])
                    top2_ratio = int(listNumtemp[3]) * 1.0 / int(listNumtemp[1])
                    other_ratio = other * 1.0 / int(listNumtemp[1])
                    listNumtemp.append(str(other))
                    listNumtemp.append(str(top1_ratio))
                    listNumtemp.append(str(top2_ratio))
                    listNumtemp.append(str(other_ratio))
                    listAlleletemp.reverse()
                    listAlleletemp.append('other')
                    listAlleletemp.append('top1_ratio')
                    listAlleletemp.append('top2_ratio')
                    listAlleletemp.append('other_ratio')
                    sampleName_allele = infile.split(".")[1] + "_" + right_allele
                    print("{}\t{}\t{}\t{}".format(sampleName_allele, "\t".join(listNumtemp), "\t".join(listAlleletemp), error_ratio))