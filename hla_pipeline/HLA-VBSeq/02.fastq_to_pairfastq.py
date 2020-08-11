#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : 02.fastq_to_pairfastq.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2018-10-15 21:26:16
version     : 1.0
Function    : The author is too lazy to write nothing
"""
import sys
import os
from collections import Counter


outputdir = os.path.dirname(os.path.abspath(sys.argv[1]))
inputfqname = os.path.basename(sys.argv[1])
prefix = inputfqname.replace("_1.fastq", "")
outputfq1 = os.path.join(outputdir, prefix + "_pair_1.fastq")
outputfq2 = os.path.join(outputdir, prefix + "_pair_2.fastq")


def dict_fastq(pathfastq):
    dict = {}
    numberline = 0
    with open(pathfastq, 'r') as f:
        for line in f:
            line = line.rstrip()
            if numberline % 4 == 0:
                key = line
                value = ""
            else:
                value = value + "|" + line
            dict[key] = value
            numberline += 1
    return dict


def strread(number, str):
    lista = []
    for i in range(0, number):
        lista.append(str)
    strread = "".join(lista)
    return strread


dictfq1 = dict_fastq(sys.argv[1])
dictfq2 = dict_fastq(sys.argv[2])
# print dict1
# for key in dictfq1:
#     print key, dictfq2[key]
listfq1 = [x[0:-2] for x in dictfq1]
listfq2 = [x[0:-2] for x in dictfq2]
listfq_12 = listfq1 + listfq2
dict_fa_12 = dict(Counter(listfq_12))
# print dict_fa_12

list = []
for i in dict_fa_12:
    if dict_fa_12[i] == 2:
        list.append(i)
# print list

with open(outputfq1, 'w') as fq1:
    with open(outputfq2, 'w') as fq2:
        for i in list:
            fq1line1 = i + "/1"
            fq1_list = dictfq1[fq1line1].split("|")
            fq1line2 = fq1_list[1]
            fq1line3 = fq1_list[2]
            fq1line4 = fq1_list[3]
            fq1.write("{}\n{}\n{}\n{}\n".format(
                fq1line1, fq1line2, fq1line3, fq1line4))
            fq2line1 = i + '/2'
            fq2_list = dictfq2[fq2line1].split("|")
            fq2line2 = fq2_list[1]
            fq2line3 = fq2_list[2]
            fq2line4 = fq2_list[3]
            fq2.write("{}\n{}\n{}\n{}\n".format(
                fq2line1, fq2line2, fq2line3, fq2line4))
list2 = []
for i in dict_fa_12:
    if dict_fa_12[i] == 1:
        list2.append(i)

with open(outputfq1, 'aw') as fq1:
    with open(outputfq2, 'aw') as fq2:
        for i in list2:
            fq1line1 = i + "/1"
            fq2line1 = i + "/2"
            if fq1line1 in dictfq1:
                fq1_list = dictfq1[fq1line1].split("|")
                fq1line2 = fq1_list[1]
                fq1line3 = fq1_list[2]
                fq1line4 = fq1_list[3]
                fq1.write("{}\n{}\n{}\n{}\n".format(
                    fq1line1, fq1line2, fq1line3, fq1line4))
                fq2line2 = strread(len(fq1line2), 'N')
                fq2.write("{}\n{}\n{}\n{}\n".format(
                    fq2line1, fq1line2, fq1line3, fq1line4))
            if fq2line1 in dictfq2:
                fq2_list = dictfq2[fq2line1].split("|")
                fq2line2 = fq2_list[1]
                fq2line3 = fq2_list[2]
                fq2line4 = fq2_list[3]
                fq2.write("{}\n{}\n{}\n{}\n".format(
                    fq2line1, fq2line2, fq2line3, fq2line4))
                fq1line2 = strread(len(fq2line2), 'N')
                fq1.write("{}\n{}\n{}\n{}\n".format(
                    fq1line1, fq1line2, fq2line3, fq2line4))
