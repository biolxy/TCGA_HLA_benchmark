#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : get_top2_table.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2020-07-19 22:17:58
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import sys


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


def get_top2_table(infile, inlist):
    # softwareName
    name = infile.split(".")[1]
    listo = []
    ff = open(infile, 'r')
    lines = ff.readlines()
    ff.close()
    for line in lines[1:]:
        line = line.strip()
        if line.startswith("\t"):
            headline = line
            continue
        tag = line.split("\t")
        allele = tag[0]
        error_num_total, allele_all_num, error_ratio = int(
            tag[-3]), int(tag[-2]), float(tag[-1])
        other_list = [ int(x) for x in tag[1:-3]]
        other_list.sort()
        other = sum(other_list)
        tmp = [allele_all_num, ]
    return listo, indict, name


def get_allele_list(infile):
    list_1 = []
    with open(infile, 'r') as ff:
        for line in ff:
            line = line.strip()
            tag = line.split("\t")
            if tag[0] != "":
                list_1.append(tag[0])
    return list_1


if __name__ == "__main__":
    tableS3 = 'class2.table3.txt'
    allele_list = get_allele_list(tableS3)
