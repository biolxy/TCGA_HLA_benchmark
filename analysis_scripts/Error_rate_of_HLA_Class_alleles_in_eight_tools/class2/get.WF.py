#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : get.WF.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-05-23 16:35:57
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


def get_soft_result(inputfile, id_d):
    dict1 = MagicDict()
    list1 = []
    with open(inputfile, 'r') as ff:
        for line in ff:
            line = line.rstrip()
            linelist = line.split("\t")
            sampleid = linelist[0]
            if sampleid in id_d:
                list1.append(sampleid)
                allelelist = []
                for i in linelist[1:]:
                    allele = i.split(",")[0]
                    allelelist.append(allele)
                dict1[sampleid] = allelelist
    return dict1, list1


def compareAlleltList(list1, list2, dict1):
    # list1 is result
    for x, item in enumerate(list1):
        if item != list2[x]:
            error_allele = list2[x]
            if not dict1[item][error_allele]:
                dict1[item][error_allele] = 0
            dict1[item][error_allele] += 1
    return dict1


def get_first_key(dict1):
    list1 = []
    for key in dict1.keys():
        list1.append(key)
    list1 = list(set(list1))
    list1.sort()
    return list1


def get_second_key(dict1):
    list1 = get_first_key(dict1)
    list2 = []
    for item in list1:
        dict2 = dict1[item]
        list3 = get_first_key(dict2)
        list2 = list2 + list3
    list2 = list(set(list2))
    list2.sort()
    return list2



def getIDlist(inIDfile):
    """ dna_normal_bam_for_hla.cl.id
    """
    id_d = {}
    with open(inIDfile, 'r') as ff:
        for line in ff:
            line = line.strip()
            id_d[line] = 1
    return id_d

# 获取 所有 10479 个 sampleID
inIDfile='dna_normal_bam_for_hla.cl.id'
id_d = getIDlist(inIDfile)

# 获得标准答案的dict
inputresultfile = 'result.t'
# 091477b4-9f99-48d0-84b3-aaaaaa362246	DRB111:01,5	DRB115:03,6	DPA101:03,3	DPA101:03,2	DPB102:01,4	DPB118:01,3	DQA101:02,5	DQA101:02,4	DQB106:02,5	DQB106:02,5
dict_result = MagicDict()
dict_result_alleleNum = MagicDict()
with open(inputresultfile, 'r') as ff:
    for line in ff:
        line = line.rstrip()
        linelist = line.split("\t")
        sampleid = linelist[0]
        allelelist = []
        if sampleid in id_d:
            for i in linelist[1:]:
                allele = i.split(",")[0]
                num = int(i.split(",")[1])
                if num >= 2:
                    allelelist.append(allele)
                else:
                    allelelist.append("None")
                allelelist2 = allele.split("|")
                for allele2 in allelelist2:
                    if not dict_result_alleleNum[allele2]:
                        dict_result_alleleNum[allele2] = 0
                    dict_result_alleleNum[allele2] += 1
            dict_result[sampleid] = allelelist


for infile in sys.argv[1:]:
    print(infile)
    soft_result, sampleidlist = get_soft_result(infile, id_d)
    soft_error_dict = MagicDict()
    for sampleid in sampleidlist:
        if sampleid in dict_result:
            result_alleltlist = dict_result[sampleid]
            soft_alleltlist = soft_result[sampleid]
            # print(soft_alleltlist)
            # print(result_alleltlist)
            soft_error_dict = compareAlleltList(result_alleltlist, soft_alleltlist, soft_error_dict)
        else:
            pass
            # print("{} not in result".format(sampleid))

    listline = get_first_key(soft_error_dict)
    if "None" in listline:
        listline.remove("None")
    listcolumn = get_second_key(soft_error_dict)
    if "None" in listcolumn:
        listcolumn.remove("None")
    print(len(listline))
    print(len(listcolumn))
    outfile = str(infile) + ".right_error.xls"
    out = open(outfile, 'w')
    out.write("{}\t{}\t{}\t{}\t{}\n".format("", "\t".join(listcolumn), "error_num_total", 'allele_all_num', "error_ratio"))
    for allele in listline:
        line_list = []
        error_num_total = 0
        for error_allele in listcolumn:
            if not soft_error_dict[allele][error_allele]:
                errorNum = 0
            else:
                errorNum = soft_error_dict[allele][error_allele]
            line_list.append(str(errorNum))
            error_num_total += errorNum
        line1 = "\t".join(line_list)
        # 统计总数
        if allele in dict_result_alleleNum:
            allele_all_num = dict_result_alleleNum[allele]
            error_ratio = error_num_total * 1.0 / int(allele_all_num)
            if error_num_total != 0:
                out.write("{}\t{}\t{}\t{}\t{}\n".format(allele, line1, error_num_total, allele_all_num, error_ratio))
    out.close()
