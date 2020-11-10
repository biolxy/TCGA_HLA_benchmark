#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : get_result2_by_bamuuid.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2018-12-25 18:02:33
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import sys


def get_str_from_list(inputlist):
    outputlist = []
    for item in inputlist:
        item2 = item.split(",")[0]
        outputlist.append(item2)
    return outputlist


# def get_score_by2list(listA, standardlist):
#     result_score = 0.0
#     result_score2 = 0.0
#     result_score3 = 0.0
#     result_score4 = 0.0
#     # rights, wrongs, nones, faileds
#     None_list = []
#     for x, item in enumerate(standardlist):
#         if listA[x] == "None":
#             result_score3 += 1
#             None_list.append(listA[x])
#         else:
#             if listA[x] == item:
#                 result_score += 1
#             else:
#                 if "|" in item:
#                     listitem = item.split("|")
#                     if listA[x] in listitem:
#                         result_score += 1
#                     else:
#                         result_score2 += 1
#                 else:
#                     result_score2 += 1
#     if len(None_list) == len(listA):
#         result_score4 = 1.0
#     return result_score, result_score2, result_score3, result_score4


def get_score_by2list(listA, standardlist):
    result_score = 0.0
    result_score2 = 0.0
    result_score3 = 0.0
    result_score4 = 0.0
    # rights, wrongs, singlenones, six_none
    None_list = []
    for x, item in enumerate(standardlist):
        if item != "None":
            # 如果标准答案非 None，只有三种情况，相等，wrongs or singlenones
            if listA[x] == item:
                # 相等
                result_score += 1
            else:
                if listA[x] == "None":
                    # singlenones
                    result_score3 += 1
                else:
                    # wrongs
                    result_score2 += 1
        else:
            # 如果标准答案是 None，只有两种情况，wrongs or singlenones
            if listA[x] == "None":
                # singlenones
                result_score3 += 1
            else:
                # wrongs
                result_score2 += 1


    if len(None_list) == len(listA):
        result_score4 = 1.0
    return result_score, result_score2, result_score3, result_score4
    # rights, wrongs, singlenones, allnone


def get_dict_from_file(file):
    dict1 = {}
    with open(file, 'r') as ff:
        for line in ff:
            line = line.rstrip()
            linelist = line.split("\t")
            key = linelist[0]
            values = get_str_from_list(linelist[1:])
            dict1[key] = values
    return dict1


def get_dict_from_file(file):
    dict1 = {}
    with open(file, 'r') as ff:
        for line in ff:
            line = line.rstrip()
            linelist = line.split("\t")
            key = linelist[0]
            values = get_str_from_list(linelist[1:])
            dict1[key] = values
    return dict1


def  get_newlist_byNum(inlist):
    list1 = []
    for i in range(0, len(inlist)):
        list1.append(inlist[i])
    return list1


standard_answer = 'result.t'
standard_result_dict = get_dict_from_file(standard_answer)

print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("softwareName", "rightAlleles_score","wrongAlleles_score", "uncalledAlleles_num", "Scuess", "Accuracy", "allnone_num", "Failed"))

list_sample_name = []
for key in standard_result_dict:
    list_sample_name.append(key)


list_10483 = []
with open('dna_normal_bam_for_hla.cl.id', 'r') as ff:
    for line in ff:
        line = line.rstrip()
        list_10483.append(line)

# 经过修改后 dna_normal_bam_for_hla.cl.id 共 10479 行，class1 2 都是如此


for file in sys.argv[1:]:
    rightAlleles_score = 0.0
    wrongAlleles_score = 0.0
    none_num = 0.0
    allnone_num = 0.0
    fail_num = 0.0
    tmpdict = get_dict_from_file(file)
    for sample in list_sample_name:
        if sample in list_10483:
            if not sample in tmpdict:
                fail_num += 1.0
                print(sample)
            else:
                standardList = standard_result_dict[sample]
                candidateList = tmpdict[sample]
                # 此处比较 软件结果和标准答案，并给出比较结果
                standardList = get_newlist_byNum(standardList)
                candidateList = get_newlist_byNum(candidateList)
                # print(standardList)
                # print(candidateList)
                rights, wrongs, singlenones, allnone = get_score_by2list(candidateList, standardList)
                rightAlleles_score += rights
                wrongAlleles_score += wrongs
                none_num += singlenones
                allnone_num += allnone

    softwareName = str(file).split(".")[1]
    if rightAlleles_score + wrongAlleles_score != 0.0:
        Scuess = rightAlleles_score / \
            (rightAlleles_score + wrongAlleles_score)
    else:
        Scuess = 0.0
    if rightAlleles_score + wrongAlleles_score + none_num != 0.0:
        Accuracy = rightAlleles_score / \
            (rightAlleles_score + wrongAlleles_score + none_num)
    else:
        Accuracy = 0.0
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(softwareName, rightAlleles_score, wrongAlleles_score, none_num, Scuess, Accuracy, allnone_num, fail_num))
