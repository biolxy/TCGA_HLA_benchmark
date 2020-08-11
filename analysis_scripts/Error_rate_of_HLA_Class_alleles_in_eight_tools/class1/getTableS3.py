#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : getTableS3.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2020-07-12 21:35:05
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import sys
import os

# error_num_total	allele_all_num	error_ratio

# 1.hlagenotyper.o.standard.class1.right_error.xls


class MagicDict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

    def walk(self):
        for key, value in self.items():
            if isinstance(value, MagicDict):
                for tup in value.walk():
                    yield (key,) + tup
            else:
                yield key, value


def get_ratio(infile, indict):
    # softwareName
    name = infile.split(".")[1]
    listo = []
    ff = open(infile, 'r')
    lines = ff.readlines()
    ff.close()
    for line in lines[1:]:
        line = line.strip()
        tag = line.split("\t")
        allele = tag[0]
        error_num_total, allele_all_num, error_ratio = int(
            tag[-3]), int(tag[-2]), float(tag[-1])
        indict[name][allele]['error_ratio'] = str(error_ratio)
        indict[name][allele]['allele_all_num'] = allele_all_num
        if error_num_total >= ERROR_ALLELE_NUM_THRESHOLD and error_ratio >= 0.05:
            listo.append(allele)
    return listo, indict, name


def get_allsoft_error_ratio(top_allele_error_ratio, softwareName, allele):
    list_1 = []
    for name in softwareName:
        if top_allele_error_ratio[name][allele] != {}:
            list_1.append(top_allele_error_ratio[name][allele]['error_ratio'])
        else:
            list_1.append('0')
    # allele 总的出现次数
    allele_all_num = top_allele_error_ratio[name][allele]['allele_all_num']
    return "\t".join(list_1), allele_all_num


if __name__ == "__main__":
    """
    polysolver    POLYSOLVER
    OptiType    OptiType
    xHLA    xHLA
    hlagenotyper    hla-genotyper
    hlahd    HLA-HD
    kourami    Kourami
    soaphla    SOAP-HLA
    HLA-VBSeq    HLA-VBSeq

    10848  rename "polysolver" "POLYSOLVER" * 
    10851  rename "hlagenotyper"    "hla-genotyper" *
    10852  rename "hlahd"    "HLA-HD" *
    10853  rename "kourami" "Kourami" *
    10855  rename "soaphla" "SOAP-HLA" *
    """
    ERROR_ALLELE_NUM_THRESHOLD = 200
    softwareName = []
    top_allele_list = []
    top_allele_error_ratio = MagicDict()
    for infile in sys.argv[1:]:
        listo, top_allele_error_ratio, name = get_ratio(
            infile, top_allele_error_ratio)
        top_allele_list.extend(listo)
        softwareName.append(name)
    top_allele_list = list(set(top_allele_list))
    top_allele_list.sort()
    # print(top_allele_list, len(top_allele_list))
    softwareName = [
        "POLYSOLVER",
        "OptiType",
        "xHLA",
        "hla-genotyper",
        "HLA-HD",
        "Kourami",
        "SOAP-HLA",
        "HLA-VBSeq"
    ]
    print("\t%s" % "\t".join(softwareName))
    for allele in top_allele_list:
        allsoft_error_ratio, allele_all_num = get_allsoft_error_ratio(
            top_allele_error_ratio, softwareName, allele)
        print("{}({})\t{}".format(allele, allele_all_num, allsoft_error_ratio))
