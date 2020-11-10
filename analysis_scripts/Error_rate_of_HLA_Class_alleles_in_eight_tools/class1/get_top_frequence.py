#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : get_top_frequence.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-07-06 19:18:11
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""

import json
import sys
from collections import Counter

infile = sys.argv[1] # class1.result.json
inlist = sys.argv[2] # class1.cl.10479.id
with open(infile, 'r') as ff:
    dict1 = json.load(ff)

listall = []
with open(inlist, 'r') as ff:
    for line in ff:
        line = line.rstrip()
        for i in ['A-1', 'A-2', 'B-1', 'B-2', 'C-1', 'C-2']:
            listall.extend(dict1[line][i])
            dict2 = Counter(listall)
print(len(listall))
print(dict2.most_common(20))