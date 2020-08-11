#!/usr/bin/env python

"""
HLA genotyping through HLAminer.

#usage: python hlagenotyper.py -w sample.conf -s /home/yeh/workspace/hengrui/2018-04-11/bam -d /home/yeh/workspace/hengrui/2018-04-11/HLA/hlagenotyper
Created by Hao Ye on Jan 29th 2018
"""

import argparse
import sys
import os
import re
import time
# import multiprocessing

data = '/storage1/TCGA/TCGA_HLA/TCGA_DATA/linkbam/5.hlagenotyper.race/bamid_2_race.d'
# Caucasian, Black, Asian or Unknown


def dict_race(file):
    dict_race = {}
    with open(data, 'r') as f:
        for line in f:
            line = line.rstrip()
            list = line.split("\t")
            sampleid = list[0]
            race = list[1]
            dict_race[sampleid] = race
    return dict_race


def get_race_fromdict(sampleid, dict):
    race = "UNK"
    if sampleid in dict:
        race = dict[sampleid]
    return race


dict_race = dict_race(data)


def loadsampleid(infile):
    li = []
    with open(infile) as fin:
        for line in fin:
            if re.findall(r'^#', line):
                continue
            li.append(line.rstrip())
    fin.close()
    return li


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--sampleconf', action='store', dest='sampleidConfFile',
                        help='REQUIRED, WGC ID config file', required=True)
    parser.add_argument('-s', '--source', action='store', dest='src',
                        help='REQUIRED, fastq path, no "/" at the end', required=True)
    parser.add_argument('-d', '--destination', action='store', dest='dest', default='.',
                        help='destination path, default is current path, no "/" at the end')
    args = parser.parse_args()
    # get sampleid
    sampleids = loadsampleid(args.sampleidConfFile)
    # pool=multiprocessing.Pool(processes=10)
    for tag in sampleids:
        # sample = tag.split('-')[1].strip()
        sample = tag.strip()
        bam = os.path.join(args.src, sample + r'.bam')
        out_dir = os.path.join(args.dest, sample)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        logfile = os.path.join(out_dir, sample + ".log")
        fLog = open(logfile, 'a')
        print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(
            time.time()))+'] python '+' '.join(sys.argv)
        # get race 
        id = tag.split("_")[0]
        print id
        race = dict_race[id]
        unmapbam = os.path.join(out_dir, sample + r'.unmapped.bam')
        command = '/mnt/pipeline-programs/samtools/samtools-1.3.1/samtools view -u -f 4 ' + \
            bam + ' -o ' + unmapbam
        command2 = "hla-genotyper  {bam} -u  {unmapbam} -e {race} -r 37 --exome -o {out_dir}".format(bam=bam, unmapbam=unmapbam, race=race, out_dir=out_dir)
        command3 = '/usr/bin/rm -f ' + unmapbam
        print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(time.time()))+'] ' + command
        os.system(command)
        print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(time.time()))+'] ' + command
        print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(time.time()))+'] ' + command2
        os.system(command2)
        print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(time.time()))+'] ' + command2
        os.system(command3)
        print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(time.time()))+'] ' + command3
    print >> fLog, '['+time.strftime("%Y-%m-%d %H:%M:%S",
                                     time.localtime(time.time()))+'] finished'
    fLog.close()
