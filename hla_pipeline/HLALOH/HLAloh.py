#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : HLAloh.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-09-27 14:53:45
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import os
import re
import sys
import argparse
from collections import Counter
from configparser import ConfigParser
__VERSION__ = 'v1.0'


def mkoutdir(indir):
    if not os.path.isdir(indir):
        os.mkdir(indir)


def buildname2col(line):
    line = line.rstrip()
    colnames = line.split("\t")
    colname2colnum = {}
    colnum2colname = {}
    for i, colname in enumerate(colnames):
        colname2colnum[colname] = i
        colnum2colname[i] = colname
    return colname2colnum, colnum2colname


def getHeterozygousList(inlist):
    # inlist = ['A*11:01', 'B*38:01', 'C*07:02', 'A*01:02']
    # If A/B/C appears only once, the modified alleles is homozygous
    ABC_list = [x.split("*")[0] for x in inlist]
    dict_alleles = dict(Counter(ABC_list))
    heterozygous_ABC = []
    for item in dict_alleles:
        if dict_alleles[item] == 2:
            heterozygous_ABC.append(item)
    heterozygous = []
    for alleles in inlist:
        if alleles.split("*")[0] in heterozygous_ABC:
            heterozygous.append(alleles)
        else:
            print("# {} is homozygous, lohhla can not support".format(alleles))
    return heterozygous


def interceptionChr6(SCTIPR_FOLDER, samtools, inBam, coordinates, outBamName):
    # hg19_coordinates='chr6:29580000-33618227'
    # hg38_coordinates='chr6:29602228-33650450'
    command = "bash {SCTIPR_FOLDER}/interceptionChr6.sh {samtools} {inBam} {coordinates} {outBamName}".format(
        SCTIPR_FOLDER=SCTIPR_FOLDER,
        samtools=samtools,
        inBam=inBam,
        coordinates=coordinates,
        outBamName=outBamName
    )
    os.system(command)


def main():
    # =================
    # creat hlas file
    # =================
    hlas = os.path.join(outputdir, 'hlas')
    SUPPORTSOFTNUM = 2
    Class1ReStr = 'A|B|C|E|F|G'  # class 1
    Class1List = []
    Class2List = []
    with open(inhlaXls, 'r') as ff1:
        for line in ff1:
            line = line.rstrip()
            listLine = line.split("\t")
            if listLine[0] == "HLA":
                continue
            else:
                if re.match(Class1ReStr, listLine[0]):
                    lista = listLine[1].split("|")
                    if len(set(lista)) >= SUPPORTSOFTNUM:
                        Class1List.append(listLine[0])
                else:
                    lista = listLine[1].split("|")
                    if len(set(lista)) >= SUPPORTSOFTNUM:
                        Class2List2.append(listLine[0])
    # print(Class1List, Class2List)
    # sed file
    datafile = os.path.join(SCRIPT_FOLDER, 'data', 'hla_gen.infor.maxlen')
    dictChange = {}
    with open(datafile, 'r') as ff:
        for line in ff:
            line = line.strip()
            linelist = line.split("\t")
            allele = linelist[-1]
            allele_4digit = linelist[2]
            dictChange[allele] = allele_4digit

    Class1List = getHeterozygousList(Class1List)
    numOfHla = 0
    with open(hlas, 'w') as ff:
        for allele in Class1List:
            if allele in dictChange:
                ff.write("{}\n".format(dictChange[allele]))
                numOfHla += 1
            else:
                print("# Alleles {} don\'t have 4digit name in {}".format(
                    allele, datafile))
    prefix = ""
    tumorbam = os.path.realpath(args.tumorbam)
    normalbam = os.path.realpath(args.normalbam)
    Tumor_SampleId = os.path.basename(tumorbam).split("_")[0]
    Normla_SampleId = os.path.basename(normalbam).split("_")[0]
    prefix = Tumor_SampleId + "_" + Normla_SampleId
    if numOfHla == 0:
        print("# All Alleles is homozygous, Hla loss cannot be calculated")
        statufile = os.path.join(outputdir, prefix + ".all-homozygous" + ".fail")
        if not os.path.exists(statufile):
            os.mknod("{}".format(statufile))
        sys.exit()
    # =================
    # chamge input bam Name
    # =================

    # the suffix of  normal bam file must be "_GL_sorted.bam"
    # the suffix of tumor bam must be "_tumor_sorted.bam"
    # don't support link file, because will used bam file

    inTbam = os.path.join(outputdir, prefix + "_tumor_sorted.bam")
    inNbam = os.path.join(outputdir, prefix + "_BS_GL_sorted.bam")

    # os.system("cp {} {}".format(tumorbam, inTbam))
    # os.system("cp {} {}".format(normalbam, inNbam))
    # os.system("cp {} {}".format(tumorbam + '.bai', inTbam + '.bai'))
    # os.system("cp {} {}".format(normalbam + '.bai', inNbam + '.bai'))
    coordinates = "chr6:29580000-33618227"
    interceptionChr6(SCRIPT_FOLDER, samtools, tumorbam, coordinates, inTbam)
    interceptionChr6(SCRIPT_FOLDER, samtools, normalbam, coordinates, inNbam)
    # =================
    # creat copyNumsolution file
    # =================
    tumorbamName2 = os.path.splitext(os.path.basename(inTbam))[0]
    copyNumsolution = os.path.join(outputdir, 'copyNumsolution.txt')
    purity = ""
    ploidy = ""
    with open(inPurity, 'r') as ff:
        for line in ff:
            line = line.strip()
            linelist = line.split("\t")
            if line.startswith("sample_name"):
                continue
            purity = linelist[1]
            ploidy = linelist[2]
            break
    with open(copyNumsolution, 'w') as ff2:
        ff2.write("Ploidy\ttumorPurity\ttumorPloidy\t\n")
        ff2.write("{tumorbamName}\t{ploidy}\t{purity}\t{ploidy}\t\n".format(
            purity=purity, ploidy=ploidy, tumorbamName=tumorbamName2))
    print(tumorbamName2, inTbam)
    # =================
    # run hla-loshls docker for calcu
    # =================
    command = "bash {SCRIPT_FOLDER}/lohhla_docker2.sh \
        {inputTumorBam} \
        {inputNormalBam} \
        {copyNumsolution} \
        {hlas} ".format(SCRIPT_FOLDER=SCRIPT_FOLDER,
                             inputTumorBam=inTbam,
                             inputNormalBam=inNbam,
                             hlas=hlas,
                             copyNumsolution=copyNumsolution)
    os.system(command)
    # end


if __name__ == '__main__':
    try:
        SCRIPT_FOLDER = os.path.abspath(os.path.dirname(__file__))
        configfile = os.path.join(SCRIPT_FOLDER, "config.ini")
        parser = argparse.ArgumentParser(
            prog="HLAloh".format(__VERSION__),
            description="get input hla for lohhla soft")
        parser.add_argument(
            '-l',
            '--inhlaXls',
            type=str,
            help="input file, such as BJ-B_resulting_number_of_software.xls",
            required=True)
        parser.add_argument('-p',
                            '--inPurity',
                            type=str,
                            help="input Purity , such as BJ-B_cn_summary.txt",
                            required=True)
        parser.add_argument('-t',
                            '--tumorbam',
                            type=str,
                            help="input tumor bam",
                            required=True)
        parser.add_argument('-n',
                            '--normalbam',
                            type=str,
                            help="input normal bam",
                            required=True)
        parser.add_argument('-o',
                            '--outputdir',
                            type=str,
                            help="specify output directory",
                            required=True)
        parser.add_argument(
            '-c',
            '--config',
            type=str,
            help="the config.ini file, default $SCRIPT_FOLDER/config.ini ",
            default=os.path.join(SCRIPT_FOLDER, 'config.ini'),
            metavar='')
        args = parser.parse_args()
        cfg = ConfigParser()
        cfg.read(configfile)
        samtools = cfg.get("soft", "SAMTOOLS")
        inhlaXls = os.path.abspath(args.inhlaXls)
        outputdir = os.path.abspath(args.outputdir)
        inPurity = os.path.join(args.inPurity)
        mkoutdir(outputdir)
        main()
    except KeyboardInterrupt:
        pass
    except IOError as e:
        raise
