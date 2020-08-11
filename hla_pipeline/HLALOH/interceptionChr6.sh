##################################################
# File Name   : interceptionChr6.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2019年09月27日 星期五 16时46分30秒
##################################################
#!/bin/bash
samtools=$1
inputbam=$2
coordinates=$3
outputbam=$4

${samtools} view -t 8 -h -b ${inputbam} ${coordinates} > ${outputbam}
${samtools} index ${outputbam}
