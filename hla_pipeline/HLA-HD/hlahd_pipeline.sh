##################################################
# File Name     : batch.hlahd.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年09月26日 星期三 16时51分22秒
##################################################
#!/bin/bash
# data.list3 中是样本名字
# ls *.fastq | grep R1.fastq > list
# docker pull biolxy/hlahd 
shell_script_dir=$(cd `dirname $0`; realpath .)
hla_ref=${shell_script_dir}/database/hla_all.fasta
inputbam=$1
bamfilename=`basename ${inputbam}`
inputbampath=`realpath ${inputbam}`
inputdir=`dirname ${inputbampath}`
outputdir=`realpath $2`
prefix=${bamfilename%%.*}
echo "start at time `date +%F' '%H:%M`"
echo "inputbam file_path is :$inputbampath"
echo "output dirname is     :$outputdir"
if [ ! -d ${outputdir}/${prefix} ];then
    mkdir -p ${outputdir}/${prefix}
fi;
python ${shell_script_dir}/01.bam_to_fastq.py $inputbam -1 ${outputdir}/${prefix}/${prefix}_1.fastq -2 ${outputdir}/${prefix}/${prefix}_2.fastq
python ${shell_script_dir}/02.fastq_to_pairfastq.py ${outputdir}/${prefix}/${prefix}_1.fastq ${outputdir}/${prefix}/${prefix}_2.fastq
fq1=${prefix}_pair_1.fastq
fq2=${prefix}_pair_2.fastq
cmd1=" source ~/.bashrc ; hlahd.sh -t 8 -m 30 -c 0.95 -f /root/hlahd.1.2.0.1/freq_data/ /input/${fq1} /input/${fq2}  /root/hlahd.1.2.0.1/HLA_gene.split.3.32.0.txt  /root/hlahd.1.2.0.1/dictionary/  ${prefix}  /output; mv /output/${prefix}/result /output/;rm -rf /output/${prefix} "
cmd="docker run -iv ${outputdir}/${prefix}:/input -v ${outputdir}/${prefix}:/output  -w ${inputdir} docker.io/biolxy/hlahd:latest  /bin/bash -c  \" ${cmd1} \" "
echo ${cmd} > .tmp.sh
bash .tmp.sh
rm -f .tmp.sh
rm -r ${outputdir}/${prefix}/*.fastq


