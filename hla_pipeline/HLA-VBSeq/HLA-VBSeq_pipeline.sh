##################################################
# File Name     : aa.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年10月15日 星期一 16时17分39秒
##################################################
#!/bin/bash
shell_script_dir=$(cd `dirname $0`; realpath .)
hla_ref=${shell_script_dir}/database/hla_all.fasta
allele=${shell_script_dir}/database/Allelelist.txt
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
fi
python ${shell_script_dir}/01.bam_to_fastq.py ${inputbampath} -1 ${outputdir}/${prefix}/${prefix}_1.fastq -2 ${outputdir}/${prefix}/${prefix}_2.fastq
python ${shell_script_dir}/02.fastq_to_pairfastq.py ${outputdir}/${prefix}/${prefix}_1.fastq ${outputdir}/${prefix}/${prefix}_2.fastq
fq1=${outputdir}/${prefix}/${prefix}_pair_1.fastq
fq2=${outputdir}/${prefix}/${prefix}_pair_2.fastq
bwa mem -t 8 -P -L 10000 -a $hla_ref $fq1 $fq2 > ${outputdir}/${prefix}/${prefix}_mapped.sam
java -jar $shell_script_dir/HLAVBSeq.jar $hla_ref ${outputdir}/${prefix}/${prefix}_mapped.sam ${outputdir}/${prefix}/${prefix}_result.txt --alpha_zero 0.01 --is_paired
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^A\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_A.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^B\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_B.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^C\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_C.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DPA1\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DPA1.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DPB1\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DPB1.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DQA1\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DQA1.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DQB1\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DQB1.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRA\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRA.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB1\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB1.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB2\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB2.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB3\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB3.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB4\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB4.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB5\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB5.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB6\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB6.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB7\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB7.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB8\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB8.txt
$shell_script_dir/parse_result.pl ${allele} ${outputdir}/${prefix}/${prefix}_result.txt | grep "^DRB9\*" | sort -k2 -n -r > ${outputdir}/${prefix}/HLA_DRB9.txt
#rm -rf ${outputdir}/${prefix}/$sample.sam
#rm -rf ${outputdir}/${prefix}/$FQ1 ${outputdir}/${prefix}/$FQ2
for i in `ls ${outputdir}/${prefix}/* | grep -v ".txt"`
do
    rm -rf ${i}
done
echo "end at time  `date +%F' '%H:%M` "
