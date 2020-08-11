##################################################
# File Name     : 3.7.batch.soaphla_pipeline.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年10月23日 星期二 13时16分19秒
##################################################
#!/bin/bash
conffile=`realpath $1`
bamdir=`realpath $2`
outputdir=`realpath $3`
logfile="$0_$(date +"%Y-%m-%d-%H-%M").log"
exec > $logfile 2>&1
# echo "The jobs number is ${num} :"
# num=$((num + 1))
# echo "[ `date +%F\ %H:%M` ] command is : bash /home/lixy/production/2018-10-12_11soft/xHLA_pipeline/xHLA_pipeline.sh ${bam} ${outputdir} > ${outputdir}/${bamName}.log 2>&1"
python /home/lixy/production/2018-10-23_other_4_pipeline/util/soaphla.py -w ${conffile} -s ${bamdir} -d ${outputdir} 
