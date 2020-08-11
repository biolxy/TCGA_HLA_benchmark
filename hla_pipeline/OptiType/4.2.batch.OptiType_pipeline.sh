##################################################
# File Name     : 3.2.batch.OptiType_pipeline.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年10月23日 星期二 13时16分19秒
##################################################
#!/bin/bash
inputdir=`realpath $1`
outputdir=`realpath $2`
logfile="$0.$(date +"%Y%m%d%H%M").log"
exec >${inputdir}/$logfile 2>&1
num=1
for bam in $(cat list)
do
    bamName=`basename ${bam}`
    echo "The jobs number is ${num} :"
    num=$((num + 1))
    echo "command is : bash /home/lixy/production/2018-10-12_11soft/OptiType_pipeline/OptiType_pipeline.sh  ${bam} ${outputdir} > ${outputdir}/${bamName}.log 2>&1"
    if [ -f ${outputdir}/${bamName}.log ];then
        logsize=$(ls -l ${outputdir}/${bamName}.log  | awk '{ print $5}')
        if [ $logsize == 568 -o $logsize == 750 ] ;then
            bash /home/lixy/production/2018-10-12_11soft/OptiType_pipeline/OptiType_pipeline2.sh ${bam} ${outputdir} > ${outputdir}/${bamName}.log 2>&1
            echo "[ `date +%F\ %H:%M` ] docker rm :"
            # docker rm $(docker ps -a | grep "fred2/optitype" | awk '{ print $1 }' )
        fi
    else
        bash /home/lixy/production/2018-10-12_11soft/OptiType_pipeline/OptiType_pipeline2.sh ${bam} ${outputdir} > ${outputdir}/${bamName}.log 2>&1
        echo "[ `date +%F\ %H:%M`  ] docker rm :"
        docker rm $(docker ps -a | grep "fred2/optitype" | awk '{ print $1  }' )
    fi
done
