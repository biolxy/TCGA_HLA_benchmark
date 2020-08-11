##################################################
# File Name   : lohhla_docker.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2019年03月05日 星期二 16时13分32秒
##################################################
#!/bin/bash

# usage： bash lohhla.sh  example /opt/lohhla/tt  /opt/lohhla/example-file/bam/example_BS_GL_sorted.bam \
# /opt/lohhla/example-file/bam /opt/lohhla/example-file/solutions.txt /opt/lohhla/example-file/hlas
start=$(date +'%s')
startTimestamp=$(date +"%y:%m:%d_%H:%M:%S")
echo "#start time stamp ${startTimestamp}" 
shell_script_dir=$(cd `dirname $0`; realpath .)
# /opt/lohhla/lohhla.sh
docker_images_name='hla-lohlos:1.0'

inputTumorbam=`realpath ${1}`
inputdir1=`dirname $inputTumorbam`
inputNormalbam=`realpath ${2}`
baseNormalbam=`basename ${inputNormalbam}`
inputdir2=`dirname $inputNormalbam`
# inputdir1=inputdir2
outputdir=$5
bambasename=`basename $inputTumorbam`
bambasename3=${bambasename%%_*}
echo $bambasename3
bambasename=$bambasename3
bambasename2=`basename $inputNormalbam`
logfile="$bambasename.$(date +"%Y%m%d%H%M").log"
prefix=${bambasename}
pid=$$
echo "Started at `date`"

# 创建 copyNumSolutions 文件
copyNumSolutions=`basename $3`
input3=`dirname $3`
# 生产 hlatype 文件

hlatype=`basename $4`
input4=`dirname $4`

echo "inputTumorbam file_path is            : $inputTumorbam"
echo "inputNormalbam file_path is           : $inputNormalbam"
echo "output dirname is                     : $outputdir"
echo "this copyNumSolutions file is         : ${input3}/${copyNumSolutions}"
echo "this hlatype file is                  : ${input4}/${hlatype}"
# exec >$outputdir/$logfile 2>&1
pid=$$
echo "This task pid is                      : $pid"

# tumor bam 改名字,格式为   {example}_tumor_sorted.bam
if [ ! -f  ${inputdir1}/${bambasename}_tumor_sorted.bam ];then
    mv ${inputTumorbam} ${inputdir1}/${bambasename}_tumor_sorted.bam
    mv ${inputTumorbam}.bai ${inputdir1}/${bambasename}_tumor_sorted.bam.bai
fi

docker run --rm -v ${inputdir1}:/input1 -v ${inputdir2}:/input2 -v ${input3}:/input3 -v ${input4}:/input4 -v ${outputdir}:/output ${docker_images_name} \
bash -c "source ~/.bashrc; bash /opt/lohhla/lohhla.sh \
    ${bambasename} \
    /output \
    /input1/${baseNormalbam} \
    /input1 \
    /input3/${copyNumSolutions} \
    /input4/${hlatype}"

end=$(date +'%s')
endTimestamp=$(date +"%y:%m:%d_%H:%M:%S")
echo "#end time stamp ${endTimestamp}"
echo "#It took $((($end - $start)/60)) minutes"

