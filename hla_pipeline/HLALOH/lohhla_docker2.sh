##################################################
# File Name   : lohhla_docker.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2019年03月05日 星期二 16时13分32秒
##################################################
#!/bin/bash
# usage： bash lohhla.sh  example /opt/lohhla/tt  /opt/lohhla/example-file/bam/example_BS_GL_sorted.bam \
# /opt/lohhla/example-file/bam /opt/lohhla/example-file/solutions.txt /opt/lohhla/example-file/hlas
function usage(){
    echo "Usage:"
    echo "      " bash $0 {inputTumorbam} {inputNormalbam} {copyNumSolutions} {hlas}
    echo "tips --> this four file should in same folder, and output dir is also this folder."
    echo ""
}
if [ $# != 4 ];then
    usage
    exit
fi


start=$(date +'%s')
startTimestamp=$(date +"%y:%m:%d_%H:%M:%S")
echo "#start time stamp ${startTimestamp}" 
shell_script_dir=$(cd `dirname $0` && pwd)
docker_images_name='hla-lohlos:1.1'

inputTumorbam=`realpath ${1}`
inputdir1=$(dirname $inputTumorbam)
BaseTumorBam=$(basename ${inputTumorbam})
inputNormalbam=`realpath ${2}`
BaseNormalBam=`basename ${inputNormalbam}`
# 创建 copyNumSolutions 文件
copyNumSolutions=$(basename $3)
# 生产 hlatype 文件
hlatype=$(basename $4)

SampleID=`basename $inputTumorbam`
SampleID=${SampleID%%_*}

outputdir=${inputdir1}


echo "inputTumorbam file_path is            : $inputTumorbam"
echo "inputNormalbam file_path is           : $inputNormalbam"
echo "output dirname is                     : $outputdir"
echo "this copyNumSolutions file is         : ${copyNumSolutions}"
echo "this hlatype file is                  : ${hlatype}"
echo "this project SampleID is              : ${SampleID}"
docker run --rm -v ${inputdir1}:/input1 ${docker_images_name} \
bash -c "source ~/.bashrc;cd /input1 ; bash /opt/lohhla/lohhla.sh \
    ${SampleID} \
    /input1 \
    ${BaseNormalBam} \
    ${BaseTumorBam} \
    ${copyNumSolutions} \
    ${hlatype} "



end=$(date +'%s')
endTimestamp=$(date +"%y:%m:%d_%H:%M:%S")
echo "#end time stamp ${endTimestamp}"
echo "#It took $((($end - $start)/60)) minutes"
