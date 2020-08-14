##################################################
# File Name   : test.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2019年07月16日 星期二 17时19分20秒
##################################################
#!/bin/bash
SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
inputTumorBam=$1
inputNormalBam=$2
copyNumsolution=`realpath $3`
hlas=`realpath $4`
output_dir=`realpath $5`
if [ -d ${output_dir} ];then
    mkdir ${output_dir}
fi

bash ${SHELL_FOLDER}/lohhla_docker.sh \
        ${inputTumorBam} \
        ${inputNormalBam} \
        ${copyNumsolution} \
        ${hlas} \
        ${output_dir}

