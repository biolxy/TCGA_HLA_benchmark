##################################################
# File Name     : run_docker.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年10月16日 星期二 15时29分06秒
##################################################
#!/bin/bash
inputbam=$1
bamfilename=`basename ${inputbam}`
inputbampath=`realpath ${inputbam}`
inputdir=`dirname ${inputbampath}`
outputdir=`realpath $2`
prefix=${bamfilename%%.*}
shell_script_dir=$(cd `dirname $0`; realpath .)

cmd1="source ~/.bashrc; bash /root/HLA-VBSeq_pipeline/HLA-VBSeq_pipeline.sh /input/${bamfilename} /output "
cmd="docker run -iv ${inputdir}:/input -v ${outputdir}:/output hla-vbseq_pipeline-v0.4  /bin/bash -c \" ${cmd1} \" "
echo ${cmd} > .tmp.sh
bash .tmp.sh
rm .tmp.sh
