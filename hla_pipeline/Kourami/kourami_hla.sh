##################################################
# File Name   : kourami_hla.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年11月28日 星期三 13时08分16秒
##################################################
#!/bin/bash
docker_images_name='kourami_hla'
inputbam=`realpath ${1}`
inputdir=`dirname ${inputbam}`
outputdir=`realpath ${2}`
bamname=`basename ${inputbam}`
output_prefix=${bamname/.bam/}
# namepre_on_KouramiPanel.bam
cmd1="/usr/local/kourami/scripts/alignAndExtract_hs38DH.sh -d /usr/local/kourami/db /output/${output_prefix} /input/${bamname}"
cmd2="java -jar /usr/local/kourami/target/Kourami.jar -d /usr/local/kourami/db/ -o /output/${output_prefix} /output/${output_prefix}_on_KouramiPanel.bam"
echo ${cmd1}
echo ${cmd2}
echo " docker run -iv ${inputdir}:/input -v ${outputdir}:/output -w ${inputdir} ${docker_images_name} /bin/bash -c \" ${cmd1}; ${cmd2} \" "
docker run -iv ${inputdir}:/input -v ${outputdir}:/output -w ${outputdir} ${docker_images_name} /bin/bash -c " ${cmd1}; ${cmd2}"