##################################################
# File Name     : OptiType_pipeline.sh
# Author      : biolxy
# E-mail      : biolxy@aliyun.com
# Created Time: 2018年10月15日 星期一 09时59分40秒
##################################################
#!/bin/bash
pid=$$
echo $pid
shell_script_dir=$(cd `dirname $0`; realpath .)
inputbam=$1
bamfilename=`basename ${inputbam}`
dockername=${bamfilename}.optitype
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
chmod +777 ${outputdir}/${prefix}
# check bam is single or paired
source ${shell_script_dir}/check_is_single.sh
type=$(check_is_single ${inputbampath})
echo $type
if [[ ${type} = 'paired' ]];then
    python ${shell_script_dir}/01.bam_to_fastq.py ${inputbampath} -1 ${outputdir}/${prefix}/${prefix}_1.fastq -2 ${outputdir}/${prefix}/${prefix}_2.fastq
    python ${shell_script_dir}/02.fastq_to_pairfastq.py ${outputdir}/${prefix}/${prefix}_1.fastq ${outputdir}/${prefix}/${prefix}_2.fastq
    fq1=${prefix}_pair_1.fastq
    fq2=${prefix}_pair_2.fastq
    docker run --name ${dockername} -i -v ${outputdir}/${prefix}:/output fred2/optitype  -i  /output/${fq1} /output/${fq2} --dna -v -o /output -p ${prefix}
    # docker run -v `pwd`:/input fred2/optitype  -i /input/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice_1.fastq  /input/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice_2.fastq --dna -v -o /input/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice
    # rm  -f ${outputdir}/${prefix}/*.fastq 
else
    bedtools bamtofastq -i ${inputbampath} -fq ${outputdir}/${prefix}/${prefix}.fq
    docker run --name ${dockername} -i -v ${outputdir}/${prefix}:/output fred2/optitype  -i /output/${prefix}.fq --dna -v -o /output -p ${prefix}
fi
docker rm ${dockername}
