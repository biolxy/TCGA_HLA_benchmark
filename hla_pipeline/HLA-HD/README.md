# hlahd pipeline

## 1、hlahd 

- <https://www.genome.med.kyoto-u.ac.jp/HLA-HD/>  

## 2、download hlahd docker 

```shell
docker pull biolxy/hlahd 
```

- https://hub.docker.com/r/biolxy/hlahd/

## 3、usage

install python package
```
pip install pysam
```
run in terminal
```shell
bash 02_print_hlahd_CMD.sh  882f41ab-a595-4328-9e05-310303b63a8c_regions_slice.bam /home/lixy/production/2018-10-12_11soft
```

argv:

- inputfile 882f41ab-a595-4328-9e05-310303b63a8c_regions_slice.bam
- outpuffolder /home/lixy/production/2018-10-12_11soft


`/home/lixy/production/2018-10-12_11soft/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice/result/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice_final.result.txt `

## 3、usage2

```shell
bash 3.3.batch.hlahd_pipeline.sh inputbam outputdir
```



