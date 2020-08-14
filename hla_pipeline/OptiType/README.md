# OptiType

## 1. download docker images

```
docker pull fred2/optitype
```

## 2. install pysam

```shell
pip install pysam
```

## 3. usage

- input1: `/storage1/TCGA/TCGA-CHOL/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice.bam`
- input2: `/home/lixy/production/2018-10-12_11soft/OptiType_pipeline`  

**run in  terminal: **

```shell
bash OptiType_pipeline.sh 882f41ab-a595-4328-9e05-310303b63a8c_regions_slice.bam
```

- output: 

```
$ tree 882f41ab-a595-4328-9e05-310303b63a8c_regions_slice
882f41ab-a595-4328-9e05-310303b63a8c_regions_slice
├── 882f41ab-a595-4328-9e05-310303b63a8c_regions_slice_coverage_plot.pdf
└── 882f41ab-a595-4328-9e05-310303b63a8c_regions_slice_result.tsv
```

`.tsv` hla type:

```shell
$ less 882f41ab-a595-4328-9e05-310303b63a8c_regions_slice/882f41ab-a595-4328-9e05-310303b63a8c_regions_slice_result.tsv
        A1      A2      B1      B2      C1      C2      Reads   Objective
0       A*31:01 A*02:01 B*51:02 B*39:05 C*07:02 C*08:01 2837.0  2760.3909999999946
```

## 4. batch run

```shell
for i in `ls /input_bam_dir/*.bam`
do
	bash OptiType_pipeline.sh $i
done
```

`/input_bam_dir`  is `bam` folder

## 5. contact

- author:  biolxy@aliyun.com


