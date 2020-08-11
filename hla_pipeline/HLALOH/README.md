# HLALOH

**train of thought:**

Tumor purity and ploidy were calculated using Sclust, hla type was calculated using HLAtype, and tumor hla deletion was calculated using HLALOH.



```shell
docker tag  124767e7abcc hla-lohlos:1.0
```



## shell usage

```shell
bash ../../lohhla_docker2.sh 60290007_tumor_sorted.bam 60290007_BS_GL_sorted.bam copyNumsolution.txt hlas 
```
