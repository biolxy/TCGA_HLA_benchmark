# TCGA_HLA_benchmark

Due to TCGA privacy concerns, the **password** for HLA allele data is provided only to researchers who have access to the primary sequence data from the TCGA. Please contact us at qiliu@tongji.edu.cn

## software version

Eight hla genotyping tools with default parameter were used to process the TCGA WES datasets. 

| Tools         | URL                                              | version             |
| ------------- | ------------------------------------------------ | ------------------- |
| POLYSOLVER    | https://hub.docker.com/r/sachet/polysolver/tags  | v4                  |
| OptiType      | https://hub.docker.com/r/fred2/optitype/tags     | v1.3.1              |
| xHLA          | https://hub.docker.com/r/humanlongevity/hla/tags | DIGEST:425487b52034 |
| HLA-HD        | https://www.genome.med.kyoto-u.ac.jp/HLA-HD/     | v1.2.0.1            |
| hla-genotyper | https://pypi.org/project/hla-genotyper/          | v0.4.2b1            |
| SOAP-HLA      | http://soap.genomics.org.cn/SOAP-HLA.html        | v2.2                |
| HLA-VBSeq     | http://nagasakilab.csml.org/hla/                 | v2                  |
| Kourami       | https://github.com/Kingsford-Group/kourami       | v0.9.6              |
| LOHHLA        | https://hub.docker.com/r/mleventhal/lohhla/tags  | DIGEST:60253161bf8e |

## Citation
**Li, X., Zhou, C., Chen, K., Huang, B., Liu, Q. and Ye, H. (2021), Benchmarking HLA genotyping and clarifying HLA impact on survival in tumor immunotherapy. Mol Oncol. https://doi.org/10.1002/1878-0261.12895 PMID:33411982**
