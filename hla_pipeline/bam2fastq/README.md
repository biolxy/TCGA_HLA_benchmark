# usage

```python

python 01.bam_to_fastq.py $inputbam -1 ${outputdir}/${prefix}/${prefix}_1.fastq -2 ${outputdir}/${prefix}/${prefix}_2.fastq
python 02.fastq_to_pairfastq.py ${outputdir}/${prefix}/${prefix}_1.fastq ${outputdir}/${prefix}/${prefix}_2.fastq

```