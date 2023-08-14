# Data analysis pipeline for RNA-seq experiments: From raw data to differenetial expression


## step 1: Quality check on the raw reads.
```
fastqc Fastq/*.fastq –o FastQC
```

## step 2:  Groom raw reads.

Remove sequences with low quality to get better alignment in the later steps. Based on the FastQC reports from above commands for this example,
the qualities of the reads are fine except random distribution of sequence content at the 5’ end of reads.Thus, we trim 10bp from the beginning of each read: 
```
awk -v s=10 -v e=0 ‘ {if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; } ‘ Fastq/MT1_1.fastq > Fastq/Trim_MT1_1.fastq
awk -v s=10 -v e=0 ‘ {if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; } ‘ Fastq/MT1_2.fastq > Fastq/Trim_MT1_2.fastq
```
Use used the awk command native to Unix system here. The s=10 indicates that the first 10bp (the 5’end) will be trimmed, where e=0 indicates none will be trimmed 
from the end (3’ end). Users should change these according to the FastQC reports. File name for input fastq reads is specified before the > sign, where the file 
name for the trimmed reads should be followed after the > sign

## step 3: map clean reads to reference genome