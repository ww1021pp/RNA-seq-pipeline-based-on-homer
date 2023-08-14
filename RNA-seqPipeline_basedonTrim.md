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
Map the trimmed sequence reads to reference genome using approciat softwares(STAR, tophat or bowtie2 et.al). Specifically, we first create a directory with the sample name to store the output then calling tophat2 and output results to the directory’s sub-directory named Tophat_Out. Since this is a mouse sample, UCSC mm10 is used as the reference genome and the reference transcriptome file and bowtie index files are stored under Indexes directory in the home directory. Bowtie2 index files contain the genome sequences to be aligned to in bowtie2 format. Tophat2 uses bowtie2 as the base sequence aligner.
```
STAR --runThreadN 12 --genomeDir genomeIndex --readFilesCommand zcat --readFilesIn R1.fastq.gz R2.fastq.gz --outFileNamePrefix $dir/$package/$filename --outSAMunmapped Within --outFilterType BySJout --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outWigType bedGraph
```
## step 4: Assemble gene expression from aligned reads

Use homer(HTSeq) to quantify the number of reads mapped to each gene and make a genome annotation with different parameters.

```
analysisRepeat.pl -h  ## use homer to do gene expression quantification
samtools view MT1/Tophat_Out/accepted_hits.sorted.bam | python -m
HTSeq.scripts.count -q -s no - ~/Indexes/Mus_musculus/UCSC/mm10/Genes/genes.gtf >MT1/MT1.count.txt
```
### After you have installed HTSeq (see Prequisites and installation), you can run htseq-count from the command line:
```
htseq-count [options] <alignment_files> <gff_file>
```
### If the file htseq-count is not in your path, you can, alternatively, call the script with

```
python -m HTSeq.scripts.count [options] <alignment_files> <gff_file>
```
## step 5: Differential gene expression analysis from assembled gene expression (DEseq2 or EdgeR)

