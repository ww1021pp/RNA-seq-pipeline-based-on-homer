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
## step 5: Differential gene expression analysis from assembled gene expression for samples with repeat(DEseq2 or EdgeR), if there is no replicants, just do odd ration for each gene

### 5.1. Launch RStudio and load necessary library
```
library(“DESeq2”)
library(reshape2)
library(stringr)
library(ggpubr)
library(DESeq2)
library(data.table)
```
### 5.2(1) Create necessary data object from HTSeqCount results.

```
sample.names <- sort(paste(c(“MT”, “WT”), rep(1:3, each=2), sep=““))
file.names <- paste(“../”, sample.names, “/”, sample.names, “.count.txt”, sep=““)
conditions <- factor(c(rep(“MT”, 3), rep(“WT”, 3)))
sampleTable <- data.frame(sampleName=sample.names,
fileName=file.names,
condition=conditions)
# read in the HTSeq count data 
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=“.”, 
design=~ condition )
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 10, ]
```
In this sub-step, we specify the sample identifiers, name of files with gene counts for each sample, and experiment condition(s) for each sample; then pass this information to DESeqDataSetFromHTSeqCount function to make a DESeqDataSet object for following analysis.

### 5.2(2) Create necessary data object from expression count matrix
```
###setwd################
rm(list = ls())
setwd("/workspace/rsrch2/panpanliu/project/RNA-seq/Jiahua_0104/hs_Liver_DE_DeSeq2")

SampleType = OutputPrefix = "DE"
genecols<-c("Transcript_ID","chr", "start", "end", "strand", "Length", "Copies", "Symbol","Annotation.Divergence")
FC = 0.586
pvalue = 0.05

library(ggplot2)
library(reshape2)
library(stringr)
library(ggpubr)
library(DESeq2)
library(data.table)
raw <- read.delim("../4_homer/rawCount.txt", header = TRUE,stringsAsFactors = FALSE)

colnames(raw)<-sapply(colnames(raw),function(x) gsub("\\.\\..*","",as.character(x)))
colnames(raw)[1]="Transcript_ID"
rownames(raw)=raw$Transcript_ID
#######################the same as the command###################################
raw$Symbol<-sapply(strsplit(raw$Annotation.Divergence,"\\|"),getElement, 1)

##########################

all_tag=c("Norm","NAT","NASH","HCC")

#######################remove the transcripts without gene annotation(As we focused on the known protein)#####
raw<-raw[which(raw$Symbol != "0.000"),]
print(dim(raw))
raw<-raw[order(raw$Symbol,raw$Length,decreasing=TRUE),]
raw_uniq <- raw[match(unique(raw$Symbol), raw$Symbol),]
print(dim(raw_uniq))



coldata=DataFrame(
  batch=factor(sapply(strsplit(colnames(raw_uniq)[9:(dim(raw_uniq)[2]-1)],"\\_"),getElement, 1),levels = c("Norm","NAT","NASH","HCC")),
 # strain=factor(sapply(strsplit(colnames(raw_uniq)[9:(dim(raw_uniq)[2]-1)],"\\."),getElement, 2),levels = c("B6","129")),
  condition=factor(sapply(strsplit(colnames(raw_uniq)[9:(dim(raw_uniq)[2]-1)],"\\_"),getElement, 1),levels = c("Norm","NAT","NASH","HCC") )
)

head(coldata);dim(coldata)
rownames(coldata)=colnames(raw_uniq)[9:(dim(raw_uniq)[2]-1)]

coldata$cnd<-sapply(coldata$condition,function(x) gsub(all_tag[1],"0",x))
coldata$cnd<-sapply(coldata$cnd,function(x) gsub(all_tag[2],"1",x))
coldata$cnd<-sapply(coldata$cnd,function(x) gsub(all_tag[3],"2",x))
coldata$cnd<-sapply(coldata$cnd,function(x) gsub(all_tag[4],"3",x))

coldata$group<-coldata$cnd
#coldata$group<-relevel(coldata$group,ref="1")

as.data.frame(coldata)
Filter_value=10
row_mean=data.frame(
Nrom = rowMeans(raw_uniq[,grepl(all_tag[1],colnames(raw_uniq))]),
NAT = rowMeans(raw_uniq[,grepl(all_tag[2],colnames(raw_uniq))]),
NASH = rowMeans(raw_uniq[,grepl(all_tag[3],colnames(raw_uniq))]),
HCC = rowMeans(raw_uniq[,grepl(all_tag[4],colnames(raw_uniq))])
)

head(row_mean)
row_mean$count<-rowSums(row_mean[,1:4]>3)
rownames(row_mean)=rownames(raw_uniq)
raw_uniq_filter<-raw_uniq[row_mean$count>2,]
head(raw_uniq_filter);dim(raw_uniq_filter)

count<-raw_uniq_filter[,9:(dim(raw_uniq_filter)[2]-1)]
rownames(count)<-raw_uniq_filter$Transcript_ID
dds <- DESeqDataSetFromMatrix(countData=round(count), colData=coldata, design= ~group)

```



5.3. Run differential gene analysis.
```
 dds <-DESeq(ddsHTSeq)



