# RNA-seq pipeline from rawdata to expression table 

## softwares used in th pipeline
1. ** fastQC ** for reads sequencing quanlity. FastQC provides a quick view on the quality of the raw sequence reads from multiple analyses, ranging from the sequence quality, GC content, to library complexity. 
2. ** fastp ** trim adapter and filter low quanlity reads
3. ** STAR ** (Spliced Transcripts Alignment to a Reference) (maybe use samtools to filter low quanlity alignments)
4. use homer to make tagDirectory based on alignment results
5. get the expression table and annotation with homer

## I construct a bash script for the pipeline from raw reads to expression matrix include Pair-ends and sigle-end reads sequencing data 
