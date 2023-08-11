# The bash script for RNA-seq pipeline for all samples in a for loop.

## This script is easy to use, just input two parameters 
1. **dir** the work directory
2. ** i ** the log file number


```
#!/usr/bin/bash
dir=$1
i=$2

PEsplit="_" 
FQsuffix=".fq.gz"    ## fastq file type
depth="sample_depth.txt" 

genomeVerstion="hg38"   ##reference genome version
Qc=1_Fastp                         
package=2_star_mapping_GCR38 
HomerOut=4_homer
index="/workspace/rsrch2/common_data/Align2Index/StarIndex/hg38" ## STAR alignment reference genome index
gtf_file="/workspace/rsrch2/common_data/Refgenome/hg38/chr.Homo_sapiens.GRCh38.94.gtf"  ## gtf file for reference genes to do quantity gene expression
chr_size="/workspace/rsrch2/common_data/Align2Index/StarIndex/hg38/chrNameLength.txt" ## the size of chromosome for make bigwig file 
 ##the MetaData download from GEOi
################################################################################

if [ ! -d $dir/$Qc/ ]    ## make directory 
then
  mkdir -p $dir/$Qc/     
fi

if [ ! -d $dir/$package ]; then
  mkdir -p $dir/$package
fi

if [ ! -d $dir/$HomerOut/tagDirectory/ ]; then
  mkdir -p $dir/$HomerOut/tagDirectory/
fi

if [ ! -d $dir/$HomerOut/bigwig/ ];  then
  mkdir -p $dir/$HomerOut/bigwig/
fi  


if [ ! -d "$dir"/logs/ ]
then 
 mkdir $dir/logs/
fi

TagDir=$dir/$HomerOut/tagDirectory 
#wigfile="Aligned_Wig_files="

##change the parameters
SampleInfo=$(cat $dir/rawdata/Gsm_Name.txt|cut -f 2|sed '1d') ## the sample information 

for Info in $SampleInfo  
  do
    
    filename=$(cat $dir/rawdata/Gsm_Name.txt| grep -w $Info| awk -F '\t' '{print $1}'|sort|uniq)
    echo -e "Now process with $Info and its sample name is $filename at $(date)\n" >>$dir/logs/StarBigwik"$i".log
   
  
    cd $dir/rawdata
    echo -e "$Info\t" >> $dir/filesize.txt
    du -cm --max-depth=1 $dir/rawdata/${Info}/${Info}_*.gz >>$dir/filesize.txt
    
    InputTagDir+=" $filename"

    #step 2 Use fastp to  reads to do QC
     echo -e "####################Now do fastp######################\n" >>$dir/logs/StarBigwik"$i".log
  fastp -i $dir/rawdata/${Info}/${Info}_*_1*${FQsuffix} -I $dir/rawdata/${Info}/${Info}_*_2*${FQsuffix}  -o $dir/$Qc/$filename"_1"$FQsuffix -O $dir/$Qc/$filename"_2"$FQsuffix -h $dir/$Qc/$filename.html

     # Step 3 use star to do mapping 
   echo -e "now is mapping $Info and its sample name is $filename\n" >>$dir/logs/StarBigwik"$i".log

  STAR --runThreadN 12 --genomeDir $index --readFilesCommand zcat --readFilesIn $dir/$Qc/"$filename"_1"$FQsuffix" $dir/$Qc/"$filename"_2"$FQsuffix" --outFileNamePrefix $dir/$package/$filename --outSAMunmapped Within --outFilterType BySJout --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outWigType bedGraph
   
   if [ -s $dir/$package/"$filename"Aligned*.bam ]
    then
    rm $dir/$Qc/"$filename"*.gz
    rm $dir/$package/${filename}*.bg
    rm -rf $dir/$package/${filename}*STAR*
    rm $dir/$package/${filename}*.tab
   fi    
  
  

    ##step 5 use homer to creat TagDirectory
    echo -e "now make TagDir for $Info and Its sample is $filename\n" >> $dir/logs/StarBigwik"$i".log 
makeTagDirectory $TagDir/$filename $dir/$package/${filename}Aligned.sortedByCoord.out.bam -sspe
    ##step 6 use homer to creat bigwik file
makeUCSCfile $TagDir/$filename -bigWig $chr_size -fsize 1e20 -style rnaseq -strand both -o $dir/$HomerOut/bigwig/"$filename".bigwig

#stringtie "$dir"/3_Unique/"$filename".sorted.unique.bam -G $gtf_file -e -A $dir/7_Stringtie/$filename.tab -o $dir/7_Stringtie/$filename.gtf
  echo -e "The process of $Info whose sample name is $filename is done at $(date) \n" >> $dir/logs/StarBigwik"$i".log


done

 
#analyzeRepeats.pl rna $genomeVerstion -strand both -count exons -d$InputTagDir -noadj > $dir/$HomerOut/rawCount.txt
#analyzeRepeats.pl rna $genomeVerstion -strand both -count exons -d$InputTagDir -quantile > $dir/$HomerOut/quatile.txt
#analyzeRepeats.pl rna $genomeVerstion -strand both -count exons -d$InputTagDir -rpkm > $dir/$HomerOut/rpkm.txt
   ## Example for Mouse genome: mm10/GRCm38, annotation version: gencode.vM25
```
