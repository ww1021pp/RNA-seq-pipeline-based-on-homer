# this is the scripts used to do RAA-seq rutine analysis from raw sequencing data

```
#!/usr/bin/bash
################################################################################
# Help                                                                         #
################################################################################
if [ "$1" == "-h" | "$1" == "--help"] 
then
    echo "Usage 1: bash `basename $0` dirtory sampleID"
    echo "This Help File: sh `basename $0` -h"
    exit 0
fi

dir=$1
PEsplit="_"
Srr_dir=$dir/rawdata/SRR ##SRR_accssion download directory
SRR_metafile=$dir/rawdata/SraRunTable.txt
FQsubDIR="rawdata/Fastq" ## Fastq file store
FQsuffix=".fastq.gz"    ## may need to change baed on the fasta file suffix
depth="sample_depth.txt"  ## to save the file size of sequenceing reads number

index="/project/guan/panpan/AlignIndex/starIndex/mm10/GRCm38.primary_assembly.genome.M25"  ## mm10 star Index
annotation="/project/guan/panpan/ReferGonome/mm10/encode.vM25.primary_assembly.annotation.gff3"  ## refgene annotation file for gene annotation and gene expression quantity
gtf_file="/project/guan/chunjie/ReferenceGenome2/mm10/chr.Mus_musculus.GRCm38.99.gtf" ## refgene annotation file for gene annotation and gene expression quantity
genomeVerstion="mm10"   ## reference genome version   ## genome version for homer annotation and gene expression quantification
Qc=1_Fastp                                 
package=2_star_mapping_mm10v25
HomerOut=4_homer
homer="/project/guan/chunjie/Software/homer"  ##specific the homer directory
chr_size="/project/guan/panpan/AlignIndex/starIndex/mm10/GRCm38.primary_assembly.genome.M25/chrNameLength.txt"  ## for the bigwig file chromosome size file
 ##the MetaData download from GEO
###############################################################################


if [ ! -d "$dir"/"$Qc"/ ] 
then
  mkdir -p "$dir"/"$Qc"/""
fi

if [ ! -d "$dir"/"$package" ]; then
  mkdir -p "$dir"/"$package"
fi

if [ ! -d "$dir"/3_Unique ]
then
   mkdir -p $dir/3_Unique
fi

if [ ! -d "$dir"/"$HomerOut"/tagDirectory/ ]; then
  mkdir -p "$dir"/"$HomerOut"/tagDirectory/
fi

if [ ! -d "$dir"/"$HomerOut"/bigwig/ ];  then
  mkdir -p "$dir"/"$HomerOut"/bigwig/
fi

if [ ! -d "$dir"/wig/ ]
then 
 mkdir $dir/wig/
fi

if [ ! -d "$dir"/logs/ ]
then
 mkdir -p $dir/logs/
fi

if [ ! -d "$dir"/7_Stringtie/ ]
then
 mkdir -p $dir/7_Stringtie/
fi


TagDir="$dir"/"$HomerOut"/tagDirectory
wigfile="Aligned_Wig_files="

##change the parameters
SampleInfo=$2
i=$3
    SampleName=$(cat $SRR_metafile|grep $Info| awk -F ',' '{print $1}')
    SrrInfo=$(cat $dir/rawdata/Gsm_Name.txt| grep -w $Info| awk -F '\t' '{print $1}')
    filename=$(cat $dir/rawdata/Gsm_Name.txt| grep -w $Info| awk -F '\t' '{print $NF}'|sort|uniq)
    echo -e "Now process with $Info and its sample name is $filename at $(date)\n" >>$dir/logs/StarBigwik"$i".log
    echo -e "Now process with $Info and its SRR name is $SrrInfo at $(date)\n" >>$dir/logs/StarBigwik"$i".log
  
if [ ! -e $dir/$FQsubDIR/"$filename"_1.fastq.gz ] 
then 
   for srr in $SrrInfo
   do 
   echo "$srr"
   cd $Srr_dir
   if [ ! -e $srr.sra ]
    then
    wget -q -t 5 "https://sra-pub-run-odp.s3.amazonaws.com/sra/"$srr/$srr -O $srr.sra 
   fi
    
       if [ -e $srr.sra ] 
   then 
       if [ ! -s "$srr.sra" ]
       then 
          echo -e "$srr.sra is not download successfully\n" >>$dir/logs/StarBigwik"$i".log 
          wget -q -t 5 "https://sra-pub-sars-cov2.s3.amazonaws.com/run/"$srr/$srr -O $srr.sra
           if [ -s $srr.sra ]
           then
          echo -e "$srr.sra is download successly from second http\n" >>$dir/logs/StarBigwik"$i".log
          fi
       else
          echo "download successly" >>$dir/logs/StarBigwik"$i".log
          echo -e "###########Fastq-dump $srr begine at $(date)#################\n" >> $dir/logs/StarBigwik"$i".log
        fi
    fi   

 
   if [ ! -e $dir/$FQsubDIR/"$srr"_1.fastq.gz ] 
   then
      fastq-dump --skip-technical  --split-files --gzip --outdir "$dir"/"$FQsubDIR" "$Srr_dir"/"$srr".sra
    fi
    
    echo -e "################$Info Fatsq-dump compeleted at $(date)#######################\n" >> $dir/logs/StarBigwik"$i".log    
     
    if [ -e  "$dir"/"$FQsubDIR"/"$srr"_1.fastq.gz ] 
    then 
        rm $Srr_dir/$srr.sra
    fi 
    
    cd "$dir"/"$FQsubDIR"
    echo -e "$srr\t" >> $dir/filesizeOf.txt
    du -cm --max-depth=1 $dir/$FQsubDIR/"$srr"*.fastq.gz >>$dir/filesizeOf.txt
    cat $dir/$FQsubDIR/"$srr"_1.fastq.gz >>$dir/$FQsubDIR/"$filename"_1.fastq.gz
    cat $dir/$FQsubDIR/"$srr"_2.fastq.gz >>$dir/$FQsubDIR/"$filename"_2.fastq.gz
    rm $dir/$FQsubDIR/"$srr"*
   
 done
fi    
    
  
  if [ -s $dir/$FQsubDIR/"$filename"_1.fastq.gz ] || [ -f "$dir"/$Qc/"$filename"_1"$FQsuffix" ]
  then 
    echo -e "$Info\t$sfilename" >>$dir/filesizeOf.txt
    du -cm --max-depth=1 $dir/$FQsubDIR/"$filename"*.gz >>$dir/filesizeOf.txt
    
    InputTagDir+=" $filename"
   
   cd $dir
    #step 2 Use fastp to  reads to do QC
     echo -e "####################Now do fastp######################\n" >>$dir/logs/StarBigwik"$i""".log
     if [ -s $dir/$FQsubDIR/"$filename"_1.fastq.gz ]
     then 
  fastp -i $dir/$FQsubDIR/"$filename"_1.fastq.gz -I $dir/$FQsubDIR/"$filename"_2.fastq.gz  -o $dir/$Qc/$filename"_1"$FQsuffix -O $dir/$Qc/$filename"_2"$FQsuffix -h $dir/$Qc/$filename.html
    fi

  if [ -f "$dir"/$Qc/"$filename"_1"$FQsuffix" ]
   then 
        rm $dir/$FQsubDIR/"$filename"*.gz
  fi

     # Step 3 use star to do mapping 
   echo -e "now is mapping $Info and its sample name is $filename\n" >>$dir/logs/StarBigwik"$i".log
  STAR --runThreadN 8 --limitBAMsortRAM 40309881964 --genomeDir $index  --readFilesCommand zcat --readFilesIn $dir/$Qc/"$filename"_1"$FQsuffix" $dir/$Qc/"$filename"_2"$FQsuffix" --outFileNamePrefix $dir/$package/$filename --outSAMunmapped Within --outFilterType BySJout --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outWigType bedGraph
    
   if [ -f "$dir"/"$package"/"$filename"*.bam ]
    then
    rm "$dir"/$Qc/"$filename"*.gz
    rm $dir/$package/"$filename"[a-zA-Z]*.out.bg
    rm -rf $dir/$package/"$filename"*STAR*
    rm $dir/$package/"$filename"*SJ.out.tab
 fi    
  
   # Step 4 now lets get the unique mapped reads and the mapped reads count
  echo -e "get the unique reads from the .bam fil\e\n" >>$dir/logs/StarBigwik"$i".log
  samtools view -bq 30 -o $dir/3_Unique/$filename.sorted.unique.bam "$dir"/"$package"/"$filename"Align*.bam
 
  
  readDepth=`samtools view -c -F 260 $dir/3_Unique/"$filename"*.bam`
  wig="wig/$filename"
  wig+=".wig"
  wigfile+="${wig},"
  echo -e "$wig\t$readDepth" >>$dir/$depth 
  echo -e "get the wiggle file from unique .bam file" >>$dir/logs/StarBigwik"$i".log
  bedtools genomecov -ibam "$dir"/3_Unique/"$filename"*.bam -bga -split -trackline > $dir/wig/$filename.wig
  

    ##step 5 use homer to creat TagDirectory
    echo -e "now make TagDir for $Info and Its sample is $filename\n" >> $dir/logs/StarBigwik"$i".log 
$homer/bin/makeTagDirectory "$TagDir"/"$filename"/ "$dir"/3_Unique/"$filename".sorted.unique.bam
rm $dir/$package/"$filename"Align*.bam   
 ##step 6 use homer to creat bigwik file
$homer/bin/makeUCSCfile "$TagDir"/"$filename" -bigWig $chr_size -fsize 1e20 -style rnaseq -strand both -o "$dir"/"$HomerOut"/bigwig/"$filename".bigwig
   
   echo -e "now running Stringtie for $Info and Its sample is $filename\n" >> $dir/logs/StarBigwik"$i".log
  /project/guan/chunjie/Software/stringtie-1.3.4d.Linux_x86_64/stringtie "$dir"/3_Unique/"$filename".sorted.unique.bam -G $gtf_file -e -A $dir/7_Stringtie/$filename.tab -o $dir/7_Stringtie/$filename.gtf 
 echo -e "The process of $Info whose sample name is $filename is done at $(date) \n" >> $dir/logs/StarBigwik"$i".log
else
continue 1
fi

```


