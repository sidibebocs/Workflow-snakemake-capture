#!/bin/bash

#WARNING : this script MUST be used through statistical_analyses_bam_files.pl.

#This script will get some stats from reads data in order to build a graph that will able us to see the impact of the different steps of the pipeline.
#1 : Total amount of raw reads
#2 : Total amount of reads after trimming
#3 : Total amount of mapped reads
#4 : Total amount of mapped reads without diplicates
#5 : Total amount of reads properly mapped


unset finalBamPath
unset bedFile
graphe=0


echo -ne "Name\tTrimmed_reads\tUnmapped_reads\tMapped" > Mapping_stats.csv

while getopts "a:b:c:d:e" options
do
  case $options in
    a) #raw read
      rawReadsPath=$OPTARG
      ;;
    b) #raw bam
      rawBamPath=$OPTARG
      ;;
    c) #finalBamPath
      finalBamPath=$OPTARG
      echo -ne "\tOff target" >> Mapping_stats.csv
      ;;
    d) #bed
      bedFile=$OPTARG
      echo -ne "\tOn target" >> Mapping_stats.csv
      ;;
    e) #plot
      graphe=1
      ;;
  esac
done

echo -ne "\n" >> Mapping_stats.csv

echo -e "\nNow starting the analyses...\n"


for file in `ls $rawBamPath/*.bam` #être plus précis sur les nomsn petite doc à faire.
do
  name=${file##*/}
  name=${name%%.bam}
  name=${name%%_*}
  echo -e "Analysing sample : $name...\n"
  totalReads=$(zcat $rawReadsPath/$name\_*R1*.fastq.gz | wc -l) #FIXME : work only with zip files FIXME : Fail when $name is the prefixe of another name
  let "totalReads=$totalReads / 2" #FIXME : Add an option to handle not paired data.
  afterTrim=$(samtools view -c -F 0x900 $file)
  properlyMapped=$(samtools view -c -F 0x904 $file)
  echo -ne "$name\t$totalReads\t$afterTrim\t$properlyMapped" >> Mapping_stats.csv
  if [ -n "$finalBamPath" ]
  then
    withoutDuplicate=$(samtools view -c -F 0x904 $finalBamPath/$name.bam)
    echo -ne "\t$withoutDuplicate" >> Mapping_stats.csv
  fi
  if [ -n "$bedFile" ] && [ -n "$finalBamPath" ]
  then
    onTarget=$(samtools view -b -F 0x904 $finalBamPath/$name.bam | samtools view -c -L $bedFile)
    echo -ne "\t$onTarget" >> Mapping_stats.csv
  elif [ -n "$bedFile" ]
  then
    onTarget=$(samtools view -b -F 0x904 $file | samtools view -c -L $bedFile)
    echo -ne "\t$onTarget" >> Mapping_stats.csv
  fi
  echo -ne "\n" >> Mapping_stats.csv
done

position=${0%/*}

if [ $graphe == 1 ]
then
  Rscript $position/autoGraph.R Mapping_stats.csv $position/graphe_aire_sum_auto.R
fi
