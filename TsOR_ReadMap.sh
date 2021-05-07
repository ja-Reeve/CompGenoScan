#!/bin/bash

### Bash script for mapping trimmed reads from the Convergent Evolution Project to a reference genome. It's assumed that the FastQ file has been trimmed before mapping.
### James Reeve - University of Calgary
### 28/11/2018
### Population: TsOR

# Filepaths
DATA_TRIM=/data/not_backup/james/Con_evo_project/trim-data
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

BWA=/data/programs/bwa-0.7.12/bwa
SAMTOOLS=/data/programs/samtools-1.3.1/samtools
PICARD=/data/programs/picard.jar
Trimmomatic=/data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar

## Read trimming
java -jar $Trimmomatic PE -threads 18 -phred33 \
 $DATA_FASTQ/TsOR_R1.fastq.gz $DATA_FASTQ/TsOR_R2.fastq.gz \
 -baseout $DATA_TRIM/TsOR \
 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:120

## Indexing Reference Genome
#$BWA index -a bwtsw $DATA_REF/threespine_genome.fa >> $LOG_PATH/TsOR.bwa_index.log 2>&1

## Read mapping with BWA-MEM
#$BWA mem -t 4 -M $DATA_REF/threespine_genome.fa \
#	$DATA_TRIM/TsOR_1P \
#	$DATA_TRIM/TsOR_2P \
#	| $SAMTOOLS sort -o TsOR.bam >> $LOG_PATH/TsOR.bwa_mem.log 2>&1

## Sort and add read group information (needed for GATK to work)
#java -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD AddOrReplaceReadGroups \
#  I=$DATA_BAM/TsOR.bam \
#  O=$DATA_BAM/TsOR_sort.bam \
#  SORT_ORDER=coordinate \
#  RGID=IDT_64 \
#  RGLB=MPS12343449-F07 \
#  RGPL=illumina \
#  RGPU=HI.4883.003 \
#  RGSM=TsOR \
#  RGCN=Genome_Quebec \
#  RGDS=Pooled_sample_52_fish \
#  RGDT=2018-11-12 \
#  RGPI=150 >> $LOG_PATH/TsOR.sorting.log 2>&1

## Mark duplicates
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD MarkDuplicates \
  I=$DATA_BAM/TsOR_sort.bam \
  O=$DATA_BAM/TsOR_DupMark.bam \
  M=$DATA_HOME/TsOR_duplicate_metrics.txt \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  CREATE_INDEX=true >> $LOG_PATH/TsOR.MarkDuplicates.log 2>&1

## Error Check
java -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
      I=$DATA_BAM/TsOR_DupMark.bam \
      MODE=SUMMARY >> $LOG_PATH/TsOR.ValidateSamFile.log 2>&1
