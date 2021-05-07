#!/bin/bash

### Bash script for mapping trimmed reads from the Convergent Evolution Project to a reference genome. It's assumed that the FastQ file has been trimmed before mapping.
### James Reeve - University of Calgary
### 28/11/2018
### Population: TsAK

# Filepaths
DATA_FASTQ=/data/fastq/james/Con_evo_project
DATA_TRIM=/data/not_backup/james/Con_evo_project/trim-data
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

BWA=/data/programs/bwa-0.7.12/bwa
SAMTOOLS=/data/programs/samtools-1.9/samtools
PICARD=/data/programs/picard.jar
Trimmomatic=/data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar

## Read trimming
java -jar $Trimmomatic PE -threads 18 -phred33 \
 $DATA_FASTQ/TsAK_R1.fastq.gz $DATA_FASTQ/TsAK_R2.fastq.gz \
 -baseout $DATA_TRIM/TsAK \
 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:120

## Indexing Reference Genome
#$BWA index -a bwtsw $DATA_REF/threespine_genome.fa >> $LOG_PATH/TsAK.bwa_index.log 2>&1

## Read mapping with BWA-MEM
#$BWA mem -t 4 -M $DATA_REF/threespine_genome.fa \
#	$DATA_TRIM/TsAK_1P $DATA_TRIM/TsAK_2P \
#	| $SAMTOOLS sort -o $DATA_BAM/TsAK.bam >> $LOG_PATH/TsAK.bwa_mem.log 2>&1

## Sort and add read group information (needed for GATK to work)
java -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD AddOrReplaceReadGroups \
  I=$DATA_BAM/TsAK.bam \
  O=$DATA_BAM/TsAK_sort.bam \
  RGID=IDT_76 \
  RGLB=MPS12343449-G07 \
  RGPL=illumina \
  RGPU=HI.4883.003 \
  RGSM=TsAK \
  RGCN=Genome_Quebec \
  RGDS=Pooled_sample_51_fish \
  RGDT=2018-11-12 \
  RGPI=150 >> $LOG_PATH/TsAK.ReadInfo.log 2>&1

## Mark duplicates
java -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD MarkDuplicates \
  I=$DATA_BAM/TsAK_sort.bam \
  O=$DATA_BAM/TsAK_DupMark.bam \
  M=$DATA_HOME/TsAK_duplicate_metrics.txt \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  CREATE_INDEX=true >> $LOG_PATH/TsAK.MarkDuplicates.log 2>&1

## Error Check
java -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
      I=$DATA_BAM/TsAK_DupMark.bam \
      MODE=SUMMARY >> $LOG_PATH/TsAK.ValidateSamFile.log 2>&1
