#!/bin/bash

### Bash script for mapping trimmed reads from the Convergent Evolution Project to a reference genome. It's assumed that the FastQ file has been trimmed before mapping.
### James Reeve - University of Calgary
### 28/11/2018
### Population: NsNUn Sequencing lane 1

# Filepaths
DATA_TRIM=/data/not_backup/james/Con_evo_project/trim-data
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

BWA=/data/programs/bwa-0.7.12/bwa
SAMTOOLS=/data/programs/samtools-1.9/samtools
PICARD=/data/programs/picard.jar

## Indexing Reference Genome
#$BWA index -a bwtsw $DATA_REF/ninespine_genome_draft >> $LOG_PATH/Ns.bwa_index.log 2>&1

## Read mapping with BWA-MEM
$BWA mem -t 4 -M $DATA_REF/ninespine_genome_draft \
	$DATA_TRIM/HI.4883.001.NsNUn_1P \
	$DATA_TRIM/HI.4883.001.NsNUn_2P \
	| $SAMTOOLS sort -o $DATA_BAM/HI.4883.001.NsNUn.bam >> $LOG_PATH/NsNUn.1.bwa_mem.log 2>&1

## Sort and add read group information (needed for GATK to work)
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD AddOrReplaceReadGroups \
  I=$DATA_BAM/HI.4883.001.NsNUn.bam \
  O=$DATA_BAM/HI.4883.001.NsNUn_sort.bam \
  RGID=Index_28 \
  RGLB=MPS12343449-C07 \
  RGPL=illumina \
  RGPU=HI.4883.001 \
  RGSM=NsNUn \
  RGCN=Genome_Quebec \
  RGDS=Pooled_sample_46_fish \
  RGDT=2018-11-12 \
  RGPI=150 >> $LOG_PATH/NsNUn.1.ReadInfo.log 2>&1

## Mark duplicates
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD MarkDuplicates \
  I=$DATA_BAM/HI.4883.001.NsNUn_sort.bam \
  O=$DATA_BAM/HI.4883.001.NsNUn_DupMark.bam \
  M=$DATA_HOME/HI.4883.001.NsNUn_duplicate_metrics.txt \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
  CREATE_INDEX=true >> $LOG_PATH/NsNUn.1.MarkDuplicates.log 2>&1

## Error Check
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
  I=$DATA_BAM/HI.4883.001.NsNUn_DupMark.bam \
  MAX_OPEN_TEMP_FILES=1000 \
  MODE=SUMMARY >> $LOG_PATH/NsNUn.1.ValidateSamFile.log 2>&1
