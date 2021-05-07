#!/bin/bash

### Bash script for mapping trimmed reads from the Convergent Evolution Project to a reference genome. It's assumed that the FastQ file has been trimmed before mapping.
### James Reeve - University of Calgary
### 23/4/2018
### Population: NsNUd Sequencing lane 1

# Filepaths
DATA_FASTQ=/data/fastq/james/Con_evo_project
DATA_TRIM=/data/not_backup/james/Con_evo_project/trim-data
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

BWA=/data/programs/bwa-0.7.12/bwa
SAMTOOLS=/data/programs/samtools-1.3.1/samtools
PICARD=/data/programs/picard.jar
Trimmomatic=/data/programs/Trimmomatic-0.36/trimmomatic-0.36.jar

### Read trimming
java -jar $Trimmomatic PE -threads 18 -phred33 \
 $DATA_FASTQ/HI.4883.001.NsNUd_R1.fastq.gz $DATA_FASTQ/HI.4883.001.NsNUd_R2.fastq.gz \
 -baseout $DATA_TRIM/HI.4883.001.NsNUd \
 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:120

## Indexing Reference Genome
#$BWA index -a bwtsw $DATA_REF/ninespine_genome_draft >> $LOG_PATH/Ns.bwa_index.log 2>&1

## Read mapping with BWA-MEM
$BWA mem -t 4 -M $DATA_REF/threespine_genome.fa \
        $DATA_TRIM/HI.4883.001.NsNUd_1P \
        $DATA_TRIM/HI.4883.001.NsNUd_2P \
        | $SAMTOOLS sort -o $DATA_BAM/HI.4883.001.NsNUd.bam >> $LOG_PATH/NsNUd.1.bwa_mem.log 2>&1
		
java -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD AddOrReplaceReadGroups \
  I=$DATA_BAM/HI.4883.001.NsNUd.bam \
  O=$DATA_BAM/HI.4883.001.NsNUd_sort.bam \
  SORT_ORDER=coordinate \
  RGID=Index_16 \
  RGLB=MPS12343449-B07 \
  RGPL=illumina \
  RGPU=HI.4883.001 \
  RGSM=NsNUd \
  RGCN=Genome_Quebec \
  RGDS=Pooled_sample_42_fish \
  RGDT=2018-11-12 \
  RGPI=150 >> $LOG_PATH/NsNUd.1.sorting.log 2>&1

## Mark duplicates
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD MarkDuplicates \
  I=$DATA_BAM/HI.4883.001.NsNUd_sort.bam \
  O=$DATA_BAM/HI.4883.001.NsNUd_DupMark.bam \
  M=$DATA_HOME/HI.4883.001.NsNUd_duplicate_metrics.txt \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  CREATE_INDEX=true >> $LOG_PATH/NsNUd.1.MarkDuplicates.log 2>&1

## Error Check
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
      I=$DATA_BAM/HI.4883.001.NsNUd_DupMark.bam \
      MODE=SUMMARY >> $LOG_PATH/NsNUd.1.ValidateSamFile.log 2>&1
