#!/bin/bash

### Bash script for mapping trimmed reads from the Convergent Evolution Project to a reference genome. It's assumed that the FastQ file has been trimmed before mapping.
### James Reeve - University of Calgary
### 23/4/2018
### Population: NsABm Sequencing lane 1

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

## Read trimming
java -jar $Trimmomatic PE -threads 18 -phred33 \
 $DATA_FASTQ/HI.4883.001.NsABm_R1.fastq.gz $DATA_FASTQ/HI.4883.001.NsABm_R2.fastq.gz \
 -baseout $DATA_TRIM/HI.4883.001.NsABm \
 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:120

## Read mapping with BWA-MEM
$BWA mem -t 4 -M $DATA_REF/threespine_genome.fa \
        $DATA_TRIM/HI.4883.001.NsABm_1P \
        $DATA_TRIM/HI.4883.001.NsABm_2P \
        | $SAMTOOLS sort -o $DATA_BAM/HI.4883.001.NsABm.bam >> $LOG_PATH/NsABm.1.bwa_mem.log 2>&1
		
java -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD AddOrReplaceReadGroups \
  I=$DATA_BAM/HI.4883.001.NsABm.bam \
  O=$DATA_BAM/HI.4883.001.NsABm_sort.bam \
  RGID=Index_36 \
  RGLB=MPS12343462-C03 \
  RGPL=illumina \
  RGPU=HI.4883.001 \
  RGSM=NsABm \
  RGCN=Genome_Quebec \
  RGDS=Pooled_sample_41_fish \
  RGDT=2018-11-12 \
  RGPI=150 >> $LOG_PATH/NsAB.1.sorting.log 2>&1

## Mark duplicates
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD MarkDuplicates \
  I=$DATA_BAM/HI.4883.001.NsABm_sort.bam \
  O=$DATA_BAM/HI.4883.001.NsABm_DupMark.bam \
  M=$DATA_HOME/HI.4883.001.NsABm_duplicate_metrics.txt \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  CREATE_INDEX=true >> $LOG_PATH/NsABm.1.MarkDuplicates.log 2>&1

## Error Check
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
      I=$DATA_BAM/HI.4883.001.NsABm_DupMark.bam \
      MODE=SUMMARY >> $LOG_PATH/NsABm.1.ValidateSamFile.log 2>&1
