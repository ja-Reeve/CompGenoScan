#!/bin/bash

### Bash script for realigning mapped reads around indels. The input bam is produced by the ReadMap.sh script. This is part of the Convergent Evolution Project.
### James Reeve - University of Calgary
### 21/01/2019
### Population: NsNUd

# Filepaths
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

SAMTOOLS=/data/programs/samtools-1.9/samtools
PICARD=/data/programs/picard.jar
BCFTOOLS=/data/programs/bcftools-1.9/bcftoools
GATK3=/data/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar

## Call Indels from UnifiedGenotyper
#java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $GATK3 \
#-T UnifiedGenotyper \
#-R $DATA_REF/ninespine_genome_draft.fasta \
#-I $DATA_BAM/NsNUd_DupMark.bam \
#-o $DATA_HOME/NsNUd_Indel.vcf \
#-out_mode EMIT_ALL_SITES \
#-ploidy 84 \
#| $BCFTOOLS view -v indels $DATA_HOME/NsNUd_Indel.vcf >> $LOG_PATH/NsNUd.IndelCall.log 2>&1

## Create Target Sites for Realignment
java -XX:ParallelGCThreads=16 -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $GATK3 \
-T RealignerTargetCreator \
-R $DATA_REF/threespine_genome.fa \
-I $DATA_BAM/Old/NsNUd_DupMark.bam \
-o $DATA_HOME/NsNUd_realignment_sites2.list >> $LOG_PATH/NsNUd.RealigmentSet2.log 2>&1

## Realign Reads
java -XX:ParallelGCThreads=16 -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $GATK3 \
-T IndelRealigner \
-R $DATA_REF/threespine_genome.fa \
-targetIntervals $DATA_HOME/NsNUd_realignment_sites2.list \
-I $DATA_BAM/Old/NsNUd_DupMark.bam \
-o $DATA_BAM/NsNUd_Realign2.bam >> $LOG_PATH/NsNUd.IndelRealigner2.log 2>&1

## Error Check
java -XX:ParallelGCThreads=16 -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
      I=$DATA_BAM/NsNUd_Realign2.bam \
      MAX_OPEN_TEMP_FILES=1000 \
      MODE=SUMMARY >> $LOG_PATH/NsNUd.ValidateRealignBam2.log 2>&1
