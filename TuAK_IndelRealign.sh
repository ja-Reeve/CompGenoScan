#!/bin/bash

### Bash script for realigning mapped reads around indels. The input bam is produced by the ReadMap.sh script. This is part of the Convergent Evolution Project.
### James Reeve - University of Calgary
### 21/01/2019
### Population: TuAK

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
#-R $DATA_REF/tubesnout_genome.v2.fasta \
#-I $DATA_BAM/TuAK_DupMark.bam \
#-o $DATA_HOME/TuAK_Indel.vcf \
#-out_mode EMIT_ALL_SITES \
#-ploidy 88 \
#| $BCFTOOLS view -v indels $DATA_HOME/TuAK_Indel.vcf >> $LOG_PATH/TuAK.IndelCall.log 2>&1

## Create Target Sites for Realignment
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $GATK3 \
-T RealignerTargetCreator \
-R $DATA_REF/tubesnout_genome.v2.fasta \
-I $DATA_BAM/TuAK_DupMark.bam \
-o $DATA_HOME/TuAK_realignment_sites.list >> $LOG_PATH/TuAK.RealigmentSet.log 2>&1

## Realign Reads
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $GATK3 \
-T IndelRealigner \
-R $DATA_REF/tubesnout_genome.v2.fasta \
-targetIntervals $DATA_HOME/TuAK_realignment_sites.list \
-I $DATA_BAM/TuAK_DupMark.bam \
-o $DATA_BAM/TuAK_Realign.bam >> $LOG_PATH/TuAK.IndelRealigner.log 2>&1

## Error Check
java -XX:ParallelGCThreads=4 -Xmx8G -Djava.io.tmpdir=/data/home/james/tmp -jar $PICARD ValidateSamFile \
      I=$DATA_BAM/TuAK_Realign.bam \
      MODE=SUMMARY >> $LOG_PATH/TuAK.ValidateRealignBam.log 2>&1
