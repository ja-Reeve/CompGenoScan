#!/bin/bash

### Bash script for calling SNPs for the Convergent Evolution Project. It's assumed that reads have been mapped and realigned around indels before running this script.
### James Reeve - University of Calgary
### 17/02/2019
### Species: Threespine Stickleback (Gasterosteus aculeatus)

# Filepaths
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

SAMTOOLS=/data/programs/samtools-1.9/samtools
VARSCAN=/data/programs/VarScan.v2.3.9.jar
POPOOLATION=/data/programs/popoolation2_1201

## Create mpileup file
#$SAMTOOLS mpileup -B -f $DATA_REF/threespine_genome.fa \
#$DATA_BAM/TsAK_Realign.bam $DATA_BAM/TsOR_Realign.bam > $DATA_HOME/mpileups/Ts.mpileup

## SNP Call VarScan
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $VARSCAN \
mpileup2snp $DATA_HOME/mpileups/Ts.mpileup --output-vcf 1 --min-coverage 50 --min-avg-qual 20 \
> $DATA_HOME/vcf-files/Ts.VarScan.vcf

## Create PoPoolation sync file
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $POPOOLATION/mpileup2sync.jar \
--input $DATA_HOME/mpileups/Ts.mpileup --output $DATA_HOME/vcf-files/Ts.sync \
--fastq-type sanger --min-qual 20

## Remove mpileup
#rm $DATA_HOME/mpileups/Ts.mpileup
