#!/bin/bash

### Bash script for calling SNPs for the Convergent Evolution Project. It's assumed that reads have been mapped and realigned around indels before running this script.
### James Reeve - University of Calgary
### 18/04/2019
### Species: Ninespine Stickleback (Pungitius pungitius)

# Filepaths
DATA_BAM=/data/aligned/james/Con_evo_project
DATA_REF=/data/home/james/genome
DATA_HOME=/data/home/james/Con_evo_project
LOG_PATH=/data/home/james/Con_evo_project/log-files

SAMTOOLS=/data/programs/samtools-1.9/samtools
VARSCAN=/data/programs/VarScan.v2.3.9.jar

## North-South 1: NUn vs ABk ##
## Create mpileup file
$SAMTOOLS mpileup -B -f $DATA_REF/ninespine_genome_draft \
$DATA_BAM/NsNUn_Realign.bam $DATA_BAM/NsABk_Realign.bam \
> $DATA_HOME/mpileups/Ns_NS1.mpileup

## SNP Call
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $VARSCAN \
mpileup2snp $DATA_HOME/mpileups/Ns_NS1.mpileup --output-vcf 1 --mincoverage 50 --min-avg-qual 20 \
> $DATA_HOME/vcf-files/Ns_NS1.VarScan.vcf

## North-South 2: NUn vs ABm ##
## Create mpileup file
$SAMTOOLS mpileup -B -f $DATA_REF/ninespine_genome_draft \
$DATA_BAM/NsNUn_Realign.bam $DATA_BAM/NsABm_Realign.bam \
> $DATA_HOME/mpileups/Ns_NS2.mpileup

## SNP Call
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $VARSCAN \
mpileup2snp $DATA_HOME/mpileups/Ns_NS2.mpileup --output-vcf 1 --mincoverage 50 --min-avg-qual 20 \
> $DATA_HOME/vcf-files/Ns_NS2.VarScan.vcf

## North-South 3: NUd vs ABk ##
## Create mpileup file
$SAMTOOLS mpileup -B -f $DATA_REF/ninespine_genome_draft \
$DATA_BAM/NsNUd_Realign.bam $DATA_BAM/NsABk_Realign.bam \
> $DATA_HOME/mpileups/Ns_NS3.mpileup

## SNP Call
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $VARSCAN \
mpileup2snp $DATA_HOME/mpileups/Ns_NS3.mpileup --output-vcf 1 --mincoverage 50 --min-avg-qual 20 \
> $DATA_HOME/vcf-files/Ns_NS3.VarScan.vcf

## North-South 4: NUd vs ABm ##
## Create mpileup file
$SAMTOOLS mpileup -B -f $DATA_REF/ninespine_genome_draft \
$DATA_BAM/NsNUd_Realign.bam $DATA_BAM/NsABm_Realign.bam \
> $DATA_HOME/mpileups/Ns_NS4.mpileup

## SNP Call
java -Xmx32G -Djava.io.tmpdir=/data/home/james/tmp -jar $VARSCAN \
mpileup2snp $DATA_HOME/mpileups/Ns_NS4.mpileup --output-vcf 1 --mincoverage 50 --min-avg-qual 20 \
> $DATA_HOME/vcf-files/Ns_NS4.VarScan.vcf
