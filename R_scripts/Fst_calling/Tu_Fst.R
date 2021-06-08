############## Calculate Fst for Convergent Evolution Project ###########################
### The following code uses PoolFstat (Hivert et al. 2018) to calculate Fst between two populations of tubesnouts
### All Fst values are calculated as the intesect of SNP calls made by PoPoolation and VarScan
### James Reeve - University of Calgary
### 20/02/2019

### Packages
library(poolfstat)
library(ggplot2)

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)

# Calcualte Fst from PoPoolation sync file
Tu.PP_pooldata <- popsync2pooldata(sync.file = "/Volumes/Con_evo_backup/vcf-files/Tu.sync", poolsizes = c(88, 100), 
                                   poolnames = c("TuAK", "TuBC"), min.cov.per.pool = 50, min.maf = 0.05)
Tu.PP_Fst <- computeFST(Tu.PP_pooldata)

# Calculate Fst from VarScan vcf
Tu.VS_pooldata <- vcf2pooldata(vcf.file = "/Volumes/Con_evo_backup/vcf-files/Tu.VarScan.vcf", poolsizes = c(88, 100), 
                               poolnames = c("TuAK", "TuBC"), min.cov.per.pool = 50, min.maf = 0.05)
Tu.VS_Fst<- computeFST(Tu.VS_pooldata)

# Extract Fst, CHROM and POS into a dataframe
Tu.PP <- data.frame("CHROM" = Tu.PP_pooldata@snp.info[,1], 
                    "POS" = as.numeric(Tu.PP_pooldata@snp.info[,2]), 
                    "Fst" = as.numeric(Tu.PP_Fst[["snp.FST"]]))

Tu.VS <- data.frame("CHROM" = Tu.VS_pooldata@snp.info[,1], 
                    "POS" = as.numeric(Tu.VS_pooldata@snp.info[,2]), 
                    "Fst" = as.numeric(Tu.VS_Fst[["snp.FST"]]))

# Rename Chromosomes
chrom_list <- read.table("/Volumes/Con_evo_backup/chrom_list/Tu_chr_list.txt")$V1
new_chrom_list <- paste0("chr", 1:23)
for (i in 1:length(new_chrom_list)) {
  Tu.PP[Tu.PP$CHROM == chrom_list[i], 1] <- new_chrom_list[i]
  Tu.VS[Tu.VS$CHROM == chrom_list[i], 1] <- new_chrom_list[i]
}; rm(i)

# Order by Chromosome
Tu.PP <- Tu.PP[order(Tu.PP$CHROM), ]
Tu.VS <- Tu.VS[order(Tu.VS$CHROM), ]

# Replace negative Fst with 0
Tu.PP[Tu.PP$Fst < 0, "Fst"] <- 0
Tu.VS[Tu.VS$Fst < 0, "Fst"] <- 0

# Add column that combines CHORM and POS
Tu.PP$LOC <- paste(Tu.PP$CHROM, Tu.PP$POS, sep = ":")
Tu.VS$LOC <- paste(Tu.VS$CHROM, Tu.VS$POS, sep = ":")

# Find Intersect of PoPoolation and VarScan
Tu.concensus <- merge(Tu.PP, Tu.VS, by = "LOC")

Tu.concensus <- Tu.concensus[, -c(5:6)]
colnames(Tu.concensus) <- c("LOC", "CHROM", "POS", "Fst.PP", "Fst.VS")

Tu.concensus$deltaFst <- Tu.concensus$Fst.PP - Tu.concensus$Fst.VS
Tu.concensus <- Tu.concensus[with(Tu.concensus, order(CHROM, POS)),]

# Save Fst datasets
write.table(Tu.PP, "/Volumes/Con_evo_backup/Fst-table/Tu.PoPoolation.Fst.txt", row.names = FALSE)
write.table(Tu.VS, "/Volumes/Con_evo_backup/Fst-table/Tu.VarScan.Fst.txt", row.names = FALSE)
write.table(Tu.concensus, "/Volumes/Con_evo_backup/Fst-table/Tu.Fst.txt", row.names = FALSE)
