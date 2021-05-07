############## Calculate Fst for Convergent Evolution Project ###########################
### The following code uses PoolFstat (Hivert et al. 2018) to calculate Fst between two populations of threespine stickleback
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
Ts.PP_pooldata <- popsync2pooldata(sync.file = "/Volumes/Con_evo_backup/vcf-files/Ts.sync", poolsizes = c(102, 104), 
                                   poolnames = c("TsAK", "TsOR"), min.cov.per.pool = 50, min.maf = 0.05)
Ts.PP_Fst <- computeFST(Ts.PP_pooldata)

# Calculate Fst from VarScan vcf
Ts.VS_pooldata <- vcf2pooldata(vcf.file = "/Volumes/Con_evo_backup/vcf-files/Ts.VarScan.vcf", poolsizes = c(102, 104), 
                               poolnames = c("TsAK", "TsOR"), min.cov.per.pool = 50, min.maf = 0.05)
Ts.VS_Fst <- computeFST(Ts.VS_pooldata)

# Extract Fst, CHROM and POS into a dataframe
Ts.PP <- data.frame("CHROM" = Ts.PP_pooldata@snp.info[,1], 
                    "POS" = as.numeric(Ts.PP_pooldata@snp.info[,2]), 
                    "Fst" = as.numeric(Ts.PP_Fst[["snp.FST"]]))

Ts.VS <- data.frame("CHROM" = Ts.VS_pooldata@snp.info[,1], 
                    "POS" = as.numeric(Ts.VS_pooldata@snp.info[,2]), 
                    "Fst" = as.numeric(Ts.VS_Fst[["snp.FST"]]))

# Rename Chromosomes to avoid issues with roman numeral sorting
Chrom_List <- c("chr1", "chr2", "chr3", "chr4", "chr9",
                "chrM", "chrUn", "chr5", "chr6", "chr7",
                "chr8", "chr10", "chr11", "chr12", "chr13",
                "chr14", "chr19", "chr15", "chr16", "chr17",
                "chr18", "chr20", "chr21")

chromosome.rename <- function(data){
  data$CHROM <- as.factor(data$CHROM) #convert to factor (more efficent)
  levels(data$CHROM) <- Chrom_List #reasign factor levels
  data$CHROM <- as.character(data$CHROM) #convert back to character
  return(data)
}

Ts.PP <- chromosome.rename(Ts.PP)
Ts.VS <- chromosome.rename(Ts.VS)

# Remove mitochondira and unmapped SNPs
Ts.PP <- Ts.PP[Ts.PP$CHROM != "chrM" & Ts.PP$CHROM != "chrUn", ]
Ts.VS <- Ts.VS[Ts.VS$CHROM != "chrM" & Ts.VS$CHROM != "chrUn", ]

# Order by Chromosome
Ts.PP <- Ts.PP[order(Ts.PP$CHROM), ]
Ts.VS <- Ts.VS[order(Ts.VS$CHROM), ]

# Replace negative Fst with 0
Ts.PP[Ts.PP$Fst < 0, "Fst"] <- 0
Ts.VS[Ts.VS$Fst < 0, "Fst"] <- 0

# Add column that combines CHORM and POS
Ts.PP$LOC <- paste(Ts.PP$CHROM, Ts.PP$POS, sep = ":")
Ts.VS$LOC <- paste(Ts.VS$CHROM, Ts.VS$POS, sep = ":")

# Find Intersect of PoPoolation and VarScan
Ts.concensus <- merge(Ts.PP, Ts.VS, by = "LOC")

Ts.concensus <- Ts.concensus[, -c(5:6)]
colnames(Ts.concensus) <- c("LOC", "CHROM", "POS", "Fst.PP", "Fst.VS")

Ts.concensus$deltaFst <- Ts.concensus$Fst.PP - Ts.concensus$Fst.VS
Ts.concensus <- Ts.concensus[with(Ts.concensus, order(CHROM, POS)),]

# Save Fst datasets
write.table(Ts.PP, "/Volumes/Con_evo_backup/Fst-table/Ts.PoPoolation.Fst.txt", row.names = FALSE, quote = FALSE)
write.table(Ts.VS, "/Volumes/Con_evo_backup/Fst-table/Ts.VarScan.Fst.txt", row.names = FALSE, quote = FALSE)
write.table(Ts.concensus, "/Volumes/Con_evo_backup/Fst-table/Ts.Fst.txt", row.names = FALSE, quote = FALSE)
