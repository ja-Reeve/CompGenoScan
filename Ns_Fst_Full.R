############## Calculate Fst for Convergent Evolution Project ###########################
### The following code uses PoolFstat (Hivert et al. 2018) to calculate Fst between four populations of ninespine stickleback
### All Fst values are calculated as the intesect of SNP calls made by PoPoolation and VarScan
### James Reeve - University of Calgary
### 05/03/2019

### Note: when using VarScan make sure to use 'min-coverage' not 'mincoverage' as a parameter ###

### Packages
library(poolfstat)
library(ggplot2)

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)

# Calcualte Fst from PoPoolation sync file
Ns.PP_pooldata <- popsync2pooldata(sync.file = "/Volumes/Con_evo_backup/vcf-files/Ns.sync", poolsizes = c(92, 84, 60, 82), 
                                   poolnames = c("NsNUn", "NsNUd", "NsABk", "NsABm"), min.cov.per.pool = 50, min.maf = 0.05)
Ns.PP_Fst <- computePairwiseFSTmatrix(Ns.PP_pooldata, output.snp.values = TRUE)

Ns.VS_pooldata <- vcf2pooldata(vcf.file = "/Volumes/Con_evo_backup/vcf-files/Ns_FULL.VarScan.vcf", poolsizes = c(92, 84, 60, 82), 
                            poolnames = c("NsNUn", "NsNUd", "NsABk", "NsABm"), min.cov.per.pool = 50, min.maf = 0.05)
Ns.VS_Fst <- computePairwiseFSTmatrix(Ns.VS_pooldata, output.snp.values = TRUE)

# Extract Fst, CHROM and POS into a dataframe
Ns.VS_pooldata <- Ns_pooldata

Ns.VS <- data.frame("CHROM" = Ns.VS_pooldata@snp.info[,1], 
                 "POS" = as.numeric(Ns.VS_pooldata@snp.info[,2]),
                 "NUn_vs_NUd" = as.numeric(Ns.VS_Fst$PairwiseSnpFST[,1]),
                 "NUn_vs_ABk" = as.numeric(Ns.VS_Fst$PairwiseSnpFST[,2]),
                 "NUn_vs_ABm" = as.numeric(Ns.VS_Fst$PairwiseSnpFST[,3]),
                 "NUd_vs_ABk" = as.numeric(Ns.VS_Fst$PairwiseSnpFST[,4]),
                 "NUd_vs_ABm" = as.numeric(Ns.VS_Fst$PairwiseSnpFST[,5])
)

# Seperate rows with NAs
Ns.PP2 <- na.omit(Ns.PP)
Ns.VS2 <- na.omit(Ns.VS)

# Remove Fst = 1
Ns.PP3 <- Ns.PP2[Ns.PP2$NUn_vs_NUd != 1 & Ns.PP2$NUn_vs_ABk != 1 & Ns.PP2$NUn_vs_ABm != 1 & Ns.PP2$NUd_vs_ABk != 1 & Ns.PP2$NUd_vs_ABm != 1,]
Ns.VS3 <- Ns.VS2[Ns.VS2$NUn_vs_NUd != 1 & Ns.VS2$NUn_vs_ABk != 1 & Ns.VS2$NUn_vs_ABm != 1 & Ns.VS2$NUd_vs_ABk != 1 & Ns.VS2$NUd_vs_ABm != 1,]

# Change negative Fst to 0s
Ns.PP3[Ns.PP3$NUn_vs_NUd < 0, 3] <- 0
Ns.PP3[Ns.PP3$NUn_vs_ABk < 0, 4] <- 0
Ns.PP3[Ns.PP3$NUn_vs_ABm < 0, 5] <- 0
Ns.PP3[Ns.PP3$NUd_vs_ABk < 0, 6] <- 0
Ns.PP3[Ns.PP3$NUd_vs_ABm < 0, 7] <- 0

Ns.VS3[Ns.VS3$NUn_vs_NUd < 0, 3] <- 0
Ns.VS3[Ns.VS3$NUn_vs_ABk < 0, 4] <- 0
Ns.VS3[Ns.VS3$NUn_vs_ABm < 0, 5] <- 0
Ns.VS3[Ns.VS3$NUd_vs_ABk < 0, 6] <- 0
Ns.VS3[Ns.VS3$NUd_vs_ABm < 0, 7] <- 0

# Add column that combines CHORM and POS
Ns.PP3$LOC <- paste(Ns.PP3$CHROM, Ns.PP3$POS, sep = ":")
Ns.VS3$LOC <- paste(Ns.VS3$CHROM, Ns.VS3$POS, sep = ":")

# Find Intersect of PoPoolation and VarScan
Ns.concensus <- merge(Ns.PP3, Ns.VS3, by = "LOC")

Ns.concensus <- Ns.concensus[, -c(9:10)]
colnames(Ns.concensus) <- c("LOC", "CHROM", "POS", 
                            "NUn_vs_NUd.PP", "NUn_vs_ABk.PP", "NUn_vs_ABm.PP", "NUd_vs_ABk.PP", "NUd_vs_ABm.PP", 
                            "NUn_vs_NUd.VS", "NUn_vs_ABk.VS", "NUn_vs_ABm.VS", "NUd_vs_ABk.VS", "NUd_vs_ABm.VS") #edit#

Ns.concensus$deltaNUn_vs_NUd <- Ns.concensus$NUn_vs_NUd.PP - Ns.concensus$NUn_vs_NUd.VS
Ns.concensus$deltaNUn_vs_ABk <- Ns.concensus$NUn_vs_ABk.PP - Ns.concensus$NUn_vs_ABk.VS
Ns.concensus$deltaNUn_vs_ABm <- Ns.concensus$NUn_vs_ABm.PP - Ns.concensus$NUn_vs_ABm.VS
Ns.concensus$deltaNUd_vs_ABk <- Ns.concensus$NUd_vs_ABk.PP - Ns.concensus$NUd_vs_ABk.VS
Ns.concensus$deltaNUd_vs_ABm <- Ns.concensus$NUd_vs_ABm.PP - Ns.concensus$NUd_vs_ABm.VS

Ns.concensus <- Ns.concensus[with(Ns.concensus, order(CHROM, POS)),]

# Save Fst datasets
write.table(Ns.PP3, "/Volumes/Con_evo_backup/Fst-table/Ns.PoPoolation.Fst.txt", row.names = FALSE)
write.table(Ns.VS3, "/Volumes/Con_evo_backup/Fst-table/Ns.VarScan.Fst.txt", row.names = FALSE)
write.table(Ns.concensus, "/Volumes/Con_evo_backup/Fst-table/Ns.Fst.txt", row.names = FALSE)

# Plot callers against each other
Fst.comp.plot <- function(pop_pair.VS, pop_pair.PP, pop_pair){
  ggplot(Ns.concensus, aes(pop_pair.VS, pop_pair.PP))+
    geom_point(size = 0.1, colour = "navy")+
    geom_abline(slope = 1, intercept = 0, lty = 1)+
    geom_vline(xintercept = quantile(pop_pair.VS, 0.999), lty = 3)+
    geom_hline(yintercept = quantile(pop_pair.PP, 0.999), lty = 3)+
    labs(x = expression(VarScan~F[ST]), y = expression(PoPoolation~F[ST]), 
         title = paste("Ninespine stickleback:", pop_pair))
}
pt.NUn_vs_NUd <- Fst.comp.plot(Ns.concensus$NUn_vs_NUd.VS, Ns.concensus$NUn_vs_NUd.PP, "NUn_vs_NUd")
pt.NUn_vs_ABk <- Fst.comp.plot(Ns.concensus$NUn_vs_ABk.VS, Ns.concensus$NUn_vs_ABk.PP, "NUn_vs_ABk")
pt.NUn_vs_ABm <- Fst.comp.plot(Ns.concensus$NUn_vs_ABm.VS, Ns.concensus$NUn_vs_ABm.PP, "NUn_vs_ABm")
pt.NUd_vs_ABk <- Fst.comp.plot(Ns.concensus$NUd_vs_ABk.VS, Ns.concensus$NUd_vs_ABk.PP, "NUd_vs_ABk")
pt.NUd_vs_ABm <- Fst.comp.plot(Ns.concensus$NUd_vs_ABm.VS, Ns.concensus$NUd_vs_ABm.PP, "NUd_vs_ABm")

library("gridExtra")
grid.arrange(pt.NUn_vs_ABk, pt.NUn_vs_ABm, pt.NUd_vs_ABk, pt.NUd_vs_ABm)

################ Testing area: Create Concensus for Fst ##################################
### Test 1. Average Fst (Treat N-S as replicates)
Ns_consensus <- rowMeans(Ns3[,3:6])
#ggplot()+
#  geom_density(aes(x = Ns_consensus))

Ns4 <- data.frame("CHROM" = Ns3$CHROM,
                  "POS" = Ns3$POS, 
                  "Fst" = rowMeans(Ns3[,3:6]), 
                  "LOC" = Ns3$LOC)

manhatten.plot(Ns4, Ns4$Fst, "Ninespine: Concensus (Mean)")

### Re calculate Fst independently for each pairing, and find mean
Ns_pooldata2 <- new("pooldata")

Ns_pooldata2@npools=2

# Add modified refallele.readcount
data.Y <- array(c(rowSums(Ns_pooldata@refallele.readcount[,1:2]),
                  rowSums(Ns_pooldata@refallele.readcount[,3:4])), 
                dim = c(Ns_pooldata@nsnp, 2))

tmp1<-as.data.frame(data.Y)
tmp2<-rownames(tmp1[(rowSums(tmp1) - tmp1[,1]) != 0,])

data.Y2<-data.Y[as.numeric(tmp2),]

Ns_pooldata2@refallele.readcount=data.Y2
rm(data.Y, tmp1, data.Y2)

# Add modified readcoverage
data.N <- array(c(rowSums(Ns_pooldata@readcoverage[,1:2]),
                  rowSums(Ns_pooldata@readcoverage[,3:4])), 
                dim = c(Ns_pooldata@nsnp, 2))

data.N2<-data.N[as.numeric(tmp2),]
Ns_pooldata2@readcoverage=data.N2
rm(data.N, data.N2)


Ns_pooldata2@snp.info=Ns_pooldata@snp.info[as.numeric(tmp2),]
rm(tmp2)

Ns_pooldata2@nsnp=nrow(Ns_pooldata2@snp.info)

Ns_pooldata2@poolsizes=as.numeric(c(176, 142))
Ns_pooldata2@poolnames=as.character(c("NsNU", "NsAB"))

Ns_Fst2 <- computeFST(Ns_pooldata2)

# Extract Fst, CHROM and POS into a dataframe
Ns5 <- data.frame("CHROM" = Ns_pooldata2@snp.info[,1], 
                  "POS" = as.numeric(Ns_pooldata2@snp.info[,2]), 
                  "Fst" = as.numeric(Ns_Fst2$snp.FST)
)

# Seperate rows with NAs
Ns6 <- na.omit(Ns5)

# Remove Fst = 1
Ns7 <- Ns6[Ns6$Fst != 1,]

# Change negative Fst to 0s
Ns7[Ns7$Fst < 0, 3] <- 0

# Add column that combines CHORM and POS
Ns7$LOC <- paste(Ns7$CHROM, Ns7$POS, sep = ":")

# Manhatten Plot
manhatten.plot(Ns7, Ns7$Fst, "Ninespine: Mean(north-south) Fst")