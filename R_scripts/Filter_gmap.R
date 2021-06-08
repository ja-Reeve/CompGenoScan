############## Filtering GMAP results ###########################
### The following code filters the results of a GMAP alignment between threespine genes and the
### draft ninespine genome (Nelson and Cresko 2017)
### James Reeve - University of Calgary
### 09/07/2019

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)
library("dplyr")

###1: Remove header and lines between genes  (bash)
#grep -v '^#' ga_pungitius_unfiltered.gff3 > tmp_ga.gff3

###2: Upload data
GFF3 <- read.table("/Volumes/Con_evo_backup/GFF/tmp_ga.gff3", sep = "\t")
colnames(GFF3) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

###3: Filter out non-mRNA enteries (only ones containing quality info)
GFF3_mRNA <- GFF3[GFF3$type =="mRNA",]

###4: Split atribute column
tmp <- strsplit(GFF3_mRNA$attributes, ";")
GFF3_mRNA$ID <- sapply(tmp, `[[`, 1)
GFF3_mRNA$cov <- as.numeric(sapply(strsplit(sapply(tmp, `[[`, 4), "[=]"),`[[`, 2))
GFF3_mRNA$iden <- as.numeric(sapply(strsplit(sapply(tmp, `[[`, 5), "[=]"),`[[`, 2))

GFF3_mRNA <- GFF3_mRNA[,-9]

###5: Remove path < 1
tmp2 <- GFF3_mRNA[sapply(strsplit(GFF3_mRNA$ID, ".mrna"), `[[`, 2) == 1, ]

###6: Filter out coverate < 80 and % identity < 90%
tmp3 <- tmp2[tmp2$cov >= 80 & tmp2$iden >= 90.0, ] ### Note: changing these parapmeters changes number of genes identified

####7: Compare filtered positions to GFF
tmp4 <- GFF3[GFF3$seqid %in% tmp2$seqid, ]
tmp5 <- tmp4[tmp4$type == "gene",]
tmp6 <- tmp5[sapply(strsplit(
  sapply(strsplit(tmp5$attributes, ";"), `[[`, 1), 
  ".path"), `[[`, 2) == 1,]
tmp6$attributes <- sapply(strsplit(tmp6$attributes, ".path1;"), `[[`, 1)

###8: Identify 1:many or many:1
genes <- unique(tmp6$attributes)

#For loop counting number of rows with the same gene ID
one2one.tmp <- list()
one2many.tmp <- list()

for (i in 1:length(genes)) {
  dat <- tmp6[tmp6$attributes == genes[i],]
  
  if(nrow(dat) == 1) {
    one2one.tmp[[i]] <- dat #single entry = 1:1
  } else {
    one2many.tmp[[i]] <- dat #multiple entry = 1:many
  }
}

one2one <- do.call("rbind", one2one.tmp)
one2many <- do.call("rbind", one2many.tmp)

# No one2many found

one2one$t2 <- sapply(strsplit(one2one$attributes, "[.]"), `[[`, 2)
one2one <- one2one[one2one$t2 == "t1", ]
one2one$t2 <- NULL

###9: Save filtered data
write.table(one2one, "/Volumes/Con_evo_backup/GFF/ga_pungitius_filtered.gff3", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################ Appendix: Exploring how filters effect data #######################
#I: Test retention for different coverage thresholds
CovThresh <- seq(0, 100, 0.1)
tmp <- sapply(CovThresh, function(th){
  nrow(one2many[one2many$coverage >= th,])
})

CovDat <- data.frame("coverage.threshold" = CovThresh, "Count" = tmp)
ggplot(CovDat, aes(coverage.threshold, Count/1000))+
  geom_point(size = 0.2, colour = "blue")+
  labs(x = "Coverage threshold", y = "Count (1000 genes)")+
  lims(x = c(0, 100),
       y = c(0,25))+
  theme_classic()

#II: Test retention for different percentage Identifty
IDThresh <- seq(0, 100, 0.1)
tmp <- sapply(IDThresh, function(th){
  nrow(one2many[one2many$identity >= th,])
})

IDDat <- data.frame("identity.threshold" = CovThresh, "Count" = tmp)
ggplot(IDDat, aes(identity.threshold, Count/1000))+
  geom_point(size = 0.2, colour = "darkgreen")+
  labs(x = "Identity threshold", y = "Count (1000 genes)")+
  lims(x = c(0, 100),
       y = c(0,25))+
  theme_classic()
