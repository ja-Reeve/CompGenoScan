#### Testing comparisons among Ts.He and Ns.He using 50Kb windows ####
### James Reeve
### 23/12/2019

rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)

### Download He data
Ts.He <- read.table("file:///G:/nucleotide_diversity/He_with_0s/Ts_Avg_Hexp_50Kwind.txt", header = TRUE)
Ns.He <- read.table("file:///G:/nucleotide_diversity/He_with_0s/Ns_Avg_Hexp_50Kwind.txt", header = TRUE)

Ts.He <- Ts.He[complete.cases(Ts.He),]
Ns.He <- Ns.He[complete.cases(Ns.He),]

### Download Fst data
Ts.Fst <- read.table("file:///G:/Fst-table/Ts.Fst.50K.txt", header = TRUE)
Ns.Fst <- read.table("file:///G:/Fst-table/Ns.Fst.50K.txt", header = TRUE)
# assign top candidate flag
Ts.TC <- read.table("file:///G:/top-candidate/Ts.TopCandidates.50Kwind.txt", header = TRUE)
Ts.TC$CAND <- Ts.TC$Observed_Outliers > Ts.TC$Expected_Outliers

Ts.Fst <- merge(Ts.Fst, Ts.TC, by.x = "WIND", by.y = "Window_Name")[,c(2:5,12,1)] 

rm(Ts.TC)

### Download gene annotations
Ts.gff <- read.table("file:///G:/GFF/Ts_gene_anno.txt", header = TRUE)
Ts.gff <- Ts.gff[Ts.gff$CHROM != "chrUn", ]
Ns.gff <- read.table("file:///G:/GFF/Ns_gene_anno.txt", header = TRUE)
# Rename gene header
Ts.gff$Gene <- paste("GAC", sapply(strsplit(Ts.gff$Gene, "[=]"), `[[`, 2), sep = "_")
Ns.gff$Gene <- paste("GAC", sapply(strsplit(Ns.gff$Gene, "[=]"), `[[`, 2), sep = "_")

### Update chromosome names
chrom_list_Ts <- read.table("file:///G:/chrom_list/Ts_chr_list.txt")$V1[-c(6:7)]
new_chrom_list_Ts <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8",
                       "chr10", "chr11", "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17",
                       "chr18", "chr20", "chr21")
for(i in 1:length(chrom_list_Ts)){
  Ts.gff[Ts.gff$CHROM == chrom_list_Ts[i], 1] <- new_chrom_list_Ts[i]
}; rm(i, chrom_list_Ts, new_chrom_list_Ts)

### Merge He and Fst data
Ts.Fst.He <- merge(Ts.He, Ts.Fst, by = "WIND")[,c(2:7,11:12,1)]
colnames(Ts.Fst.He) <- c("CHROM", "START_POS", "END_POS", "TsAK_He", "TsOR_He",
                         "Avg_He", "Avg_Fst", "CAND", "WINDOW")
Ns.Fst.He <- merge(Ns.He, Ns.Fst, by = "WIND")[,c(2:9,13:14,1)]
colnames(Ns.Fst.He) <- c("CHROM", "START_POS", "END_POS", "NsNUn_He", "NsNUd_He", "NsABk_He", "NsABm_He",
                         "Avg_He", "Avg_Fst_NS1", "Avg_Fst_NS2", "WINDOW")

# Remove NAs
Ts.Fst.He <- Ts.Fst.He[complete.cases(Ts.Fst.He),]
Ns.Fst.He <- Ns.Fst.He[complete.cases(Ns.Fst.He),]

### Assing window He to genes within windows ####

GeneAsign <- function(data, gene.annotation){
  tmp <- list()
  
  for(i in 1:nrow(data)){
    Chrom <- data[i, "CHROM"]
    Spos <- data[i, "START_POS"]
    Epos <- data[i, "END_POS"]
    
    dat <- gene.annotation[gene.annotation$CHROM == Chrom &
                             gene.annotation$Start_Pos >= Spos &
                             gene.annotation$End_Pos <= Epos, ]
    if(nrow(dat) == 0)next
    dat$North_He <- rep(data[i, 4], nrow(dat))
    dat$South_He <- rep(data[i, 5], nrow(dat))
    dat$Avg_He <- rep(data[i, "Avg_He"], nrow(dat))
    dat$Avg_Fst <- rep(data[i, "Avg_Fst"], nrow(dat))
    dat$CAND <- rep(data[i, "CAND"], nrow(dat))
    dat$WINDOW <- rep(data[i, "WINDOW"], nrow(dat))
    
    tmp[[i]] <- dat
    
    rm(Chrom, Spos, Epos, dat)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  return(tmp2)
}

Ts.gene <- GeneAsign(Ts.Fst.He, Ts.gff)

### Assign genes in between windows to window with longer section of gene

GeneAsign2 <- function(data, gene.annotation, assinged.genes){
  # Find genes between windows
  btw <- gene.annotation[!(gene.annotation[,"Gene"] %in% assinged.genes[,"Gene"]),]
  
  # Empty dataset
  tmp <- list()
  
  # For loop assinging candidate info to 'btw'
  for(i in 1:nrow(btw)){
    Chrom <- btw[i, "CHROM"]
    Spos <- btw[i, "Start_Pos"]
    Epos <- btw[i, "End_Pos"]
    
    Pos <- Spos:Epos
    
    # Find windows covering the overlapped gene
    dat <- data[data$CHROM == Chrom & data$END_POS %in% Pos | data$CHROM == Chrom & data$START_POS %in% Pos, ]
    
    # Skip any windows without genes inbetween them
    if(nrow(dat) != 2)next
    
    wind1 <- dat$END_POS[1] - Spos
    wind2 <- Epos - dat$START_POS[2]
    
    if(wind1 > wind2){
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "North_He" = dat[1,4],
                             "South_He" = dat[1,5],
                             "Avg_He" = dat[1, "Avg_He"],
                             "Avg_Fst" = dat[1, "Avg_Fst"],
                             "CAND" = dat[1, "CAND"],
                             "WINDOW" = dat[1,"WINDOW"])
    } else {
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "North_He" = dat[2,4],
                             "South_He" = dat[2,5],
                             "Avg_He" = dat[2, "Avg_He"],
                             "Avg_Fst" = dat[2, "Avg_Fst"],
                             "CAND" = dat[2, "CAND"],
                             "WINDOW" = dat[2,"WINDOW"])
    }
    rm(Chrom, Spos, Epos, Pos, dat, wind1, wind2)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  # Join genes to end of previously assinged genes
  tmp3 <- rbind(assinged.genes, tmp2)
  
  # sort
  tmp4 <- tmp3[with(tmp3, order(CHROM, Start_Pos)), ]
  
  return(tmp4) 
}

Ts.gene2 <- GeneAsign2(Ts.Fst.He, Ts.gff, Ts.gene)

# Asign genes to ninespien windows
tmp <- list()
  
for(i in 1:nrow(Ns.Fst.He)){
  Chrom <- Ns.Fst.He[i, "CHROM"]
  Spos <- Ns.Fst.He[i, "START_POS"]
  Epos <- Ns.Fst.He[i, "END_POS"]
    
  dat <- Ns.gff[Ns.gff$CHROM == Chrom &
                    Ns.gff$Start_Pos >= Spos &
                    Ns.gff$End_Pos <= Epos, ]
  if(nrow(dat) == 0)next
  dat$NsNUn_He <- rep(Ns.Fst.He[i, "NsNUn_He"], nrow(dat))
  dat$NsNUd_He <- rep(Ns.Fst.He[i, "NsNUd_He"], nrow(dat))
  dat$NsABk_He <- rep(Ns.Fst.He[i, "NsABk_He"], nrow(dat))
  dat$NsABm_He <- rep(Ns.Fst.He[i, "NsABm_He"], nrow(dat))
  dat$Avg_He <- rep(Ns.Fst.He[i, "Avg_He"], nrow(dat))
  dat$Avg_Fst_NS1 <- rep(Ns.Fst.He[i, "Avg_Fst_NS1"], nrow(dat))
  dat$Avg_Fst_NS2 <- rep(Ns.Fst.He[i, "Avg_Fst_NS2"], nrow(dat))
  dat$WINDOW <- rep(Ns.Fst.He[i, "WINDOW"], nrow(dat))
    
  tmp[[i]] <- dat
    
  rm(Chrom, Spos, Epos, dat)
}; rm(i)
  
Ns.gene <- do.call(rbind.data.frame, tmp)

rm(tmp)

GeneAsign2.Ns <- function(data, gene.annotation, assinged.genes){
  # Find genes between windows
  btw <- gene.annotation[!(gene.annotation[,4] %in% assinged.genes[,4]),]
  
  # Empty dataset
  tmp <- list()
  
  # For loop assinging candidate info to 'btw'
  for(i in 1:nrow(btw)){
    Chrom <- btw[i, "CHROM"]
    Spos <-btw[i, "Start_Pos"]
    Epos <- btw[i, "End_Pos"]
    
    Pos <- Spos:Epos
    
    # Find windows covering the overlapped gene
    dat <- data[data$CHROM == Chrom & data$END_POS %in% Pos | data$CHROM == Chrom & data$START_POS %in% Pos, ]
    
    # Skip any windows without genes inbetween them
    if(nrow(dat) != 2)next
    
    wind1 <- dat$END_POS[1] - Spos
    wind2 <- Epos - dat$START_POS[2]
    
    if(wind1 > wind2){
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "NsNUn_He" = dat[1,"NsNUn_He"],
                             "NsNUd_He" = dat[1,"NsNUd_He"],
                             "NsABk_He" = dat[1,"NsABk_He"],
                             "NsABm_He" = dat[1,"NsABm_He"],
                             "Avg_He" = dat[1, "Avg_He"],
                             "Avg_Fst_NS1" = dat[1, "Avg_Fst_NS1"],
                             "Avg_Fst_NS2" = dat[1, "Avg_Fst_NS2"],
                             "WINDOW" = dat[1,"WINDOW"])
    } else {
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "NsNUn_He" = dat[2,"NsNUn_He"],
                             "NsNUd_He" = dat[2,"NsNUd_He"],
                             "NsABk_He" = dat[2,"NsABk_He"],
                             "NsABm_He" = dat[2,"NsABm_He"],
                             "Avg_He" = dat[2, "Avg_He"],
                             "Avg_Fst_NS1" = dat[2, "Avg_Fst_NS1"],
                             "Avg_Fst_NS2" = dat[2, "Avg_Fst_NS2"],
                             "WINDOW" = dat[2,"WINDOW"])
    }
    rm(Chrom, Spos, Epos, Pos, dat, wind1, wind2)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  # Join genes to end of previously assinged genes
  tmp3 <- rbind(assinged.genes, tmp2)
  
  # sort
  tmp4 <- tmp3[with(tmp3, order(CHROM, Start_Pos)), ]
  
  return(tmp4) 
}

Ns.gene2 <- GeneAsign2.Ns(Ns.Fst.He, Ns.gff, Ns.gene)

###3. Merge Threespine and Ninespine
tmp1 <- merge(Ts.gene2, Ns.gene2, by = "Gene")
Ts_Ns <- tmp1[, c(2:21,1)]
colnames(Ts_Ns) <- c("CHROM.Ts", "Start_Pos.Ts", "End_Pos.Ts", "TsAK_He", "TsOR_He",
                     "Avg_He.Ts", "Avg_Fst.Ts", "CAND.Ts", "Window.Ts",
                     "CHROM.Ns", "Start_Pos.Ns", "End_Pos.Ns", "NsNUn_He", "NsNUd_He","NsABk_He", "NsABm_He",
                     "Avg_He.Ns", "Avg_Fst_NS1.Ns", "Avg_Fst_NS2.Ns", "Window.Ns", "Gene.Ts")
rm(tmp1)

###4. Find outluiers 
assign.outliers <- function(data, full.genome.data, upper.threshold, lower.threshold){
  
  logi <- array(NA, dim = c(length(data), 3))
  
  logi[,1] <- data > quantile(full.genome.data, upper.threshold)
  logi[,2] <- data < quantile(full.genome.data, lower.threshold)
  
  logi[1:nrow(logi),3] <- "non-outlier" 
  logi[which(logi[,1] == TRUE), 3] <-"upper"
  logi[which(logi[,2] == TRUE), 3] <-"lower"
  
  return(logi[,3])
}

# Threespine
Ts_Ns$Outliers.Ts <- assign.outliers(Ts_Ns$Avg_He.Ts, Ts.He$Avg_He, 0.99, 0.01)

# Ninespine
Ts_Ns$Outliers.Ns <- assign.outliers(Ts_Ns$Avg_He.Ns, Ns.He$Avg_He, 0.99, 0.01)


# both
outlier.overlaps <- function(x, y, label.x, label.y){
  
  logi <- rep("non-outlier", length(x))
  
  logi[which(x != "non-outlier")] <- label.x
  logi[which(y != "non-outlier")] <- label.y
  logi[which(x != "non-outlier" & y != "non-outlier")] <- "Both"
  
  return(logi)
}

Ts_Ns$Outliers.Ts_Ns <- outlier.overlaps(Ts_Ns$Outliers.Ts, Ts_Ns$Outliers.Ns, "Threespine", "Ninespine")


###5. Correlations
cor.test(Ts_Ns$TsAK_He, Ts_Ns$TsOR_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsAK_He, Ts_Ns$NsNUn_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsAK_He, Ts_Ns$NsNUd_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsAK_He, Ts_Ns$NsABk_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsAK_He, Ts_Ns$NsABm_He, alternative = "two.sided", method = "spearman")

cor.test(Ts_Ns$TsOR_He, Ts_Ns$NsNUn_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsOR_He, Ts_Ns$NsNUd_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsOR_He, Ts_Ns$NsABk_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$TsOR_He, Ts_Ns$NsABm_He, alternative = "two.sided", method = "spearman")

cor.test(Ts_Ns$NsNUn_He, Ts_Ns$NsNUd_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$NsNUn_He, Ts_Ns$NsABk_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$NsNUn_He, Ts_Ns$NsABm_He, alternative = "two.sided", method = "spearman")

cor.test(Ts_Ns$NsNUd_He, Ts_Ns$NsABk_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$NsNUd_He, Ts_Ns$NsABm_He, alternative = "two.sided", method = "spearman")

cor.test(Ts_Ns$NsABk_He, Ts_Ns$NsABm_He, alternative = "two.sided", method = "spearman")

cor.test(Ts_Ns$Avg_He.Ts, Ts_Ns$Avg_He.Ns, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$Avg_Fst.Ts, Ts_Ns$NS1_Fst.Ns, alternative = "two.sided", method = "spearman")
cor.test(Ts_Ns$Avg_Fst.Ts, Ts_Ns$NS2_Fst.Ns, alternative = "two.sided", method = "spearman")

###6. Plot Genetic diversity overlaps
ggplot(Ts_Ns)+
  geom_point(aes(x = Avg_He.Ts, y = Avg_He.Ns, colour = Outliers.Ts_Ns))+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.99, na.rm = TRUE), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.01, na.rm = TRUE), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He$Avg_He, 0.99, na.rm = TRUE), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He$Avg_He, 0.01, na.rm = TRUE), lty = 2, colour = "firebrick")+
  labs(x = expression(Threespine~bar(H)[E]), y = expression(Ninespine~bar(H)[E]),
       colour = "Outlier")+
  lims(y = c(0, 0.012), x = c(0, 0.012))+
  scale_colour_manual(values = c("purple", "firebrick", "grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

#####################################################
# NS1
ggplot(Ts_Ns)+
  geom_point(aes(x = Avg_Fst.Ts, y = Avg_Fst_NS1.Ns, colour = CAND.Ts, alpha = CAND.Ts))+
  labs(x = expression(Threespine~F[ST]), y = expression(Ninespine~F[ST]~NS1),
       colour = "Outlier")+
  lims(y = c(0, 1), x = c(0, 0.75))+
  scale_colour_manual(values = c("grey", "navy"))+
  scale_alpha_manual(values = c(0.1,1,1))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

# NS2
ggplot(Ts_Ns)+
  geom_point(aes(x = Avg_Fst.Ts, y = Avg_Fst_NS2.Ns, colour = CAND.Ts, alpha = CAND.Ts))+
  labs(x = expression(Threespine~F[ST]), y = expression(Ninespine~F[ST]~NS2),
       colour = "Outlier")+
  lims(y = c(0, 1), x = c(0, 0.75))+
  scale_colour_manual(values = c("grey", "navy"))+
  scale_alpha_manual(values = c(0.1,1,1))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

###7. Save data
write.table(Ts_Ns, "file:///G:/Ts_Ns_50Kdata.txt")
