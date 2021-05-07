#### Testing comparisons among Ts.He and Tu.He using 50Kb windows ####
### James Reeve
### 19/09/2019

rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)

### Download He data
Ts.He <- read.table("file:///G:/nucleotide_diversity/Ts_Avg_Hexp_50K.txt", header = TRUE)
Tu.He <- read.table("file:///G:/nucleotide_diversity/Tu_Avg_Hexp_50K.txt", header = TRUE)

Ts.He <- Ts.He[complete.cases(Ts.He),]
Tu.He <- Tu.He[complete.cases(Tu.He),]

### Download Fst data
Ts.Fst <- read.table("file:///G:/Fst-table/Ts.Fst.50K.txt", header = TRUE)
Tu.Fst <- read.table("file:///G:/Fst-table/Tu.Fst.50K.txt", header = TRUE)
# assign top candidate flag
Ts.TC <- read.table("file:///G:/top-candidate/Ts.TopCandidates.50Kwind.txt", header = TRUE)
Ts.TC$CAND <- Ts.TC$Observed_Outliers > Ts.TC$Expected_Outliers
Tu.TC <- read.table("file:///G:/top-candidate/Tu.TopCandidates.50Kwind.txt", header = TRUE)
Tu.TC$CAND <- Tu.TC$Observed_Outliers > Tu.TC$Expected_Outliers

Ts.Fst <- merge(Ts.Fst, Ts.TC, by.x = "WIND", by.y = "Window_Name")[,c(2:5,12,1)] 
Tu.Fst <- merge(Tu.Fst, Tu.TC, by.x = "WIND", by.y = "Window_Name")[,c(2:5,12,1)] 

rm(Ts.TC, Tu.TC)

### Download gene annotations
Ts.gff <- read.table("file:///G:/GFF/Ts_gene_anno.txt", header = TRUE)
Ts.gff <- Ts.gff[Ts.gff$CHROM != "chrUn", ]
Tu.gff <- read.table("file:///G:/GFF/Tu_gene_anno.txt", header = TRUE)
# Rename gene header
Ts.gff$Gene <- paste("GAC", sapply(strsplit(Ts.gff$Gene, "[=]"), `[[`, 2), sep = "_")
Tu.gff$Gene <- paste("AFL", sapply(strsplit(Tu.gff$Gene, "[=]"), `[[`, 2), sep = "_")

### Download ortholog table
ortho <- read.table("file:///G:/GFF/GAC-AFL.txt", sep = "\t")[,-c(1:2,6)]
colnames(ortho) <- c("threespine_ID", "tubesnout_ID", "type")
ortho <- ortho[ortho$type == "1:1",]

### Update chromosome names
chrom_list_Ts <- read.table("file:///G:/chrom_list/Ts_chr_list.txt")$V1[-c(6:7)]
new_chrom_list_Ts <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8",
                    "chr10", "chr11", "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17",
                    "chr18", "chr20", "chr21")
for(i in 1:length(chrom_list_Ts)){
  Ts.gff[Ts.gff$CHROM == chrom_list_Ts[i], 1] <- new_chrom_list_Ts[i]
}; rm(i, chrom_list_Ts, new_chrom_list_Ts)

chrom_list_Tu <- read.table("file:///G:/chrom_list/Tu_chr_list.txt")$V1
new_chrom_list_Tu <- paste0("chr",1:23)
for(i in 1:length(chrom_list_Tu)){
  Tu.gff[Tu.gff$CHROM == chrom_list_Tu[i], 1] <- new_chrom_list_Tu[i]
}; rm(i, chrom_list_Tu, new_chrom_list_Tu)


### Merge He and Fst data
Ts.Fst.He <- merge(Ts.He, Ts.Fst, by = "WIND")[,c(2:7,11:12,1)]
colnames(Ts.Fst.He) <- c("CHROM", "START_POS", "END_POS", "TsAK_He", "TsOR_He",
                         "Avg_He", "Avg_Fst", "CAND", "WINDOW")
Tu.Fst.He <- merge(Tu.He, Tu.Fst, by = "WIND")[,c(2:7,11:12,1)]
colnames(Tu.Fst.He) <- c("CHROM", "START_POS", "END_POS", "TuAK_He", "TuBC_He",
                         "Avg_He", "Avg_Fst", "CAND", "WINDOW")

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
Tu.gene <- GeneAsign(Tu.Fst.He, Tu.gff)

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
Tu.gene2 <- GeneAsign2(Tu.Fst.He, Tu.gff, Tu.gene)

###3. Merge Threespine and tubesnouts
tmp1 <- merge(Ts.gene2, ortho, by.x = "Gene", by.y = "threespine_ID")
tmp2 <- merge(tmp1, Tu.gene2, by.x = "tubesnout_ID", by.y = "Gene")
Ts_Tu <- tmp2[, c(3:11,2,13:21,1)]
colnames(Ts_Tu) <- c("CHROM.Ts", "Start_Pos.Ts", "End_Pos.Ts", "TsAK_He", "TsOR_He",
                    "Avg_He.Ts", "Avg_Fst.Ts", "CAND.Ts", "Window.Ts", "Gene.Ts",
                     "CHROM.Tu", "Start_Pos.Tu", "End_Pos.Tu", "TuAK_He", "TuBC_He",
                     "Avg_He.Tu", "Avg_Fst.Tu", "CAND.Tu", "Window.Tu", "Gene.Tu")
rm(tmp1, tmp2)

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
Ts_Tu$Outliers.Ts <- assign.outliers(Ts_Tu$Avg_He.Ts, Ts.He$Avg_He, 0.99, 0.01)

# Ninespine
Ts_Tu$Outliers.Tu <- assign.outliers(Ts_Tu$Avg_He.Tu, Tu.He$Avg_He, 0.99, 0.01)


# both
outlier.overlaps <- function(x, y, label.x, label.y){
  
  logi <- rep("non-outlier", length(x))
  
  logi[which(x != "non-outlier")] <- label.x
  logi[which(y != "non-outlier")] <- label.y
  logi[which(x != "non-outlier" & y != "non-outlier")] <- "Both"
  
  return(logi)
}

Ts_Tu$Outliers.Ts_Tu <- outlier.overlaps(Ts_Tu$Outliers.Ts, Ts_Tu$Outliers.Tu, "Threespine", "Tubesnout")

###5. Correlations
cor.test(Ts_Tu$TsAK_He, Ts_Tu$TsOR_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Tu$TsAK_He, Ts_Tu$TuAK_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Tu$TsAK_He, Ts_Tu$TuBC_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Tu$TsOR_He, Ts_Tu$TuAK_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Tu$TsOR_He, Ts_Tu$TuBC_He, alternative = "two.sided", method = "spearman")
cor.test(Ts_Tu$TuAK_He, Ts_Tu$TuBC_He, alternative = "two.sided", method = "spearman")

cor.test(Ts_Tu$Avg_He.Ts, Ts_Tu$Avg_He.Tu, alternative = "two.sided", method = "spearman")
cor.test(Ts_Tu$Avg_Fst.Ts, Ts_Tu$Avg_Fst.Tu, alternative = "two.sided", method = "spearman")

###6. Plot Genetic diversity overlaps
ggplot(Ts_Tu)+
  geom_point(aes(x = Avg_He.Ts, y = Avg_He.Tu, colour = Outliers.Ts_Tu))+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.01), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Tu.He$Avg_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Tu.He$Avg_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(Threespine~bar(H)[E]), y = expression(Tubesnout~bar(H)[E]),
       colour = "Outlier")+
  lims(y = c(0, 0.45), x = c(0, 0.45))+
  scale_colour_manual(values = c("purple", "grey", "firebrick", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),legend.position = "NULL")
##############################################
logi <- array(data = FALSE, dim = c(nrow(Ts_Tu), 3))
logi[,1] <- Ts_Tu$CAND.Ts != "FALSE"
logi[,2] <- Ts_Tu$CAND.Tu != "FALSE"

logi[which(logi[,1] == TRUE), 3] <-"Threespine"
logi[which(logi[,2] == TRUE), 3] <-"Tubesnout"
logi[which(logi[,1] == TRUE & logi[,2] == TRUE), 3] <-"Both"
logi[which(logi[,3] == FALSE), 3] <- "Non-cand"

Ts_Tu$CAND.both <-  logi[,3]
rm(logi)
#####################################################
ggplot(Ts_Tu)+
  geom_point(aes(x = Avg_Fst.Ts, y = Avg_Fst.Tu, colour = CAND.both, alpha = CAND.both))+
  labs(x = expression(Threespine~F[ST]), y = expression(Tubesnout~F[ST]),
       colour = "Outlier")+
  lims(y = c(0, 1), x = c(0, 1))+
  scale_colour_manual(values = c("purple", "grey", "firebrick", "navy"))+
  scale_alpha_manual(values = c(1,0.1,1,1))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),legend.position = "NULL")



###7. Save data
write.table(Ts_Tu, "file:///G:/Ts_Tu_50Kdata.txt")

