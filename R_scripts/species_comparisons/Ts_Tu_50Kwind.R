#### Testing comparisons among Ts.He and Tu.He using 50Kb windows ####
### James Reeve
### 19/09/2019
### Edit: replace candiates detected in 50Kb window with candidates per gene (26/08/2021)

rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)

### Download He data
Ts.He <- read.table("file:///F:/nucleotide_diversity/Ts_Avg_Hexp_50Kwind.txt", header = TRUE)
Tu.He <- read.table("file:///F:/nucleotide_diversity/Tu_Avg_Hexp_50Kwind.txt", header = TRUE)

Ts.He <- Ts.He[complete.cases(Ts.He),]
Tu.He <- Tu.He[complete.cases(Tu.He),]

### Download Fst data
Ts.Fst <- read.table("file:///F:/Fst-table/Ts.Fst.50K.txt", header = TRUE)
Tu.Fst <- read.table("file:///F:/Fst-table/Tu.Fst.50K.txt", header = TRUE)

### Download top canadidate data
Ts.TC <- read.table("file:///F:/top-candidate/Ts.VarScan.TopCandidates.txt", header = TRUE)
Tu.TC <- read.table("file:///F:/top-candidate/Tu.VarScan.TopCandidates.txt", header = TRUE)
# Remove inter-geneic regions
Ts.TC <- Ts.TC[Ts.TC$Gene_Name != 'itgr',]
Tu.TC <- Tu.TC[Tu.TC$Gene_Name != 'itgr',]
# Create flag for top candidates
Ts.TC$CAND <- Ts.TC$Observed_Outliers > Ts.TC$Expected_Outliers
Tu.TC$CAND <- Tu.TC$Observed_Outliers > Tu.TC$Expected_Outliers
# Rename gene prefix
Ts.TC$Gene_Name <- gsub("ID=", "GAC_", Ts.TC$Gene_Name)
Tu.TC$Gene_Name <- gsub("ID=", "AFL_", Tu.TC$Gene_Name)

### Download ortholog table
ortho <- read.table("file:///F:/GFF/GAC-AFL.txt", sep = "\t")[,-c(1:2,6)]
colnames(ortho) <- c("threespine_ID", "tubesnout_ID", "type")
ortho <- ortho[ortho$type == "1:1",]

### Merge He and Fst data
Ts.Fst.He <- merge(Ts.He, Ts.Fst, by = "WIND")[,c(2:7,11,1)]
colnames(Ts.Fst.He) <- c("CHROM", "START_POS", "END_POS", "TsAK_He", "TsOR_He",
                         "Avg_He", "Avg_Fst", "WINDOW")
Tu.Fst.He <- merge(Tu.He, Tu.Fst, by = "WIND")[,c(2:7,11,1)]
colnames(Tu.Fst.He) <- c("CHROM", "START_POS", "END_POS", "TuAK_He", "TuBC_He",
                         "Avg_He", "Avg_Fst", "WINDOW")

### Assing window He to genes within windows ####

GeneAsign <- function(data, gene.annotation){
  tmp <- list()
  
  for(i in 1:nrow(data)){
    Chrom <- data[i, "CHROM"]
    Spos <- data[i, "START_POS"]
    Epos <- data[i, "END_POS"]
    
    dat <- gene.annotation[gene.annotation$Chromosome == Chrom &
                             gene.annotation$Gene_Start >= Spos &
                             gene.annotation$Gene_End <= Epos, ]
    if(nrow(dat) == 0)next

    tmp[[i]] <- data.frame("CHROM" = dat$Chromosome, 
                           "Start_Pos" = dat$Gene_Start, 
                           "End_Pos" = dat$Gene_End, 
                           "Gene" = dat$Gene_Name,
                           "North_He" = data[i, 4],
                           "South_He" = data[i, 5],
                           "Avg_He" = data[i, "Avg_He"],
                           "Avg_Fst" = data[i, "Avg_Fst"],
                           "CAND" = dat$CAND,
                           "WINDOW" = data[i,"WINDOW"])
    
    rm(Chrom, Spos, Epos, dat)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  return(tmp2)
}

Ts.gene <- GeneAsign(Ts.Fst.He, Ts.TC)
Tu.gene <- GeneAsign(Tu.Fst.He, Tu.TC)

### Assign genes in between windows to window with longer section of gene

GeneAsign2 <- function(data, gene.annotation, assigned.genes){
  # Find genes between windows
  btw <- gene.annotation[!(gene.annotation[,"Gene_Name"] %in% assigned.genes[,"Gene"]),]
  
  # Empty dataset
  tmp <- list()
  
  # For loop assinging candidate info to 'btw'
  for(i in 1:nrow(btw)){
    Chrom <- btw[i, "Chromosome"]
    Spos <- btw[i, "Gene_Start"]
    Epos <- btw[i, "Gene_End"]
    
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
                             "Gene" = btw$Gene_Name[i],
                             "North_He" = dat[1,4],
                             "South_He" = dat[1,5],
                             "Avg_He" = dat[1, "Avg_He"],
                             "Avg_Fst" = dat[1, "Avg_Fst"],
                             "CAND" = btw$CAND[i],
                             "WINDOW" = dat[1,"WINDOW"])
    } else {
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene_Name[i],
                             "North_He" = dat[2,4],
                             "South_He" = dat[2,5],
                             "Avg_He" = dat[2, "Avg_He"],
                             "Avg_Fst" = dat[2, "Avg_Fst"],
                             "CAND" = btw$CAND[i],
                             "WINDOW" = dat[2,"WINDOW"])
    }
    rm(Chrom, Spos, Epos, Pos, dat, wind1, wind2)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
 
  # Join genes to end of previously assinged genes
  tmp3 <- rbind(assigned.genes, tmp2)
  
  # sort
  tmp4 <- tmp3[with(tmp3, order(CHROM, Start_Pos)), ]
  
  return(tmp4) 
}

Ts.gene2 <- GeneAsign2(Ts.Fst.He, Ts.TC, Ts.gene)
Tu.gene2 <- GeneAsign2(Tu.Fst.He, Tu.TC, Tu.gene)

###3. Merge Threespine and tubesnouts
tmp1 <- merge(Ts.gene2, ortho, by.x = "Gene", by.y = "threespine_ID")
tmp2 <- merge(tmp1, Tu.gene2, by.x = "tubesnout_ID", by.y = "Gene")
Ts_Tu <- tmp2[, c(3:9, 11, 2, 10, 13:19, 21, 1, 20)]
colnames(Ts_Tu) <- c("CHROM.Ts", "Start_Pos.Ts", "End_Pos.Ts",
                     "TsAK_He", "TsOR_He","Avg_He.Ts", "Avg_Fst.Ts",  
                     "Window.Ts", "Gene.Ts", "CAND.Ts",
                     "CHROM.Tu", "Start_Pos.Tu", "End_Pos.Tu", 
                     "TuAK_He", "TuBC_He", "Avg_He.Tu", "Avg_Fst.Tu", 
                     "Window.Tu", "Gene.Tu", "CAND.Tu")
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
write.table(Ts_Tu, "file:///F:/Ts_Tu_50Kdata.txt")

