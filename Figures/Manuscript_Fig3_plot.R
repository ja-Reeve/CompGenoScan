####### Manuscript figure 3 ################
### This script creates the plot used in figure 3 of my manuscript.
### Fig. 3 show the correlations of genetic diversity and Fst among species
### James Reeve
### G?teborgs Universitet
### 17/02/2020

#### Preperation ####

### clear workingspace
rm(list = ls())
dev.off()
options(stringsAsFactors = FALSE)

### Packages
library("ggplot2")
library("gridExtra")
library("cowplot")

#### Fig. 3a: genetic diversity correlations ####

### Download data
Ts.He <- read.table("file:///F:/nucleotide_diversity/Ts_Avg_Hexp_50Kwind.txt", header = TRUE)
Tu.He <- read.table("file:///F:/nucleotide_diversity/Tu_Avg_Hexp_50Kwind.txt", header = TRUE)
Ns.He <- read.table("file:///F:/nucleotide_diversity/Ns_Avg_Hexp_50Kwind.txt", header = TRUE)

### Remove NAs
Ts.He <- Ts.He[complete.cases(Ts.He),]
Tu.He <- Tu.He[complete.cases(Tu.He),]
Ns.He <- Ns.He[complete.cases(Ns.He),]

### Download gene annotations
Ts.gff <- read.table("file:///F:/GFF/Ts_gene_anno.txt", header = TRUE)
Ts.gff <- Ts.gff[Ts.gff$CHROM != "chrUn", ]
Tu.gff <- read.table("file:///F:/GFF/Tu_gene_anno.txt", header = TRUE)
Ns.gff <- read.table("file:///F:/GFF/Ns_gene_anno.txt", header = TRUE)

# Rename gene header
Ts.gff$Gene <- paste("GAC", sapply(strsplit(Ts.gff$Gene, "[=]"), `[[`, 2), sep = "_")
Tu.gff$Gene <- paste("AFL", sapply(strsplit(Tu.gff$Gene, "[=]"), `[[`, 2), sep = "_")
Ns.gff$Gene <- paste("GAC", sapply(strsplit(Ns.gff$Gene, "[=]"), `[[`, 2), sep = "_")

### Update chromosome names
# threespine
chrom_list_Ts <- read.table("file:///F:/chrom_list/Ts_chr_list.txt")$V1[-c(6:7)]
new_chrom_list_Ts <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8",
                       "chr10", "chr11", "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17",
                       "chr18", "chr20", "chr21")
for(i in 1:length(chrom_list_Ts)){
  Ts.gff[Ts.gff$CHROM == chrom_list_Ts[i], 1] <- new_chrom_list_Ts[i]
}; rm(i, chrom_list_Ts, new_chrom_list_Ts)

# tubesnout
chrom_list_Tu <- read.table("file:///F:/chrom_list/Tu_chr_list.txt")$V1
new_chrom_list_Tu <- paste0("chr", 1:23)
for(i in 1:length(chrom_list_Tu)){
  Tu.gff[Tu.gff$CHROM == chrom_list_Tu[i], 1] <- new_chrom_list_Tu[i]
}; rm(i, chrom_list_Tu, new_chrom_list_Tu)

### Assing window He to genes within windows
### Note: this leads to multicolinearity among genes in same window

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
    dat$North_He <- rep(data[i, 5], nrow(dat))
    dat$South_He <- rep(data[i, 6], nrow(dat))
    dat$Avg_He <- rep(data[i, "Avg_He"], nrow(dat))
    dat$WINDOW <- rep(data[i, "WIND"], nrow(dat))
    
    tmp[[i]] <- dat
    
    rm(Chrom, Spos, Epos, dat)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  return(tmp2)
}

Ts.gene <- GeneAsign(Ts.He, Ts.gff)
Tu.gene <- GeneAsign(Tu.He, Tu.gff)

# seperate function for ninespines
GeneAsign.Ns <- function(data, gene.annotation){
  tmp <- list()
  
  for(i in 1:nrow(data)){
    Chrom <- data[i, "CHROM"]
    Spos <- data[i, "START_POS"]
    Epos <- data[i, "END_POS"]
    
    dat <- gene.annotation[gene.annotation$CHROM == Chrom &
                             gene.annotation$Start_Pos >= Spos &
                             gene.annotation$End_Pos <= Epos, ]
    
    if(nrow(dat) == 0)next
    dat$NsNUn_He <- rep(data[i, 5], nrow(dat))
    dat$NsNUd_He <- rep(data[i, 6], nrow(dat))
    dat$NsABk_He <- rep(data[i, 7], nrow(dat))
    dat$NsABm_He <- rep(data[i, 8], nrow(dat))
    dat$Avg_He <- rep(data[i, "Avg_He"], nrow(dat))
    dat$WINDOW <- rep(data[i, "WIND"], nrow(dat))
    
    tmp[[i]] <- dat
    
    
    rm(Chrom, Spos, Epos, dat)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  return(tmp2)
}

Ns.gene <- GeneAsign.Ns(Ns.He, Ns.gff)

### Assing He for genes that straddel two windows
### The score of the window with the majority of the gene's length is used
GeneAsign2 <- function(data, gene.annotation, assinged.genes){
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
    #if(nrow(dat) != 2)next
    
    wind1 <- dat$END_POS[1] - Spos
    wind2 <- Epos - dat$START_POS[2]
    
    if(wind1 > wind2){
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "North_He" = dat[1,5],
                             "South_He" = dat[1,6],
                             "Avg_He" = dat[1, "Avg_He"],
                             "WINDOW" = dat[1,"WIND"])
    } else {
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "North_He" = dat[2,5],
                             "South_He" = dat[2,6],
                             "Avg_He" = dat[2, "Avg_He"],
                             "WINDOW" = dat[2,"WIND"])
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

Ts.gene2 <- GeneAsign2(Ts.He, Ts.gff, Ts.gene)
Tu.gene2 <- GeneAsign2(Tu.He, Tu.gff, Tu.gene)

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
                             "NsNUn_He" = dat[1,5],
                             "NsNUd_He" = dat[1,6],
                             "NsABk_He" = dat[1,7],
                             "NsABm_He" = dat[1,8],
                             "Avg_He" = dat[1, "Avg_He"],
                             "WINDOW" = dat[1,"WIND"])
    } else {
      tmp[[i]] <- data.frame("CHROM" = Chrom, 
                             "Start_Pos" = Spos, 
                             "End_Pos" = Epos, 
                             "Gene" = btw$Gene[i],
                             "NsNUn_He" = dat[2,5],
                             "NsNUd_He" = dat[2,6],
                             "NsABk_He" = dat[2,7],
                             "NsABm_He" = dat[2,8],
                             "Avg_He" = dat[2, "Avg_He"],
                             "WINDOW" = dat[2,"WIND"])
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

Ns.gene2 <- GeneAsign2.Ns(Ns.He, Ns.gff, Ns.gene)



### Join threespine and ninespine data
Ts_Ns <- merge(Ts.gene2, Ns.gene2, by = "Gene")[c(2:8,1,9:17)]
names(Ts_Ns) <- c("CHROM.Ts", "START_POS.Ts", "END_POS.Ts", "TsAK_He", "TsOR_He", "Avg_He.Ts", "WIND.Ts", "GENE.Ts",
                  "CHROM.Ns", "START_POS.Ns", "END_POS.Ns", "NsNUn_He", "NsNUd_He", "NsABk_He", "NsABm_He", "Avg_He.Ns", "WIND.Ns")

### Join Ts_Ns with Tubesnout data

ortho <- read.table("file:///F:/GFF/GAC-AFL.txt", skip = 4, sep = "\t")[, 3:5]

#filter out one to many orthologs
ortho <- ortho[ortho$V5 == "1:1", 1:2]

tmp <- merge(Ts_Ns, ortho, by.x = "GENE.Ts", by.y = "V3")
Ts_Ns_Tu <- merge(tmp, Tu.gene2, by.x = "V4", by.y = "Gene")[,c(3:9, 2, 10:25, 1)]
names(Ts_Ns_Tu)[18:25] <- c("CHROM.Tu", "START_POS.Tu", "END_POS.Tu", "TuAK_He", "TuBC_He", "Avg_He.Tu", "WIND.Tu", "GENE.Tu")

### Find outliers

# Function
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
Ts_Ns_Tu$Outliers.Ts <- assign.outliers(Ts_Ns_Tu$Avg_He.Ts, Ts.He$Avg_He, 0.99, 0.01)

# Ninespine
Ts_Ns_Tu$Outliers.Ns <- assign.outliers(Ts_Ns_Tu$Avg_He.Ns, Ns.He$Avg_He, 0.99, 0.01)

# Threespine
Ts_Ns_Tu$Outliers.Tu <- assign.outliers(Ts_Ns_Tu$Avg_He.Tu, Tu.He$Avg_He, 0.99, 0.01)

###  Find Outlier overlaps
outlier.overlaps <- function(x, y, label.x, label.y){
  
  logi <- rep("non-outlier", length(x))
  
  logi[which(x != "non-outlier")] <- label.x
  logi[which(y != "non-outlier")] <- label.y
  logi[which(x != "non-outlier" & y != "non-outlier")] <- "Both"
  
  return(logi)
}

# Threespine vs Ninespine
Ts_Ns_Tu$Outliers.Ts_Ns <- outlier.overlaps(Ts_Ns_Tu$Outliers.Ts, Ts_Ns_Tu$Outliers.Ns, "Threespine", "Ninespine")

# Threespine vs Tubesnout
Ts_Ns_Tu$Outliers.Ts_Tu <- outlier.overlaps(Ts_Ns_Tu$Outliers.Ts, Ts_Ns_Tu$Outliers.Tu, "Threespine", "Tubesnout")

# Ninespine vs Tubesnout
Ts_Ns_Tu$Outliers.Ns_Tu <- outlier.overlaps(Ts_Ns_Tu$Outliers.Ns, Ts_Ns_Tu$Outliers.Tu, "Ninespine", "Tubesnout")


### Create plots
txt_size <- 6 # Label text size
pnt_size <- 5 # Axis tick size

Ts_Ns.plot.He <- ggplot(Ts_Ns_Tu)+
  geom_point(aes(x = Avg_He.Ts, y = Avg_He.Ns, colour = Outliers.Ts_Ns), size = 0.5)+
  geom_vline(xintercept = c(quantile(Ts.He$Avg_He, 0.01), quantile(Ts.He$Avg_He, 0.99)), colour = "red", lty = 3)+
  geom_hline(yintercept = c(quantile(Ns.He$Avg_He, 0.01), quantile(Ns.He$Avg_He, 0.99)), colour = "red", lty = 3)+
  labs(x = expression(Threespine~bar(H)[E]), y = expression(Ninespine~bar(H)[E]))+
  scale_colour_manual(values = c("navy", "grey", "navy"))+
  theme_classic()+ theme(axis.title = element_text(size = txt_size), axis.text = element_text(size = pnt_size -1), legend.position = "none")

Ts_Tu.plot.He <- ggplot(Ts_Ns_Tu)+
  geom_point(aes(x = Avg_He.Ts, y = Avg_He.Tu, colour = Outliers.Ts_Tu), size = 0.5)+
  geom_vline(xintercept = c(quantile(Ts.He$Avg_He, 0.01), quantile(Ts.He$Avg_He, 0.99)), colour = "red", lty = 3)+
  geom_hline(yintercept = c(quantile(Tu.He$Avg_He, 0.01), quantile(Tu.He$Avg_He, 0.99)), colour = "red", lty = 3)+
  labs(x = expression(Threespine~bar(H)[E]), y = expression(Tubesnout~bar(H)[E]))+
  scale_colour_manual(values = c("firebrick", "grey", "navy", "navy"))+
  theme_classic()+ theme(axis.title = element_text(size = txt_size), axis.text = element_text(size = pnt_size - 1), legend.position = "none")

Ns_Tu.plot.He <- ggplot(Ts_Ns_Tu)+
  geom_point(aes(x = Avg_He.Ns, y = Avg_He.Tu, colour = Outliers.Ns_Tu), size = 0.5)+
  geom_vline(xintercept = c(quantile(Ns.He$Avg_He, 0.01), quantile(Ns.He$Avg_He, 0.99)), colour = "red", lty = 3)+
  geom_hline(yintercept = c(quantile(Tu.He$Avg_He, 0.01), quantile(Tu.He$Avg_He, 0.99)), colour = "red", lty = 3)+
  labs(x = expression(Ninespine~bar(H)[E]), y = expression(Tubesnout~bar(H)[E]))+
  scale_colour_manual(values = c("navy", "grey", "navy"))+
  theme_classic()+ theme(axis.title = element_text(size = txt_size), axis.text = element_text(size = pnt_size - 1), legend.position = "none")

Fig3a <- grid.arrange(grobs = list(Ts_Ns.plot.He, Ts_Tu.plot.He, Ns_Tu.plot.He), ncol = 1)

#### Figure 3B: Among population correlations heatmap ####

### Determine correlations

# Within threespines
t1 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$TsOR_He, method = "spearman")

# Within ninespine
t2 <- cor.test(x = Ts_Ns_Tu$NsNUn_He, y = Ts_Ns_Tu$NsNUd_He, method = "spearman")
t3 <- cor.test(x = Ts_Ns_Tu$NsNUn_He, y = Ts_Ns_Tu$NsABk_He, method = "spearman")
t4 <- cor.test(x = Ts_Ns_Tu$NsNUn_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

t5 <- cor.test(x = Ts_Ns_Tu$NsNUd_He, y = Ts_Ns_Tu$NsABk_He, method = "spearman")
t6 <- cor.test(x = Ts_Ns_Tu$NsNUd_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

t7 <- cor.test(x = Ts_Ns_Tu$NsABk_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

# Within tubesnout
t8 <- cor.test(x = Ts_Ns_Tu$TuAK_He, y = Ts_Ns_Tu$TuBC_He, method = "spearman")

# Among threespine and tubesnout populations
t9 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$TuAK_He, method = "spearman")
t10 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$TuBC_He, method = "spearman")

t11 <- cor.test(x = Ts_Ns_Tu$TsOR_He, y = Ts_Ns_Tu$TuAK_He, method = "spearman")
t12 <- cor.test(x = Ts_Ns_Tu$TsOR_He, y = Ts_Ns_Tu$TuBC_He, method = "spearman")

# Among threespine and ninespine populations
t13 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$NsNUn_He, method = "spearman")
t14 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$NsNUd_He, method = "spearman")
t15 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$NsABk_He, method = "spearman")
t16 <- cor.test(x = Ts_Ns_Tu$TsAK_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

t17 <- cor.test(x = Ts_Ns_Tu$TsOR_He, y = Ts_Ns_Tu$NsNUn_He, method = "spearman")
t18 <- cor.test(x = Ts_Ns_Tu$TsOR_He, y = Ts_Ns_Tu$NsNUd_He, method = "spearman")
t19 <- cor.test(x = Ts_Ns_Tu$TsOR_He, y = Ts_Ns_Tu$NsABk_He, method = "spearman")
t20 <- cor.test(x = Ts_Ns_Tu$TsOR_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

# Among threespine and ninespine populations
t21 <- cor.test(x = Ts_Ns_Tu$TuAK_He, y = Ts_Ns_Tu$NsNUn_He, method = "spearman")
t22 <- cor.test(x = Ts_Ns_Tu$TuAK_He, y = Ts_Ns_Tu$NsNUd_He, method = "spearman")
t23 <- cor.test(x = Ts_Ns_Tu$TuAK_He, y = Ts_Ns_Tu$NsABk_He, method = "spearman")
t24 <- cor.test(x = Ts_Ns_Tu$TuAK_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

t25 <- cor.test(x = Ts_Ns_Tu$TuBC_He, y = Ts_Ns_Tu$NsNUn_He, method = "spearman")
t26 <- cor.test(x = Ts_Ns_Tu$TuBC_He, y = Ts_Ns_Tu$NsNUd_He, method = "spearman")
t27 <- cor.test(x = Ts_Ns_Tu$TuBC_He, y = Ts_Ns_Tu$NsABk_He, method = "spearman")
t28 <- cor.test(x = Ts_Ns_Tu$TuBC_He, y = Ts_Ns_Tu$NsABm_He, method = "spearman")

### Join correlations into a single dataset
spearman_correlations <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,
                              t21,t22,t23,t24,t25,t26,t27,t28)

Res <- list() 
for(i in 1:length(spearman_correlations)){
  dat <- spearman_correlations[[i]]
  
  tmp.name <- unlist(strsplit(dat$data.name, " and "))
  pop1 <- sapply(strsplit(tmp.name[1], "[$]"), `[[`, 2)
  pop2 <- sapply(strsplit(tmp.name[2], "[$]"), `[[`, 2)
  
  rho <- dat$estimate
  p.value <- dat$p.value
  
  Res[[i]] <- data.frame(pop1, pop2, rho, p.value)
  
  rm(dat, tmp.name, pop1, pop2, rho, p.value)
};rm(i)

correlations <- do.call(rbind.data.frame, Res)

correlations2 <- correlations[,c(2,1,3:4)] # reverse order of populations
names(correlations2) <- c("pop1", "pop2", "rho", "p.value")
correlations2 <- rbind.data.frame(correlations, correlations2)
correlations2$pop1 <- unlist(strsplit(correlations2$pop1, "_He"))
correlations2$pop2 <- unlist(strsplit(correlations2$pop2, "_He"))

correlations2$sig <- cut(correlations2$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

rm(Res, t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,
   t21,t22,t23,t24,t25,t26,t27,t28, spearman_correlations)

### Create heatmap plot
fig3b <- ggplot(correlations2, aes(x = pop1, y = pop2))+
  geom_tile(aes(fill = abs(rho)))+
  geom_text(aes(label = sig), size = pnt_size)+
  geom_hline(yintercept = c(2.5, 6.5), lty = 2, colour = "grey")+
  geom_vline(xintercept = c(2.5, 6.5), lty = 2, colour = "grey")+
  labs(fill = expression(Spearman~rho))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "orange", midpoint = 0.2,
                       breaks = seq(-0.2, 1.0, 0.1), labels = seq(-0.2, 1.0, 0.1))+
  scale_x_discrete(limits = c("TsAK", "TsOR", "NsNUn", "NsNUd", "NsABk", "NsABm", "TuAK", "TuBC"), position = "top")+
  scale_y_discrete(limits = rev(c("TsAK", "TsOR", "NsNUn", "NsNUd", "NsABk", "NsABm", "TuAK", "TuBC")))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(colour = "grey", fill=NA),
        legend.position = "bottom",
        legend.title = element_text(size = txt_size, hjust = 1, vjust = 0.7),
        legend.text = element_text(size = txt_size-2)#,
        )#legend.key.width = unit(3, "cm"))



#### Figure 3c: correlation of gene averaged Fst between species ####

### Dowload gene-averaged Fst
Ts_Ns.Fst <- read.table("file:///F:/Ts_Ns_50Kdata.txt", header = TRUE)
Ts_Tu.Fst <- read.table("file:///F:/Ts_Tu_50Kdata.txt", header = TRUE)

### Merge data
Ts_Ns_Tu.Fst <- merge(Ts_Ns.Fst, Ts_Tu.Fst, by = "Gene.Ts")
Ts_Ns_Tu.Fst$Avg_Fst.Ns <- rowMeans(Ts_Ns_Tu.Fst[,19:20])
### Plots
# Ts vs Ns
Ts_Ns.Fstplot <- ggplot(Ts_Ns_Tu.Fst)+
  geom_point(aes(x = Avg_Fst.Ts.x, y = Avg_Fst.Ns, colour = CAND.both, alpha = CAND.both), size = 0.5)+
  labs(x = expression(Threespine~F[ST]), y = expression(Ninespine~F[ST]),
       colour = "Outlier")+
  lims(y = c(0, 1), x = c(0, 1))+
  coord_fixed()+
  scale_colour_manual(values = c("grey", "firebrick", "navy"))+
  scale_alpha_manual(values = c(0.1,1,1))+
  theme(axis.title = element_text(size = txt_size), 
        axis.text = element_text(size = pnt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),legend.position = "NULL")

# Ns vs Tu
Ns_Tu.Fstplot <- ggplot(Ts_Ns_Tu.Fst)+
  geom_point(aes(x = Avg_Fst.Ns, y = Avg_Fst.Tu, colour = CAND.both, alpha = CAND.both), size = 0.5)+
  labs(x = expression(Ninespine~F[ST]), y = expression(Tubesnout~F[ST]),
       colour = "Outlier")+
  lims(y = c(0, 1), x = c(0, 1))+
  coord_fixed()+
  scale_colour_manual(values = c("grey", "firebrick", "navy"))+
  scale_alpha_manual(values = c(0.1,1,1))+
  theme(axis.title = element_text(size = txt_size), 
        axis.text = element_text(size = pnt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),legend.position = "NULL")

# Ts vs Tu
#tmp.5 <- cor.test(Ts_Ns_Tu.Fst$Avg_Fst.Ts.x, Ts_Ns_Tu.Fst$Avg_Fst.Tu, method = "spearman")

Ts_Tu.Fstplot <- ggplot(Ts_Ns_Tu.Fst)+
  geom_point(aes(x = Avg_Fst.Ts.x, y = Avg_Fst.Tu, colour = CAND.both, alpha = CAND.both), size = 0.5)+
  #geom_text(aes(label = paste('rho', "==", round(tmp.5$estimate, 4)), x = 0.85, y = 0.85), parse = TRUE, size = 8)+
  labs(x = expression(Threespine~F[ST]), y = expression(Tubesnout~F[ST]),
       colour = "Outlier")+
  lims(y = c(0, 1), x = c(0, 1))+
  coord_fixed()+
  scale_colour_manual(values = c("grey", "firebrick", "navy"))+
  scale_alpha_manual(values = c(0.1,1,1))+
  theme(axis.title = element_text(size = txt_size), 
        axis.text = element_text(size = pnt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),legend.position = "NULL")


# Multi-panel plot
Fig3c <- grid.arrange(grobs = list(Ts_Ns.Fstplot, Ts_Tu.Fstplot, Ns_Tu.Fstplot), ncol = 1)
#rm(tmp.1, tmp.2, tmp.3, tmp.4, tmp.5)

# Save plot
tiff("/Users/james/Dropbox (Personal)/Comp_geno_2020/Figure_3.tiff", width = 18, height = 10,
     units = "cm", res = 300)
plot_grid(Fig3a, fig3b, Fig3c, nrow = 1, rel_widths = c(1,2,1), labels = "AUTO", label_size = txt_size * 2, scale = 0.98)
dev.off()
