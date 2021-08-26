### Script to create Manhattan plots of Fst and He for each species ###
### Fst scores are downsamples to every 100th SNP
### Tubesnout and Ninespine SNPs are plotted on the Threespine genome
### James Reeve - University of Gothenburg
### 05/04/2021

### Prepreation
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()
setwd("C:/Users/James/Desktop")

### Packages
library("tidyverse")
library("cowplot")
library("gridExtra")

### A: Access Fst data ####

### Fst per SNP data
# threespine
Ts.Fst <- read.table("file:///F:/Fst-table/Ts.VarScan.Fst.txt", header = T)
# tubesnout
Tu.Fst <- read.table("file:///F:/Fst-table/Tu.VarScan.Fst.txt", header = T)
Tu.Fst$CHROM <- gsub("chr", "HiC_chromosome_", Tu.Fst$CHROM)
# ninespine
Ns.Fst <- read.table("file:///F:/Fst-table/Ns.VarScan.Fst.txt", header = T)


### B: Access He data ####

### Genetic diversity (He) data
# threespine
Ts.He <- read.table("file:///F:/nucleotide_diversity/Ts_Avg_Hexp_genes.txt", header = TRUE)
Ts.He <- Ts.He[Ts.He$GENE != "itgr", ] # Remove He for windows in intergeneic windows
Ts.He$MID_POS <- round(apply(Ts.He[,2:3], 1, median)) # Add median position of each gene - helps ploting
# tubesnout
Tu.He <- read.table("file:///F:/nucleotide_diversity/Tu_Avg_Hexp_genes.txt", header = TRUE)
Tu.He <- Tu.He[Tu.He$GENE != "itgr", ]
Tu.He$MID_POS <- round(apply(Tu.He[,2:3], 1, median))
# ninespine
Ns.He <- read.table("file:///F:/nucleotide_diversity/Ns_Avg_Hexp_genes.txt", header = TRUE)
Ns.He <- Ns.He[Ns.He$GENE != "itgr", ]
Ns.He$MID_POS <- round(apply(Ns.He[,2:3], 1, median))

### C: Map all scores onto threespine genome ####

### Load gene annotations
Ts.gff <- read.table("file:///F:/GFF/Ts_gene_anno.txt", header = TRUE)
Ns.gff <- read.table("file:///F:/GFF/Ns_gene_anno.txt", header = TRUE)
Tu.gff <- read.table("file:///F:/GFF/Tu_gene_anno.txt", header = TRUE)

### Identify gene matches between threespine and tubesnout
# Ortholog table
ortho <- read.table("file:///F:/GFF/GAC-AFL.txt", sep = "\t")[,-c(1:2,6)]
colnames(ortho) <- c("threespine_ID", "tubesnout_ID", "type")
ortho <- ortho[ortho$type == "1:1",]
# Replace gene column in Tu.GFF with threespine gene
Tu.gff.2 <- Tu.gff
Tu.gff.2$Gene <- gsub("ID=", "AFL_", Tu.gff.2$Gene)
Tu.gff.2 <- merge(Tu.gff.2, ortho, by.x = "Gene", by.y = "tubesnout_ID")[,2:5]
colnames(Tu.gff.2) <- c(colnames(Tu.gff.2)[1:3], "Gene")
Tu.gff.2$Gene <- gsub("GAC_", "ID=", Tu.gff.2$Gene)

### Function to map species positions onto threepine genome
# I've set this up so you have the option of running either Fst or He scores with the same function
Map_to_threespine <- function(GFF_x, GFF_y, genetic.score, data){
  if(genetic.score != "Fst" & genetic.score != "He")stop('Error: genetic.score must be either "Fst" or "He".')
  ### Merge GFFs
  tmp <- merge(GFF_x, GFF_y, by = "Gene")
  tmp2 <- tmp[,c(2:7,1)]
  colnames(tmp2) <- c("chrom_x", "start_x", "end_x", "chrom_y", "start_y", "end_y", "gene_ID")
  ### Fst function
  if(genetic.score == "Fst"){
    tmp3 <- apply(tmp2, 1, function(gene){
      # Key vectors
      gene_chrom <- gene[4]
      gene_start <- as.numeric(gene[5])
      gene_end <- as.numeric(gene[6])
      
      # Get Fst scores for this gene
      dat <- data[data$CHROM == gene_chrom & data$POS >= gene_start & data$POS <= gene_end,]
      
      # Skip genes without SNPs
      if(nrow(dat) > 0){
        
        # Convert species 1's positions to species 2
        t1 <- (dat$POS - gene_start) / (gene_end - gene_start)      # quantile of gene position
        t2 <- as.numeric(gene[3]) - as.numeric(gene[2])             # length gene in species 2
        t3 <- round(t1*t2)                                          # find position (bp) along genome
        t4 <- as.numeric(gene[2]) + t3                              # Add start postion to t3
        
        # Write output
        Res <- as.data.frame(cbind(dat, gene[1], t4))
        
        return(Res)}
    })
    
    tmp4 <- do.call(rbind.data.frame, tmp3)
    colnames(tmp4) <- c(names(data), "CHROM_Ts", "POS_Ts")
    return(tmp4)
  }
  ### He function
  if(genetic.score == "He"){
    tmp3 <- merge(tmp2, data, by.x = "gene_ID", by.y = "GENE")
    tmp4 <- tmp3[,c(8:(ncol(data)+6), 1:4)]
    colnames(tmp4) <- c(colnames(data)[c(1:3,5:ncol(data),4)], "CHROM_Ts", "START_Ts", "END_Ts")
    tmp4$MID_Ts <- round(apply(tmp4[,c("START_Ts", "END_Ts")], 1, median))
    return(tmp4)
  }
}

### Map scores to threespine genome
# Ninespine Fst
Ns.Fst.2 <- Map_to_threespine(Ts.gff, Ns.gff, genetic.score = "Fst", Ns.Fst)
Ns.Fst.2 <- Ns.Fst.2[Ns.Fst.2$CHROM_Ts != "chrUn",]
# Ninespine He
Ns.He.2 <- Map_to_threespine(Ts.gff, Ns.gff, genetic.score = "He", Ns.He)
Ns.He.2 <- Ns.He.2[Ns.He.2$CHROM_Ts != "chrUn",]
# Tubesnout Fst
# Not used in the final figure, as few SNPs map to threespine chr6
#Tu.Fst.2 <- Map_to_threespine(Ts.gff, Tu.gff.2, genetic.score = "Fst", Tu.Fst)
# Tubesnout He
#Tu.He.2 <- Map_to_threespine(Ts.gff, Tu.gff.2, genetic.score = "He", Tu.He)

### Convert chromosome names
# Threespine chr to arabic numerals
chrom_list <- read.table("file:///F:/chrom_list/Ts_chr_list.txt")$V1[-c(6:7)]

new_chrom_list <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8", "chr10",
                       "chr11", "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17", "chr18",
                       "chr20", "chr21")
for (i in 1:length(new_chrom_list)) {
  Ts.He[Ts.He$CHROM == chrom_list[i],1] <- new_chrom_list[i]
  #Tu.He.2[Tu.He.2$CHROM_Ts == chrom_list[i], "CHROM_Ts"] <- new_chrom_list[i]
  Ns.He.2[Ns.He.2$CHROM_Ts == chrom_list[i], "CHROM_Ts"] <- new_chrom_list[i]
  
  #Tu.Fst.2[Tu.Fst.2$CHROM_Ts == chrom_list[i], "CHROM_Ts"] <- new_chrom_list[i]
  Ns.Fst.2[Ns.Fst.2$CHROM_Ts == chrom_list[i], "CHROM_Ts"] <- new_chrom_list[i]
}; rm(i, chrom_list, new_chrom_list)

# Chage "HiC_chromosome_" to "chr" for Tu
Tu.He$CHROM <- gsub("HiC_chromosome_", "chr", Tu.He$CHROM)
Tu.Fst$CHROM <- gsub("HiC_chromosome_", "chr", Tu.Fst$CHROM)

### D: Downsample Fst ####
### Downsample Fst scores - this saves on processing power while still showing genetic patterns
Ts.Fst.2 <- Ts.Fst[seq(1, nrow(Ts.Fst), 100),]
Tu.Fst.2 <- Tu.Fst[seq(1, nrow(Tu.Fst), 100),]
# Not nessecary for Ninespine. There are too few SNPs after matching species

### D: Get cumulative postion on each chromosome ####

# Threespine Fst
Ts.Fst.2$CHROM_test <- factor(Ts.Fst.2$CHROM, levels = paste0("chr",1:21))

Ts.Fst.3 <- Ts.Fst.2 %>% group_by(CHROM_test) %>% summarise(chrom_len=max(POS)) %>%
  mutate(tot=cumsum(chrom_len)-chrom_len) %>% select(-chrom_len) %>%
  left_join(Ts.Fst.2, ., by = c("CHROM_test"="CHROM_test")) %>%
  arrange(CHROM_test, POS) %>% mutate(POScumm=POS+tot)

# Tubesnout Fst
Tu.Fst.2$CHROM_test <- factor(Tu.Fst.2$CHROM, levels = paste0("chr",1:23))

Tu.Fst.3 <- Tu.Fst.2 %>% group_by(CHROM_test) %>% summarise(chrom_len=max(POS)) %>%
  mutate(tot=cumsum(chrom_len)-chrom_len) %>% select(-chrom_len) %>%
  left_join(Tu.Fst.2, ., by = c("CHROM_test"="CHROM_test")) %>%
  arrange(CHROM_test, POS) %>% mutate(POScumm=POS+tot)

# Ninespine Fst
Ns.Fst.2$CHROM_test <- factor(Ns.Fst.2$CHROM_Ts, levels = paste0("chr",1:21))

Ns.Fst.3 <- Ns.Fst.2 %>% group_by(CHROM_test) %>% summarise(chrom_len=max(POS_Ts)) %>%
  mutate(tot=cumsum(chrom_len)-chrom_len) %>% select(-chrom_len) %>%
  left_join(Ns.Fst.2, ., by = c("CHROM_test"="CHROM_test")) %>%
  arrange(CHROM_test, POS_Ts) %>% mutate(POScumm=POS_Ts+tot)


# Threespine He
Ts.He$CHROM_test <- factor(Ts.He$CHROM, levels = paste0("chr",1:21))

Ts.He.2 <- Ts.He %>% group_by(CHROM_test) %>% summarise(chrom_len=max(MID_POS)) %>%
  mutate(tot=cumsum(chrom_len)-chrom_len) %>% select(-chrom_len) %>%
  left_join(Ts.He, ., by = c("CHROM_test"="CHROM_test")) %>%
  arrange(CHROM_test, MID_POS) %>% mutate(POScumm=MID_POS+tot)

# Tubesnout He
Tu.He$CHROM_test <- factor(Tu.He$CHROM, levels = paste0("chr",1:23))

Tu.He.2 <- Tu.He %>% group_by(CHROM_test) %>% summarise(chrom_len=max(MID_POS)) %>%
  mutate(tot=cumsum(chrom_len)-chrom_len) %>% select(-chrom_len) %>%
  left_join(Tu.He, ., by = c("CHROM_test"="CHROM_test")) %>%
  arrange(CHROM_test, MID_POS) %>% mutate(POScumm=MID_POS+tot)

# Ninespine He
Ns.He.2$CHROM_test <- factor(Ns.He.2$CHROM_Ts, levels = paste0("chr",1:21))

Ns.He.3 <- Ns.He.2 %>% group_by(CHROM_test) %>% summarise(chrom_len=max(MID_Ts)) %>%
  mutate(tot=cumsum(chrom_len)-chrom_len) %>% select(-chrom_len) %>%
  left_join(Ns.He.2, ., by = c("CHROM_test"="CHROM_test")) %>%
  arrange(CHROM_test, MID_Ts) %>% mutate(POScumm=MID_Ts+tot)


### E: Manhattan plots ####

txt_size <- 12 # Label text size
pnt_size <- 5 # Axis tick size

### Threespine
# Find centre point of each chromosome
Ts.tmp <- Ts.Fst.3 %>% group_by(CHROM) %>% summarize(center=( max(POScumm) + min(POScumm) ) / 2 )

Ts.Manhat.Fst <- ggplot(Ts.Fst.3, aes(x=POScumm, y=Fst)) +
  
  # quality threshold
  geom_hline(aes(yintercept = quantile(Ts.Fst$Fst, 0.999)), colour = "red", lty = 2)+
  
  # Show all points
  geom_point( aes(color=CHROM_test), alpha=0.1, size=0.01) +
  scale_color_manual(values = rep(c("navy", "orange"), 11 ))+
  
  # custom X axis:
  scale_x_continuous( label = sapply(strsplit(Ts.tmp$CHROM, "r"), `[[`, 2), breaks = Ts.tmp$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=seq(0,1,0.5)) +  # remove space between plot area and x axis
  
  # Label and axis limits
  labs(title = "Threespine Stickleback (Ts)")+
  
  # Customize the theme:
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "None")

Ts.Manhat.He <- ggplot(Ts.He.2, aes(x=POScumm, y=Avg_PI)) +
  geom_hline(yintercept = quantile(Ts.He$Avg_PI, 0.99), colour = "red", lty = 2)+
  geom_point( aes(color=CHROM_test), alpha=0.1, size=0.01) +
  scale_color_manual(values = rep(c("purple", "black"), 11 ))+
  scale_x_continuous( label = sapply(strsplit(Ts.tmp$CHROM, "r"), `[[`, 2), breaks = Ts.tmp$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.02), breaks=seq(0,0.02,0.01)) +     # remove space between plot area and x axis
  labs(title = "Threespine Stickleback (Ts)")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "None")

# 86 windows with He > 0.02


### Tubesnout

# Find centre point of each chromosome
Tu.tmp <- Tu.Fst.3 %>% group_by(CHROM) %>% summarize(center=( max(POScumm) + min(POScumm) ) / 2 )

Tu.Manhat.Fst <- ggplot(Tu.Fst.3, aes(x=POScumm, y=Fst)) +
  geom_hline(aes(yintercept = quantile(Tu.Fst$Fst, 0.999)), colour = "red", lty = 2)+
  geom_point( aes(color=CHROM_test), alpha=0.1, size=0.01) +
  scale_color_manual(values = rep(c("navy", "orange"), 12 ))+
  scale_x_continuous( label = sapply(strsplit(Tu.tmp$CHROM, "r"), `[[`, 2), breaks = Tu.tmp$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=seq(0,1,0.5)) +
  labs(title = "Tubesnout (Tu)")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "None")

Tu.Manhat.He <- ggplot(Tu.He.2, aes(x=POScumm, y=Avg_PI)) +
  geom_hline(yintercept = quantile(Tu.He$Avg_PI, 0.99), colour = "red", lty = 2)+
  geom_point( aes(color=CHROM_test), alpha=0.1, size=0.01) +
  scale_color_manual(values = rep(c("purple", "black"), 12 ))+
  scale_x_continuous( label = sapply(strsplit(Tu.tmp$CHROM, "r"), `[[`, 2), breaks = Tu.tmp$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.02), breaks=seq(0,0.02,0.01)) +     # remove space between plot area and x axis
  labs(title = "Tubesnout (Tu)")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "None")

# 78 windows with He > 0.02


### Ninespine 
Ns.Fst.3$Avg_Fst <- rowSums(Ns.Fst.3[,4:7])
  
Ns.Manhat.Fst <- ggplot(Ns.Fst.3, aes(x=POScumm, y=Avg_Fst))+
  geom_hline(aes(yintercept = quantile(Ns.Fst$Avg_Fst, 0.999)), colour = "red", lty = 2)+
  geom_point( aes(color=CHROM_test), alpha=0.1, size=0.01)+
  scale_color_manual(values = rep(c("navy", "orange"), 11))+
  scale_x_continuous( label = sapply(strsplit(Ts.tmp$CHROM, "r"), `[[`, 2), breaks = Ts.tmp$center )+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=seq(0,1,0.5))+
  labs(title = "Ninespine Stickleback (Ns)")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "None")

Ns.Manhat.He <- ggplot(Ns.He.3, aes(x=POScumm, y=Avg_PI)) +
  geom_hline(yintercept = quantile(Ns.He$Avg_PI, 0.99), colour = "red", lty = 2)+
  geom_point( aes(color=CHROM_test), alpha=0.1, size=0.01) +
  scale_color_manual(values = rep(c("purple", "black"), 11 ))+
  scale_x_continuous( label = sapply(strsplit(Ts.tmp$CHROM, "r"), `[[`, 2), breaks = Ts.tmp$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.02), breaks=seq(0,02,0.01)) +     # remove space between plot area and x axis
  labs(title = "Ninespine Stickleback (Ns)")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = txt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "None")

# 57 genes He > 0.02

### Multi-panel plot
labsY <- textGrob(expression(F[ST]), rot = 90, gp = gpar(fontsize = txt_size+2))
p1 <- grid.arrange(Ts.Manhat.Fst, Ns.Manhat.Fst, Tu.Manhat.Fst, ncol = 1, left = labsY)

labsY <- textGrob(expression(bar(H[E])), rot = 90, gp = gpar(fontsize = txt_size+2))
labsX <- textGrob("Position", gp = gpar(fontsize = txt_size))
p2 <- grid.arrange(Ts.Manhat.He, Ns.Manhat.He, Tu.Manhat.He, ncol = 1, left = labsY, bottom = labsX)

Fig2 <- plot_grid(p1, p2, labels = "AUTO", label_size = txt_size*2, ncol = 1)

# Save plot
tiff("/Users/james/Dropbox (Personal)/Comp_geno_2020/Figure_2.tiff", width = 18, height = 16,
     units = "cm", res = 200)
grid.draw(Fig2)
dev.off()
