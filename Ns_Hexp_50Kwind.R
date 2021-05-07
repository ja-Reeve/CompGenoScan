############## Calculate Genetic Diversity 50K windows ###########################
### The following code uses a custom function to calculate the average genetic 
### diversity of two populations of tubesnout.
### Genetic diversity within each population is estimated as the expected 
### heterozygosity across a 50,000bp windows, assuming each non-SNP = 0
### James Reeve - University of Calgary
### 17/09/2019

### Packages
library("dplyr")
library("ggplot2")
library("gridExtra")

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

### Download VCF
Ns.vcf <- read.table("file:///G:/vcf-files/Ns_FULL.VarScan.vcf", skip =  23) # Skip 23 lines of header
colnames(Ns.vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", 
                      "NsNUn", "NsNUd", "NsABk", "NsABm")

### Filter multi-allelic sites
tmp1 <- duplicated(Ns.vcf[,1:2])
sum(tmp1)

### Read in the gene annotation data
gene_anno <- read.table("file:///G:/GFF/Ns_gene_anno.txt", header = TRUE)

### Add rows for intragenetic regions
# Maximum length per chromosome
chrom_list <- read.table("file:///G:/chrom_list/Ns_chr_list.txt")$V1
max_length <- read.table("file:///G:/chrom_list/Ns_chr_length.txt")$V1

# Remove contigs without genes
tmp <- chrom_list %in% gene_anno$CHROM
chrom_list2 <- chrom_list[tmp]
max_length2 <- max_length[tmp]
rm(tmp)

### Function to calcualte expected heterozygosity over 50Kb window (each population)
nucleotide.diversity.window <- function(data, sample.name, window.size, list.chromosomes, chromosome.max.length){
  ### Error messages
  if(!(is.character(sample.name))){stop("sample.name is not a string")}
  if(data$FORMAT[1] != "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"){
    stop("VCF is not in VarScan format; check the FORMAT column = GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR")}
  if(length(list.chromosomes) != length(chromosome.max.length)){
    stop("chromosome.list and chromosome.max.length are unequal lengths")}
  
  ### Calculate Expected heterozygosity per SNP
  # Get vcf entry for a single population
  tmp <- strsplit(data[,sample.name], ":")
  # Filter out any fixed SNP (i.e. ./.)
  logi <- lapply(tmp, length) != 3
  tmp2 <- tmp[logi]
  data2 <- data[logi, ]
  
  # Get alternative allele frequency from VCF
  AF.q <- as.numeric(unlist(strsplit(sapply(tmp2, `[[`, 7), "%"))) / 100
  # Expected heterozygosity
  Hexp <- 2*AF.q*(1-AF.q)
  # Make dataframe
  tmp3 <- data.frame("CHROM" = data2$CHROM,
                     "POS" = data2$POS,
                     "Heter_Exp" = Hexp)
  
  res_name <- c("CHROM", "START_POS", "END_POS", "He", "WIND")
  ### Average expected heterozygotisty across a window
  ### Non-SNPs: Hexp = 0
  Res<-list()
  for (i in 1:length(list.chromosomes)) {
    sub1 <- tmp3[tmp3$CHROM == list.chromosomes[i],]
    
    Nwinds <- ceiling(chromosome.max.length[i] / window.size)
    results<- matrix(0, nrow = Nwinds, ncol = 5)
    
    # Loop by window
    for(j in 1:Nwinds){
      wind_end <- window.size * j
      wind_start <- 1 + window.size * (j - 1)
      sub2 <- sub1[sub1$POS <= wind_end & sub1$POS >= wind_start,]
      
      N_fixed <- window.size - nrow(sub2) #number of non-SNPs in gene
      zeros <- rep(0, times = N_fixed) #creates list of 0s for non-SNPs
      Hexp_list <- mean(c(sub2$Heter_Exp, zeros))#append(sub2$Heter_Exp, zeros) #joins list of 0s to end of Hexp scores
      
      results[j,1] <- list.chromosomes[i] # Chromosome
      results[j,2] <- wind_start # Window start position (bp)
      results[j,3] <- wind_end # Window end position (bp)
      results[j,4] <- Hexp_list # Average expected hterozygosity
      results[j,5] <- paste0("Ns_", list.chromosomes[i], ":", wind_start, "-", wind_end) # Window ID
      
      rm(wind_start, wind_end, N_fixed, zeros, sub2, Hexp_list)
    }
    colnames(results) <- res_name
    Res[[i]] <- results
    
    rm(sub1, Nwinds, results)
  }; rm(i,j)
  
  Res<-do.call(rbind.data.frame, Res) #Unlisting "Res"
  Res <- Res[Res$CHROM != 0,]
  
  # Changing select columns to numeric data
  Res[,2:4] <- sapply(Res[,2:4], as.numeric)
  
  return(Res) 
}

### Calculate average expected heterozygosity
# Northen population 1
NsNUn.He <- nucleotide.diversity.window(data = Ns.vcf, sample.name = "NsNUn", window.size = 50000,
                                      list.chromosomes = chrom_list2, chromosome.max.length = max_length2)
# Northen population 2
NsNUd.He <- nucleotide.diversity.window(data = Ns.vcf, sample.name = "NsNUd", window.size = 50000,
                                      list.chromosomes = chrom_list2, chromosome.max.length = max_length2)
# Southern population 1
NsABk.He <- nucleotide.diversity.window(data = Ns.vcf, sample.name = "NsABk", window.size = 50000,
                                      list.chromosomes = chrom_list2, chromosome.max.length = max_length2)
# Southern population 2
NsABm.He <- nucleotide.diversity.window(data = Ns.vcf, sample.name = "NsABm", window.size = 50000,
                                      list.chromosomes = chrom_list2, chromosome.max.length = max_length2)
# Merge into single dataset
Ns.He <- data.frame("CHROM" = NsNUn.He$CHROM, 
                    "START_POS" = NsNUn.He$START_POS, "END_POS" = NsNUn.He$END_POS,
                    "WIND" = NsNUn.He$WIND, 
                    "NsNUn_He" = NsNUn.He$He, "NsNUd_He" = NsNUd.He$He,
                    "NsABk_He" = NsABk.He$He, "NsABm_He" = NsABm.He$He, 
                    "Avg_He" = rowMeans(cbind(NsNUn.He$He, NsNUd.He$He, NsABk.He$He, NsABm.He$He)))

# Remove windows without SNPs
Ns.He2 <- Ns.He[!(Ns.He$Avg_He == 0 & Ns.He$START_POS == 1),]

###### Note: last window for each chromosome is underestimate becasue of incomplete window ###
##### Note 2: Too many contigs have low He ######

### Save data
write.table(Ns.He, "file:///G:/nucleotide_diversity/Ns_Avg_Hexp_50Kwind.txt", row.names = FALSE, quote = FALSE)

### Plots
contigs <- unique(Ns.He2$CHROM)
tmp <- list()

spos <- 1
for( i in 1: length(contigs)){
  dat <- Ns.He2[Ns.He2$CHROM == contigs[i],]
  
  dat$MID_POS <- dat$END_POS - 25000 + spos
  
  tmp[[i]] <- dat
  
  mpos <- tail(dat$MID_POS, 1)
  spos <- mpos + 25000
  
}; rm(i, spos, mpos)

Ns.He3 <- do.call(rbind.data.frame, tmp)

#plot function
He.plot.chrom <- function(data){
  ggplot(data, aes(x = MID_POS / 1000000, y = Avg_He))+
    geom_line(size = 0.2, colour = "darkgrey")+
    geom_point(data = data[data$Avg_He == 0,],
               aes(x = MID_POS /1000000, y = Avg_He), pch = 18, colour = "black")+
    geom_hline(yintercept = quantile(data$Avg_He, 0.99), colour = "red", lty = 2)+
    geom_hline(yintercept = quantile(data$Avg_He, 0.01), colour = "red", lty = 2)+
    #labs(title = title)+
    lims(y = c(0,0.01))+
    theme(axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line())
}

# Loop through chromosomes
Plots <- list()  
tmp2 <- seq(1, nrow(Ns.He3), 3500)[-1]
for(i in 1:length(tmp2)){
  Plots[[i]] <- He.plot.chrom(data = Ns.He3[(tmp2[i]-3500):tmp2[i],])
}; rm(i)

# Multi panel plot
Ns.pts <- grid.arrange(grobs = Plots, ncol = 1, nrow = 4,
                       top = "Ninespine stickleback", left = "Genetic diversity (He)", bottom = "Genome position (Mb)")
#ggsave("~/Documents/plots/Ns_Hexp_chromosome.png", Ts.pts, width = 15.30, height = 7.20)


#Plot TsAK vs TsOR
Ns.He2$Outlier <- Ns.He2$Avg_He < quantile(Ns.He2$Avg_He, 0.01) | Ns.He2$Avg_He > quantile(Ns.He2$Avg_He, 0.99)

#NsNUn vs NsNUd
ggplot(Ns.He2)+
  geom_point(aes(x = NsNUn_He, y = NsNUd_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Ns.He2$NsNUd_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He2$NsNUd_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUn_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUn_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(NsNUn~bar(H)[E]), y = expression(NsNUd~bar(H)[E]))+
  lims(y = c(0, 0.015), x = c(0, 0.015))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")
#NsNUn vs NsABm
ggplot(Ns.He2)+
  geom_point(aes(x = NsNUn_He, y = NsABm_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Ns.He2$NsABm_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He2$NsABm_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUn_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUn_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(NsNUn~bar(H)[E]), y = expression(NsABm~bar(H)[E]))+
  lims(y = c(0, 0.015), x = c(0, 0.015))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")
#NsNUn vs NsABk
ggplot(Ns.He2)+
  geom_point(aes(x = NsNUn_He, y = NsABk_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Ns.He2$NsABk_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He2$NsABk_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUn_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUn_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(NsNUn~bar(H)[E]), y = expression(NsABk~bar(H)[E]))+
  lims(y = c(0, 0.015), x = c(0, 0.015))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")
#NsNUd vs NsABm
ggplot(Ns.He2)+
  geom_point(aes(x = NsNUd_He, y = NsABm_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Ns.He2$NsABm_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He2$NsABm_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUd_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUd_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(NsNUd~bar(H)[E]), y = expression(NsABm~bar(H)[E]))+
  lims(y = c(0, 0.015), x = c(0, 0.015))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")
#NsNUd vs NsABk
ggplot(Ns.He2)+
  geom_point(aes(x = NsNUd_He, y = NsABk_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Ns.He2$NsABk_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He2$NsABk_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUd_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsNUd_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(NsNUd~bar(H)[E]), y = expression(NsABk~bar(H)[E]))+
  lims(y = c(0, 0.015), x = c(0, 0.015))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")
#NsABm vs NsABk
ggplot(Ns.He2)+
  geom_point(aes(x = NsABm_He, y = NsABk_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Ns.He2$NsABk_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Ns.He2$NsABk_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsABm_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Ns.He2$NsABm_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(NsABm~bar(H)[E]), y = expression(NsABk~bar(H)[E]))+
  lims(y = c(0, 0.015), x = c(0, 0.015))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")


