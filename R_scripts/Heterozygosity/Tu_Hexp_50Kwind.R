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
Tu.vcf <- read.table("file:///F:/vcf-files/Tu.VarScan.vcf", skip =  23) # Skip 23 lines of header
colnames(Tu.vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TuAK", "TuBC")
Tu.vcf <- Tu.vcf[Tu.vcf$CHROM != "chrM" & Tu.vcf$CHROM != "chrUn",]

### Update chromosome names
# Old names and chromosome lengths
chrom_list <- read.table("file:///F:/chrom_list/Tu_chr_list.txt")$V1
max_length <- read.table("file:///F:/chrom_list/Tu_chr_length.txt")$V1

# Updating names
new_chrom_list <- paste0("chr", 1:23)
for(i in 1:length(chrom_list)){
  Tu.vcf[Tu.vcf$CHROM == chrom_list[i], "CHROM"] <- new_chrom_list[i]
}
chrom_list <- new_chrom_list
rm(new_chrom_list, i)

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
      results[j,5] <- paste0("Tu_", list.chromosomes[i], ":", wind_start, "-", wind_end) # Window ID
      
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
# Northen population
TuAK.He <- nucleotide.diversity.window(data = Tu.vcf, sample.name = "TuAK", window.size = 50000,
                                     list.chromosomes = chrom_list, chromosome.max.length = max_length)
# Southern population
TuBC.He <- nucleotide.diversity.window(data = Tu.vcf, sample.name = "TuBC", window.size = 50000,
                                     list.chromosomes = chrom_list, chromosome.max.length = max_length)
# Merge into single dataset
Tu.He <- data.frame("CHROM" = TuAK.He$CHROM, 
                    "START_POS" = TuAK.He$START_POS, "END_POS" = TuAK.He$END_POS,
                    "WIND" = TuAK.He$WIND, 
                    "TuAK_He" = TuAK.He$He, "TuBC_He" = TuBC.He$He,
                    "Avg_He" = rowMeans(cbind(TuAK.He$He, TuBC.He$He)))

###### Note: last window for each chromosome is underestimate becasue of incomplete window ###
### Save data
write.table(Tu.He, "file:///F:/nucleotide_diversity/Tu_Avg_Hexp_50Kwind.txt", row.names = FALSE, quote = FALSE)


### Plots
#Plot TsAK vs TsOR
Tu.He$Outlier <- Tu.He$Avg_He < quantile(Tu.He$Avg_He, 0.01) | Tu.He$Avg_He > quantile(Tu.He$Avg_He, 0.99)

ggplot(Tu.He)+
  geom_point(aes(x = TuAK_He, y = TuBC_He, colour = Outlier), size = 0.1)+
  geom_hline(yintercept = quantile(Tu.He$TuBC_He, 0.99), lty = 2, colour = "firebrick")+
  geom_hline(yintercept = quantile(Tu.He$TuBC_He, 0.01), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Tu.He$TuAK_He, 0.99), lty = 2, colour = "firebrick")+
  geom_vline(xintercept = quantile(Tu.He$TuAK_He, 0.01), lty = 2, colour = "firebrick")+
  labs(x = expression(TuAK~bar(H)[E]), y = expression(TuBC~bar(H)[E]))+
  lims(y = c(0, 0.013), x = c(0, 0.013))+
  scale_colour_manual(values = c("grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

# Plot function
He.plot.chrom <- function(data, title){
  ggplot(data, aes(x = (START_POS + 25000)/1000000, y = Avg_He))+
    geom_line(size = 0.2, colour = "darkgrey")+
    geom_point(data = data[data$Avg_He == 0,],
               aes(x = (START_POS + 25000)/1000000, y = Avg_He), pch = 18, colour = "black")+
    geom_hline(yintercept = quantile(data$Avg_He, 0.99), colour = "red", lty = 2)+
    geom_hline(yintercept = quantile(data$Avg_He, 0.01), colour = "red", lty = 2)+
    labs(title = title)+
    lims(y = c(0,0.01), x = c(0,32.2))+
    theme(axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line())
}


# Loop through chromosomes
Plots <- list()
for(i in 1:length(chrom_list)){
  Plots[[i]] <- He.plot.chrom(data = Tu.He[Tu.He$CHROM == chrom_list[i],], title = chrom_list[i])
}

# Multi panel plot
Tu.pts <- grid.arrange(grobs = Plots, ncol = 6, nrow = 4,
                       top = "Tubesnout", left = "Genetic diversity (He)", bottom = "Mid point of window (Mb)")
#ggsave("~/Documents/plots/Ts_Hexp_chromosome.png", Ts.pts, width = 15.30, height = 7.20)
