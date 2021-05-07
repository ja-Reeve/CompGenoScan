############## Identify Top Candidates for 50Kb windows ###########################
### The following code uses a custom function to identify top candidate loci related 
### latitute adaptation between two populations of threespine stickleback.
### The calculations are based of Fst values calculated by PoolFstat (Hivert et al. 2018).
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
Ts.Fst <- read.table("file:///F:/Fst-table/Ts.VarScan.Fst.txt", header = TRUE)

### Read in the gene annotation data
gene_anno <- read.table("file:///F:/GFF/Ts_gene_anno.txt", header = TRUE)

### Add rows for intragenetic regions
# Maximum length per chromosome
chrom_list <- read.table("file:///F:/chrom_list/Ts_chr_list.txt")$V1[-c(6:7)]
max_length <- read.table("file:///F:/chrom_list/Ts_chr_length.txt")$V1[-c(6:7)]

### Update chromosome names
new_chrom_list <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8",
                    "chr10", "chr11", "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17",
                    "chr18", "chr20", "chr21")
for(i in 1:length(chrom_list)){
  gene_anno[gene_anno$CHROM == chrom_list[i], 1] <- new_chrom_list[i]
}
chrom_list <- new_chrom_list
rm(new_chrom_list, i)

### Top Candidate Function
top_candidate_test<-function(data, window_size, binomial_cutoff, threshold, chromosome_max){
  # Create row names
  res_name<-c("Chromosome","Window_Start", "Window_End", "Window_Name", "SNPs_in_Window","Observed_Outliers")
  
  # Calculate the Fst value at the quantile specified by the threshold parameter
  Fst_threshold<-quantile(data$Fst, threshold)

  #Loop by chromosome
  Res<-list()
  Count <- 0
  for (i in 1:length(chrom_list)) {
    sub1 <- data[data[, "CHROM"] == chrom_list[i],]
    Nwinds <- ceiling(chromosome_max[i] / window_size)
    
    results<- matrix(NA, nrow = Nwinds, ncol = 6)
    
    #loop by gene
    for(t in 1:Nwinds){
      wind_end <- window_size * t
      wind_start <- 1 + window_size * (t - 1)
      sub2 <- sub1[sub1[,"POS"] <= wind_end & sub1[,"POS"] >= wind_start,]
      SNP_in_Wind <- nrow(sub2)
      
      results[t,1] <-chrom_list[i] #chromosome number
      results[t,2] <-wind_start #window start position (bp) 
      results[t,3] <-wind_end #window end position (bp)
      results[t,4] <-paste0("Ts_", chrom_list[i], ":", wind_start, "-", wind_end) #Window ID
      results[t,5] <-SNP_in_Wind #number of SNPs in window
      results[t,6] <-sum(sub2[,"Fst"] > Fst_threshold) # the number of outliers in the window
    }
    colnames(results)<-res_name
    Res[[i]]<- results
  }
  Res<-lapply(Res,na.omit) #Removing rows with NAs
  Res<-lapply(Res,as.data.frame) #Converting to a data.frame
  Res<-do.call(rbind.data.frame, Res) #Unlisting "Res"
  
  # Changing select columns to numeric data
  Res[,c(2:3, 5:6)] <- sapply(Res[,c(2:3, 5:6)],  as.numeric)
  
  # calculate expected number of ouliers using bionomial distribution
  Pout <- sum(Res[,6]) / sum(Res[,5]) #probability of outliers
  Res$Expected_Outliers <- qbinom(binomial_cutoff, Res[,5], Pout)
  
  return(Res)
}

### Run top candidate test
Cands <- top_candidate_test(Ts.Fst, 50000, 0.9999, 0.999, max_length)

### Save the data
write.table(Cands, "file:///F:/top-candidate/Ts.TopCandidates.50Kwind.txt", row.names = FALSE, quote = FALSE)

### Plot top candidates
ggplot(Cands %>% arrange(SNPs_in_Window))+
  geom_point(aes(SNPs_in_Window, Observed_Outliers, 
                 colour = Observed_Outliers > Expected_Outliers), size = 0.1)+
  geom_line(aes(SNPs_in_Window, Expected_Outliers), colour="firebrick")+
  labs(x="SNPs in Gene", y="Observed Outliers")+
  ggtitle("Threespine Stickleback: B = 0.9999 Q = 999th")+
  scale_colour_manual(values = c("TRUE" = "orange", "FALSE" = "navy"))+
  theme_bw()+
  ylim(c(0,500))+
  xlim(c(0, 1500))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
