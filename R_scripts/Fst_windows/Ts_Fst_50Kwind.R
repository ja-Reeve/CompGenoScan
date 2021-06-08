############## Average Fst over 50K windows ###########################
### The following code uses a custom function to calculate the average Fst 
### between two populations of threespine stickleback.
### James Reeve - University of Calgary
### 25/09/2019

### Packages
library("dplyr")
library("ggplot2")
library("gridExtra")

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

### Download Fst per SNP
Ts.Fst <- read.table("file:///F:/Fst-table/Ts.VarScan.Fst.txt", header = TRUE)

### Download chromosome information
chrom_list <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8",
                "chr10", "chr11", "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17",
                "chr18", "chr20", "chr21")
max_length <- read.table("file:///F:/chrom_list/Ts_chr_length.txt")$V1[-c(6:7)]

### Averaging Fst across a window
# Create a function
average_Fst <-function(data, window_size, chromosome_max, chromosome_list){
  # Create row names
  res_name<-c("CHROM","START_POS", "END_POS", "Avg_Fst", "WIND")
  
  #Loop by chromosome
  Res<-list()
  Count <- 0
  for (i in chromosome_list) {
    Count <- Count + 1
    sub1 <- data[data$CHROM == chromosome_list[Count],]
    Nwinds <- ceiling(chromosome_max[Count] / window_size)
    
    results<- matrix(NA, nrow = Nwinds, ncol = 5)
    
    #loop by window
    for(t in 1:Nwinds){
      wind_end <- window_size * t
      wind_start <- 1 + window_size * (t - 1)
      sub2 <- sub1[sub1[,"POS"] <= wind_end & sub1[,"POS"] >= wind_start,]
      SNP_in_Wind <- nrow(sub2)
      
      results[t,1] <-i #chromosome number
      results[t,2] <-wind_start #window start position (bp) 
      results[t,3] <-wind_end #window end position (bp)
      results[t,4] <-mean(sub2$Fst) #avergae Fst
      results[t,5] <-paste0("Ts_", chromosome_list[Count], ":", wind_start, "-", wind_end) #Window ID
    }
    colnames(results)<-res_name
    Res[[i]]<- results
  }
  Res<-lapply(Res,na.omit) #Removing rows with NAs
  Res<-lapply(Res,as.data.frame) #Converting to a data.frame
  Res<-do.call(rbind.data.frame, Res) #Unlisting "Res"
  
  # Changing select columns to numeric data
  Res[,2:4] <- sapply(Res[,2:4],  as.numeric)
  
  return(Res)
}

# Run function
Ts.Fst.50k <- average_Fst(Ts.Fst, 50000, max_length, chrom_list)

# Set NAs to 0s (windows without SNPs)
Ts.Fst.50k[is.na(Ts.Fst.50k$Avg_Fst),"Avg_Fst"] <- 0

# Save data
write.table(Ts.Fst.50k, "file:///F:/Fst-table/Ts.Fst.50Kwind.txt", row.names = FALSE, quote = FALSE)

### Plot average Fst
# Add top candidates
TC <- read.table("file:///F:/top-candidate/Ts.TopCandidates.50Kwind.txt", header = TRUE)
TC$CAND <- TC$Observed_Outliers > TC$Expected_Outliers
TC$LOC <- paste(TC$Chromosome, TC$Window_Start, sep = ":")
Ts.Fst.50k$LOC <- paste(Ts.Fst.50k$CHROM, Ts.Fst.50k$START_POS, sep = ":")
tmp <- merge(TC, Ts.Fst.50k, by = "LOC")
tmp2 <- tmp[,c(10:13,9)]
Ts.Fst.50k <- tmp2
rm(TC,tmp,tmp2)

# Plot function
Fst.plot.chrom <- function(data, title){
  ggplot(data, aes(x = (START_POS + 25000)/1000000, y = Avg_Fst))+
    geom_line(size = 0.2, colour = "darkgrey")+
    geom_point(data = data[data$CAND == TRUE,],
               aes(x = (START_POS + 25000)/1000000, y = Avg_Fst), colour = "navy")+
    geom_hline(yintercept = quantile(Ts.Fst.50k$Avg_Fst, 0.999), colour = "red", lty = 2)+
    labs(title = title)+
    lims(y = c(0,1), x = c(0,35))+
    theme(axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line())
}

# Loop through chromosomes
Plots <- list()  
for(i in 1:length(chrom_list)){
  Plots[[i]] <- Fst.plot.chrom(data = Ts.Fst.50k[Ts.Fst.50k$CHROM == chrom_list[i],], title = chrom_list[i])
}

# Multi panel plot
Ts.pts <- grid.arrange(grobs = Plots[c(1:4,6:9,5,10:14,16:19,15,20:21)], ncol = 7, nrow = 3,
                       top = "Threespine stickleback", left = expression(Average~F[ST]), bottom = "Mid point of window (Mb)")
#ggsave("~/Documents/plots/Ts_Hexp_chromosome.png", Ts.pts, width = 15.30, height = 7.20)
