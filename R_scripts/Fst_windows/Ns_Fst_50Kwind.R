############## Average Fst over 50K windows ###########################
### The following code uses a custom function to calculate the average Fst 
### among four populations of ninespine stickleback.
### James Reeve - University of Calgary
### 25/11/2019

### Packages
library("dplyr")
library("ggplot2")
library("gridExtra")

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

### Download Fst per SNP
Ns.Fst <- read.table("file:///G:/Fst-table/Ns.VarScan.Fst.txt", header = TRUE)

### Download chromosome information

# Maximum length per chromosome
chrom_list <- read.table("file:///G:/chrom_list/Ns_chr_list.txt")$V1
max_length <- read.table("file:///G:/chrom_list/Ns_chr_length.txt")$V1

# Remove contigs without genes
gene_anno <- read.table("file:///G:/GFF/Ns_gene_anno.txt", header = TRUE)
tmp <- chrom_list %in% gene_anno$CHROM
chrom_list2 <- chrom_list[tmp]
max_length2 <- max_length[tmp]
rm(tmp, gene_anno)

### Averaging Fst across a window
# Create a function
average_Fst <-function(data, window_size, chromosome_max, chromosome_list){
  # Create row names
  res_name<-c("CHROM","START_POS", "END_POS", "NS1_Fst", "NS2_Fst", "WIND")
  
  #Loop by chromosome
  Res<-list()
  Count <- 0
  for (i in chromosome_list) {
    Count <- Count + 1
    sub1 <- data[data$CHROM == chromosome_list[Count],]
    Nwinds <- ceiling(chromosome_max[Count] / window_size)
    
    results<- matrix(NA, nrow = Nwinds, ncol = 6)
    
    #loop by window
    for(t in 1:Nwinds){
      wind_end <- window_size * t
      wind_start <- 1 + window_size * (t - 1)
      sub2 <- sub1[sub1[,"POS"] <= wind_end & sub1[,"POS"] >= wind_start,]
      SNP_in_Wind <- nrow(sub2)
      
      results[t,1] <-i #chromosome number
      results[t,2] <-wind_start #window start position (bp) 
      results[t,3] <-wind_end #window end position (bp)
      results[t,4] <-mean(rowMeans(sub2[,4:5])) #avergae Fst NsNUn vs Southern Ns
      results[t,5] <-mean(rowMeans(sub2[,6:7])) #avergae Fst NsNUd vs Southern Ns
      results[t,6] <-paste0("Ns_", chromosome_list[Count], ":", wind_start, "-", wind_end) #Window ID
    }
    colnames(results)<-res_name
    Res[[i]]<- results
  }
  Res<-lapply(Res,na.omit) #Removing rows with NAs
  Res<-lapply(Res,as.data.frame) #Converting to a data.frame
  Res<-do.call(rbind.data.frame, Res) #Unlisting "Res"
  
  # Changing select columns to numeric data
  Res[,2:5] <- sapply(Res[,2:5],  as.numeric)
  
  return(Res)
}

# Run function
Ns.Fst.50k <- average_Fst(Ns.Fst, 50000, max_length2, chrom_list2)

# Remove
Ns.Fst.50k <- Ns.Fst.50k[complete.cases(Ns.Fst.50k),]

# Save data
write.table(Ns.Fst.50k, "file:///G:/Fst-table/Ns.Fst.50Kwind.txt", row.names = FALSE, quote = FALSE)

### Plot average Fst
  
contigs <- unique(Ns.Fst.50k$CHROM)
tmp <- list()

spos <- 1
for( i in 1: length(contigs)){
  dat <- Ns.Fst.50k[Ns.Fst.50k$CHROM == contigs[i],]
  
  dat$MID_POS <- dat$END_POS - 25000 + spos
  
  tmp[[i]] <- dat
  
  mpos <- tail(dat$MID_POS, 1)
  spos <- mpos + 25000
  
}; rm(i, spos, mpos)

Ns.Fst.50k.2 <- do.call(rbind.data.frame, tmp)

#plot function
Fst.plot.chrom <- function(data){
  ggplot(data)+
    geom_line(aes(x = MID_POS / 1000000, y = NS1_Fst), size = 0.2, colour = "black")+
    geom_line(aes(x = MID_POS / 1000000, y = NS2_Fst), size = 0.2, colour = "darkgrey")+
    geom_hline(yintercept = quantile(Ns.Fst.50k$NS1_Fst, 0.999), colour = "red", lty = 2)+
    geom_hline(yintercept = quantile(Ns.Fst.50k$NS2_Fst, 0.999), colour = "red", lty = 3)+
    #labs(title = title)+
    lims(y = c(0,1))+
    theme(axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line())
}

# Loop through chromosomes
Plots <- list()  
tmp2 <- seq(1, nrow(Ns.Fst.50k.2), 3000)[-1]
for(i in 1:length(tmp2)){
  Plots[[i]] <- Fst.plot.chrom(data = Ns.Fst.50k.2[(tmp2[i]-3000):tmp2[i],])
}; rm(i)

# Multi panel plot
Ns.pts <- grid.arrange(grobs = Plots, ncol = 1, nrow = 4,
                       top = "Ninespine stickleback", left = expression(Average~F[ST]), bottom = "Genome position (Mb)")
#ggsave("~/Documents/plots/Ns_Hexp_chromosome.png", Ts.pts, width = 15.30, height = 7.20)

