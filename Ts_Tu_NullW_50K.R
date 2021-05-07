############## Perform Null-W test between species (50Kb windows) ###########################
### The following code runs the Null-W test among genes that are top-candidates of local adaptation in two species
### This involves both identifying orthologs of top-candidates and statisically testing their strength of selection
### James Reeve - University of Calgary
### 02/11/2019

### Packages
library(ggplot2)
#library(qvalue) - Defunct package!!!

### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)

### Download data
#1 Fst
Ts.Fst <- read.table("file:///F:/Fst-table/Ts.VarScan.Fst.txt", header = TRUE)
Tu.Fst <- read.table("file:///F:/Fst-table/Tu.VarScan.Fst.txt", header = TRUE)
#2 Top-candidates
Ts.cands <- read.table("file:///F:/top-candidate/Ts.TopCandidates.50Kwind.txt", header = TRUE)
Tu.cands <- read.table("file:///F:/top-candidate/Tu.TopCandidates.50Kwind.txt", header = TRUE)
#3 Ortholog table
ortho <- read.table("file:///F:/GFF/GAC-AFL.txt", sep = "\t")[,-c(1:2,6)]
colnames(ortho) <- c("threespine_ID", "tubesnout_ID", "type")
ortho <- ortho[ortho$type == "1:1",]
#4 Gene annotations
Ts.anno <- read.table("file:///F:/GFF/Ts_gene_anno.txt", header = TRUE)
Tu.anno <- read.table("file:///F:/GFF/Tu_gene_anno.txt", header = TRUE)

### Filter annotations
Ts.anno$Gene <- paste("GAC", sapply(strsplit(Ts.anno$Gene, "[=]"), `[[`, 2), sep = "_")
Tu.anno$Gene <- paste("AFL", sapply(strsplit(Tu.anno$Gene, "[=]"), `[[`, 2), sep = "_")

colnames(Ts.anno) <- c("chromosome", "start_pos", "end_pos", "gene")
colnames(Tu.anno) <- c("chromosome", "start_pos", "end_pos", "gene")

### Rename chromosomes
old_chrom_list_Ts <- read.table("file:///F:/chrom_list/Ts_chr_list.txt")$V1[-c(6:7)]
new_chrom_list_Ts <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", "chr6", "chr7", "chr8", "chr10", "chr11",
                       "chr12", "chr13", "chr14", "chr19", "chr15", "chr16", "chr17", "chr18", "chr20", "chr21")
old_chrom_list_Tu <- read.table("file:///F:/chrom_list/Tu_chr_list.txt")$V1
new_chrom_list_Tu <- sort(paste0("chr", 1:23))

for(i in 1:length(old_chrom_list_Ts)){
  Ts.cands[Ts.cands$Chromosome == old_chrom_list_Ts[i], 1] <- new_chrom_list_Ts[i]
  Ts.anno[Ts.anno$chromosome == old_chrom_list_Ts[i], 1] <- new_chrom_list_Ts[i]
};rm(i)

for(i in 1:length(old_chrom_list_Tu)){
  Tu.cands[Tu.cands$Chromosome == old_chrom_list_Tu[i], 1] <- new_chrom_list_Tu[i]
  Tu.anno[Tu.anno$chromosome == old_chrom_list_Tu[i], 1] <- new_chrom_list_Tu[i]
}; rm(i)

### Label genes as candidates if they are in candidate window

GeneAsign <- function(top.candidate.data, gene.annotation){
  tmp <- list()
  
  for(i in 1:nrow(top.candidate.data)){
    Chrom <- top.candidate.data[i, "Chromosome"]
    Spos <- top.candidate.data[i, "Window_Start"]
    Epos <- top.candidate.data[i, "Window_End"]
    Cand <- top.candidate.data[i, "Observed_Outliers"] > top.candidate.data[i, "Expected_Outliers"]
    
    dat <- gene.annotation[gene.annotation$chromosome == Chrom &
                             gene.annotation$start_pos >= Spos &
                             gene.annotation$end_pos <= Epos, ]
    if(nrow(dat) == 0)next
    dat$cand <- rep(Cand, nrow(dat))
    dat$window <- rep(top.candidate.data[i, "Window_Name"], nrow(dat))
    
    tmp[[i]] <- dat
    rm(Chrom, Spos, Epos, Cand, dat)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  return(tmp2)
}

Ts.gene <- GeneAsign(Ts.cands, Ts.anno)
Tu.gene <- GeneAsign(Tu.cands, Tu.anno)

### Assign genes in between windows to window with longer section of gene

GeneAsign2 <- function(top.candidate.data, gene.annotation, assinged.genes){
  # Find genes between windows
  btw <- gene.annotation[!(gene.annotation[,4] %in% assinged.genes[,4]),]
  
  # Empty dataset
  tmp <- list()
  
  # For loop assinging candidate info to 'btw'
  for(i in 1:nrow(btw)){
    Chrom <- btw[i, "chromosome"]
    Spos <- btw[i, "start_pos"]
    Epos <- btw[i, "end_pos"]
    
    Pos <- Spos:Epos
    
    dat <- top.candidate.data[top.candidate.data$Chromosome == Chrom,]
    dat2 <- dat[dat$Window_Start %in% Pos | dat$Window_End %in% Pos, ]
    if(nrow(dat2) != 2)next
    ### !!! Note: this is a temproray fix. Several genes are missing that are in the last window of each chromosome in the threespine!!! ###
    
    wind1 <- abs(Spos - dat2[1, "Window_End"])
    wind2 <- abs(Epos - dat2[2, "Window_Start"])
    
    if(wind1 > wind2){
      tmp[[i]] <- data.frame("chromosome" = Chrom, 
                             "start_pos" = Spos, 
                             "end_pos" = Epos, 
                             "gene" = btw[i, 4],
                             "cand" = dat2[2,"Observed_Outliers"] > dat2[2,"Expected_Outliers"],
                             "window" = dat2[2,"Window_Name"])
    } else {
      tmp[[i]] <- data.frame("chromosome" = Chrom, 
                             "start_pos" = Spos, 
                             "end_pos" = Epos, 
                             "gene" = btw[i, 4],
                             "cand" = dat2[1,"Observed_Outliers"] > dat2[1,"Expected_Outliers"],
                             "window" = dat2[1,"Window_Name"])
    }
    rm(Chrom, Spos, Epos, Pos, dat, dat2, wind1, wind2)
  }; rm(i)
  
  tmp2 <- do.call(rbind.data.frame, tmp)
  
  # Join genes to end of previously assinged genes
  tmp3 <- rbind(assinged.genes, tmp2)
  
  # sort
  tmp4 <- tmp3[with(tmp3, order(chromosome, start_pos)), ]
  
  return(tmp4) 
}

Ts.gene2 <- GeneAsign2(Ts.cands, Ts.anno, Ts.gene)
Tu.gene2 <- GeneAsign2(Tu.cands, Tu.anno, Tu.gene)

### Check for how many genes are missing!!!
1 - (nrow(Ts.gene2) / nrow(Ts.anno))
1- (nrow(Tu.gene2) / nrow(Tu.anno))

#################### Null- W test ####################

### Create function

Null_W.test <- function(top_candidate, ortholog_annotation, ortholog_Fst, ortholog_table, species_identifier){
  if(length(species_identifier) != 2) {print("Error: please specify two strings in species identifiers, ortholog 2nd")}
  if(class(top_candidate) != "data.frame") {print("Error: top_candidate must be a data.frame")}
  if(class(ortholog_annotation) != "data.frame") {print("Error: ortholog_annotation must be a data.frame")}
  if(class(ortholog_Fst) != "data.frame") {print("Error: ortholog_Fst must be a data.frame")}
  if(class(ortholog_table) != "data.frame") {print("Error: ortholog_table must be a data.frame")}
  
  #1. Add ortholog positions to top-canidate da
  tmp <- merge(top_candidate, ortholog_table, by.x = "gene", by.y = species_identifier[1], all = FALSE)
  tmp2 <- merge(tmp, ortholog_annotation, by.x = species_identifier[2], by.y = "gene")[, -c(6,8,12)]
  colnames(tmp2) <- c(species_identifier[2], species_identifier[1], "cand_chromosome", "cand_start", "candidate_end", 
                      "cand_window", "ortho_chromosome", "ortho_start", "ortho_end", "ortho_window")
  
  #2. Split data into canidates and non-candidate sets
  cand_ortho <- tmp[tmp$cand == TRUE, c(1:7)]
  non_cand_ortho <- tmp[tmp$cand == FALSE, c(1:7)]
  
  #3. Create subset of 10,000 randomly drawn Fst scores
  set.seed(1000)
  SNP10k <- ortholog_Fst[sample(nrow(ortholog_Fst), 10000), ]
  
  #4. Loop through each gene to calculate difference to SNP10k (Wilcoxon Ranked Sum Test)
  Res <- matrix(0, nrow = nrow(tmp2), ncol = 7)
  
  rm(gene)
  for(gene in 1:nrow(tmp2)){
    #Extract ortholog information
    dat <- tmp2[gene,]
    
    gene_start <- dat[,9] #store start position
    gene_end <- dat[,10] #store end position
    gene_chrom <- dat[,8] #store chromosome
    
    #Find Fst scores for ortholog
    dat2 <- ortholog_Fst[ortholog_Fst$CHROM == gene_chrom &
                           ortholog_Fst$POS >= gene_start & 
                           ortholog_Fst$POS <= gene_end,]
    if(nrow(dat2) == 0)next #Skip any genes without SNPs
    
    #Wilcoxon Ranked Sum Test
    Wscores <- wilcox.test(dat2$Fst, SNP10k$Fst)
    
    #Adjust W-scores for large sample size of test
    Zscore <- (2*Wscores$statistic - nrow(dat2)*10000) / sqrt(nrow(dat2)*10000*(nrow(dat2)+10000+1)/3)
    
    #Save results
    Res[gene,1] <- dat[,2] #Candidate ID
    Res[gene,2] <- dat[,1] #Ortholog ID
    Res[gene,3] <- gene_chrom #ortholog chromosome
    Res[gene,4] <- gene_start #ortholog start
    Res[gene,5] <- gene_end #ortholog end
    Res[gene,6] <- Zscore #Z statistic
    Res[gene,7] <- Wscores$p.value #Raw P-value
    
    #Clean up
    rm(dat, dat2, gene_chrom, gene_end, gene_start, Wscores, Zscore)
  }
  
  #4.5. Convert results into dataframe
  Res <- as.data.frame(Res)
  Res[, 4:7] <- sapply(Res[,4:7], as.numeric)
  colnames(Res) <- c("candidate_ID", "ortholog_ID", "ortholog_chrom", "ortholog_start", "ortholog_end",
                     "Zscore", "p_value")
  
  #5. Remove any rows missing information (i.e. no SNPs in gene)
  Res <- Res[Res$candidate_ID != "0",]
  
  #6. Split cands and non-cands 
  #####(ask Sam; will it effect sample size of test!) ####
  Wilcox_cand <- Res[Res$candidate_ID %in% cand_ortho$Gene_Name, ]
  Wilcox_control <- Res[Res$candidate_ID %in% non_cand_ortho$Gene_Name, ]
  
  #7. Calculate emperical P and correct FDR
  Wilcox_cand$emperical_Pvalue <- empPvals(Wilcox_cand$Zscore, Wilcox_control$Zscore)
  #Wilcox_cand$q_value <- qvalue(Wilcox_cand$emperical_Pvalue)$qvalue
  
  #Clean up and return results
  NullW <- list(Wilcox_cand, Wilcox_control)
  return(NullW)
  
  rm(Res, SNP10k, cand_ortho, non_cand_ortho, tmp, tmp2)
}

### Run Null-W test threespine stickleback to tubesnout  
Ts_Tu <- Null_W.test(top_candidate = Ts.gene2, ortholog_annotation = Tu.gene2, ortholog_Fst = Tu.Fst,
                     ortholog_table = ortho, species_identifier = c("threespine_ID", "tubesnout_ID"))
write.table(Ts_Tu[[1]], "file:///F:/Null-W/Ts_Tu_NullW_50K_cands.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Ts_Tu[[2]], "file:///F:/Null-W/Ts_Tu_NullW_50K_control.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

### Run Null-W test tubesnout to threespine stickleback
Tu_Ts <- Null_W.test(top_candidate = Tu.gene2, ortholog_annotation = Ts.gene2, ortholog_Fst = Ts.Fst,
                     ortholog_table = ortho, species_identifier = c("tubesnout_ID", "threespine_ID"))

write.table(Tu_Ts[[1]], "file:///F:/Null-W/Tu_Ts_NullW_50K_cands.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Tu_Ts[[2]], "file:///F:/Null-W/Tu_Ts_NullW_50K_control.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")