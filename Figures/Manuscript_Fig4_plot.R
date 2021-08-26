################# Fig. 4 Con Evo Project ##########################
### The following code was used to generate Figure 4 in my manuscript
### This figure shows a two panel plot of the Null-W test A) between threespine top candidates and
### tubesnout orthologs (Ts_Tu) and B) between tubesnout top candidates and threespine orthologs
### James Reeve - University of Calgary
### 08/10/2019

### A: Preperation ####
rm(list = ls())
options(stringsAsFactors =  FASLE)
dev.off()

### Packages
library("ggplot2")
library("grid")
library("gridExtra")
library("cowplot")

### Plot size parameters
txt_size <- 6 # Label text size
pnt_size <- 5 # Axis tick size

### B: Access data ####

### Top candidate data
Cands.Ts <- read.table("file:///F:/top-candidate/Ts.VarScan.TopCandidates.txt", header = TRUE)
Cands.Tu <- read.table("file:///F:/top-candidate/Tu.VarScan.TopCandidates.txt", header = TRUE)
# Remove non-genes
Cands.Ts <- Cands.Ts[Cands.Ts$Gene_Name != "itgr",]
Cands.Tu <- Cands.Tu[Cands.Tu$Gene_Name != "itgr",]

### Null-W results
Ts_Tu_cands <- read.table("file:///F:/Null-W/Ts_Tu_NullW_cands.txt", header = TRUE, sep = "\t")
Ts_Tu_control <- read.table("file:///F:/Null-W/Ts_Tu_NullW_control.txt", header = TRUE, sep = "\t")
Tu_Ts_cands <- read.table("file:///F:/Null-W/Tu_Ts_NullW_cands.txt", header = TRUE, sep = "\t")
Tu_Ts_control <- read.table("file:///F:/Null-W/Tu_Ts_NullW_control.txt", header = TRUE, sep = "\t")

### C: Plot top candidate plots
### Plot
cands_plot <- function(data, species){
  ggplot(data %>% arrange(SNPs_in_Gene))+
    geom_point(aes(SNPs_in_Gene, Observed_Outliers, 
                   colour = Observed_Outliers > Expected_Outliers), 
               size = 0.4)+
    geom_line(aes(SNPs_in_Gene, Expected_Outliers), colour="firebrick")+
    ggtitle(paste0(species, ": B = 0.9999 Q = 999th"))+
    scale_colour_manual(values = c("TRUE" = "orange", "FALSE" = "navy"))+
    theme_bw()+
    ylim(c(0,120))+
    xlim(c(0, 1800))+
    theme(axis.title = element_blank(),
          plot.title = element_text(size = txt_size-2),
          axis.text = element_text(size = pnt_size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}

Plot.Ts <- cands_plot(Cands.Ts, "Threespine stickleback")
Plot.Tu <- cands_plot(Cands.Tu, "Tubesnout")

### Combine into multi-panel plot
tmp_plot <- plot_grid(Plot.Ts, Plot.Tu, labels = c("A", "B"), label_size = txt_size*1.5)
labsX <- textGrob("SNPs in Gene", gp=gpar(fontface = "bold", fontsize = txt_size))
labsY <- textGrob("Observed Outliers", gp=gpar(fontface="bold", fontsize = txt_size), rot = 90)

Fig4a <- grid.arrange(arrangeGrob(tmp_plot), left = labsY, bottom = labsX)

### D: Plot Null-W test ####
Ts.Tu.plot <- ggplot()+
  geom_density(data = Ts_Tu_control, aes(x = Zscore), fill = "lightgrey")+
  geom_jitter(data = Ts_Tu_cands, aes(x = Zscore, y = 0.025), 
              size = 0.4, colour = "navy", height = 0.01, width = 0)+
  geom_vline(xintercept = quantile(Ts_Tu_control$Zscore, 0.95), colour = "red", lty = 2)+
  labs(title = "Tubesnout orthologs to threespine candidates")+
  lims(x = c(-30, 40), y = c(0, 0.17))+
  theme(axis.title = element_blank(),
        plot.title = element_text(size = txt_size-2),
        axis.text = element_text(size = pnt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line())

Tu.Ts.plot <- ggplot()+
  geom_density(data = Tu_Ts_control, aes(x = Zscore), fill = "lightgrey")+
  geom_jitter(data = Tu_Ts_cands, aes(x = Zscore, y = 0.025), 
              size = 0.4, colour = "navy", height = 0.01, width = 0)+
  geom_vline(xintercept = quantile(Tu_Ts_control$Zscore, 0.95), colour = "red", lty = 2)+
  labs(title = "Threespine orthologs to tubesnout candidates")+
  lims(x = c(-30, 40), y = c(0, 0.17))+
  theme(axis.title = element_blank(),
        plot.title = element_text(size = txt_size-2),
        axis.text = element_text(size = pnt_size),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line())

### Combine into multi-panel plot
tmp_plot <- plot_grid(Ts.Tu.plot, Tu.Ts.plot, labels = c("C", "D"), label_size = txt_size*1.5)
labsX <- textGrob("Z-scores", gp=gpar(fontface = "bold", fontsize = txt_size))
labsY <- textGrob("Density", gp=gpar(fontface="bold", fontsize = txt_size), rot = 90)

Fig4b <- grid.arrange(arrangeGrob(tmp_plot), left = labsY, bottom = labsX)


### E: Full plot ####
tiff("/Users/james/Dropbox (Personal)/Comp_geno_2020/Figure_4.tiff", width = 8, height = 6,
     units = "cm", res = 100)
plot_grid(Fig4a, Fig4b, nrow = 2)
dev.off()
