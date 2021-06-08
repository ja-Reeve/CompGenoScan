### Preperation
rm(list = ls())
options(stringsAsFactors = FALSE)

### Packages
library(ggplot2)
library(gridExtra)

### Nucleotide diversity
Ts.He <- read.table("file:///G:/nucleotide_diversity/Ts_Avg_Hexp_50K.txt", header = TRUE)
Tu.He <- read.table("file:///G:/nucleotide_diversity/Tu_Avg_Hexp_50K.txt", header = TRUE)
Ns.He <- read.table("file:///G:/nucleotide_diversity/Ns_Avg_Hexp_50K.txt", header = TRUE)

Ts.He <- Ts.He[complete.cases(Ts.He),]
Tu.He <- Tu.He[complete.cases(Tu.He),]
Ns.He <- Ns.He[complete.cases(Ns.He),]

### Load comparison data
Ts_Tu <- read.table("file:///G:/Ts_Tu_50Kdata.txt", header = TRUE)
Ts_Ns <- read.table("file:///G:/Ts_Ns_50Kdata.txt", header = TRUE)

all50K <- merge(Ts_Tu, Ts_Ns, by = "Gene.Ts")[,-c(24:33, 45)]
colnames(all50K) <- c("Gene.Ts", "CHROM.Ts", "Start_Pos.Ts", "End_Pos.Ts", "TsAK_He", "TsOR_He", "Avg_He.Ts", "Avg_Fst.Ts", "CAND.Ts", "Window.Ts",   
                      "CHROM.Tu", "Start_Pos.Tu", "End_Pos.Tu", "TuAK_He", "TuBC_He", "Avg_He.Tu", "Avg_Fst.Tu", "CAND.Tu", "Window.Tu", "Gene.Tu", 
                      "Outliers.Ts", "Outliers.Tu", "Outliers.Ts_Tu",
                      "CHROM.Ns", "Start_Pos.Ns", "End_Pos.Ns", "NsNUn_He", "NsNUd_He", "NsABk_He", "NsABm_He", "Avg_He.Ns", "NS1_Fst.Ns", "NS2_Fst.Ns", "Window.Ns",
                      "Outliers.Ns", "Outliers.Ts_Ns")
### Find shared outlier between ninespine and tubesnout
outlier.overlaps <- function(x, y, label.x, label.y){
  
  logi <- rep("non-outlier", length(x))
  
  logi[which(x != "non-outlier")] <- label.x
  logi[which(y != "non-outlier")] <- label.y
  logi[which(x != "non-outlier" & y != "non-outlier")] <- "Both"
  
  return(logi)
}

all50K$Outliers.Ns_Tu <- outlier.overlaps(all50K$Outliers.Ns, all50K$Outliers.Tu, "Ninespine", "Tubesnout")

### Correlations among populations
cor.test(x = all50K$TsAK_He, y = all50K$TsOR_He, method = "spearman")

cor.test(x = all50K$TsAK_He, y = all50K$NsNUn_He, method = "spearman")
cor.test(x = all50K$TsAK_He, y = all50K$NsNUd_He, method = "spearman")
cor.test(x = all50K$TsAK_He, y = all50K$NsABk_He, method = "spearman")
cor.test(x = all50K$TsAK_He, y = all50K$NsABm_He, method = "spearman")

cor.test(x = all50K$TsAK_He, y = all50K$TuAK_He, method = "spearman")
cor.test(x = all50K$TsAK_He, y = all50K$TuBC_He, method = "spearman")

cor.test(x = all50K$TsOR_He, y = all50K$NsNUn_He, method = "spearman")
cor.test(x = all50K$TsOR_He, y = all50K$NsNUd_He, method = "spearman")
cor.test(x = all50K$TsOR_He, y = all50K$NsABk_He, method = "spearman")
cor.test(x = all50K$TsOR_He, y = all50K$NsABm_He, method = "spearman")

cor.test(x = all50K$TsOR_He, y = all50K$TuAK_He, method = "spearman")
cor.test(x = all50K$TsOR_He, y = all50K$TuBC_He, method = "spearman")

cor.test(x = all50K$NsNUn_He, y = all50K$NsNUd_He, method = "spearman")
cor.test(x = all50K$NsNUn_He, y = all50K$NsABk_He, method = "spearman")
cor.test(x = all50K$NsNUn_He, y = all50K$NsABm_He, method = "spearman")

cor.test(x = all50K$NsNUn_He, y = all50K$TuAK_He, method = "spearman")
cor.test(x = all50K$NsNUn_He, y = all50K$TuBC_He, method = "spearman")

cor.test(x = all50K$NsNUd_He, y = all50K$NsABk_He, method = "spearman")
cor.test(x = all50K$NsNUd_He, y = all50K$NsABm_He, method = "spearman")

cor.test(x = all50K$NsNUd_He, y = all50K$TuAK_He, method = "spearman")
cor.test(x = all50K$NsNUd_He, y = all50K$TuBC_He, method = "spearman")

cor.test(x = all50K$NsABk_He, y = all50K$NsABm_He, method = "spearman")

cor.test(x = all50K$NsABk_He, y = all50K$TuAK_He, method = "spearman")
cor.test(x = all50K$NsABk_He, y = all50K$TuBC_He, method = "spearman")

cor.test(x = all50K$NsABm_He, y = all50K$TuAK_He, method = "spearman")
cor.test(x = all50K$NsABm_He, y = all50K$TuBC_He, method = "spearman")

cor.test(x = all50K$TuAK_He, y = all50K$TuBC_He, method = "spearman")

### Among species correlations
# Genetic diversity
cor.test(x = all50K$Avg_He.Ts, all50K$Avg_He.Ns, method = "spearman")
cor.test(x = all50K$Avg_He.Ts, all50K$Avg_He.Tu, method = "spearman")
cor.test(x = all50K$Avg_He.Ns, all50K$Avg_He.Tu, method = "spearman")

# Fst
cor.test(x = all50K$Avg_Fst.Ts, y = all50K$NS1_Fst.Ns, method = "spearman")
cor.test(x = all50K$Avg_Fst.Ts, y = all50K$NS2_Fst.Ns, method = "spearman")
cor.test(x = all50K$Avg_Fst.Ts, y = all50K$Avg_Fst.Tu, method = "spearman")

cor.test(x = all50K$NS1_Fst.Ns, y = all50K$NS2_Fst.Ns, method = "spearman")

cor.test(x = all50K$NS1_Fst.Ns, y = all50K$Avg_Fst.Tu, method = "spearman")
cor.test(x = all50K$NS2_Fst.Ns, y = all50K$Avg_Fst.Tu, method = "spearman")

### Correlation plots
Ts_Ns_Plot <- ggplot(all50K)+
  geom_point(aes(x = Avg_He.Ts, y = Avg_He.Ns, colour = Outliers.Ts_Ns), size = 0.75)+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.99), colour = "red", lty = 2)+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.01), colour = "red", lty = 2)+
  geom_hline(yintercept = quantile(Ns.He$Avg_He, 0.99), colour = "red", lty = 2)+
  geom_hline(yintercept = quantile(Ns.He$Avg_He, 0.01), colour = "red", lty = 2)+
  labs(x = expression(Threespine~bar(H)[E]), y = expression(Ninespine~bar(H)[E]))+
  lims(y = c(0, 0.45), x = c(0, 0.45))+
  scale_colour_manual(values = c("orange", "navy", "grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

Ts_Tu_Plot <- ggplot(all50K)+
  geom_point(aes(x = Avg_He.Ts, y = Avg_He.Tu, colour = Outliers.Ts_Tu), size = 0.75)+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.99), colour = "red", lty = 2)+
  geom_vline(xintercept = quantile(Ts.He$Avg_He, 0.01), colour = "red", lty = 2)+
  geom_hline(yintercept = quantile(Tu.He$Avg_He, 0.99), colour = "red", lty = 2)+
  geom_hline(yintercept = quantile(Tu.He$Avg_He, 0.01), colour = "red", lty = 2)+
  labs(x = expression(Threespine~bar(H)[E]), y = expression(Tubesnout~bar(H)[E]))+
  lims(y = c(0, 0.45), x = c(0, 0.45))+
  scale_colour_manual(values = c("orange", "grey", "navy", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

Ns_Tu_Plot <- ggplot(all50K)+
  geom_point(aes(x = Avg_He.Ns, y = Avg_He.Tu, colour = Outliers.Ns_Tu), size = 0.75)+
  geom_vline(xintercept = quantile(Ns.He$Avg_He, 0.99), colour = "red", lty = 2)+
  geom_vline(xintercept = quantile(Ns.He$Avg_He, 0.01), colour = "red", lty = 2)+
  geom_hline(yintercept = quantile(Tu.He$Avg_He, 0.99), colour = "red", lty = 2)+
  geom_hline(yintercept = quantile(Tu.He$Avg_He, 0.01), colour = "red", lty = 2)+
  labs(x = expression(Ninespine~bar(H)[E]), y = expression(Tubesnout~bar(H)[E]))+
  lims(y = c(0, 0.45), x = c(0, 0.45))+
  scale_colour_manual(values = c("orange", "navy", "grey", "navy"))+
  theme(axis.title = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        legend.position = "NULL")

Plots.50K <- list(Ts_Ns_Plot, Ts_Tu_Plot, Ns_Tu_Plot)
grid.arrange(grobs = Plots.50K, ncol = 1)
