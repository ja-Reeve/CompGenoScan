####### Manuscript figure 4 ################
### This script creates the plot used in figure 4 of my manuscript.
### Fig. 4 show the relationship between interspecies correlations and window size
### James Reeve
### Göteborgs Universitet
### 20/03/2020

### Prepreation
rm(list = ls())
dev.off()

### Packages
library("ggplot2"); library("cowplot")

### Access data
# Genetic diversity
Corr.dat.He <- read.csv("file:///C:/Users/James/Dropbox (Personal)/Master's/correlation_table.csv", header = TRUE)
Corr.dat.He$window.size <- factor(Corr.dat.He$window.size, levels = c("GENE", "10K", "25K", "50K", "75K", "100K"))

# Fst
Corr.dat.Fst <- read.csv("file:///C:/Users/James/Dropbox (Personal)/Master's/correlation_table_Fst.csv", header = TRUE)
Corr.dat.Fst$window.size <- factor(Corr.dat.Fst$window.size, levels = c("GENE", "10K", "25K", "50K", "75K", "100K"))

### Plot correlations within species
# genetic diversity
tmp <- Corr.dat.He[Corr.dat.He$contrast == "TsAK_TuAK" | Corr.dat.He$contrast == "TuAK_TuBC" | Corr.dat.He$contrast == "NsNUn_NsNUd" |
                     Corr.dat.He$contrast == "NsNUn_NsABm" | Corr.dat.He$contrast == "NsNUn_NsABk" | Corr.dat.He$contrast == "NsNUd_NsABk" |
                     Corr.dat.He$contrast == "NsNUd_NsABm" | Corr.dat.He$contrast == "NsABk_NsABm", ]

tmp$species <- substring(tmp$contrast, 1, 2)

P1 <- ggplot(tmp, aes(x = window.size, y = spearman, colour = species, group = contrast))+
  geom_hline(yintercept = 0.0, lty = 2, colour = "lightgrey")+
  geom_line(aes(lty = contrast))+
  geom_point(aes(pch = contrast), show.legend = F)+
  scale_linetype_manual(values = c(2:6, 2, 1, 1))+
  scale_shape_manual(values = c(2:6, 2, 16, 16))+
  labs(y = expression("Within-species Spearman's"~rho), lty = "Population contrast")+
  ylim(c(-0.2, 1.0))+
  theme_classic()+ theme(axis.title.x = element_blank(), 
                         axis.title.y = element_text(size = 14))

### Plot correlations among species
# Fst
P2 <- ggplot(Corr.dat.Fst, aes(x = window.size, y = spearman, group = contrast, colour = contrast))+
  geom_hline(yintercept = 0.0, lty = 2, colour = "lightgrey")+
  geom_line(aes(lty = contrast))+
  geom_point(aes(pch = contrast))+
  scale_colour_manual(values = c("black", "firebrick", "firebrick", "navy", "navy", "orange"))+
  labs(y = expression(~"Spearman's"~rho), colour = "Population contrast",
       title = expression(F[ST]))+
  ylim(c(-0.05, 0.1))+
  theme_classic()+ theme(axis.title.x = element_blank(), 
                         axis.title.y = element_text(size = 14))

#genetic diversity
tmp <- Corr.dat.He <- Corr.dat.He[Corr.dat.He$contrast == "AvgTs_AvgNs" | Corr.dat.He$contrast == "AvgTs_AvgTu" | 
                                    Corr.dat.He$contrast == "AvgNs_AvgTu", ]

P3 <- ggplot(tmp, aes(x = window.size, y = spearman, group = contrast, colour = contrast))+
  geom_hline(yintercept = 0.0, lty = 2, colour = "lightgrey")+
  geom_line(aes(lty = contrast))+
  geom_point(aes(pch = contrast))+
  scale_colour_manual(values = c("firebrick", "navy", "orange"))+
  labs(x = "Window size", y = expression("Spearman's"~rho), colour = "Population contrast",
       title = expression(bar(H)[E]))+
  ylim(c(-0.1, 0.2))+
  theme_classic()+ theme(axis.title = element_text(size = 14))

### Full figure

plot_grid(P2, P3, labels = c("A", "B"), label_size = 18, nrow = 2, rel_heights = c(1, 2))
