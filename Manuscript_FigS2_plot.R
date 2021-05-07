rm(list = ls())
options(stringsAsFactors = FALSE)
dev.off()

Ts.cov <- read.table("file:///F:/coverage/Ts_cov_50Kwind.txt", header = TRUE)
Tu.cov <- read.table("file:///F:/coverage/Tu_cov_50Kwind.txt", header = TRUE)

Ts.He <- read.table("file:///F:/nucleotide_diversity/Ts_Avg_Hexp_50Kwind.txt", header = TRUE)
Tu.He <- read.table("file:///F:/nucleotide_diversity/Tu_Avg_Hexp_50Kwind.txt", header = TRUE)

#Ns.He <- read.table("file:///F:/nucleotide_diversity/Ns_Avg_Hexp_50Kwind.txt", header = TRUE)

Ts.cov$WIND <- paste0("Ts_", Ts.cov$CHROM, ":", Ts.cov$START_POS, "-", Ts.cov$END_POS)
Tu.cov$WIND <- paste0("Tu_", Tu.cov$CHROM, ":", Tu.cov$START_POS, "-", Tu.cov$END_POS)

Ts.cov_He <- merge(Ts.cov, Ts.He, by = "WIND")
Tu.cov_He <- merge(Tu.cov, Tu.He, by = "WIND")

p1.1 <- ggplot(Ts.cov)+
  geom_histogram(aes(x = COV_TsAK), binwidth = 1, fill = "darkgreen", alpha = 0.4)+
  geom_histogram(aes(x = COV_TsOR), binwidth = 1, fill = "green", alpha = 0.4)+
  xlim(c(0, 200))+
  labs(y = "Count", title = "Threespine stickleback (Ts)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8))

p1.2 <- ggplot(Ts.cov_He)+
  geom_point(aes(x = COV_TsAK, y = TsAK_He), alpha = 0.6, col = "darkgreen")+
  labs(y = expression(bar(H)[E]), title = "TsAK")+
  lims(x = c(0, 200), y = c(0,0.012))+
  theme(axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line())
p1.3 <- ggplot(Ts.cov_He)+
  geom_point(aes(x = COV_TsOR, y = TsOR_He), alpha = 0.6, col = "green")+
  labs(x = "Coverage per window", y = expression(bar(H)[E]), title = "TsOR")+
  lims(x = c(0, 200), y = c(0,0.012))+
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 8),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line())

labX <- textGrob("Coverage per window", gp=gpar(fontface="bold", fontsize=12))
p1 <- plot_grid(p1.1, plot_grid(p1.2, p1.3, nrow = 1, rel_widths = c(1, 0.9), rel_heights = c(1,1)), labX, 
          ncol = 1, rel_heights = c(1.5,1,0.2))


p2.1 <- ggplot(Tu.cov)+
  geom_histogram(aes(x = COV_TuAK), binwidth = 1, fill = "firebrick", alpha = 0.4)+
  geom_histogram(aes(x = COV_TuBC), binwidth = 1, fill = "orange", alpha = 0.4)+
  xlim(c(0, 200))+
  labs(y = "Count", title = "Tubesnout (Tu)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8))

p2.2 <- ggplot(Tu.cov_He)+
  geom_point(aes(x = COV_TuAK, y = TuAK_He), alpha = 0.6, col = "firebrick")+
  labs(x = "Coverage per window", y = expression(bar(H)[E]), title = "TuAK")+
  lims(x = c(0, 200), y = c(0,0.012))+
  theme(axis.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line())

p2.3 <- ggplot(Tu.cov_He)+
  geom_point(aes(x = COV_TuBC, y = TuBC_He), alpha = 0.6, col = "firebrick")+
  labs(x = "Coverage per window", y = expression(bar(H)[E]), title = "TuBC")+
  lims(x = c(0, 200), y = c(0,0.012))+
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 8),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line())

labX <- textGrob("Coverage per window", gp=gpar(fontface="bold", fontsize=12))
p2 <- plot_grid(p2.1, plot_grid(p2.2, p2.3, nrow = 1, rel_widths = c(1, 0.9), rel_heights = c(1,1)), labX, 
                ncol = 1, rel_heights = c(1.5,1,0.2))  

Fig.S2 <- plot_grid(p1, p2, nrow = 1)
