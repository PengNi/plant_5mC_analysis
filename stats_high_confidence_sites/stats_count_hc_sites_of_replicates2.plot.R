library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(gridExtra)
library(grid)


sitesnum_a_CG <- read.table("stats_high_confidence_sites/ninanjie-2.bs_3replicates.hc_sites_count.main_genome.CG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_a_CHG <- read.table("stats_high_confidence_sites/ninanjie-2.bs_3replicates.hc_sites_count.main_genome.CHG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_a_CHH <- read.table("stats_high_confidence_sites/ninanjie-2.bs_3replicates.hc_sites_count.main_genome.CHH.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_a <- rbind.data.frame(sitesnum_a_CG, sitesnum_a_CHG, sitesnum_a_CHH)

cbPalette <- c("#e41a1c", "#377eb8", "#4daf4a")
p_sitesnum_a <- ggplot() +
  geom_bar(data = sitesnum_a, aes(y=count, x=range, fill = replicate),
           colour="black", stat="identity", position="dodge") +
  facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)), 
                     limits = c(1, 10e7)) +
  scale_fill_manual(values=cbPalette,
                    breaks=c("rep1", "rep2", "rep3"),
                    labels=c("rep1     ", "rep2     ", "rep3")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=18, family = "serif"),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=25, hjust = 0.5, face = "italic", family="serif"),
        strip.background = element_rect(colour="white", fill="white",
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 18),
        panel.border = element_blank(),
        text=element_text(size=18,  family="serif"), 
        axis.text.x = element_text(size=18,  family="serif"), 
        axis.text.y = element_text(size=20,  family="serif")) +
  labs(x="Methylation frequency", y="Number of sites (log10)", title = "A. thaliana")
p_sitesnum_a


sitesnum_o_CG <- read.table("stats_high_confidence_sites/shuidao.bs_rep2-1_rep1-2.hc_sites_count.main_genome.CG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_o_CHG <- read.table("stats_high_confidence_sites/shuidao.bs_rep2-1_rep1-2.hc_sites_count.main_genome.CHG.txt", 
                             header = T, sep = "\t", stringsAsFactors = F)
sitesnum_o_CHH <- read.table("stats_high_confidence_sites/shuidao.bs_rep2-1_rep1-2.hc_sites_count.main_genome.CHH.txt", 
                             header = T, sep = "\t", stringsAsFactors = F)
sitesnum_o <- rbind.data.frame(sitesnum_o_CG, sitesnum_o_CHG, sitesnum_o_CHH)

cbPalette <- c("#e41a1c", "#377eb8", "#4daf4a")
p_sitesnum_o <- ggplot() +
  geom_bar(data = sitesnum_o, aes(y=count, x=range, fill = replicate),
           colour="black", stat="identity", position="dodge") +
  facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)), 
                     limits = c(1, 10e7)) +
  scale_fill_manual(values=cbPalette,
                    breaks=c("rep1", "rep2"),
                    labels=c("rep1     ", "rep2")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=18, family = "serif"),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=25, hjust = 0.5, face = "italic", family="serif"),
        strip.background = element_rect(colour="white", fill="white",
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 18),
        panel.border = element_blank(),
        text=element_text(size=18,  family="serif"), 
        axis.text.x = element_text(size=13,  family="serif"), 
        axis.text.y = element_text(size=20,  family="serif")) +
  labs(x="Methylation frequency", y="Number of sites (log10)", title = "O. sativa")
p_sitesnum_o


ppi= 300
png("stats_high_confidence_sites/stats_count_hc_sites_of_replicates2.fig.raw.png", 
    width = 31, 
    height = 16, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_sitesnum_a, 
                         grid.rect(gp=gpar(col="white")),
                         p_sitesnum_o,
                         nrow = 1,
                         ncol = 3,
                         widths=c(15, 1, 11)))
dev.off()









