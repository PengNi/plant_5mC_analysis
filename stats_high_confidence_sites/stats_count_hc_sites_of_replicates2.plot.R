library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(gridExtra)
library(grid)

#font_import()
#loadfonts(device = "win")
#fonts()


sitesnum_a_CG <- read.table("stats_high_confidence_sites/ninanjie-2.bs_3replicates.hc_sites_count.main_genome.CG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_a_CG$motif = "CpG"
sitesnum_a_CHG <- read.table("stats_high_confidence_sites/ninanjie-2.bs_3replicates.hc_sites_count.main_genome.CHG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_a_CHH <- read.table("stats_high_confidence_sites/ninanjie-2.bs_3replicates.hc_sites_count.main_genome.CHH.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_a <- rbind.data.frame(sitesnum_a_CG, sitesnum_a_CHG, sitesnum_a_CHH)
sitesnum_a$motif <- factor(sitesnum_a$motif, levels = c("CpG", "CHG", "CHH"))
sitesnum_a$range <- factor(sitesnum_a$range, levels = c("=0.0", "=1.0", ">=0.9"))

cbPalette <- c("#1b9e77", "#d95f02", "#7570b3")
p_sitesnum_a <- ggplot() +
  geom_bar(data = sitesnum_a, aes(y=count, x=range, fill = replicate),
           colour="black", size=0.3, stat="identity", position="dodge") +
  facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)), 
                     limits = c(1, 10e7)) +
  scale_fill_manual(values=cbPalette,
                    breaks=c("rep1", "rep2", "rep3"),
                    labels=c(" replicate1      ", " replicate2      ", " replicate3")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, family = "Arial"),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=15, hjust = 0.5, face = "italic", family="Arial"),
        strip.background = element_rect(colour="white", fill="white",
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(),
        text=element_text(size=12,  family="Arial"), 
        axis.text.x = element_text(size=9,  family="Arial"), 
        axis.text.y = element_text(size=12,  family="Arial")) +
  labs(x="Methylation frequency", y="Number of sites (log10)", title = "A. thaliana")
# p_sitesnum_a


sitesnum_o_CG <- read.table("stats_high_confidence_sites/shuidao.bs_rep2-1_rep1-2.hc_sites_count.main_genome.CG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
sitesnum_o_CG$motif = "CpG"
sitesnum_o_CHG <- read.table("stats_high_confidence_sites/shuidao.bs_rep2-1_rep1-2.hc_sites_count.main_genome.CHG.txt", 
                             header = T, sep = "\t", stringsAsFactors = F)
sitesnum_o_CHH <- read.table("stats_high_confidence_sites/shuidao.bs_rep2-1_rep1-2.hc_sites_count.main_genome.CHH.txt", 
                             header = T, sep = "\t", stringsAsFactors = F)
sitesnum_o <- rbind.data.frame(sitesnum_o_CG, sitesnum_o_CHG, sitesnum_o_CHH)
sitesnum_o$motif <- factor(sitesnum_o$motif, levels = c("CpG", "CHG", "CHH"))
sitesnum_o$range <- factor(sitesnum_o$range, levels = c("=0.0", "=1.0", ">=0.9"))

cbPalette <- c("#1b9e77", "#d95f02", "#7570b3")
p_sitesnum_o <- ggplot() +
  geom_bar(data = sitesnum_o, aes(y=count, x=range, fill = replicate),
           colour="black", size=0.3, stat="identity", position="dodge") +
  facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)), 
                     limits = c(1, 10e7)) +
  scale_fill_manual(values=cbPalette,
                    breaks=c("rep1", "rep2"),
                    labels=c(" sample1      ", " sample2")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, family = "Arial"),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=15, hjust = 0.5, face = "italic", family="Arial"),
        strip.background = element_rect(colour="white", fill="white",
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(),
        text=element_text(size=12,  family="Arial"), 
        axis.text.x = element_text(size=6,  family="Arial"), 
        axis.text.y = element_text(size=12,  family="Arial")) +
  labs(x="Methylation frequency", y="Number of sites (log10)", title = "O. sativa")
# p_sitesnum_o


ppi= 300
# 31:16
png("stats_high_confidence_sites/stats_count_hc_sites_of_replicates2.fig.raw.png", 
    width = 18, 
    height = 9.3, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_sitesnum_a, 
                         grid.rect(gp=gpar(col="white")),
                         p_sitesnum_o,
                         nrow = 1,
                         ncol = 3,
                         widths=c(29, 1, 20)))
dev.off()
svg("stats_high_confidence_sites/stats_count_hc_sites_of_replicates2.fig.raw.svg",
    width = 18.0/2.54,
    height = 9.3/2.54)
grid.arrange(arrangeGrob(p_sitesnum_a,
                         grid.rect(gp=gpar(col="white")),
                         p_sitesnum_o,
                         nrow = 1,
                         ncol = 3,
                         widths=c(29, 1, 20)))
dev.off()
# pdf("stats_high_confidence_sites/stats_count_hc_sites_of_replicates2.fig.raw.pdf",
#     width = 18.0/2.54, 
#     height = 9.3/2.54)
# grid.arrange(arrangeGrob(p_sitesnum_a, 
#                          grid.rect(gp=gpar(col="white")),
#                          p_sitesnum_o,
#                          nrow = 1,
#                          ncol = 3,
#                          widths=c(29, 1, 20)))
# dev.off()







