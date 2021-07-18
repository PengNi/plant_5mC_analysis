library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

library(grid)
library(gridExtra)



rmet_diffnano_arab <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.c_count.bs_3reps_vs_nano50x_dp2_cnn_p0.8.only_nanopore_detected_Cs.freq.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rmet_diffnano_arab) <- c("chrom", "pos", "strand", "coverage", "rmet", "motif")

rmet_diffnano_rice <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.c_count.bs_rep2-1_vs_nano50x_dp2_cnn_p0.8.only_nanopore_detected_Cs.freq.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rmet_diffnano_rice) <- c("chrom", "pos", "strand", "coverage", "rmet", "motif")

rmet_diffnano_rice2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.c_count.bs_rep1-2_vs_nano50x_dp2_cnn_p0.8.only_nanopore_detected_Cs.freq.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rmet_diffnano_rice2) <- c("chrom", "pos", "strand", "coverage", "rmet", "motif")


rmet_arab <- data.frame(rmet=rmet_diffnano_arab$rmet, 
                        motif=rmet_diffnano_arab$motif, 
                        species="A. thaliana")
rmet_rice <- data.frame(rmet=rmet_diffnano_rice$rmet, 
                        motif=rmet_diffnano_rice$motif, 
                        species="O. sativa")
rmet_all <- rbind(rmet_arab, rmet_rice)

rmet_rice2 <- data.frame(rmet=rmet_diffnano_rice2$rmet, 
                        motif=rmet_diffnano_rice2$motif, 
                        species="O. sativa")

himethy_arab = sum(rmet_diffnano_arab$rmet>=0.7)
himethy_arab
himethy_arab/nrow(rmet_diffnano_arab)
lomethy_arab = sum(rmet_diffnano_arab$rmet<=0.3)
lomethy_arab
lomethy_arab/nrow(rmet_diffnano_arab)

himethy_rice = sum(rmet_diffnano_rice$rmet>=0.7)
himethy_rice
himethy_rice/nrow(rmet_diffnano_rice)
lomethy_rice = sum(rmet_diffnano_rice$rmet<=0.3)
lomethy_rice
lomethy_rice/nrow(rmet_diffnano_rice)

himethy_rice2 = sum(rmet_diffnano_rice2$rmet>=0.7)
himethy_rice2
himethy_rice2/nrow(rmet_diffnano_rice2)
lomethy_rice2 = sum(rmet_diffnano_rice2$rmet<=0.3)
lomethy_rice2
lomethy_rice2/nrow(rmet_diffnano_rice2)


# ==================== cytosine
cbPalette <- c("#1f78b4", "#33a02c")

p_rmet_c_a <-  ggplot(rmet_arab, aes(x=rmet)) + 
  geom_density(colour="black", alpha=0.8, size=0.4, 
               fill=cbPalette[1]) + 
  geom_vline(xintercept=0.3, 
             colour="black", 
             linetype="dashed", size=0.8) +
  geom_vline(xintercept=0.7, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        plot.title = element_text(size=18, hjust = 0.5), 
        axis.title = element_text(size=15, family = "Arial"),
        axis.text = element_text(size=15, family = "Arial")) + 
  scale_x_continuous(limits=c(0, 1), 
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 8), 
                     breaks = seq(0, 8, 1)) + 
  labs(x="Methylation frequency", 
       y="Density", title = "A.thaliana")

p_rmet_c_rice <-  ggplot(rmet_rice, aes(x=rmet)) + 
  geom_density(colour="black", alpha=0.8, size=0.4, 
               fill=cbPalette[2]) + 
  geom_vline(xintercept=0.3, 
             colour="black", 
             linetype="dashed", size=0.8) +
  geom_vline(xintercept=0.7, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        plot.title = element_text(size=18, hjust = 0.5), 
        axis.title = element_text(size=15, family = "Arial"),
        axis.text = element_text(size=15, family = "Arial")) + 
  scale_x_continuous(limits=c(0, 1), 
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 8), 
                     breaks = seq(0, 8, 1)) + 
  labs(x="Methylation frequency", 
       y="Density", title = expression(paste(italic("O. sativa"),  " (sample1)", sep = "")))

# # fontsize: 17, 15, 38
# p_rmet_c <- ggplot(rmet_all, aes(x=rmet, fill=species)) + 
#   geom_density(colour="black", alpha=0.75, size=0.5) + 
#   geom_vline(xintercept=0.3, 
#              colour="black", 
#              linetype="dashed", size=0.8) +
#   geom_vline(xintercept=0.7, 
#              colour="black", 
#              linetype="dashed", size=0.8) +
#   theme_bw() + 
#   theme(text=element_text(size=17, family = "Arial"),
#         legend.position = "bottom", 
#         legend.title = element_blank(), 
#         legend.text = element_text(size=15, face = "italic"), 
#         legend.key.size = unit(0.35, "cm"), 
#         plot.title = element_text(size=38, face="bold", hjust = -0.1)) + 
#   scale_fill_manual(values=cbPalette, 
#                     breaks = c("A. thaliana", "O. sativa"), 
#                     labels = c("A. thaliana    ", "O. sativa")) + 
#   scale_x_continuous(limits=c(0, 1), 
#                      breaks = seq(0, 1, 0.1)) +
#   scale_y_continuous(limits = c(0, 8), 
#                      breaks = seq(0, 8, 1)) + 
#   labs(x="Methylation frequency", 
#        y="Density", title = "")
# # p_rmet_c

# # combine with venn a, b (compare_detected_cytosines.plot.R)
# ppi= 300
# png("stats_nanopore_only_cytosines_profiling/fig_nanopore_only_cs.fig_abc.notitle.png", 
#      width = 37, 
#      height = 10, units = "cm", res=ppi)
# grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c, 
#              grid.rect(gp=gpar(col="white")), p_rmet_c,
#              widths = c(12, 1, 12, 1, 13))
# dev.off()
# svg("stats_nanopore_only_cytosines_profiling/fig_nanopore_only_cs.fig_abc.notitle.svg", 
#     width = 37/2.54, 
#     height = 10/2.54)
# grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c, 
#              grid.rect(gp=gpar(col="white")), p_rmet_c,
#              widths = c(12, 1, 12, 1, 13))
# dev.off()

# =============== rice rep2
# cbPalette <- c("#fc8d62")
# fontsize: 17, 15, 38
p_rmet_c_rice2 <-  ggplot(rmet_rice2, aes(x=rmet)) + 
  geom_density(colour="black", alpha=0.8, size=0.4, 
               fill=cbPalette[2]) + 
  geom_vline(xintercept=0.3, 
             colour="black", 
             linetype="dashed", size=0.8) +
  geom_vline(xintercept=0.7, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        plot.title = element_text(size=18, hjust = 0.5), 
        axis.title = element_text(size=15, family = "Arial"),
        axis.text = element_text(size=15, family = "Arial")) + 
  scale_x_continuous(limits=c(0, 1), 
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 8), 
                     breaks = seq(0, 8, 1)) + 
  labs(x="Methylation frequency", 
       y="Density", title = expression(paste(italic("O. sativa"),  " (sample2)", sep = "")))

p_rmet_c_rice2

# # combine with venn rice rep2 and rmet rice rep2 (compare_detected_cytosines.plot.R)
# ppi= 300
# png("stats_nanopore_only_cytosines_profiling/fig_nanopore_only_cs.rice1-1.ab.notitle.png", 
#     width = 20, 
#     height = 8, units = "cm", res=ppi)
# grid.arrange(p_rice_c2, 
#              grid.rect(gp=gpar(col="white")), p_rmet_c_rice2,
#              widths = c(12, 1, 13))
# dev.off()
# svg("stats_nanopore_only_cytosines_profiling/fig_nanopore_only_cs.rice1-1.ab.notitle.svg", 
#     width = 20/2.54, 
#     height = 8/2.54)
# grid.arrange(p_rice_c2, 
#              grid.rect(gp=gpar(col="white")), p_rmet_c_rice2,
#              widths = c(12, 1, 13))
# dev.off()


# methy freq distribution arab + rice1, rice2
ppi= 300
png("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c.raw.png", 
    width = 36, 
    height = 8, units = "cm", res=ppi)
grid.arrange(p_rmet_c_a, grid.rect(gp=gpar(col="white")), 
             p_rmet_c_rice, 
             grid.rect(gp=gpar(col="white")), 
             p_rmet_c_rice2,
             widths = c(11.3, 1, 11.3, 1, 11.3))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c.raw.svg", 
    width = 36/2.54, 
    height = 8/2.54)
grid.arrange(p_rmet_c_a, grid.rect(gp=gpar(col="white")), 
             p_rmet_c_rice, 
             grid.rect(gp=gpar(col="white")), 
             p_rmet_c_rice2,
             widths = c(11.3, 1, 11.3, 1, 11.3))
dev.off()






# ==================== cytosine motif
print_highlow_methy <- function(rmetdf){
  motifs = unique(rmetdf$motif)
  for(motif in motifs){
    rmet_tmp <- rmetdf[rmetdf$motif == motif, ]
    himethy_cnt = sum(rmet_tmp$rmet>=0.7)
    message(sprintf("motif %s, %d, highly methylated, %d, %f", 
                    motif, nrow(rmet_tmp), himethy_cnt, himethy_cnt/nrow(rmet_tmp)))
    lomethy_cnt = sum(rmet_tmp$rmet<=0.3)
    message(sprintf("motif %s, %d, lowly methylated, %d, %f", 
                    motif, nrow(rmet_tmp), lomethy_cnt, lomethy_cnt/nrow(rmet_tmp)))
  }
}

# arab =
print_highlow_methy(rmet_arab)

cbPalette <- c("#1f78b4")
p_arab_cg <- ggplot(rmet_arab[rmet_arab$motif == "CG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, family = "Arial", face="bold", 
                                  hjust = -0.3)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="a", 
       x="Methylation frequency", 
       y="Count") + 
  annotation_custom(grobTree(textGrob("CpG", 
                                      x=0.4,  y=0.55, hjust=0, 
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_arab_cg

p_arab_chg <- ggplot(rmet_arab[rmet_arab$motif == "CHG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.26)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="b", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CHG", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_arab_chg

p_arab_chh <- ggplot(rmet_arab[rmet_arab$motif == "CHH", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, 
                            family = "Arial"), 
        plot.title = element_text(size=28,
                                  family = "Arial", face="bold",
                                  hjust = -0.26)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="c", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CHH", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_arab_chh


# rice =
print_highlow_methy(rmet_rice)

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
cbPalette <- c("#33a02c")
p_rice_cg <- ggplot(rmet_rice[rmet_rice$motif == "CG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.22)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="d", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CpG", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_rice_cg

p_rice_chg <- ggplot(rmet_rice[rmet_rice$motif == "CHG", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.24)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="e", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CHG", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_rice_chg

p_rice_chh <- ggplot(rmet_rice[rmet_rice$motif == "CHH", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.21)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="f", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CHH", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_rice_chh


# rice rep2 (shuidao1-1)
cbPalette <- c("#33a02c")
p_rice_cg2 <- ggplot(rmet_rice2[rmet_rice2$motif == "CG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.21)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="g", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CpG", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_rice_cg2

p_rice_chg2 <- ggplot(rmet_rice2[rmet_rice2$motif == "CHG", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.23)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="h", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CHG", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_rice_chg2

p_rice_chh2 <- ggplot(rmet_rice2[rmet_rice2$motif == "CHH", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=28, 
                                  family = "Arial", face="bold",
                                  hjust = -0.2)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="i", 
       x="Methylation frequency", 
       y="Count")+ 
  annotation_custom(grobTree(textGrob("CHH", 
                                      x=0.4,  y=0.55, hjust=0,
                                      gp=gpar(col="black", fontsize=24, fontfamily="Arial"))))
# p_rice_chh2


ppi= 300
#45, 22
png("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c_motif.png", 
     width = 36, 
     height = 26, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chg, 
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chh, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg2,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg2, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh2, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             heights=c(12, 1, 12, 1, 12))
dev.off()

svg("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c_motif.svg", 
    width = 36/2.54, 
    height = 26/2.54)
grid.arrange(arrangeGrob(p_arab_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chg, 
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chh, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg2,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg2, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh2, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             heights=c(12, 1, 12, 1, 12))
dev.off()













































