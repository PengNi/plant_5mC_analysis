library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)




rmet_diffnano_arab <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.c_count.bs_3reps_vs_nano_guppy_part2_50x.only_nanopore_detected_Cs.freq.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rmet_diffnano_arab) <- c("chrom", "pos", "strand", "coverage", "rmet", "motif")

rmet_diffnano_rice <- read.table("stats_nanopore_only_cytosines_profiling/shuidao.c_count.bs_2-1n1-2_vs_nano_guppy_shuidao1-1_50x.only_nanopore_detected_Cs.freq.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rmet_diffnano_rice) <- c("chrom", "pos", "strand", "coverage", "rmet", "motif")


rmet_arab <- data.frame(rmet=rmet_diffnano_arab$rmet, 
                        motif=rmet_diffnano_arab$motif, 
                        species="A. thaliana")
rmet_rice <- data.frame(rmet=rmet_diffnano_rice$rmet, 
                        motif=rmet_diffnano_rice$motif, 
                        species="O. sativa")
rmet_all <- rbind(rmet_arab, rmet_rice)

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


# ==================== cytosine
cbPalette <- c("#66c2a5", "#fc8d62")
# fontsize: 17, 15, 38
p_rmet_c <- ggplot(rmet_all, aes(x=rmet, fill=species)) + 
  geom_density(colour="black", alpha=0.75, size=0.5) + 
  geom_vline(xintercept=0.3, 
             colour="black", 
             linetype="dashed", size=0.8) +
  geom_vline(xintercept=0.7, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  theme(text=element_text(size=17, family = "serif"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15, face = "italic"), 
        legend.key.size = unit(0.35, "cm"), 
        plot.title = element_text(size=38, face="bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(limits=c(0, 1), 
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 8), 
                     breaks = seq(0, 8, 1)) + 
  labs(x="Methylation Frequency", 
       y="Density", title = "c")
# p_rmet_c

# ppi= 300
# jpeg("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c.jpg", 
#      width = 13, 
#      height = 10, units = "cm", res=ppi)
# p_rmet_c
# dev.off()

# combine with venn a, b (compare_detected_cytosines.plot.R)
ppi= 300
png("stats_nanopore_only_cytosines_profiling/fig_nanopore_only_cs.fig_abc.png", 
     width = 37, 
     height = 10, units = "cm", res=ppi)
grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c, 
             grid.rect(gp=gpar(col="white")), p_rmet_c,
             widths = c(12, 1, 12, 1, 13))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/fig_nanopore_only_cs.fig_abc.svg", 
    width = 37/2.54, 
    height = 10/2.54)
grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c, 
             grid.rect(gp=gpar(col="white")), p_rmet_c,
             widths = c(12, 1, 12, 1, 13))
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

cbPalette <- c("#66c2a5")
p_arab_cg <- ggplot(rmet_arab[rmet_arab$motif == "CG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "serif"), 
        plot.title = element_text(size=32, family = "serif", face="bold", 
                                  hjust = -0.2)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="a", 
       x="Methylation Frequency", 
       y="Count")
# p_arab_cg

p_arab_chg <- ggplot(rmet_arab[rmet_arab$motif == "CHG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "serif"), 
        plot.title = element_text(size=32, 
                                  family = "serif", face="bold",
                                  hjust = -0.25)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="b", 
       x="Methylation Frequency", 
       y="Count")
# p_arab_chg

p_arab_chh <- ggplot(rmet_arab[rmet_arab$motif == "CHH", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, 
                            family = "serif"), 
        plot.title = element_text(size=32,
                                  family = "serif", face="bold",
                                  hjust = -0.25)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="c", 
       x="Methylation Frequency", 
       y="Count")
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
cbPalette <- c("#fc8d62")
p_rice_cg <- ggplot(rmet_rice[rmet_rice$motif == "CG", ], 
                    aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "serif"), 
        plot.title = element_text(size=32, 
                                  family = "serif", face="bold",
                                  hjust = -0.2)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="d", 
       x="Methylation Frequency", 
       y="Count")
# p_rice_cg

p_rice_chg <- ggplot(rmet_rice[rmet_rice$motif == "CHG", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "serif"), 
        plot.title = element_text(size=32, 
                                  family = "serif", face="bold",
                                  hjust = -0.25)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="e", 
       x="Methylation Frequency", 
       y="Count")
# p_rice_chg

p_rice_chh <- ggplot(rmet_rice[rmet_rice$motif == "CHH", ], 
                     aes(x=rmet, fill=motif)) + 
  geom_histogram(binwidth=.01, position="dodge", 
                 colour="black", size=0.3) + 
  theme_bw() +
  theme(legend.position = "None", 
        legend.title = element_blank(), 
        text = element_text(size=15, family = "serif"), 
        plot.title = element_text(size=32, 
                                  family = "serif", face="bold",
                                  hjust = -0.25)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(label=scientific_format()) +
  labs(title="f", 
       x="Methylation Frequency", 
       y="Count")
# p_rice_chh

ppi= 300
#45, 22
png("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c_motif.png", 
     width = 36, 
     height = 17, units = "cm", res=ppi)
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
             heights=c(12, 1, 12))
dev.off()

svg("stats_nanopore_only_cytosines_profiling/generate_nanopore_only_cytosines_info.methyfreq.c_motif.svg", 
    width = 36/2.54, 
    height = 17/2.54)
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
             heights=c(12, 1, 12))
dev.off()













































