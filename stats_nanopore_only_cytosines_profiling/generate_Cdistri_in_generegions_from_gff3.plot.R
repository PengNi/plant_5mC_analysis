library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

# arab gene
arab_gene <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.c_count.bs_3reps_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_gene.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)
arab_gene$motif_ratio <- arab_gene$moitf_num / arab_gene$motif_total
arab_gene$region <- factor(arab_gene$region, 
                           levels = c("protein-coding genes", "non-coding genes", 
                                         "pseudogenes", "transposons"))

p_arab_gene <- ggplot(data=arab_gene, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=28, face="bold", hjust = -0.15, vjust = -0.01), 
        axis.text.x  = element_text(size=7.5)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c(" CpG    ", " CHG    ", " CHH")) +
  scale_x_discrete(breaks=c("protein-coding genes", "non-coding genes", 
                            "pseudogenes", "transposons"), 
                   labels=c("protein-coding genes", "non-coding genes", 
                            "pseudogenes", "transposons")) +
  scale_y_continuous(limits=c(0, 0.2), breaks = seq(0, 0.2, 0.05), 
                     labels = seq(0, 0.2, 0.05) * 100) +
  labs(x="", y="Ratio (%)", title="b")
# p



# arab transcript
arab_rna <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.c_count.bs_3reps_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_transcript.txt", 
                        header = T, sep = "\t", stringsAsFactors = F, quote = "")
arab_rna$motif_ratio <- arab_rna$moitf_num / arab_rna$motif_total
arab_rna$region <- factor(arab_rna$region, 
                          levels = c("5'UTR", "CDS", 
                                      "3'UTR"))
p_arab_rna <- ggplot(data=arab_rna, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=28, face="bold", hjust = -0.15, vjust = -0.01)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c(" CpG    ", " CHG    ", " CHH")) +
  scale_y_continuous(limits=c(0, 0.12), breaks = seq(0, 0.12, 0.03), 
                     labels = seq(0, 0.12, 0.03) * 100) +
  labs(x="", y="Ratio (%)", title="c")




# rice gene
rice_gene <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.c_count.bs_rep2-1_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_gene.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)
rice_gene$motif_ratio <- rice_gene$moitf_num / rice_gene$motif_total
rice_gene$region <- factor(rice_gene$region, 
                           levels = c("protein-coding genes", "non-coding genes"))

p_rice_gene <- ggplot(data=rice_gene, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=28, face="bold", hjust = -0.10, vjust = -0.01)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c(" CpG    ", " CHG    ", " CHH")) +
  scale_y_continuous(limits=c(0, 0.08), breaks = seq(0, 0.08, 0.02), 
                     labels = seq(0, 0.08, 0.02) * 100) +
  labs(x="", y="Ratio (%)", title = "e")




# rice transcript
rice_rna <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.c_count.bs_rep2-1_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_transcript.txt", 
                       header = T, sep = "\t", stringsAsFactors = F, quote = "")
rice_rna$motif_ratio <- rice_rna$moitf_num / rice_rna$motif_total
rice_rna$region <- factor(rice_rna$region, 
                          levels = c("5'UTR", "CDS", 
                                     "3'UTR"))
p_rice_rna <- ggplot(data=rice_rna, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=28, face="bold", hjust = -0.14, vjust = -0.01)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c(" CpG    ", " CHG    ", " CHH")) +
  scale_y_continuous(limits=c(0, 0.02), breaks = seq(0, 0.02, 0.005), 
                     labels = seq(0, 0.02, 0.005) * 100) +
  labs(x="", y="Ratio (%)", title="f")




# rice gene rep2
rice_gene2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.c_count.bs_rep1-2_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_gene.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)
rice_gene2$motif_ratio <- rice_gene2$moitf_num / rice_gene2$motif_total
rice_gene2$region <- factor(rice_gene2$region, 
                           levels = c("protein-coding genes", "non-coding genes"))

p_rice_gene2 <- ggplot(data=rice_gene2, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=28, face="bold", hjust = -0.10, vjust = -0.01)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c(" CpG    ", " CHG    ", " CHH")) +
  scale_y_continuous(limits=c(0, 0.08), breaks = seq(0, 0.08, 0.02), 
                     labels = seq(0, 0.08, 0.02) * 100) +
  labs(x="", y="Ratio (%)", title = "h")



# rice transcript
rice_rna2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.c_count.bs_rep1-2_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_transcript.txt", 
                       header = T, sep = "\t", stringsAsFactors = F, quote = "")
rice_rna2$motif_ratio <- rice_rna2$moitf_num / rice_rna2$motif_total
rice_rna2$region <- factor(rice_rna2$region, 
                          levels = c("5'UTR", "CDS", 
                                     "3'UTR"))
p_rice_rna2 <- ggplot(data=rice_rna2, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=15, family = "Arial"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=28, face="bold", hjust = -0.14, vjust = -0.01)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c(" CpG    ", " CHG    ", " CHH")) +
  scale_y_continuous(limits=c(0, 0.026), breaks = seq(0, 0.025, 0.005), 
                     labels = seq(0, 0.025, 0.005) * 100) +
  labs(x="", y="Ratio (%)", title="i")





# combined figs
ppi= 300
#45, 22
png("stats_nanopore_only_cytosines_profiling/generate_Cdistri_in_repeats_n_genes.png", 
     width = 36, 
     height = 32, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab_rep,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_gene, 
                         grid.rect(gp=gpar(col="white")),
                         p_arab_rna, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_rep,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_gene, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_rna, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_rep2,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_gene2, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_rna2, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             heights=c(12, 1, 12, 1, 12))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/generate_Cdistri_in_repeats_n_genes.svg", 
    width = 36/2.54, 
    height = 32/2.54)
grid.arrange(arrangeGrob(p_arab_rep,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_gene, 
                         grid.rect(gp=gpar(col="white")),
                         p_arab_rna, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_rep,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_gene, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_rna, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_rep2,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_gene2, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice_rna2, 
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             heights=c(12, 1, 12))
dev.off()









