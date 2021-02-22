library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

library(grid)
library(gridExtra)


# arab =============
# repeats ====
repeat_arab_bs <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_repeats.CG_CHG_CHH.bs_rep123.covcf_5.txt", 
                                header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_arab_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                                 "region_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                                 "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_arab_bs <- repeat_arab_bs[repeat_arab_bs$repeattype != "Repeats(WindowMasker)", ]

repeat_arab_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_repeats.CG_CHG_CHH.dp2_p0.8.covcf_5.txt", 
                                header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_arab_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                                  "region_len", "motif", "total_num", "covcf_num", "covcf_ratio", 
                                 "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_arab_dp2 <- repeat_arab_dp2[repeat_arab_dp2$repeattype != "Repeats(WindowMasker)", ]

repeat_arab <- rbind.data.frame(repeat_arab_bs, repeat_arab_dp2)
repeat_arab$repeattype <- factor(repeat_arab$repeattype, levels = c("Repeats(RepeatMasker)", 
                                                                          "Tandem repeats", 
                                                                          "Inverted repeats"))

# keys ===
repeat_arab_bs$key <- paste(repeat_arab_bs$chrom, repeat_arab_bs$start, 
                            repeat_arab_bs$end, repeat_arab_bs$strand, 
                            repeat_arab_bs$motif, sep = "||")
keys_test <- repeat_arab_bs[repeat_arab_bs$covcf_ratio<=0.9, ]$key

repeat_arab_1k <- repeat_arab[repeat_arab$region_len>=1, ]  # 1 or 1000
repeat_arab_1k$key <- paste(repeat_arab_1k$chrom, repeat_arab_1k$start, 
                            repeat_arab_1k$end, repeat_arab_1k$strand, 
                            repeat_arab_1k$motif, sep = "||")
repeat_arab_1k <- repeat_arab_1k[repeat_arab_1k$key %in% keys_test, ]
# =
# repeat_arab_1k[repeat_arab_1k$motif=="CG", ]$motif = "CpG"
# repeat_arab_1k$motif <- factor(repeat_arab_1k$motif, levels = c("CpG", "CHG", "CHH"))

cbPalette <- c('#fc8d62', '#66c2a5')
p_repeat_arab <- ggplot(repeat_arab_1k, aes(x=motif, y=covcf_ratio, fill=method)) + 
  # geom_boxplot(position=position_dodge(0.9),
  #             lwd=0.5, outlier.size=1, fatten=2) + 
  geom_violin(position=position_dodge(), scale = "width", adjust = .8) + 
  theme_bw() + 
  facet_grid(. ~ repeattype, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=15, family = "Arial"),
        strip.text.x = element_text(size = 12),
        text = element_text(size=15, family="Arial"), 
        axis.title.y = element_text(size = 13, family = "Arial"), 
        axis.text=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=18, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1.0), 
                     breaks = seq(0, 1, 0.2), 
                     labels = seq(0, 1, 0.2)*100) + 
  scale_x_discrete(breaks=c("CG", "CHG", "CHH"), 
                   labels=c("CpG", "CHG", "CHH")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "Nanopore"), 
                    labels=c(" bisulfite    ", " Nanopore")) +
  xlab("") + 
  ylab("Percent of profiled cytosines (%)")
# p_repeat_arab

# genes ===
gene_arab_bs <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_genes.CG_CHG_CHH.bs_rep123.covcf_5.txt", 
                             header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gene_arab_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                              "gene_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                              "cov1_num", "cov1_ratio", "genetype", "method")
gene_arab_bs <- gene_arab_bs[!(gene_arab_bs$genetype %in% c("5'UTR", "CDS", "3'UTR")), ]

gene_arab_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_genes.CG_CHG_CHH.dp2_p0.8.covcf_5.txt", 
                           header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gene_arab_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                            "gene_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                            "cov1_num", "cov1_ratio", "genetype", "method")
gene_arab_dp2 <- gene_arab_dp2[!(gene_arab_dp2$genetype %in% c("5'UTR", "CDS", "3'UTR")), ]

gene_arab <- rbind.data.frame(gene_arab_bs, gene_arab_dp2)
gene_arab$genetype <- factor(gene_arab$genetype, levels = c("protein-coding genes",
                                                            "non-coding genes",
                                                            "pseudogenes",
                                                            "transposons"))

# keys ===
gene_arab_bs$key <- paste(gene_arab_bs$chrom, gene_arab_bs$start, 
                          gene_arab_bs$end, gene_arab_bs$strand, 
                          gene_arab_bs$motif, sep = "||")
keys_test <- gene_arab_bs[gene_arab_bs$covcf_ratio<=0.9, ]$key


gene_arab_1k <- gene_arab[gene_arab$gene_len>=1, ]
gene_arab_1k$key <- paste(gene_arab_1k$chrom, gene_arab_1k$start, 
                          gene_arab_1k$end, gene_arab_1k$strand, 
                          gene_arab_1k$motif, sep = "||")
gene_arab_1k <- gene_arab_1k[gene_arab_1k$key %in% keys_test, ]

cbPalette <- c('#fc8d62', '#66c2a5')
p_gene_arab <- ggplot(gene_arab_1k, aes(x=motif, y=covcf_ratio, fill=method)) + 
  # geom_boxplot(position=position_dodge(0.9),
  #             lwd=0.5, outlier.size=1, fatten=2) + 
  geom_violin(position=position_dodge(), scale = "width", adjust = .8) + 
  theme_bw() + 
  facet_grid(. ~ genetype, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=15, family = "Arial"),
        strip.text.x = element_text(size = 11),
        text = element_text(size=15, family="Arial"), 
        axis.title.y = element_text(size = 13, family = "Arial"), 
        axis.text.y=element_text(size=15, family = "Arial"), 
        axis.text.x=element_text(size=13, family = "Arial"), 
        plot.title = element_text(size=18, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2)*100) + 
  scale_x_discrete(breaks=c("CG", "CHG", "CHH"), 
                   labels=c("CpG", "CHG", "CHH")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "Nanopore"), 
                    labels=c(" bisulfite    ", " Nanopore")) +
  xlab("") + 
  ylab("Percent of profiled cytosines (%)")
# p_gene_arab


# shuidao2-1 =================
# repeats ====
repeat_rice_bs <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.cytosine_coverage.in_repeats.CG_CHG_CHH.bs_rep2-1.covcf_5.txt", 
                             header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_rice_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                              "region_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                              "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_rice_bs <- repeat_rice_bs[repeat_rice_bs$repeattype != "Repeats(WindowMasker)", ]

repeat_rice_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.cytosine_coverage.in_repeats.CG_CHG_CHH.dp2_p0.8.covcf_5.txt", 
                              header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_rice_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                               "region_len", "motif", "total_num", "covcf_num", "covcf_ratio", 
                               "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_rice_dp2 <- repeat_rice_dp2[repeat_rice_dp2$repeattype != "Repeats(WindowMasker)", ]

repeat_rice <- rbind.data.frame(repeat_rice_bs, repeat_rice_dp2)
repeat_rice$repeattype <- factor(repeat_rice$repeattype, levels = c("Repeats(RepeatMasker)", 
                                                                    "Tandem repeats", 
                                                                    "Inverted repeats"))

# keys ===
repeat_rice_bs$key <- paste(repeat_rice_bs$chrom, repeat_rice_bs$start, 
                            repeat_rice_bs$end, repeat_rice_bs$strand, 
                            repeat_rice_bs$motif, sep = "||")
keys_test <- repeat_rice_bs[repeat_rice_bs$covcf_ratio<=0.9, ]$key

repeat_rice_1k <- repeat_rice[repeat_rice$region_len>=1, ]
# repeat_rice_1k <- repeat_rice_1k[repeat_rice_1k$repeattype=="Tandem repeats",]
repeat_rice_1k$key <- paste(repeat_rice_1k$chrom, repeat_rice_1k$start, 
                            repeat_rice_1k$end, repeat_rice_1k$strand, 
                            repeat_rice_1k$motif, sep = "||")
repeat_rice_1k <- repeat_rice_1k[repeat_rice_1k$key %in% keys_test, ]

cbPalette <- c('#fc8d62', '#66c2a5')
p_repeat_rice <- ggplot(repeat_rice_1k, aes(x=motif, y=covcf_ratio, fill=method)) + 
  # geom_boxplot(position=position_dodge(0.9),
  #             lwd=0.5, outlier.size=1, fatten=2) + 
  geom_violin(position=position_dodge(), scale = "width", adjust = .8) + 
  theme_bw() + 
  facet_grid(. ~ repeattype, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=15, family = "Arial"),
        strip.text.x = element_text(size = 12),
        text = element_text(size=15, family="Arial"), 
        axis.title.y = element_text(size = 13, family = "Arial"), 
        axis.text=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=18, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2)*100) + 
  scale_x_discrete(breaks=c("CG", "CHG", "CHH"), 
                   labels=c("CpG", "CHG", "CHH")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "Nanopore"), 
                    labels=c(" bisulfite    ", " Nanopore")) +
  xlab("") + 
  ylab("Percent of profiled cytosines (%)")
# p_repeat_rice


# genes ===
gene_rice_bs <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.cytosine_coverage.in_genes.CG_CHG_CHH.bs_rep2-1.covcf_5.txt", 
                           header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gene_rice_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                            "gene_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                            "cov1_num", "cov1_ratio", "genetype", "method")
gene_rice_bs <- gene_rice_bs[!(gene_rice_bs$genetype %in% c("5'UTR", "CDS", "3'UTR")), ]

gene_rice_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.cytosine_coverage.in_genes.CG_CHG_CHH.dp2_p0.8.covcf_5.txt", 
                            header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gene_rice_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                             "gene_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                             "cov1_num", "cov1_ratio", "genetype", "method")
gene_rice_dp2 <- gene_rice_dp2[!(gene_rice_dp2$genetype %in% c("5'UTR", "CDS", "3'UTR")), ]

gene_rice <- rbind.data.frame(gene_rice_bs, gene_rice_dp2)
gene_rice$genetype <- factor(gene_rice$genetype, levels = c("protein-coding genes",
                                                            "non-coding genes"))

# keys ===
gene_rice_bs$key <- paste(gene_rice_bs$chrom, gene_rice_bs$start, 
                          gene_rice_bs$end, gene_rice_bs$strand, 
                          gene_rice_bs$motif, sep = "||")
keys_test <- gene_rice_bs[gene_rice_bs$covcf_ratio<=0.9, ]$key


gene_rice_1k <- gene_rice[gene_rice$gene_len>=1, ]
gene_rice_1k$key <- paste(gene_rice_1k$chrom, gene_rice_1k$start, 
                          gene_rice_1k$end, gene_rice_1k$strand, 
                          gene_rice_1k$motif, sep = "||")
gene_rice_1k <- gene_rice_1k[gene_rice_1k$key %in% keys_test, ]

cbPalette <- c('#fc8d62', '#66c2a5')
p_gene_rice <- ggplot(gene_rice_1k, aes(x=motif, y=covcf_ratio, fill=method)) + 
  # geom_boxplot(position=position_dodge(0.9),
  #             lwd=0.5, outlier.size=1, fatten=2) + 
  geom_violin(position=position_dodge(), scale = "width", adjust = .8) + 
  theme_bw() + 
  facet_grid(. ~ genetype, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=15, family = "Arial"),
        strip.text.x = element_text(size = 12),
        text = element_text(size=15, family="Arial"), 
        axis.title.y = element_text(size = 13, family = "Arial"), 
        axis.text=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=18, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2)*100) + 
  scale_x_discrete(breaks=c("CG", "CHG", "CHH"), 
                   labels=c("CpG", "CHG", "CHH")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "Nanopore"), 
                    labels=c(" bisulfite    ", " Nanopore")) +
  xlab("") + 
  ylab("Percent of profiled cytosines (%)")
# p_gene_rice



# shuidao1-1 =================
# repeats ====
repeat_rice2_bs <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.cytosine_coverage.in_repeats.CG_CHG_CHH.bs_rep1-2.covcf_5.txt", 
                             header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_rice2_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                              "region_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                              "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_rice2_bs <- repeat_rice2_bs[repeat_rice2_bs$repeattype != "Repeats(WindowMasker)", ]

repeat_rice2_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.cytosine_coverage.in_repeats.CG_CHG_CHH.dp2_p0.8.covcf_5.txt", 
                              header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_rice2_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                               "region_len", "motif", "total_num", "covcf_num", "covcf_ratio", 
                               "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_rice2_dp2 <- repeat_rice2_dp2[repeat_rice2_dp2$repeattype != "Repeats(WindowMasker)", ]

repeat_rice2 <- rbind.data.frame(repeat_rice2_bs, repeat_rice2_dp2)
repeat_rice2$repeattype <- factor(repeat_rice2$repeattype, levels = c("Repeats(RepeatMasker)", 
                                                                    "Tandem repeats", 
                                                                    "Inverted repeats"))

# keys ===
repeat_rice2_bs$key <- paste(repeat_rice2_bs$chrom, repeat_rice2_bs$start, 
                             repeat_rice2_bs$end, repeat_rice2_bs$strand, 
                             repeat_rice2_bs$motif, sep = "||")
keys_test <- repeat_rice2_bs[repeat_rice2_bs$covcf_ratio<=0.9, ]$key

repeat_rice2_1k <- repeat_rice2[repeat_rice2$region_len>=1, ]
# repeat_rice_1k <- repeat_rice_1k[repeat_rice_1k$repeattype=="Tandem repeats",]
repeat_rice2_1k$key <- paste(repeat_rice2_1k$chrom, repeat_rice2_1k$start, 
                             repeat_rice2_1k$end, repeat_rice2_1k$strand, 
                             repeat_rice2_1k$motif, sep = "||")
repeat_rice2_1k <- repeat_rice2_1k[repeat_rice2_1k$key %in% keys_test, ]

cbPalette <- c('#fc8d62', '#66c2a5')
p_repeat_rice2 <- ggplot(repeat_rice2_1k, aes(x=motif, y=covcf_ratio, fill=method)) + 
  # geom_boxplot(position=position_dodge(0.9),
  #             lwd=0.5, outlier.size=1, fatten=2) + 
  geom_violin(position=position_dodge(), scale = "width", adjust = .8) + 
  theme_bw() + 
  facet_grid(. ~ repeattype, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=15, family = "Arial"),
        strip.text.x = element_text(size = 12),
        text = element_text(size=15, family="Arial"), 
        axis.title.y = element_text(size = 13, family = "Arial"),  
        axis.text=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=18, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2)*100) + 
  scale_x_discrete(breaks=c("CG", "CHG", "CHH"), 
                   labels=c("CpG", "CHG", "CHH")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "Nanopore"), 
                    labels=c(" bisulfite    ", " Nanopore")) +
  xlab("") + 
  ylab("Percent of profiled cytosines (%)")
# p_repeat_rice2


# genes ===
gene_rice2_bs <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.cytosine_coverage.in_genes.CG_CHG_CHH.bs_rep1-2.covcf_5.txt", 
                           header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gene_rice2_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                            "gene_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                            "cov1_num", "cov1_ratio", "genetype", "method")
gene_rice2_bs <- gene_rice2_bs[!(gene_rice2_bs$genetype %in% c("5'UTR", "CDS", "3'UTR")), ]

gene_rice2_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.cytosine_coverage.in_genes.CG_CHG_CHH.dp2_p0.8.covcf_5.txt", 
                            header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gene_rice2_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                             "gene_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                             "cov1_num", "cov1_ratio", "genetype", "method")
gene_rice2_dp2 <- gene_rice2_dp2[!(gene_rice2_dp2$genetype %in% c("5'UTR", "CDS", "3'UTR")), ]

gene_rice2 <- rbind.data.frame(gene_rice2_bs, gene_rice2_dp2)
gene_rice2$genetype <- factor(gene_rice2$genetype, levels = c("protein-coding genes",
                                                            "non-coding genes"))

# keys ===
gene_rice2_bs$key <- paste(gene_rice2_bs$chrom, gene_rice2_bs$start, 
                           gene_rice2_bs$end, gene_rice2_bs$strand, 
                           gene_rice2_bs$motif, sep = "||")
keys_test <- gene_rice2_bs[gene_rice2_bs$covcf_ratio<=0.9, ]$key


gene_rice2_1k <- gene_rice2[gene_rice2$gene_len>=1, ]
gene_rice2_1k$key <- paste(gene_rice2_1k$chrom, gene_rice2_1k$start, 
                           gene_rice2_1k$end, gene_rice2_1k$strand, 
                           gene_rice2_1k$motif, sep = "||")
gene_rice2_1k <- gene_rice2_1k[gene_rice2_1k$key %in% keys_test, ]

cbPalette <- c('#fc8d62', '#66c2a5')
p_gene_rice2 <- ggplot(gene_rice2_1k, aes(x=motif, y=covcf_ratio, fill=method)) + 
  # geom_boxplot(position=position_dodge(0.9),
  #             lwd=0.5, outlier.size=1, fatten=2) + 
  geom_violin(position=position_dodge(), scale = "width", adjust = .8) + 
  theme_bw() + 
  facet_grid(. ~ genetype, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=15, family = "Arial"),
        strip.text.x = element_text(size = 12),
        text = element_text(size=15, family="Arial"), 
        axis.title.y = element_text(size = 13, family = "Arial"), 
        axis.text=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size=18, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2)*100) + 
  scale_x_discrete(breaks=c("CG", "CHG", "CHH"), 
                   labels=c("CpG", "CHG", "CHH")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "Nanopore"), 
                    labels=c(" bisulfite    ", " Nanopore")) +
  xlab("") + 
  ylab("Percent of profiled cytosines (%)")
# p_gene_rice2


# combined figs
ppi= 300
#45, 39
png("stats_nanopore_only_cytosines_profiling/generate_Ccoverperc_in_repeats_n_genes_len1_bs0.9.png", 
    width = 36, 
    height = 32, units = "cm", res=ppi)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_repeat_arab,
                         grid.rect(gp=gpar(col="white")),
                         p_gene_arab, 
                         nrow=1, 
                         widths = c(22, 1, 22)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_repeat_rice,
                         grid.rect(gp=gpar(col="white")),
                         p_gene_rice, 
                         nrow=1, 
                         widths = c(22, 1, 22)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_repeat_rice2,
                         grid.rect(gp=gpar(col="white")),
                         p_gene_rice2, 
                         nrow=1, 
                         widths = c(22, 1, 22)), 
             heights=c(1, 12, 1, 12, 1, 12))
dev.off()

svg("stats_nanopore_only_cytosines_profiling/generate_Ccoverperc_in_repeats_n_genes_len1_bs0.9.svg", 
    width = 36/2.54, 
    height = 32/2.54)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_repeat_arab,
                         grid.rect(gp=gpar(col="white")),
                         p_gene_arab, 
                         nrow=1, 
                         widths = c(22, 1, 22)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_repeat_rice,
                         grid.rect(gp=gpar(col="white")),
                         p_gene_rice, 
                         nrow=1, 
                         widths = c(22, 1, 22)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_repeat_rice2,
                         grid.rect(gp=gpar(col="white")),
                         p_gene_rice2, 
                         nrow=1, 
                         widths = c(22, 1, 22)), 
             heights=c(1, 12, 1, 12, 1, 12))
dev.off()






# find example to plot






















