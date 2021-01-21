library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

library(grid)
library(gridExtra)

cbPalette <- c('#fc8d62', '#66c2a5')

# repeat =================================
repeat_arab_bs <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_repeats.CG.bs_rep123.covcf_5.txt", 
                                header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_arab_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                                 "region_len", "motif",  "total_num", "covcf_num", "covcf_ratio", 
                                 "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_arab_bs <- repeat_arab_bs[repeat_arab_bs$repeattype != "Repeats(WindowMasker)", ]

repeat_arab_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_repeats.CG.dp2_p0.8.covcf_5.txt", 
                                 header = F, stringsAsFactors = F, sep = "\t")
colnames(repeat_arab_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                                  "region_len", "motif", "total_num", "covcf_num", "covcf_ratio", 
                                  "cov1_num", "cov1_ratio", "repeattype", "method")
repeat_arab_dp2 <- repeat_arab_dp2[repeat_arab_dp2$repeattype != "Repeats(WindowMasker)", ]

# repeat_arab <- rbind.data.frame(repeat_arab_bs, repeat_arab_dp2)


# ======
repeat_arab_test <- cbind.data.frame(repeat_arab_bs[, c(1,2,3,4,5,6,7,8,9, 11)],
                                        repeat_arab_dp2[, c(11,14)])
colnames(repeat_arab_test)[10:11] <-c("ratio_bs", "ratio_dp2") 
repeat_arab_test <- repeat_arab_test[repeat_arab_test$ratio_dp2-repeat_arab_test$ratio_bs>=0.5 & repeat_arab_test$region_len>=5000, ]


# gene =================================
gene_arab_bs <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_genes.CG.bs_rep123.covcf_5.txt", 
                                header = F, stringsAsFactors = F, sep = "\t")
colnames(gene_arab_bs) <- c("chrom", "start", "end", "name", "score", "strand", 
                                 "motif", "total_num", "covcf_num", "covcf_ratio", 
                                 "cov1_num", "cov1_ratio", "genetype")
gene_arab_bs$method <- "bisulfite"

gene_arab_dp2 <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.cytosine_coverage.in_genes.CG.dp2_p0.8.covcf_5.txt", 
                              header = F, stringsAsFactors = F, sep = "\t")
colnames(gene_arab_dp2) <- c("chrom", "start", "end", "name", "score", "strand", 
                               "motif", "total_num", "covcf_num", "covcf_ratio", 
                               "cov1_num", "cov1_ratio", "genetype")
gene_arab_dp2$method <- "Nanopore"

# ======
gene_arab_test <- cbind.data.frame(gene_arab_bs[, c(1,2,3,4,5,6,7,8,10)],
                                      gene_arab_dp2[, c(10,13)])
colnames(gene_arab_test)[9:10] <-c("ratio_bs", "ratio_dp2") 
gene_arab_test <- gene_arab_test[gene_arab_test$ratio_dp2-gene_arab_test$ratio_bs>=0.4 & gene_arab_test$total_num>20, ]






# rice2-1 ========
# repeat =================================
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

# repeat_rice <- rbind.data.frame(repeat_rice_bs, repeat_rice_dp2)

# ======
repeat_rice_test <- cbind.data.frame(repeat_rice_bs[, c(1,2,3,4,5,6,7,8,9, 11)],
                                        repeat_rice_dp2[, c(11,14)])
repeat_rice_test <- repeat_rice_test[repeat_rice_test$motif=="CG", ]
repeat_rice_test <- repeat_rice_test[repeat_rice_test$repeattype=="Repeats(RepeatMasker)", ]
colnames(repeat_rice_test)[10:11] <-c("ratio_bs", "ratio_dp2") 
repeat_rice_test <- repeat_rice_test[repeat_rice_test$ratio_dp2-repeat_rice_test$ratio_bs>=0.4, ]


# gene ===
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

# ======
gene_rice_test <- cbind.data.frame(gene_rice_bs[, c(1,2,3,4,5,6,7,8,9, 11)],
                                   gene_rice_dp2[, c(11,14)])
gene_rice_test <- gene_rice_test[gene_rice_test$motif=="CG", ]
colnames(gene_rice_test)[10:11] <-c("ratio_bs", "ratio_dp2") 
gene_rice_test <- gene_rice_test[gene_rice_test$ratio_dp2-gene_rice_test$ratio_bs>=0.4, ]



