library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

library(grid)
library(gridExtra)
library(dplyr)

cbPalette <- c('#fc8d62', '#66c2a5')

# arab =======================
# repeats =========
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

# # keys ===
# repeat_arab_bs$key <- paste(repeat_arab_bs$chrom, repeat_arab_bs$start, 
#                             repeat_arab_bs$end, repeat_arab_bs$strand, 
#                             repeat_arab_bs$motif, sep = "||")
# keys_test <- repeat_arab_bs[repeat_arab_bs$covcf_ratio<=0.9, ]$key
# 
# repeat_arab_1k <- repeat_arab[repeat_arab$region_len>=1, ]  # 1 or 1000
# repeat_arab_1k$key <- paste(repeat_arab_1k$chrom, repeat_arab_1k$start, 
#                             repeat_arab_1k$end, repeat_arab_1k$strand, 
#                             repeat_arab_1k$motif, sep = "||")
# repeat_arab_1k <- repeat_arab_1k[repeat_arab_1k$key %in% keys_test, ]


# combine bs dp2===
repeat_arab_bs <- repeat_arab_bs[order(repeat_arab_bs$chrom, 
                                       repeat_arab_bs$start, 
                                       repeat_arab_bs$end, 
                                       repeat_arab_bs$motif), ]
repeat_arab_dp2 <- repeat_arab_dp2[order(repeat_arab_dp2$chrom, 
                                         repeat_arab_dp2$start, 
                                         repeat_arab_dp2$end, 
                                         repeat_arab_dp2$motif), ]
repeat_arab_comb <- cbind(repeat_arab_bs[, c("chrom", "start", "end", "name", "score", "strand", 
                                             "region_len", "motif",  "total_num",  
                                             "covcf_ratio", "repeattype")], 
                          repeat_arab_dp2[, c("covcf_ratio")])
colnames(repeat_arab_comb) <- c("chrom", "start", "end", "name", "score", "strand", 
                                "region_len", "motif",  "total_num",  
                                "covcf_ratio_bs", "repeattype", "covcf_ratio_dp2")

repeat_arab_comb_f <- repeat_arab_comb[repeat_arab_comb$covcf_ratio_dp2-
                                         repeat_arab_comb$covcf_ratio_bs >=0.5, ]



# genes ==============
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
# add annotation
gid2anno <- read.table("stats_nanopore_only_cytosines_profiling/Araport11_GFF3_genes_transposons.201606.gid2anno.tsv", 
                       header = F, stringsAsFactors = F, sep = "\t", quote="")
colnames(gid2anno) <- c("geneid", "annotation")



# # keys ===
# gene_arab_bs$key <- paste(gene_arab_bs$chrom, gene_arab_bs$start, 
#                           gene_arab_bs$end, gene_arab_bs$strand, 
#                           gene_arab_bs$motif, sep = "||")
# keys_test <- gene_arab_bs[gene_arab_bs$covcf_ratio<=0.9, ]$key
# 
# 
# gene_arab_1k <- gene_arab[gene_arab$gene_len>=1, ]
# gene_arab_1k$key <- paste(gene_arab_1k$chrom, gene_arab_1k$start, 
#                           gene_arab_1k$end, gene_arab_1k$strand, 
#                           gene_arab_1k$motif, sep = "||")
# gene_arab_1k <- gene_arab_1k[gene_arab_1k$key %in% keys_test, ]

# combine bs dp2===
gene_arab_bs <- gene_arab_bs[order(gene_arab_bs$chrom, 
                                   gene_arab_bs$start, 
                                   gene_arab_bs$end, 
                                   gene_arab_bs$motif), ]
gene_arab_dp2 <- gene_arab_dp2[order(gene_arab_dp2$chrom, 
                                     gene_arab_dp2$start, 
                                     gene_arab_dp2$end, 
                                     gene_arab_dp2$motif), ]
gene_arab_comb <- cbind(gene_arab_bs[, c("chrom", "start", "end", "name", "score", "strand", 
                                         "gene_len", "motif",  "total_num", "covcf_num", 
                                         "covcf_ratio")], 
                        gene_arab_dp2[, c("covcf_num", "covcf_ratio", "genetype")])
colnames(gene_arab_comb) <- c("chrom", "start", "end", "name", "score", "strand", 
                              "gene_len", "motif",  "total_num", 
                              "covcf_num_bs", "covcf_ratio_bs", 
                              "covcf_num_nano", "covcf_ratio_nano", "gene_type")

gene_arab_comb_c <- gene_arab_comb %>% 
  group_by(chrom, start, end, name, score, strand, gene_len, gene_type) %>% 
  dplyr::summarise(across(where(is.numeric), sum))
gene_arab_comb_c$covcf_ratio_bs <- gene_arab_comb_c$covcf_num_bs / gene_arab_comb_c$total_num
gene_arab_comb_c$covcf_ratio_nano <- gene_arab_comb_c$covcf_num_nano / gene_arab_comb_c$total_num
gene_arab_comb_c$motif = "C"
gene_arab_comb_c <- gene_arab_comb_c[, c("chrom", "start", "end", "name", "score", "strand",
                                     "gene_len", "gene_type", "motif", "total_num", "covcf_num_bs", 
                                     "covcf_ratio_bs", "covcf_num_nano", "covcf_ratio_nano")]

gene_arab_comb_c_f <- gene_arab_comb_c[gene_arab_comb_c$covcf_ratio_bs==0, ]
gene_arab_comb_c_f[gene_arab_comb_c_f$chrom=="NC_003070.9", ]$chrom = "chr1"
gene_arab_comb_c_f[gene_arab_comb_c_f$chrom=="NC_003071.7", ]$chrom = "chr2"
gene_arab_comb_c_f[gene_arab_comb_c_f$chrom=="NC_003074.8", ]$chrom = "chr3"
gene_arab_comb_c_f[gene_arab_comb_c_f$chrom=="NC_003075.7", ]$chrom = "chr4"
gene_arab_comb_c_f[gene_arab_comb_c_f$chrom=="NC_003076.8", ]$chrom = "chr5"

nrow(gene_arab_comb_c_f)
sum(gene_arab_comb_c_f$covcf_ratio_nano>=0.9)

write.table(gene_arab_comb_c_f, "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.arabidopsis.tsv", 
            quote = F, sep = "\t", row.names = F)
gene_arab_comb_c_f <- gene_arab_comb_c_f[gene_arab_comb_c_f$covcf_ratio_nano>=0.9, ]

# add annotation
gid2anno_c_f <- gid2anno[gid2anno$geneid %in% gene_arab_comb_c_f$name, ]
gid2anno_c_f <- gid2anno_c_f[order(gid2anno_c_f$geneid), ]
gene_arab_comb_c_f <- gene_arab_comb_c_f[order(gene_arab_comb_c_f$name), ]
gid2anno_c_f$geneid == gene_arab_comb_c_f$name
gene_arab_comb_c_f$annotation <- gid2anno_c_f$annotation
gene_arab_comb_c_f <- gene_arab_comb_c_f[order(gene_arab_comb_c_f$chrom, 
                                               gene_arab_comb_c_f$start, 
                                               gene_arab_comb_c_f$end, 
                                               gene_arab_comb_c_f$motif), ]

write.table(gene_arab_comb_c_f, "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.nano_cov90.arabidopsis.tsv", 
            quote = F, sep = "\t", row.names = F)



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

# # keys ===
# repeat_rice_bs$key <- paste(repeat_rice_bs$chrom, repeat_rice_bs$start, 
#                             repeat_rice_bs$end, repeat_rice_bs$strand, 
#                             repeat_rice_bs$motif, sep = "||")
# keys_test <- repeat_rice_bs[repeat_rice_bs$covcf_ratio<=0.9, ]$key
# 
# repeat_rice_1k <- repeat_rice[repeat_rice$region_len>=1, ]
# # repeat_rice_1k <- repeat_rice_1k[repeat_rice_1k$repeattype=="Tandem repeats",]
# repeat_rice_1k$key <- paste(repeat_rice_1k$chrom, repeat_rice_1k$start, 
#                             repeat_rice_1k$end, repeat_rice_1k$strand, 
#                             repeat_rice_1k$motif, sep = "||")
# repeat_rice_1k <- repeat_rice_1k[repeat_rice_1k$key %in% keys_test, ]

# combine bs dp2===
repeat_rice_bs <- repeat_rice_bs[order(repeat_rice_bs$chrom, 
                                       repeat_rice_bs$start, 
                                       repeat_rice_bs$end, 
                                       repeat_rice_bs$motif), ]
repeat_rice_dp2 <- repeat_rice_dp2[order(repeat_rice_dp2$chrom, 
                                         repeat_rice_dp2$start, 
                                         repeat_rice_dp2$end, 
                                         repeat_rice_dp2$motif), ]
repeat_rice_comb <- cbind(repeat_rice_bs[, c("chrom", "start", "end", "name", "score", "strand", 
                                             "region_len", "motif",  "total_num",  
                                             "covcf_ratio", "repeattype")], 
                          repeat_rice_dp2[, c("covcf_ratio")])
colnames(repeat_rice_comb) <- c("chrom", "start", "end", "name", "score", "strand", 
                                "region_len", "motif",  "total_num",  
                                "covcf_ratio_bs", "repeattype", "covcf_ratio_dp2")

repeat_rice_comb_f <- repeat_rice_comb[repeat_rice_comb$covcf_ratio_dp2-
                                         repeat_rice_comb$covcf_ratio_bs >=0.5, ]


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

# # keys ===
# gene_rice_bs$key <- paste(gene_rice_bs$chrom, gene_rice_bs$start, 
#                           gene_rice_bs$end, gene_rice_bs$strand, 
#                           gene_rice_bs$motif, sep = "||")
# keys_test <- gene_rice_bs[gene_rice_bs$covcf_ratio<=0.9, ]$key
# 
# 
# gene_rice_1k <- gene_rice[gene_rice$gene_len>=1, ]
# gene_rice_1k$key <- paste(gene_rice_1k$chrom, gene_rice_1k$start, 
#                           gene_rice_1k$end, gene_rice_1k$strand, 
#                           gene_rice_1k$motif, sep = "||")
# gene_rice_1k <- gene_rice_1k[gene_rice_1k$key %in% keys_test, ]

# combine bs dp2===
gene_rice_bs <- gene_rice_bs[order(gene_rice_bs$chrom, 
                                   gene_rice_bs$start, 
                                   gene_rice_bs$end, 
                                   gene_rice_bs$motif), ]
gene_rice_dp2 <- gene_rice_dp2[order(gene_rice_dp2$chrom, 
                                     gene_rice_dp2$start, 
                                     gene_rice_dp2$end, 
                                     gene_rice_dp2$motif), ]
gene_rice_comb <- cbind(gene_rice_bs[, c("chrom", "start", "end", "name", "score", "strand", 
                                         "gene_len", "motif",  "total_num", "covcf_num", 
                                         "covcf_ratio")], 
                        gene_rice_dp2[, c("covcf_num", "covcf_ratio", "genetype")])
colnames(gene_rice_comb) <- c("chrom", "start", "end", "name", "score", "strand", 
                              "gene_len", "motif",  "total_num", 
                              "covcf_num_bs", "covcf_ratio_bs", 
                              "covcf_num_nano", "covcf_ratio_nano", "gene_type")

gene_rice_comb_c <- gene_rice_comb %>% 
  group_by(chrom, start, end, name, score, strand, gene_len, gene_type) %>% 
  dplyr::summarise(across(where(is.numeric), sum))
gene_rice_comb_c$covcf_ratio_bs <- gene_rice_comb_c$covcf_num_bs / gene_rice_comb_c$total_num
gene_rice_comb_c$covcf_ratio_nano <- gene_rice_comb_c$covcf_num_nano / gene_rice_comb_c$total_num
gene_rice_comb_c$motif = "C"
gene_rice_comb_c <- gene_rice_comb_c[, c("chrom", "start", "end", "name", "score", "strand",
                                         "gene_len", "gene_type", "motif", "total_num", "covcf_num_bs", 
                                         "covcf_ratio_bs", "covcf_num_nano", "covcf_ratio_nano")]

gene_rice_comb_c_f <- gene_rice_comb_c[gene_rice_comb_c$covcf_ratio_bs==0, ]
gene_rice_comb_c_f$chrom <- paste("chr", gene_rice_comb_c_f$chrom, sep = "")

nrow(gene_rice_comb_c_f)
sum(gene_rice_comb_c_f$covcf_ratio_nano>=0.9)

write.table(gene_rice_comb_c_f, "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.rice_sample1.tsv", 
            quote = F, sep = "\t", row.names = F)
gene_rice_comb_c_f <- gene_rice_comb_c_f[gene_rice_comb_c_f$covcf_ratio_nano>=0.9, ]
gene_rice_comb_c_f$annotation <- "-"
write.table(gene_rice_comb_c_f, "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.nano_cov90.rice_sample1.tsv", 
            quote = F, sep = "\t", row.names = F)


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

# # keys ===
# repeat_rice2_bs$key <- paste(repeat_rice2_bs$chrom, repeat_rice2_bs$start, 
#                              repeat_rice2_bs$end, repeat_rice2_bs$strand, 
#                              repeat_rice2_bs$motif, sep = "||")
# keys_test <- repeat_rice2_bs[repeat_rice2_bs$covcf_ratio<=0.9, ]$key
# 
# repeat_rice2_1k <- repeat_rice2[repeat_rice2$region_len>=1, ]
# # repeat_rice_1k <- repeat_rice_1k[repeat_rice_1k$repeattype=="Tandem repeats",]
# repeat_rice2_1k$key <- paste(repeat_rice2_1k$chrom, repeat_rice2_1k$start, 
#                              repeat_rice2_1k$end, repeat_rice2_1k$strand, 
#                              repeat_rice2_1k$motif, sep = "||")
# repeat_rice2_1k <- repeat_rice2_1k[repeat_rice2_1k$key %in% keys_test, ]

# combine bs dp2===
repeat_rice2_bs <- repeat_rice2_bs[order(repeat_rice2_bs$chrom, 
                                         repeat_rice2_bs$start, 
                                         repeat_rice2_bs$end, 
                                         repeat_rice2_bs$motif), ]
repeat_rice2_dp2 <- repeat_rice2_dp2[order(repeat_rice2_dp2$chrom, 
                                           repeat_rice2_dp2$start, 
                                           repeat_rice2_dp2$end, 
                                           repeat_rice2_dp2$motif), ]
repeat_rice2_comb <- cbind(repeat_rice2_bs[, c("chrom", "start", "end", "name", "score", "strand", 
                                             "region_len", "motif",  "total_num",  
                                             "covcf_ratio", "repeattype")], 
                           repeat_rice2_dp2[, c("covcf_ratio")])
colnames(repeat_rice2_comb) <- c("chrom", "start", "end", "name", "score", "strand", 
                                "region_len", "motif",  "total_num",  
                                "covcf_ratio_bs", "repeattype", "covcf_ratio_dp2")

repeat_rice2_comb_f <- repeat_rice2_comb[repeat_rice2_comb$covcf_ratio_dp2-
                                         repeat_rice2_comb$covcf_ratio_bs >=0.5, ]


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

# # keys ===
# gene_rice2_bs$key <- paste(gene_rice2_bs$chrom, gene_rice2_bs$start, 
#                            gene_rice2_bs$end, gene_rice2_bs$strand, 
#                            gene_rice2_bs$motif, sep = "||")
# keys_test <- gene_rice2_bs[gene_rice2_bs$covcf_ratio<=0.9, ]$key
# 
# 
# gene_rice2_1k <- gene_rice2[gene_rice2$gene_len>=1, ]
# gene_rice2_1k$key <- paste(gene_rice2_1k$chrom, gene_rice2_1k$start, 
#                            gene_rice2_1k$end, gene_rice2_1k$strand, 
#                            gene_rice2_1k$motif, sep = "||")
# gene_rice2_1k <- gene_rice2_1k[gene_rice2_1k$key %in% keys_test, ]

# combine bs dp2===
gene_rice2_bs <- gene_rice2_bs[order(gene_rice2_bs$chrom, 
                                     gene_rice2_bs$start, 
                                     gene_rice2_bs$end, 
                                     gene_rice2_bs$motif), ]
gene_rice2_dp2 <- gene_rice2_dp2[order(gene_rice2_dp2$chrom, 
                                       gene_rice2_dp2$start, 
                                       gene_rice2_dp2$end, 
                                       gene_rice2_dp2$motif), ]
gene_rice2_comb <- cbind(gene_rice2_bs[, c("chrom", "start", "end", "name", "score", "strand", 
                                         "gene_len", "motif",  "total_num", "covcf_num", 
                                         "covcf_ratio")], 
                        gene_rice2_dp2[, c("covcf_num", "covcf_ratio", "genetype")])
colnames(gene_rice2_comb) <- c("chrom", "start", "end", "name", "score", "strand", 
                              "gene_len", "motif",  "total_num", 
                              "covcf_num_bs", "covcf_ratio_bs", 
                              "covcf_num_nano", "covcf_ratio_nano", "gene_type")

gene_rice2_comb_c <- gene_rice2_comb %>% 
  group_by(chrom, start, end, name, score, strand, gene_len, gene_type) %>% 
  dplyr::summarise(across(where(is.numeric), sum))
gene_rice2_comb_c$covcf_ratio_bs <- gene_rice2_comb_c$covcf_num_bs / gene_rice2_comb_c$total_num
gene_rice2_comb_c$covcf_ratio_nano <- gene_rice2_comb_c$covcf_num_nano / gene_rice2_comb_c$total_num
gene_rice2_comb_c$motif = "C"
gene_rice2_comb_c <- gene_rice2_comb_c[, c("chrom", "start", "end", "name", "score", "strand",
                                         "gene_len", "gene_type", "motif", "total_num", "covcf_num_bs", 
                                         "covcf_ratio_bs", "covcf_num_nano", "covcf_ratio_nano")]

gene_rice2_comb_c_f <- gene_rice2_comb_c[gene_rice2_comb_c$covcf_ratio_bs==0, ]
gene_rice2_comb_c_f$chrom <- paste("chr", gene_rice2_comb_c_f$chrom, sep = "")

nrow(gene_rice2_comb_c_f)
sum(gene_rice2_comb_c_f$covcf_ratio_nano>=0.9)

write.table(gene_rice2_comb_c_f, "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.rice_sample2.tsv", 
            quote = F, sep = "\t", row.names = F)
gene_rice2_comb_c_f <- gene_rice2_comb_c_f[gene_rice2_comb_c_f$covcf_ratio_nano>=0.9, ]
gene_rice2_comb_c_f$annotation <- "-"
write.table(gene_rice2_comb_c_f, "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.nano_cov90.rice_sample2.tsv", 
            quote = F, sep = "\t", row.names = F)


# rice rice2 intersection
gene_rice_comb_c_f$key <- paste(gene_rice_comb_c_f$chrom, gene_rice_comb_c_f$start, 
                                gene_rice_comb_c_f$end, gene_rice_comb_c_f$strand, 
                                sep = "||")
gene_rice2_comb_c_f$key <- paste(gene_rice2_comb_c_f$chrom, gene_rice2_comb_c_f$start, 
                                 gene_rice2_comb_c_f$end, gene_rice2_comb_c_f$strand, 
                                 sep = "||")
overlap_key <- intersect(gene_rice_comb_c_f$key, 
                         gene_rice2_comb_c_f$key)
gene_rice_comb_c_f <- gene_rice_comb_c_f[gene_rice_comb_c_f$key %in% overlap_key, ]
gene_rice2_comb_c_f <- gene_rice2_comb_c_f[gene_rice2_comb_c_f$key %in% overlap_key, ]
gene_ricei_comb_c_f <- cbind.data.frame(gene_rice_comb_c_f[, c("chrom", "start", "end", "name", 
                                                               "score", 
                                              "strand", "gene_len", "gene_type", "motif", 
                                              "total_num", "covcf_num_bs", "covcf_ratio_bs", 
                                              "covcf_num_nano", "covcf_ratio_nano")], 
                             gene_rice2_comb_c_f[, c("covcf_num_nano", "covcf_ratio_nano")])
colnames(gene_ricei_comb_c_f)<- c("chrom", "start", "end", "name", "score", 
                                  "strand", "gene_len", "gene_type", "motif", 
                                  "total_num", "covcf_num_bs", "covcf_ratio_bs", 
                                  "covcf_num_nano1", "covcf_ratio_nano1", 
                                  "covcf_num_nano2", "covcf_ratio_nano2")
gene_ricei_comb_c_f$annotation <- "-"
write.table(gene_ricei_comb_c_f, 
            "stats_nanopore_only_cytosines_profiling/genes_no_bs_covered.nano_cov90.rice_inter.tsv", 
            quote = F, sep = "\t", row.names = F)






