library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(eulerr)
library(grid)
library(gridExtra)
library(UpSetR)

stat_colnames <- c("rpair_name", "num_all", "num_cov", "num_diff", "ratio_cov", 
  "ratio_diff2all", "ratio_diff2cov", "repeat_len")

cbPalette <- c("#1f78b4", "#33a02c")

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

# arab ====================================================================
arab_dp2_CG_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CG.deepsignal.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_dp2_CG_stats) <- stat_colnames
arab_dp2_CG_stats <- arab_dp2_CG_stats[arab_dp2_CG_stats$repeat_len>=100, ]
arab_dp2_diffid_CG <- arab_dp2_CG_stats[arab_dp2_CG_stats$ratio_diff2all>=0.1,]$rpair_name

arab_dp2_CHG_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHG.deepsignal.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_dp2_CHG_stats) <- stat_colnames
arab_dp2_CHG_stats <- arab_dp2_CHG_stats[arab_dp2_CHG_stats$repeat_len>=100, ]
arab_dp2_diffid_CHG <- arab_dp2_CHG_stats[arab_dp2_CHG_stats$ratio_diff2all>=0.1,]$rpair_name

arab_dp2_CHH_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHH.deepsignal.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_dp2_CHH_stats) <- stat_colnames
arab_dp2_CHH_stats <- arab_dp2_CHH_stats[arab_dp2_CHH_stats$repeat_len>=100, ]
arab_dp2_diffid_CHH <- arab_dp2_CHH_stats[arab_dp2_CHH_stats$ratio_diff2all>=0.1,]$rpair_name

arab_dp2_C_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.C.deepsignal.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_dp2_C_stats) <- stat_colnames
arab_dp2_C_stats <- arab_dp2_C_stats[arab_dp2_C_stats$repeat_len>=100, ]
arab_dp2_diffid_C <- arab_dp2_C_stats[arab_dp2_C_stats$ratio_diff2all>=0.1,]$rpair_name

# ==================== cytosine histogram
p_diffratio_arab <- ggplot(arab_dp2_C_stats, aes(x=ratio_diff2all)) + 
  geom_histogram(colour="black", binwidth = 0.02, size=0.3, fill=cbPalette[1]) + 
  geom_vline(xintercept=0.1, 
             colour="black", 
             linetype="dashed", size=0.6) +
  theme_bw() + 
  theme(text=element_text(size=12, family = "Arial"),
        legend.position = "None", 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=8, family = "Arial"), 
        axis.title.y=element_text(size=7.2, family = "Arial")) + 
  scale_x_continuous(# limits=c(0, 1), 
                     breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=3),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_arab



arab_dp2_CG_stats$motif <- "CpG"
arab_dp2_CHG_stats$motif <- "CHG"
arab_dp2_CHH_stats$motif <- "CHH"
arab_dp2_motif_stats <- rbind.data.frame(arab_dp2_CG_stats, 
                                         arab_dp2_CHG_stats, 
                                         arab_dp2_CHH_stats)
arab_dp2_motif_stats$motif <- factor(arab_dp2_motif_stats$motif, levels = c("CpG", 
                                                                            "CHG", 
                                                                            "CHH"))
p_diffratio_arab_motif <- ggplot(arab_dp2_motif_stats, aes(x=ratio_diff2all, fill=motif)) + 
  geom_histogram(colour="black", binwidth = 0.02, size=0.3) + 
  geom_vline(xintercept=0.1, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  facet_grid(.~motif, scales = "free") + 
  theme(text=element_text(size=18, family = "Arial"),
        legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15), 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=15, family = "Arial"), 
        axis.title.y=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size = 18, family = "Arial", hjust = 0.5)) + 
  scale_fill_brewer(palette = "Set2", 
                    breaks = c("CG", "CHG", "CHH"), 
                    labels = c("CpG    ", "CHG    ", "CHH")) + 
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Density of repeat pairs", 
       title = expression(italic("A. thaliana"))) + 
  scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x, n=3),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_arab_motif



# ==================== intersect venn
arab_dp2_rdiff_pairs <- list(C=arab_dp2_diffid_C, 
                    CpG=arab_dp2_diffid_CG, 
                    CHG=arab_dp2_diffid_CHG,
                    CHH=arab_dp2_diffid_CHH)
# plot venn with euler =================
# combs <- 
#   unlist(lapply(1:length(rdiff_pairs), 
#                 function(j) combn(names(rdiff_pairs), j, simplify = FALSE)),
#          recursive = FALSE)
# names(combs) <- sapply(combs, function(i) paste0(i, collapse = "&"))
# # str(combs)
# elements <- 
#   lapply(combs, function(i) Setdiff(rdiff_pairs[i], rdiff_pairs[setdiff(names(rdiff_pairs), i)]))
# n.elements <- sapply(elements, length)
# # print(n.elements)
# arab_dp2_diffid_venn <- euler(n.elements)
# vennpalette <- c("#e78ac3", "#66c2a5", "#fc8d62", "#8da0cb")
# p_arab_dp2_diffid_venn <- plot(arab_dp2_diffid_venn, 
#                  quantities = list(fontsize=16, fontfamily="Arial", 
#                                    labels=unname(n.elements)), 
#                  fills = list(fill=vennpalette, alpha=0.7),
#                  edges =TRUE,
#                  labels = list(font=2, fontsize=16, fontfamily="Arial"), 
#                  main = list(label="", 
#                              fontsize=24,
#                              font=2,
#                              hjust=8,
#                              vjust=1, 
#                              fontfamily="Arial"))
# p_arab_dp2_diffid_venn
# plot venn with UpSetR =================
p_arab_dp2_diffid_venn <- upset(fromList(arab_dp2_rdiff_pairs), 
      sets = c("CHH", "CHG", "CpG", "C"), 
      order.by = "freq", 
      keep.order = TRUE,
      empty.intersections = "on",
      point.size = 4, line.size = 1.5,
      mainbar.y.label = "Intersection size", sets.x.label = "Set size", 
      matrix.color = cbPalette[1], main.bar.color = cbPalette[1],
      sets.bar.color = cbPalette[1], text.scale = c(2, 2, 2, 2, 2, 2))
p_arab_dp2_diffid_venn


# rice ====================================================================
# rice 2-1
rice_dp2_CG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CG.deepsignal_2-1.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_dp2_CG_stats) <- stat_colnames
rice_dp2_CG_stats <- rice_dp2_CG_stats[rice_dp2_CG_stats$repeat_len>=100, ]
rice_dp2_diffid_CG <- rice_dp2_CG_stats[rice_dp2_CG_stats$ratio_diff2all>=0.1,]$rpair_name

rice_dp2_CHG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHG.deepsignal_2-1.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_dp2_CHG_stats) <- stat_colnames
rice_dp2_CHG_stats <- rice_dp2_CHG_stats[rice_dp2_CHG_stats$repeat_len>=100, ]
rice_dp2_diffid_CHG <- rice_dp2_CHG_stats[rice_dp2_CHG_stats$ratio_diff2all>=0.1,]$rpair_name

rice_dp2_CHH_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHH.deepsignal_2-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_dp2_CHH_stats) <- stat_colnames
rice_dp2_CHH_stats <- rice_dp2_CHH_stats[rice_dp2_CHH_stats$repeat_len>=100, ]
rice_dp2_diffid_CHH <- rice_dp2_CHH_stats[rice_dp2_CHH_stats$ratio_diff2all>=0.1,]$rpair_name

rice_dp2_C_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.C.deepsignal_2-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_dp2_C_stats) <- stat_colnames
rice_dp2_C_stats <- rice_dp2_C_stats[rice_dp2_C_stats$repeat_len>=100, ]
rice_dp2_diffid_C <- rice_dp2_C_stats[rice_dp2_C_stats$ratio_diff2all>=0.1,]$rpair_name

# rice 1-1
rice2_dp2_CG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CG.deepsignal_1-1.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_dp2_CG_stats) <- stat_colnames
rice2_dp2_CG_stats <- rice2_dp2_CG_stats[rice2_dp2_CG_stats$repeat_len>=100, ]
rice2_dp2_diffid_CG <- rice2_dp2_CG_stats[rice2_dp2_CG_stats$ratio_diff2all>=0.1,]$rpair_name

rice2_dp2_CHG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHG.deepsignal_1-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_dp2_CHG_stats) <- stat_colnames
rice2_dp2_CHG_stats <- rice2_dp2_CHG_stats[rice2_dp2_CHG_stats$repeat_len>=100, ]
rice2_dp2_diffid_CHG <- rice2_dp2_CHG_stats[rice2_dp2_CHG_stats$ratio_diff2all>=0.1,]$rpair_name

rice2_dp2_CHH_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHH.deepsignal_1-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_dp2_CHH_stats) <- stat_colnames
rice2_dp2_CHH_stats <- rice2_dp2_CHH_stats[rice2_dp2_CHH_stats$repeat_len>=100, ]
rice2_dp2_diffid_CHH <- rice2_dp2_CHH_stats[rice2_dp2_CHH_stats$ratio_diff2all>=0.1,]$rpair_name

rice2_dp2_C_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.C.deepsignal_1-1.all2_repeat_group_stats.txt", 
                               header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_dp2_C_stats) <- stat_colnames
rice2_dp2_C_stats <- rice2_dp2_C_stats[rice2_dp2_C_stats$repeat_len>=100, ]
rice2_dp2_diffid_C <- rice2_dp2_C_stats[rice2_dp2_C_stats$ratio_diff2all>=0.1,]$rpair_name


# ==================== cytosine histogram
# rice_dp2_C_stats$replicate <- "rep1"
# rice2_dp2_C_stats$replicate <- "rep2"
# rice_dp2_C_stats_rep12 <- rbind.data.frame(rice_dp2_C_stats, rice2_dp2_C_stats)
# cbPalette_r <- c("#fdcc8a", "#d7301f")
p_diffratio_rice <- ggplot(rice_dp2_C_stats, aes(x=ratio_diff2all)) + 
  geom_histogram(colour="black", binwidth = 0.02, size=0.3, fill=cbPalette[2]) + 
  geom_vline(xintercept=0.1, 
             colour="black", 
             linetype="dashed", size=0.6) +
  theme_bw() + 
  theme(text=element_text(size=12, family = "Arial"),
        legend.position = "None", 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=8, family = "Arial"), 
        axis.title.y=element_text(size=7.2, family = "Arial")) + 
  scale_x_continuous(# limits=c(0, 1), 
    breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_rice


rice_dp2_CG_stats$motif <- "CpG"
rice_dp2_CHG_stats$motif <- "CHG"
rice_dp2_CHH_stats$motif <- "CHH"
rice_dp2_motif_stats <- rbind.data.frame(rice_dp2_CG_stats, 
                                         rice_dp2_CHG_stats, 
                                         rice_dp2_CHH_stats)
rice_dp2_motif_stats$motif <- factor(rice_dp2_motif_stats$motif, levels = c("CpG", 
                                                                            "CHG", 
                                                                            "CHH"))
p_diffratio_rice_motif <- ggplot(rice_dp2_motif_stats, aes(x=ratio_diff2all, fill=motif)) + 
  geom_histogram(colour="black", binwidth = 0.02, size=0.3) + 
  geom_vline(xintercept=0.1, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  facet_grid(.~motif, scales = "free") + 
  theme(text=element_text(size=18, family = "Arial"),
        legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15), 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=15, family = "Arial"), 
        axis.title.y=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size = 18, family = "Arial", hjust = 0.5)) + 
  scale_fill_brewer(palette = "Set2", 
                    breaks = c("CG", "CHG", "CHH"), 
                    labels = c("CpG    ", "CHG    ", "CHH")) + 
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)", 
       title = expression(paste(italic("O. sativa"), " (sample1)", sep=""))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_rice_motif





p_diffratio_rice2 <- ggplot(rice2_dp2_C_stats, aes(x=ratio_diff2all)) + 
  geom_histogram(colour="black", binwidth = 0.02, size=0.3, fill=cbPalette[2]) + 
  geom_vline(xintercept=0.1, 
             colour="black", 
             linetype="dashed", size=0.6) +
  theme_bw() + 
  theme(text=element_text(size=12, family = "Arial"),
        legend.position = "None", 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=8, family = "Arial"), 
        axis.title.y=element_text(size=7.2, family = "Arial")) + 
  scale_x_continuous(# limits=c(0, 1), 
    breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_rice2



rice2_dp2_CG_stats$motif <- "CpG"
rice2_dp2_CHG_stats$motif <- "CHG"
rice2_dp2_CHH_stats$motif <- "CHH"
rice2_dp2_motif_stats <- rbind.data.frame(rice2_dp2_CG_stats, 
                                         rice2_dp2_CHG_stats, 
                                         rice2_dp2_CHH_stats)
rice2_dp2_motif_stats$motif <- factor(rice2_dp2_motif_stats$motif, levels = c("CpG", 
                                                                            "CHG", 
                                                                            "CHH"))
p_diffratio_rice2_motif <- ggplot(rice2_dp2_motif_stats, aes(x=ratio_diff2all, fill=motif)) + 
  geom_histogram(colour="black", binwidth = 0.02, size=0.3) + 
  geom_vline(xintercept=0.1, 
             colour="black", 
             linetype="dashed", size=0.8) +
  theme_bw() + 
  facet_grid(.~motif, scales = "free") + 
  theme(text=element_text(size=18, family = "Arial"),
        legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15), 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=15, family = "Arial"), 
        axis.title.y=element_text(size=15, family = "Arial"), 
        plot.title = element_text(size = 18, family = "Arial", hjust = 0.5)) + 
  scale_fill_brewer(palette = "Set2", 
                    breaks = c("CG", "CHG", "CHH"), 
                    labels = c("CpG    ", "CHG    ", "CHH")) + 
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)", 
       title = expression(paste(italic("O. sativa"), " (sample2)", sep=""))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_rice2_motif


# ==================== intersect venn
# rice 2-1
rice_dp2_rdiff_pairs <- list(C=rice_dp2_diffid_C, 
                             CpG=rice_dp2_diffid_CG, 
                             CHG=rice_dp2_diffid_CHG,
                             CHH=rice_dp2_diffid_CHH)
# plot venn with euler =================
# combs <-
#   unlist(lapply(1:length(rice_dp2_rdiff_pairs),
#                 function(j) combn(names(rice_dp2_rdiff_pairs), j, simplify = FALSE)),
#          recursive = FALSE)
# names(combs) <- sapply(combs, function(i) paste0(i, collapse = "&"))
# # str(combs)
# elements <-
#   lapply(combs, function(i) Setdiff(rice_dp2_rdiff_pairs[i], rice_dp2_rdiff_pairs[setdiff(names(rice_dp2_rdiff_pairs), i)]))
# n.elements <- sapply(elements, length)
# # print(n.elements)
# rice_dp2_diffid_venn <- euler(n.elements)
# vennpalette <- c("#e78ac3", "#66c2a5", "#fc8d62", "#8da0cb")
# p_rice_dp2_diffid_venn <- plot(rice_dp2_diffid_venn,
#                  quantities = list(fontsize=16, fontfamily="Arial",
#                                    labels=unname(n.elements)),
#                  fills = list(fill=vennpalette, alpha=0.7),
#                  edges =TRUE,
#                  labels = list(font=2, fontsize=16, fontfamily="Arial"),
#                  main = list(label="",
#                              fontsize=24,
#                              font=2,
#                              hjust=8,
#                              vjust=1,
#                              fontfamily="Arial"))
# plot venn with UpSetR =================
p_rice_dp2_diffid_venn <- upset(fromList(rice_dp2_rdiff_pairs), 
                sets = c("CHH", "CHG", "CpG", "C"), 
                order.by = "freq", 
                keep.order = TRUE,
                empty.intersections = "on",
                point.size = 4, line.size = 1.5,
                mainbar.y.label = "Intersection size", sets.x.label = "Set size", 
                matrix.color = cbPalette[2], main.bar.color = cbPalette[2],
                sets.bar.color = cbPalette[2], text.scale = c(2, 2, 2, 1, 2, 1.5))
p_rice_dp2_diffid_venn

# rice 1-1
rice2_dp2_rdiff_pairs <- list(C=rice2_dp2_diffid_C, 
                             CpG=rice2_dp2_diffid_CG, 
                             CHG=rice2_dp2_diffid_CHG,
                             CHH=rice2_dp2_diffid_CHH)
# plot venn with UpSetR =================
p_rice2_dp2_diffid_venn <- upset(fromList(rice2_dp2_rdiff_pairs), 
                 sets = c("CHH", "CHG", "CpG", "C"), 
                 order.by = "freq", 
                 keep.order = TRUE,
                 empty.intersections = "on",
                 point.size = 4, line.size = 1.5,
                 mainbar.y.label = "Intersection size", sets.x.label = "Set size", 
                 matrix.color = cbPalette[2], main.bar.color = cbPalette[2],
                 sets.bar.color = cbPalette[2], text.scale = c(2, 2, 2, 1, 2, 1.5))
p_rice2_dp2_diffid_venn


# venn rice dp2 rep1 (2-1) vs rep2 (1-1) =======================================
rice_dp2_diffC_rep12_venn_list <- list(sample1=rice_dp2_diffid_C, 
                                         sample2=rice2_dp2_diffid_C)

# plot venn with euler =================
rice_dp2_diffC_rep12_combs <-
  unlist(lapply(1:length(rice_dp2_diffC_rep12_venn_list),
                function(j) combn(names(rice_dp2_diffC_rep12_venn_list), j, simplify = FALSE)),
         recursive = FALSE)
names(rice_dp2_diffC_rep12_combs) <- sapply(rice_dp2_diffC_rep12_combs, function(i) paste0(i, collapse = "&"))
# str(combs)
rice_dp2_diffC_rep12_elements <-
  lapply(rice_dp2_diffC_rep12_combs, function(i) Setdiff(rice_dp2_diffC_rep12_venn_list[i], rice_dp2_diffC_rep12_venn_list[setdiff(names(rice_dp2_diffC_rep12_venn_list), i)]))
rice_dp2_diffC_rep12_n.elements <- sapply(rice_dp2_diffC_rep12_elements, length)
# print(n.elements)
rice_dp2_diffC_rep12_venn <- euler(rice_dp2_diffC_rep12_n.elements)
vennpalette <- c("#78c679", "#006837")
p_rice_dp2_diffC_rep12_venn <- plot(rice_dp2_diffC_rep12_venn,
                 quantities = list(fontsize=12, fontfamily="Arial",
                                   labels=trimws(format(unname(rice_dp2_diffC_rep12_n.elements), 
                                                        big.mark = ",", 
                                                        scientific = F))),
                 fills = list(fill=vennpalette, alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=12, fontfamily="Arial"), 
                 main = list(label="O. sativa", 
                             fontsize=13,
                             font=3,
                             vjust=0.8,
                             fontfamily="Arial"))
p_rice_dp2_diffC_rep12_venn


# ppi= 300
# png("methyrep/fig_diffmethyregions_info2.raw.png", 
#     width = 48, 
#     height = 17, units = "cm", res=ppi)
# grid.arrange(arrangeGrob(p_diffratio_arab,
#                          grid.rect(gp=gpar(col="white")),
#                          p_diffratio_rice,
#                          ncol=1,
#                          nrow=3,
#                          heights = c(8, 1, 8)), 
#              p_arab_dp2_diffid_venn,
#              p_rice_dp2_diffid_venn, 
#              p_rice_dp2_diffC_rep12_venn, 
#              nrow = 1,
#              ncol = 4,
#              widths=c(12, 12, 12, 12))
# dev.off()


# fig in main text ============================================
ppi= 300
png("methyrep/fig_diffmethyregions_info2.a.raw.png", 
    width = 7, 
    height = 4.7, units = "cm", res=ppi)
p_diffratio_arab
dev.off()

png("methyrep/fig_diffmethyregions_info2.c.raw.png", 
    width = 7, 
    height = 4.7, units = "cm", res=ppi)
p_diffratio_rice
dev.off()

png("methyrep/fig_diffmethyregions_info2.b.raw.png", 
    width = 14.5, 
    height = 15, units = "cm", res=ppi)
p_arab_dp2_diffid_venn
dev.off()

png("methyrep/fig_diffmethyregions_info2.d.raw.png", 
    width = 14.5, 
    height = 15, units = "cm", res=ppi)
p_rice_dp2_diffid_venn
dev.off()

png("methyrep/fig_diffmethyregions_info2.e.raw.png", 
    width = 7, 
    height = 6, units = "cm", res=ppi)
p_rice_dp2_diffC_rep12_venn
dev.off()


# fig of rice rep2 ==================
png("methyrep/fig_diffmethyregions_info2.rice_rep2.a.raw.png", 
    width =7, 
    height = 4.7, units = "cm", res=ppi)
p_diffratio_rice2
dev.off()

png("methyrep/fig_diffmethyregions_info2.rice_rep2.b.raw.png", 
    width = 14.5, 
    height = 15, units = "cm", res=ppi)
p_rice2_dp2_diffid_venn
dev.off()


# motif histogram ===================
ppi=300
png("methyrep/fig_diffmethyregions_info2.motif_histo.raw.png", 
    width = 36, 
    height = 32, units = "cm", res=ppi)
grid.arrange(p_diffratio_arab_motif, 
             grid.rect(gp=gpar(col="white")),
             p_diffratio_rice_motif,
             grid.rect(gp=gpar(col="white")),
             p_diffratio_rice2_motif,
             nrow = 5,
             ncol = 1,
             heights=c(10, 1, 10, 1, 10))
dev.off()
svg("methyrep/fig_diffmethyregions_info2.motif_histo.raw.svg", 
    width = 36/2.54, 
    height = 32/2.54)
grid.arrange(p_diffratio_arab_motif, 
             grid.rect(gp=gpar(col="white")),
             p_diffratio_rice_motif,
             grid.rect(gp=gpar(col="white")),
             p_diffratio_rice2_motif,
             nrow = 5,
             ncol = 1,
             heights=c(10, 1, 10, 1, 10))
dev.off()













# vs bs ======================================
# function ===
make_two_sets_venn_data <- function(bsset, nanoset){
  venn_list <- list(bisulfite=bsset, Nanopore=nanoset)
  # plot venn with euler =================
  venn_combs <-
    unlist(lapply(1:length(venn_list),
                  function(j) combn(names(venn_list), j, simplify = FALSE)),
           recursive = FALSE)
  names(venn_combs) <- sapply(venn_combs, function(i) paste0(i, collapse = "&"))
  # str(combs)
  venn_elements <-
    lapply(venn_combs, function(i) Setdiff(venn_list[i], venn_list[setdiff(names(venn_list), i)]))
  venn_n.elements <- sapply(venn_elements, length)
  return(venn_n.elements)
}


# arab ====================================================================
arab_bs_CG_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CG.bedmethyl.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_bs_CG_stats) <- stat_colnames
arab_bs_CG_stats <- arab_bs_CG_stats[arab_bs_CG_stats$repeat_len>=100, ]
arab_bs_diffid_CG <- arab_bs_CG_stats[arab_bs_CG_stats$ratio_diff2all>=0.1,]$rpair_name

arab_bs_CHG_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHG.bedmethyl.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_bs_CHG_stats) <- stat_colnames
arab_bs_CHG_stats <- arab_bs_CHG_stats[arab_bs_CHG_stats$repeat_len>=100, ]
arab_bs_diffid_CHG <- arab_bs_CHG_stats[arab_bs_CHG_stats$ratio_diff2all>=0.1,]$rpair_name

arab_bs_CHH_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHH.bedmethyl.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_bs_CHH_stats) <- stat_colnames
arab_bs_CHH_stats <- arab_bs_CHH_stats[arab_bs_CHH_stats$repeat_len>=100, ]
arab_bs_diffid_CHH <- arab_bs_CHH_stats[arab_bs_CHH_stats$ratio_diff2all>=0.1,]$rpair_name

arab_bs_C_stats <- read.table("methyrep/data/GCF_000001735.4_TAIR10.1_genomic.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.C.bedmethyl.all2_repeat_group_stats.txt", 
                               header = F, sep = "\t", stringsAsFactors = F)
colnames(arab_bs_C_stats) <- stat_colnames
arab_bs_C_stats <- arab_bs_C_stats[arab_bs_C_stats$repeat_len>=100, ]
arab_bs_diffid_C <- arab_bs_C_stats[arab_bs_C_stats$ratio_diff2all>=0.1,]$rpair_name


# rice 2-1 ===================================================================
rice_bs_CG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CG.bs_rmet_2-1.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_bs_CG_stats) <- stat_colnames
rice_bs_CG_stats <- rice_bs_CG_stats[rice_bs_CG_stats$repeat_len>=100, ]
rice_bs_diffid_CG <- rice_bs_CG_stats[rice_bs_CG_stats$ratio_diff2all>=0.1,]$rpair_name

rice_bs_CHG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHG.bs_rmet_2-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_bs_CHG_stats) <- stat_colnames
rice_bs_CHG_stats <- rice_bs_CHG_stats[rice_bs_CHG_stats$repeat_len>=100, ]
rice_bs_diffid_CHG <- rice_bs_CHG_stats[rice_bs_CHG_stats$ratio_diff2all>=0.1,]$rpair_name

rice_bs_CHH_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHH.bs_rmet_2-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_bs_CHH_stats) <- stat_colnames
rice_bs_CHH_stats <- rice_bs_CHH_stats[rice_bs_CHH_stats$repeat_len>=100, ]
rice_bs_diffid_CHH <- rice_bs_CHH_stats[rice_bs_CHH_stats$ratio_diff2all>=0.1,]$rpair_name

rice_bs_C_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.C.bs_rmet_2-1.all2_repeat_group_stats.txt", 
                               header = F, sep = "\t", stringsAsFactors = F)
colnames(rice_bs_C_stats) <- stat_colnames
rice_bs_C_stats <- rice_bs_C_stats[rice_bs_C_stats$repeat_len>=100, ]
rice_bs_diffid_C <- rice_bs_C_stats[rice_bs_C_stats$ratio_diff2all>=0.1,]$rpair_name

# rice 1-1===================================================================
rice2_bs_CG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CG.bs_rmet_1-1.all2_repeat_group_stats.txt", 
                                 header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_bs_CG_stats) <- stat_colnames
rice2_bs_CG_stats <- rice2_bs_CG_stats[rice2_bs_CG_stats$repeat_len>=100, ]
rice2_bs_diffid_CG <- rice2_bs_CG_stats[rice2_bs_CG_stats$ratio_diff2all>=0.1,]$rpair_name

rice2_bs_CHG_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHG.bs_rmet_1-1.all2_repeat_group_stats.txt", 
                                  header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_bs_CHG_stats) <- stat_colnames
rice2_bs_CHG_stats <- rice2_bs_CHG_stats[rice2_bs_CHG_stats$repeat_len>=100, ]
rice2_bs_diffid_CHG <- rice2_bs_CHG_stats[rice2_bs_CHG_stats$ratio_diff2all>=0.1,]$rpair_name

rice2_bs_CHH_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.CHH.bs_rmet_1-1.all2_repeat_group_stats.txt", 
                                  header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_bs_CHH_stats) <- stat_colnames
rice2_bs_CHH_stats <- rice2_bs_CHH_stats[rice2_bs_CHH_stats$repeat_len>=100, ]
rice2_bs_diffid_CHH <- rice2_bs_CHH_stats[rice2_bs_CHH_stats$ratio_diff2all>=0.1,]$rpair_name

rice2_bs_C_stats <- read.table("methyrep/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.rmN_cut.overlap_mummer.chr_offset2absloc.coords.100_0.99.main_contig.reformat.group.sites.C.bs_rmet_1-1.all2_repeat_group_stats.txt", 
                                header = F, sep = "\t", stringsAsFactors = F)
colnames(rice2_bs_C_stats) <- stat_colnames
rice2_bs_C_stats <- rice2_bs_C_stats[rice2_bs_C_stats$repeat_len>=100, ]
rice2_bs_diffid_C <- rice2_bs_C_stats[rice2_bs_C_stats$ratio_diff2all>=0.1,]$rpair_name







# plots ====
arab_diffid_C_bsvdp2 <- make_two_sets_venn_data(arab_bs_diffid_C, arab_dp2_diffid_C)
arab_diffid_C_bsvdp2_venn <- euler(arab_diffid_C_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_arab_diffid_C_bsvdp2_venn <- plot(arab_diffid_C_bsvdp2_venn,
                                     quantities = list(fontsize=18, fontfamily="Arial",
                                                       labels=trimws(format(unname(arab_diffid_C_bsvdp2), 
                                                                            big.mark = ",", 
                                                                            scientific = F))),
                                     fills = list(fill=palette_bsdp2, alpha=0.7),
                                     edges =TRUE,
                                     labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                     main = list(label="C", 
                                                 fontsize=18,
                                                 font=1,
                                                 vjust=0.8,
                                                 fontfamily="Arial"))
p_arab_diffid_C_bsvdp2_venn

arab_diffid_CG_bsvdp2 <- make_two_sets_venn_data(arab_bs_diffid_CG, arab_dp2_diffid_CG)
arab_diffid_CG_bsvdp2_venn <- euler(arab_diffid_CG_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_arab_diffid_CG_bsvdp2_venn <- plot(arab_diffid_CG_bsvdp2_venn,
                                    quantities = list(fontsize=18, fontfamily="Arial",
                                                      labels=trimws(format(unname(arab_diffid_CG_bsvdp2), 
                                                                           big.mark = ",", 
                                                                           scientific = F))),
                                    fills = list(fill=palette_bsdp2, alpha=0.7),
                                    edges =TRUE,
                                    labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                    main = list(label="CpG", 
                                                fontsize=18,
                                                font=1,
                                                vjust=0.8,
                                                fontfamily="Arial"))
p_arab_diffid_CG_bsvdp2_venn

arab_diffid_CHG_bsvdp2 <- make_two_sets_venn_data(arab_bs_diffid_CHG, arab_dp2_diffid_CHG)
arab_diffid_CHG_bsvdp2_venn <- euler(arab_diffid_CHG_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_arab_diffid_CHG_bsvdp2_venn <- plot(arab_diffid_CHG_bsvdp2_venn,
                                     quantities = list(fontsize=18, fontfamily="Arial",
                                                       labels=trimws(format(unname(arab_diffid_CHG_bsvdp2), 
                                                                            big.mark = ",", 
                                                                            scientific = F))),
                                     fills = list(fill=palette_bsdp2, alpha=0.7),
                                     edges =TRUE,
                                     labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                     main = list(label="CHG", 
                                                 fontsize=18,
                                                 font=1,
                                                 vjust=0.8,
                                                 fontfamily="Arial"))
p_arab_diffid_CHG_bsvdp2_venn

arab_diffid_CHH_bsvdp2 <- make_two_sets_venn_data(arab_bs_diffid_CHH, arab_dp2_diffid_CHH)
arab_diffid_CHH_bsvdp2_venn <- euler(arab_diffid_CHH_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_arab_diffid_CHH_bsvdp2_venn <- plot(arab_diffid_CHH_bsvdp2_venn,
                                      quantities = list(fontsize=18, fontfamily="Arial",
                                                        labels=trimws(format(unname(arab_diffid_CHH_bsvdp2), 
                                                                             big.mark = ",", 
                                                                             scientific = F))),
                                      fills = list(fill=palette_bsdp2, alpha=0.7),
                                      edges =TRUE,
                                      labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                      main = list(label="CHH", 
                                                  fontsize=18,
                                                  font=1,
                                                  vjust=0.8,
                                                  fontfamily="Arial"))
p_arab_diffid_CHH_bsvdp2_venn

# rice 2-1 
rice_diffid_C_bsvdp2 <- make_two_sets_venn_data(rice_bs_diffid_C, rice_dp2_diffid_C)
rice_diffid_C_bsvdp2_venn <- euler(rice_diffid_C_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice_diffid_C_bsvdp2_venn <- plot(rice_diffid_C_bsvdp2_venn,
                                    quantities = list(fontsize=18, fontfamily="Arial",
                                                      labels=trimws(format(unname(rice_diffid_C_bsvdp2), 
                                                                           big.mark = ",", 
                                                                           scientific = F))),
                                    fills = list(fill=palette_bsdp2, alpha=0.7),
                                    edges =TRUE,
                                    labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                    main = list(label="C", 
                                                fontsize=18,
                                                font=1,
                                                vjust=0.8,
                                                fontfamily="Arial"))
p_rice_diffid_C_bsvdp2_venn

rice_diffid_CG_bsvdp2 <- make_two_sets_venn_data(rice_bs_diffid_CG, rice_dp2_diffid_CG)
rice_diffid_CG_bsvdp2_venn <- euler(rice_diffid_CG_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice_diffid_CG_bsvdp2_venn <- plot(rice_diffid_CG_bsvdp2_venn,
                                    quantities = list(fontsize=18, fontfamily="Arial",
                                                      labels=trimws(format(unname(rice_diffid_CG_bsvdp2), 
                                                                           big.mark = ",", 
                                                                           scientific = F))),
                                    fills = list(fill=palette_bsdp2, alpha=0.7),
                                    edges =TRUE,
                                    labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                    main = list(label="CpG", 
                                                fontsize=18,
                                                font=1,
                                                vjust=0.8,
                                                fontfamily="Arial"))
p_rice_diffid_CG_bsvdp2_venn

rice_diffid_CHG_bsvdp2 <- make_two_sets_venn_data(rice_bs_diffid_CHG, rice_dp2_diffid_CHG)
rice_diffid_CHG_bsvdp2_venn <- euler(rice_diffid_CHG_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice_diffid_CHG_bsvdp2_venn <- plot(rice_diffid_CHG_bsvdp2_venn,
                                     quantities = list(fontsize=18, fontfamily="Arial",
                                                       labels=trimws(format(unname(rice_diffid_CHG_bsvdp2), 
                                                                            big.mark = ",", 
                                                                            scientific = F))),
                                     fills = list(fill=palette_bsdp2, alpha=0.7),
                                     edges =TRUE,
                                     labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                     main = list(label="CHG", 
                                                 fontsize=18,
                                                 font=1,
                                                 vjust=0.8,
                                                 fontfamily="Arial"))
p_rice_diffid_CHG_bsvdp2_venn

rice_diffid_CHH_bsvdp2 <- make_two_sets_venn_data(rice_bs_diffid_CHH, rice_dp2_diffid_CHH)
rice_diffid_CHH_bsvdp2_venn <- euler(rice_diffid_CHH_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice_diffid_CHH_bsvdp2_venn <- plot(rice_diffid_CHH_bsvdp2_venn,
                                      quantities = list(fontsize=18, fontfamily="Arial",
                                                        labels=trimws(format(unname(rice_diffid_CHH_bsvdp2), 
                                                                             big.mark = ",", 
                                                                             scientific = F))),
                                      fills = list(fill=palette_bsdp2, alpha=0.7),
                                      edges =TRUE,
                                      labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                      main = list(label="CHH", 
                                                  fontsize=18,
                                                  font=1,
                                                  vjust=0.8,
                                                  fontfamily="Arial"))
p_rice_diffid_CHH_bsvdp2_venn


# rice 1-1 
rice2_diffid_C_bsvdp2 <- make_two_sets_venn_data(rice2_bs_diffid_C, rice2_dp2_diffid_C)
rice2_diffid_C_bsvdp2_venn <- euler(rice2_diffid_C_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice2_diffid_C_bsvdp2_venn <- plot(rice2_diffid_C_bsvdp2_venn,
                                    quantities = list(fontsize=18, fontfamily="Arial",
                                                      labels=trimws(format(unname(rice2_diffid_C_bsvdp2), 
                                                                           big.mark = ",", 
                                                                           scientific = F))),
                                    fills = list(fill=palette_bsdp2, alpha=0.7),
                                    edges =TRUE,
                                    labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                    main = list(label="C", 
                                                fontsize=18,
                                                font=1,
                                                vjust=0.8,
                                                fontfamily="Arial"))
p_rice2_diffid_C_bsvdp2_venn

rice2_diffid_CG_bsvdp2 <- make_two_sets_venn_data(rice2_bs_diffid_CG, rice2_dp2_diffid_CG)
rice2_diffid_CG_bsvdp2_venn <- euler(rice2_diffid_CG_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice2_diffid_CG_bsvdp2_venn <- plot(rice2_diffid_CG_bsvdp2_venn,
                                     quantities = list(fontsize=18, fontfamily="Arial",
                                                       labels=trimws(format(unname(rice2_diffid_CG_bsvdp2), 
                                                                     big.mark = ",", 
                                                                     scientific = F))),
                                     fills = list(fill=palette_bsdp2, alpha=0.7),
                                     edges =TRUE,
                                     labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                     main = list(label="CpG", 
                                                 fontsize=18,
                                                 font=1,
                                                 vjust=0.8,
                                                 fontfamily="Arial"))
p_rice2_diffid_CG_bsvdp2_venn

rice2_diffid_CHG_bsvdp2 <- make_two_sets_venn_data(rice2_bs_diffid_CHG, rice2_dp2_diffid_CHG)
rice2_diffid_CHG_bsvdp2_venn <- euler(rice2_diffid_CHG_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice2_diffid_CHG_bsvdp2_venn <- plot(rice2_diffid_CHG_bsvdp2_venn,
                                      quantities = list(fontsize=18, fontfamily="Arial",
                                                        labels=trimws(format(unname(rice2_diffid_CHG_bsvdp2), 
                                                                      big.mark = ",", 
                                                                      scientific = F))),
                                      fills = list(fill=palette_bsdp2, alpha=0.7),
                                      edges =TRUE,
                                      labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                      main = list(label="CHG", 
                                                  fontsize=18,
                                                  font=1,
                                                  vjust=0.8,
                                                  fontfamily="Arial"))
p_rice2_diffid_CHG_bsvdp2_venn

rice2_diffid_CHH_bsvdp2 <- make_two_sets_venn_data(rice2_bs_diffid_CHH, rice2_dp2_diffid_CHH)
rice2_diffid_CHH_bsvdp2_venn <- euler(rice2_diffid_CHH_bsvdp2)
palette_bsdp2 <- c('#fc8d62', '#66c2a5')
p_rice2_diffid_CHH_bsvdp2_venn <- plot(rice2_diffid_CHH_bsvdp2_venn,
                                      quantities = list(fontsize=18, fontfamily="Arial",
                                                        labels=trimws(format(unname(rice2_diffid_CHH_bsvdp2), 
                                                                             big.mark = ",", 
                                                                             scientific = F))),
                                      fills = list(fill=palette_bsdp2, alpha=0.7),
                                      edges =TRUE,
                                      labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                      main = list(label="CHH", 
                                                  fontsize=18,
                                                  font=1,
                                                  vjust=0.8,
                                                  fontfamily="Arial"))
p_rice2_diffid_CHH_bsvdp2_venn

# 56/35
ppi=300
png("methyrep/fig_diffmethyregions_info2.motif_venn_bsvsdp2.raw.png", 
    width = 36, 
    height = 22.5, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab_diffid_C_bsvdp2_venn, 
                         p_arab_diffid_CG_bsvdp2_venn, 
                         p_arab_diffid_CHG_bsvdp2_venn, 
                         p_arab_diffid_CHH_bsvdp2_venn,  
                         nrow=1, 
                         widths=c(14, 14, 14, 14)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_diffid_C_bsvdp2_venn,
                         p_rice_diffid_CG_bsvdp2_venn, 
                         p_rice_diffid_CHG_bsvdp2_venn, 
                         p_rice_diffid_CHH_bsvdp2_venn,  
                         nrow=1, 
                         widths=c(14, 14, 14, 14)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice2_diffid_C_bsvdp2_venn,
                         p_rice2_diffid_CG_bsvdp2_venn, 
                         p_rice2_diffid_CHG_bsvdp2_venn, 
                         p_rice2_diffid_CHH_bsvdp2_venn,  
                         nrow=1, 
                         widths=c(14, 14, 14, 14)),
             nrow = 5,
             ncol = 1,
             heights=c(11, 1, 11, 1, 11))
dev.off()
svg("methyrep/fig_diffmethyregions_info2.motif_venn_bsvsdp2.raw.svg", 
    width = 36/2.54, 
    height = 22.5/2.54)
grid.arrange(arrangeGrob(p_arab_diffid_C_bsvdp2_venn, 
                         p_arab_diffid_CG_bsvdp2_venn, 
                         p_arab_diffid_CHG_bsvdp2_venn, 
                         p_arab_diffid_CHH_bsvdp2_venn,  
                         nrow=1, 
                         widths=c(14, 14, 14, 14)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_diffid_C_bsvdp2_venn,
                         p_rice_diffid_CG_bsvdp2_venn, 
                         p_rice_diffid_CHG_bsvdp2_venn, 
                         p_rice_diffid_CHH_bsvdp2_venn,  
                         nrow=1, 
                         widths=c(14, 14, 14, 14)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice2_diffid_C_bsvdp2_venn,
                         p_rice2_diffid_CG_bsvdp2_venn, 
                         p_rice2_diffid_CHG_bsvdp2_venn, 
                         p_rice2_diffid_CHH_bsvdp2_venn,  
                         nrow=1, 
                         widths=c(14, 14, 14, 14)),
             nrow = 5,
             ncol = 1,
             heights=c(11, 1, 11, 1, 11))
dev.off()





# fig rice rep1 vs rep2 ==========================================
# function ===
make_two_sets_venn_data <- function(rep1set, rep2set){
  venn_list <- list(sample1=rep1set, sample2=rep2set)
  # plot venn with euler =================
  venn_combs <-
    unlist(lapply(1:length(venn_list),
                  function(j) combn(names(venn_list), j, simplify = FALSE)),
           recursive = FALSE)
  names(venn_combs) <- sapply(venn_combs, function(i) paste0(i, collapse = "&"))
  # str(combs)
  venn_elements <-
    lapply(venn_combs, function(i) Setdiff(venn_list[i], venn_list[setdiff(names(venn_list), i)]))
  venn_n.elements <- sapply(venn_elements, length)
  return(venn_n.elements)
}


rice_dp2_diffCG_rep12 <- make_two_sets_venn_data(rice_dp2_diffid_CG, rice2_dp2_diffid_CG)
rice_diffid_CG_bsvdp2_venn <- euler(rice_dp2_diffCG_rep12)
vennpalette <- c("#78c679", "#006837")
p_rice_diffid_CG_bsvdp2_venn <- plot(rice_diffid_CG_bsvdp2_venn,
                                    quantities = list(fontsize=18, fontfamily="Arial",
                                                      labels=trimws(format(unname(rice_dp2_diffCG_rep12), 
                                                                           big.mark = ",", 
                                                                           scientific = F))),
                                    fills = list(fill=vennpalette, alpha=0.7),
                                    edges =TRUE,
                                    labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                    main = list(label="CpG", 
                                                fontsize=18,
                                                font=1,
                                                vjust=0.8,
                                                fontfamily="Arial"))
p_rice_diffid_CG_bsvdp2_venn

rice_dp2_diffCHG_rep12 <- make_two_sets_venn_data(rice_dp2_diffid_CHG, rice2_dp2_diffid_CHG)
rice_diffid_CHG_bsvdp2_venn <- euler(rice_dp2_diffCHG_rep12)
vennpalette <- c("#78c679", "#006837")
p_rice_diffid_CHG_bsvdp2_venn <- plot(rice_diffid_CHG_bsvdp2_venn,
                                     quantities = list(fontsize=18, fontfamily="Arial",
                                                       labels=trimws(format(unname(rice_dp2_diffCHG_rep12), 
                                                                            big.mark = ",", 
                                                                            scientific = F))),
                                     fills = list(fill=vennpalette, alpha=0.7),
                                     edges =TRUE,
                                     labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                     main = list(label="CHG", 
                                                 fontsize=18,
                                                 font=1,
                                                 vjust=0.8,
                                                 fontfamily="Arial"))
p_rice_diffid_CHG_bsvdp2_venn

rice_dp2_diffCHH_rep12 <- make_two_sets_venn_data(rice_dp2_diffid_CHH, rice2_dp2_diffid_CHH)
rice_diffid_CHH_bsvdp2_venn <- euler(rice_dp2_diffCHH_rep12)
vennpalette <- c("#78c679", "#006837")
p_rice_diffid_CHH_bsvdp2_venn <- plot(rice_diffid_CHH_bsvdp2_venn,
                                      quantities = list(fontsize=18, fontfamily="Arial",
                                                        labels=trimws(format(unname(rice_dp2_diffCHH_rep12), 
                                                                             big.mark = ",", 
                                                                             scientific = F))),
                                      fills = list(fill=vennpalette, alpha=0.7),
                                      edges =TRUE,
                                      labels = list(font=2, fontsize=18, fontfamily="Arial"), 
                                      main = list(label="CHH", 
                                                  fontsize=18,
                                                  font=1,
                                                  vjust=0.8,
                                                  fontfamily="Arial"))
p_rice_diffid_CHH_bsvdp2_venn

ppi=300
# 42/11
png("methyrep/fig_diffmethyregions_info2.motif_venn_ricedp2_rep12.raw.png", 
    width = 36, 
    height = 9.45, units = "cm", res=ppi)
grid.arrange(p_rice_diffid_CG_bsvdp2_venn, 
             p_rice_diffid_CHG_bsvdp2_venn, 
             p_rice_diffid_CHH_bsvdp2_venn,
             nrow = 1,
             ncol = 3,
             widths=c(14, 14, 14))
dev.off()
svg("methyrep/fig_diffmethyregions_info2.motif_venn_ricedp2_rep12.raw.svg", 
    width = 36/2.54, 
    height = 9.45/2.54)
grid.arrange(p_rice_diffid_CG_bsvdp2_venn, 
             p_rice_diffid_CHG_bsvdp2_venn, 
             p_rice_diffid_CHH_bsvdp2_venn,
             nrow = 1,
             ncol = 3,
             widths=c(14, 14, 14))
dev.off()






