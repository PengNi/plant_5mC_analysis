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
  theme(text=element_text(size=10, family = "Arial"),
        legend.position = "None", 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=9, family = "Arial"), 
        axis.title.y=element_text(size=7.5, family = "Arial")) + 
  scale_x_continuous(# limits=c(0, 1), 
    breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=3),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_arab







# ==================== intersect venn
arab_dp2_rdiff_pairs <- list(C=arab_dp2_diffid_C, 
                             CpG=arab_dp2_diffid_CG, 
                             CHG=arab_dp2_diffid_CHG,
                             CHH=arab_dp2_diffid_CHH)

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
  theme(text=element_text(size=10, family = "Arial"),
        legend.position = "None", 
        legend.key.size = unit(0.35, "cm"), 
        axis.title.x=element_text(size=9, family = "Arial"), 
        axis.title.y=element_text(size=7.5, family = "Arial")) + 
  scale_x_continuous(# limits=c(0, 1), 
    breaks = seq(0, 1, 0.2)) +
  labs(x="Ratio of differentially methylated cytosines", 
       y="Number of repeat pairs (log10)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
p_diffratio_rice


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




# fig in main text ============================================
ppi= 500
png("methyrep/fig_diffmethyregions_info2.a.raw2.png", 
    width = 7, 
    height = 4.7, units = "cm", res=ppi)
p_diffratio_arab
dev.off()
pdf("methyrep/fig_diffmethyregions_info2.a.raw2.pdf", 
    width = 7/2.54, 
    height = 4.7/2.54)
p_diffratio_arab
dev.off()

png("methyrep/fig_diffmethyregions_info2.c.raw2.png", 
    width = 7, 
    height = 4.7, units = "cm", res=ppi)
p_diffratio_rice
dev.off()
pdf("methyrep/fig_diffmethyregions_info2.c.raw2.pdf", 
    width = 7/2.54, 
    height = 4.7/2.54)
p_diffratio_rice
dev.off()

png("methyrep/fig_diffmethyregions_info2.b.raw2.png", 
    width = 14.5, 
    height = 15, units = "cm", res=ppi)
p_arab_dp2_diffid_venn
dev.off()
pdf("methyrep/fig_diffmethyregions_info2.b.raw2.pdf", 
    width = 14.5/2.54, 
    height = 15/2.54)
p_arab_dp2_diffid_venn
dev.off()

png("methyrep/fig_diffmethyregions_info2.d.raw2.png", 
    width = 14.5, 
    height = 15, units = "cm", res=ppi)
p_rice_dp2_diffid_venn
dev.off()
pdf("methyrep/fig_diffmethyregions_info2.d.raw2.pdf", 
    width = 14.5/2.54, 
    height = 15/2.54)
p_rice_dp2_diffid_venn
dev.off()

png("methyrep/fig_diffmethyregions_info2.e.raw2.png", 
    width = 7, 
    height = 6, units = "cm", res=ppi)
p_rice_dp2_diffC_rep12_venn
dev.off()
pdf("methyrep/fig_diffmethyregions_info2.e.raw2.pdf", 
    width = 7/2.54, 
    height = 6/2.54)
p_rice_dp2_diffC_rep12_venn
dev.off()








