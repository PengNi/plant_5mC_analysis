library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

# arab, rice2-1

# === cg
# df_cg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CG.comb.tsv", 
#                     header = F, sep = "\t", stringsAsFactors = F)
df_cg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CG.comb.dp2p0.8_megbscov15.tsv", 
                      header = F, sep = "\t", stringsAsFactors = F)
colnames(df_cg_a) <- c("key", "rmet", "bins", "method")
df_cg_a <- df_cg_a[, c("rmet", "bins", "method")]
df_cg_a$binname <- "low methylation"
df_cg_a[df_cg_a$bins=="mediate",]$binname <- "intermediate methylation"
df_cg_a[df_cg_a$bins=="high",]$binname <- "high methylation"
df_cg_a <- df_cg_a[, c("rmet", "binname", "method")]


df_cg_a$binname <- factor(df_cg_a$binname, levels=c("low methylation", "intermediate methylation", 
                                                "high methylation"))
df_cg_a$method <- factor(df_cg_a$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_cg_a <- ggplot(df_cg_a, aes(x=binname, y=rmet, fill=method)) + 
  geom_boxplot(position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=17, family = "serif"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="serif"), 
        axis.title = element_text(size = 17, family = "serif"), 
        axis.text=element_text(size=17, family = "serif"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="serif", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low methylation", "intermediate methylation", 
                            "high methylation"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("A. thaliana"), "     CG", sep = "")))


# CHG
# df_chg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_chg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHG.comb.dp2p0.8_megbscov15.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chg_a) <- c("key", "rmet", "bins", "method")
df_chg_a <- df_chg_a[, c("rmet", "bins", "method")]
df_chg_a$binname <- "low methylation"
df_chg_a[df_chg_a$bins=="mediate",]$binname <- "intermediate methylation"
df_chg_a[df_chg_a$bins=="high",]$binname <- "high methylation"
df_chg_a <- df_chg_a[, c("rmet", "binname", "method")]


df_chg_a$binname <- factor(df_chg_a$binname, levels=c("low methylation", "intermediate methylation", 
                                                    "high methylation"))
df_chg_a$method <- factor(df_chg_a$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chg_a <- ggplot(df_chg_a, aes(x=binname, y=rmet, fill=method)) + 
  geom_boxplot(position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=17, family = "serif"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="serif"), 
        axis.title = element_text(size = 17, family = "serif"), 
        axis.text=element_text(size=17, family = "serif"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="serif", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low methylation", "intermediate methylation", 
                            "high methylation"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("A. thaliana"), "     CHG", sep = "")))


# CHH
# df_chh_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHH.comb.tsv", 
#                        header = F, sep = "\t", stringsAsFactors = F)
df_chh_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHH.comb.dp2p0.8_megbscov15.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chh_a) <- c("key", "rmet", "bins", "method")
df_chh_a <- df_chh_a[, c("rmet", "bins", "method")]
df_chh_a$binname <- "low methylation"
df_chh_a[df_chh_a$bins=="mediate",]$binname <- "intermediate methylation"
df_chh_a[df_chh_a$bins=="high",]$binname <- "high methylation"
df_chh_a <- df_chh_a[, c("rmet", "binname", "method")]


df_chh_a$binname <- factor(df_chh_a$binname, levels=c("low methylation", "intermediate methylation", 
                                                      "high methylation"))
df_chh_a$method <- factor(df_chh_a$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chh_a <- ggplot(df_chh_a, aes(x=binname, y=rmet, fill=method)) + 
  geom_boxplot(position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=17, family = "serif"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="serif"), 
        axis.title = element_text(size = 17, family = "serif"), 
        axis.text=element_text(size=17, family = "serif"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="serif", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low methylation", "intermediate methylation", 
                            "high methylation"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("A. thaliana"), "     CHH", sep = "")))


# === cg rice
# df_cg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_cg_o <- read.table("tools_to_cmp/rice2-1.methyfreq_bins.CG.comb.dp2p0.8_megbscov15.tsv", 
                      header = F, sep = "\t", stringsAsFactors = F)
colnames(df_cg_o) <- c("key", "rmet", "bins", "method")
df_cg_o <- df_cg_o[, c("rmet", "bins", "method")]
df_cg_o$binname <- "low methylation"
df_cg_o[df_cg_o$bins=="mediate",]$binname <- "intermediate methylation"
df_cg_o[df_cg_o$bins=="high",]$binname <- "high methylation"
df_cg_o <- df_cg_o[, c("rmet", "binname", "method")]


df_cg_o$binname <- factor(df_cg_o$binname, levels=c("low methylation", "intermediate methylation", 
                                                    "high methylation"))
df_cg_o$method <- factor(df_cg_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_cg_o <- ggplot(df_cg_o, aes(x=binname, y=rmet, fill=method)) + 
  geom_boxplot(position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=17, family = "serif"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="serif"), 
        axis.title = element_text(size = 17, family = "serif"), 
        axis.text=element_text(size=17, family = "serif"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="serif", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low methylation", "intermediate methylation", 
                            "high methylation"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("O. sativa"), " (rep1)    CG", sep = "")))


# CHG
# df_chg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_chg_o <- read.table("tools_to_cmp/rice2-1.methyfreq_bins.CHG.comb.dp2p0.8_megbscov15.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chg_o) <- c("key", "rmet", "bins", "method")
df_chg_o <- df_chg_o[, c("rmet", "bins", "method")]
df_chg_o$binname <- "low methylation"
df_chg_o[df_chg_o$bins=="mediate",]$binname <- "intermediate methylation"
df_chg_o[df_chg_o$bins=="high",]$binname <- "high methylation"
df_chg_o <- df_chg_o[, c("rmet", "binname", "method")]


df_chg_o$binname <- factor(df_chg_o$binname, levels=c("low methylation", "intermediate methylation", 
                                                    "high methylation"))
df_chg_o$method <- factor(df_chg_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chg_o <- ggplot(df_chg_o, aes(x=binname, y=rmet, fill=method)) + 
  geom_boxplot(position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=17, family = "serif"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="serif"), 
        axis.title = element_text(size = 17, family = "serif"), 
        axis.text=element_text(size=17, family = "serif"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="serif", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low methylation", "intermediate methylation", 
                            "high methylation"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("O. sativa"), " (rep1)    CHG", sep = "")))


# CHH
# df_chh_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHH.comb.tsv", 
#                        header = F, sep = "\t", stringsAsFactors = F)
df_chh_o <- read.table("tools_to_cmp/rice2-1.methyfreq_bins.CHH.comb.dp2p0.8_megbscov15.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chh_o) <- c("key", "rmet", "bins", "method")
df_chh_o <- df_chh_o[, c("rmet", "bins", "method")]
df_chh_o$binname <- "low methylation"
df_chh_o[df_chh_o$bins=="mediate",]$binname <- "intermediate methylation"
df_chh_o[df_chh_o$bins=="high",]$binname <- "high methylation"
df_chh_o <- df_chh_o[, c("rmet", "binname", "method")]


df_chh_o$binname <- factor(df_chh_o$binname, levels=c("low methylation", "intermediate methylation", 
                                                      "high methylation"))
df_chh_o$method <- factor(df_chh_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chh_o <- ggplot(df_chh_o, aes(x=binname, y=rmet, fill=method)) + 
  geom_boxplot(position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname, scales = "free") + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-2, 0, 0, 0),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=17, family = "serif"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="serif"), 
        axis.title = element_text(size = 17, family = "serif"), 
        axis.text=element_text(size=17, family = "serif"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="serif", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low methylation", "intermediate methylation", 
                            "high methylation"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("O. sativa"), " (rep1)    CHH", sep = "")))




ppi= 300
png("tools_to_cmp/comparison_highly_lowly_methylated_sites2.boxplot.arab_n_rice2-1.raw.png", 
    width = 39, 
    height = 36, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_cg_a, grid.rect(gp=gpar(col="white")), p_cg_o,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chg_a, grid.rect(gp=gpar(col="white")), p_chg_o,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh_a, grid.rect(gp=gpar(col="white")), p_chh_o,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12, 1, 12))
dev.off()

svg("tools_to_cmp/comparison_highly_lowly_methylated_sites2.boxplot.arab_n_rice2-1.raw.svg", 
    width = 39/2.54, 
    height = 36/2.54)
grid.arrange(arrangeGrob(p_cg_a, grid.rect(gp=gpar(col="white")), p_cg_o,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chg_a, grid.rect(gp=gpar(col="white")), p_chg_o,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh_a, grid.rect(gp=gpar(col="white")), p_chh_o,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12, 1, 12))
dev.off()








