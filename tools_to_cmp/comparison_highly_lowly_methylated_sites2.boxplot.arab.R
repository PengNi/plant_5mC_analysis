library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(extrafont)
library(plyr)

# arab

# === cg
# df_cg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CG.comb.tsv", 
#                     header = F, sep = "\t", stringsAsFactors = F)
df_cg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CG.comb_dp2p0.8_meg_CNN.tsv", 
                      header = F, sep = "\t", stringsAsFactors = F)
colnames(df_cg_a) <- c("key", "rmet", "bins", "method")
df_cg_a <- df_cg_a[, c("rmet", "bins", "method")]
df_cg_a$binname <- "low frequency"
df_cg_a[df_cg_a$bins=="mediate",]$binname <- "intermediate frequency"
df_cg_a[df_cg_a$bins=="high",]$binname <- "high frequency"
df_cg_a <- df_cg_a[, c("rmet", "binname", "method")]


df_cg_a$binname <- factor(df_cg_a$binname, levels=c("low frequency", "intermediate frequency", 
                                                    "high frequency"))
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
        legend.text = element_text(size=17, family = "Arial"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="Arial"), 
        axis.title = element_text(size = 17, family = "Arial"), 
        axis.text=element_text(size=17, family = "Arial"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low frequency", "intermediate frequency", 
                            "high frequency"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("A. thaliana"), "     CpG", sep = "")))


# CHG
# df_chg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_chg_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHG.comb_dp2p0.8_meg_CNN.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chg_a) <- c("key", "rmet", "bins", "method")
df_chg_a <- df_chg_a[, c("rmet", "bins", "method")]
df_chg_a$binname <- "low frequency"
df_chg_a[df_chg_a$bins=="mediate",]$binname <- "intermediate frequency"
df_chg_a[df_chg_a$bins=="high",]$binname <- "high frequency"
df_chg_a <- df_chg_a[, c("rmet", "binname", "method")]


df_chg_a$binname <- factor(df_chg_a$binname, levels=c("low frequency", "intermediate frequency", 
                                                      "high frequency"))
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
        legend.text = element_text(size=17, family = "Arial"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="Arial"), 
        axis.title = element_text(size = 17, family = "Arial"), 
        axis.text=element_text(size=17, family = "Arial"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low frequency", "intermediate frequency", 
                            "high frequency"), 
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
df_chh_a <- read.table("tools_to_cmp/athaliana.methyfreq_bins.CHH.comb_dp2p0.8_meg_CNN.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chh_a) <- c("key", "rmet", "bins", "method")
df_chh_a <- df_chh_a[, c("rmet", "bins", "method")]
df_chh_a$binname <- "low frequency"
df_chh_a[df_chh_a$bins=="mediate",]$binname <- "intermediate frequency"
df_chh_a[df_chh_a$bins=="high",]$binname <- "high frequency"
df_chh_a <- df_chh_a[, c("rmet", "binname", "method")]


df_chh_a$binname <- factor(df_chh_a$binname, levels=c("low frequency", "intermediate frequency", 
                                                      "high frequency"))
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
        legend.text = element_text(size=17, family = "Arial"),
        strip.text.x = element_text(size = 13),
        text = element_text(size=17, family="Arial"), 
        axis.title = element_text(size = 17, family = "Arial"), 
        axis.text=element_text(size=17, family = "Arial"), 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="Arial", face="italic")) + 
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low frequency", "intermediate frequency", 
                            "high frequency"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite", "deepsignal2", "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  ggtitle(expression(paste(italic("A. thaliana"), "     CHH", sep = "")))







ppi= 300
png("tools_to_cmp/comparison_highly_lowly_methylated_sites2.boxplot.arab.raw.png", 
    width = 39, 
    height = 24, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_cg_a, grid.rect(gp=gpar(col="white")), p_chg_a,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh_a, grid.rect(gp=gpar(col="white")), 
                         grid.rect(gp=gpar(col="white")),
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12))
dev.off()

svg("tools_to_cmp/comparison_highly_lowly_methylated_sites2.boxplot.arab.raw.svg", 
    width = 39/2.54, 
    height = 24/2.54)
grid.arrange(arrangeGrob(p_cg_a, grid.rect(gp=gpar(col="white")), p_chg_a,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh_a, grid.rect(gp=gpar(col="white")), 
                         grid.rect(gp=gpar(col="white")),
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12))
dev.off()








