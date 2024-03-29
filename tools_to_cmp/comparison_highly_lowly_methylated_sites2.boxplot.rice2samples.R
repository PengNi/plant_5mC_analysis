library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(extrafont)
library(plyr)

# rice2-1, rice1-1

# === cg rice2-1
# df_cg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_cg_o <- read.table("tools_to_cmp/rice2-1.methyfreq_bins.CG.comb_dp2p0.8_meg_CNN.tsv", 
                      header = F, sep = "\t", stringsAsFactors = F)
colnames(df_cg_o) <- c("key", "rmet", "bins", "method")
df_cg_o <- df_cg_o[, c("rmet", "bins", "method")]
df_cg_o$binname <- "low frequency"
df_cg_o[df_cg_o$bins=="mediate",]$binname <- "intermediate frequency"
df_cg_o[df_cg_o$bins=="high",]$binname <- "high frequency"
df_cg_o <- df_cg_o[, c("rmet", "binname", "method")]

cg_n_bs_low <- sum(df_cg_o$binname=="low frequency" & df_cg_o$method=="bisulfite")
cg_n_bs_int <- sum(df_cg_o$binname=="intermediate frequency" & df_cg_o$method=="bisulfite")
cg_n_bs_hig <- sum(df_cg_o$binname=="high frequency" & df_cg_o$method=="bisulfite")
df_cg_o$n_exp <- format(cg_n_bs_low, big.mark = ",", scientific = F)
df_cg_o[df_cg_o$binname=="intermediate frequency",]$n_exp <- format(cg_n_bs_int, 
                                                                    big.mark = ",", scientific = F)
df_cg_o[df_cg_o$binname=="high frequency",]$n_exp <- format(cg_n_bs_hig, big.mark = ",", scientific = F)


df_cg_o$binname <- factor(df_cg_o$binname, levels=c("low frequency", "intermediate frequency", 
                                                    "high frequency"))
df_cg_o$method <- factor(df_cg_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_cg_o <- ggplot(df_cg_o, aes(x=binname, y=rmet, fill=method)) + 
  stat_boxplot(aes(group=interaction(binname, method)), geom = "errorbar", 
               position=position_dodge(0.9)) +
  geom_boxplot(aes(fill=method), 
               position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname + n_exp, scales = "free", 
             labeller = label_bquote(cols = atop(.(as.character(binname)), (italic(n)==.(n_exp))))) +
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
  ggtitle(expression(paste(italic("O. sativa"), " (sample1)    CpG", sep = "")))


# CHG
# df_chg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_chg_o <- read.table("tools_to_cmp/rice2-1.methyfreq_bins.CHG.comb_dp2p0.0_meg_CNN.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chg_o) <- c("key", "rmet", "bins", "method")
df_chg_o <- df_chg_o[, c("rmet", "bins", "method")]
df_chg_o$binname <- "low frequency"
df_chg_o[df_chg_o$bins=="mediate",]$binname <- "intermediate frequency"
df_chg_o[df_chg_o$bins=="high",]$binname <- "high frequency"
df_chg_o <- df_chg_o[, c("rmet", "binname", "method")]

chg_n_bs_low <- sum(df_chg_o$binname=="low frequency" & df_chg_o$method=="bisulfite")
chg_n_bs_int <- sum(df_chg_o$binname=="intermediate frequency" & df_chg_o$method=="bisulfite")
chg_n_bs_hig <- sum(df_chg_o$binname=="high frequency" & df_chg_o$method=="bisulfite")
df_chg_o$n_exp <- format(chg_n_bs_low, big.mark = ",", scientific = F)
df_chg_o[df_chg_o$binname=="intermediate frequency",]$n_exp <- format(chg_n_bs_int, 
                                                                    big.mark = ",", scientific = F)
df_chg_o[df_chg_o$binname=="high frequency",]$n_exp <- format(chg_n_bs_hig, big.mark = ",", scientific = F)


df_chg_o$binname <- factor(df_chg_o$binname, levels=c("low frequency", "intermediate frequency", 
                                                      "high frequency"))
df_chg_o$method <- factor(df_chg_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chg_o <- ggplot(df_chg_o, aes(x=binname, y=rmet, fill=method)) + 
  stat_boxplot(aes(group=interaction(binname, method)), geom = "errorbar", 
               position=position_dodge(0.9)) +
  geom_boxplot(aes(fill=method), 
               position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname + n_exp, scales = "free", 
             labeller = label_bquote(cols = atop(.(as.character(binname)), (italic(n)==.(n_exp))))) +
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
  ggtitle(expression(paste(italic("O. sativa"), " (sample1)    CHG", sep = "")))


# CHH
# df_chh_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHH.comb.tsv", 
#                        header = F, sep = "\t", stringsAsFactors = F)
df_chh_o <- read.table("tools_to_cmp/rice2-1.methyfreq_bins.CHH.comb_dp2p0.8_meg_CNN.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chh_o) <- c("key", "rmet", "bins", "method")
df_chh_o <- df_chh_o[, c("rmet", "bins", "method")]
df_chh_o$binname <- "low frequency"
df_chh_o[df_chh_o$bins=="mediate",]$binname <- "intermediate frequency"
df_chh_o[df_chh_o$bins=="high",]$binname <- "high frequency"
df_chh_o <- df_chh_o[, c("rmet", "binname", "method")]

chh_n_bs_low <- sum(df_chh_o$binname=="low frequency" & df_chh_o$method=="bisulfite")
chh_n_bs_int <- sum(df_chh_o$binname=="intermediate frequency" & df_chh_o$method=="bisulfite")
chh_n_bs_hig <- sum(df_chh_o$binname=="high frequency" & df_chh_o$method=="bisulfite")
df_chh_o$n_exp <- format(chh_n_bs_low, big.mark = ",", scientific = F)
df_chh_o[df_chh_o$binname=="intermediate frequency",]$n_exp <- format(chh_n_bs_int, 
                                                                    big.mark = ",", scientific = F)
df_chh_o[df_chh_o$binname=="high frequency",]$n_exp <- format(chh_n_bs_hig, big.mark = ",", scientific = F)


df_chh_o$binname <- factor(df_chh_o$binname, levels=c("low frequency", "intermediate frequency", 
                                                      "high frequency"))
df_chh_o$method <- factor(df_chh_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chh_o <- ggplot(df_chh_o, aes(x=binname, y=rmet, fill=method)) + 
  stat_boxplot(aes(group=interaction(binname, method)), geom = "errorbar", 
               position=position_dodge(0.9)) +
  geom_boxplot(aes(fill=method), 
               position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname + n_exp, scales = "free", 
             labeller = label_bquote(cols = atop(.(as.character(binname)), (italic(n)==.(n_exp))))) +
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
  ggtitle(expression(paste(italic("O. sativa"), " (sample1)    CHH", sep = "")))

# === cg rice1-1
# df_cg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_cg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CG.comb_dp2p0.8_meg_CNN.tsv", 
                      header = F, sep = "\t", stringsAsFactors = F)
colnames(df_cg_o) <- c("key", "rmet", "bins", "method")
df_cg_o <- df_cg_o[, c("rmet", "bins", "method")]
df_cg_o$binname <- "low frequency"
df_cg_o[df_cg_o$bins=="mediate",]$binname <- "intermediate frequency"
df_cg_o[df_cg_o$bins=="high",]$binname <- "high frequency"
df_cg_o <- df_cg_o[, c("rmet", "binname", "method")]

cg_n_bs_low <- sum(df_cg_o$binname=="low frequency" & df_cg_o$method=="bisulfite")
cg_n_bs_int <- sum(df_cg_o$binname=="intermediate frequency" & df_cg_o$method=="bisulfite")
cg_n_bs_hig <- sum(df_cg_o$binname=="high frequency" & df_cg_o$method=="bisulfite")
df_cg_o$n_exp <- format(cg_n_bs_low, big.mark = ",", scientific = F)
df_cg_o[df_cg_o$binname=="intermediate frequency",]$n_exp <- format(cg_n_bs_int, 
                                                                    big.mark = ",", scientific = F)
df_cg_o[df_cg_o$binname=="high frequency",]$n_exp <- format(cg_n_bs_hig, big.mark = ",", scientific = F)


df_cg_o$binname <- factor(df_cg_o$binname, levels=c("low frequency", "intermediate frequency", 
                                                    "high frequency"))
df_cg_o$method <- factor(df_cg_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_cg_o2 <- ggplot(df_cg_o, aes(x=binname, y=rmet, fill=method)) + 
  stat_boxplot(aes(group=interaction(binname, method)), geom = "errorbar", 
               position=position_dodge(0.9)) +
  geom_boxplot(aes(fill=method), 
               position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname + n_exp, scales = "free", 
             labeller = label_bquote(cols = atop(.(as.character(binname)), (italic(n)==.(n_exp))))) +
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
  ggtitle(expression(paste(italic("O. sativa"), " (sample2)    CpG", sep = "")))


# CHG
# df_chg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHG.comb.tsv", 
#                       header = F, sep = "\t", stringsAsFactors = F)
df_chg_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHG.comb_dp2p0.0_meg_CNN.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chg_o) <- c("key", "rmet", "bins", "method")
df_chg_o <- df_chg_o[, c("rmet", "bins", "method")]
df_chg_o$binname <- "low frequency"
df_chg_o[df_chg_o$bins=="mediate",]$binname <- "intermediate frequency"
df_chg_o[df_chg_o$bins=="high",]$binname <- "high frequency"
df_chg_o <- df_chg_o[, c("rmet", "binname", "method")]

chg_n_bs_low <- sum(df_chg_o$binname=="low frequency" & df_chg_o$method=="bisulfite")
chg_n_bs_int <- sum(df_chg_o$binname=="intermediate frequency" & df_chg_o$method=="bisulfite")
chg_n_bs_hig <- sum(df_chg_o$binname=="high frequency" & df_chg_o$method=="bisulfite")
df_chg_o$n_exp <- format(chg_n_bs_low, big.mark = ",", scientific = F)
df_chg_o[df_chg_o$binname=="intermediate frequency",]$n_exp <- format(chg_n_bs_int, 
                                                                      big.mark = ",", scientific = F)
df_chg_o[df_chg_o$binname=="high frequency",]$n_exp <- format(chg_n_bs_hig, big.mark = ",", scientific = F)


df_chg_o$binname <- factor(df_chg_o$binname, levels=c("low frequency", "intermediate frequency", 
                                                      "high frequency"))
df_chg_o$method <- factor(df_chg_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chg_o2 <- ggplot(df_chg_o, aes(x=binname, y=rmet, fill=method)) + 
  stat_boxplot(aes(group=interaction(binname, method)), geom = "errorbar", 
               position=position_dodge(0.9)) +
  geom_boxplot(aes(fill=method), 
               position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname + n_exp, scales = "free", 
             labeller = label_bquote(cols = atop(.(as.character(binname)), (italic(n)==.(n_exp))))) +
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
  ggtitle(expression(paste(italic("O. sativa"), " (sample2)    CHG", sep = "")))


# CHH
# df_chh_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHH.comb.tsv", 
#                        header = F, sep = "\t", stringsAsFactors = F)
df_chh_o <- read.table("tools_to_cmp/rice1-1.methyfreq_bins.CHH.comb_dp2p0.8_meg_CNN.tsv", 
                       header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chh_o) <- c("key", "rmet", "bins", "method")
df_chh_o <- df_chh_o[, c("rmet", "bins", "method")]
df_chh_o$binname <- "low frequency"
df_chh_o[df_chh_o$bins=="mediate",]$binname <- "intermediate frequency"
df_chh_o[df_chh_o$bins=="high",]$binname <- "high frequency"
df_chh_o <- df_chh_o[, c("rmet", "binname", "method")]

chh_n_bs_low <- sum(df_chh_o$binname=="low frequency" & df_chh_o$method=="bisulfite")
chh_n_bs_int <- sum(df_chh_o$binname=="intermediate frequency" & df_chh_o$method=="bisulfite")
chh_n_bs_hig <- sum(df_chh_o$binname=="high frequency" & df_chh_o$method=="bisulfite")
df_chh_o$n_exp <- format(chh_n_bs_low, big.mark = ",", scientific = F)
df_chh_o[df_chh_o$binname=="intermediate frequency",]$n_exp <- format(chh_n_bs_int, 
                                                                      big.mark = ",", scientific = F)
df_chh_o[df_chh_o$binname=="high frequency",]$n_exp <- format(chh_n_bs_hig, big.mark = ",", scientific = F)


df_chh_o$binname <- factor(df_chh_o$binname, levels=c("low frequency", "intermediate frequency", 
                                                      "high frequency"))
df_chh_o$method <- factor(df_chh_o$method, levels=c("bisulfite", "deepsignal2", "megalodon"))

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chh_o2 <- ggplot(df_chh_o, aes(x=binname, y=rmet, fill=method)) + 
  stat_boxplot(aes(group=interaction(binname, method)), geom = "errorbar", 
               position=position_dodge(0.9)) +
  geom_boxplot(aes(fill=method), 
               position=position_dodge(0.9),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  facet_grid(. ~ binname + n_exp, scales = "free", 
             labeller = label_bquote(cols = atop(.(as.character(binname)), (italic(n)==.(n_exp))))) +
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
  ggtitle(expression(paste(italic("O. sativa"), " (sample2)    CHH", sep = "")))




ppi= 300
png("tools_to_cmp/comparison_highly_lowly_methylated_sites2.boxplot.rice2samples.raw2.png", 
    width = 39, 
    height = 39, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_cg_o, grid.rect(gp=gpar(col="white")), p_cg_o2,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chg_o, grid.rect(gp=gpar(col="white")), p_chg_o2,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh_o, grid.rect(gp=gpar(col="white")), p_chh_o2,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12, 1, 12))
dev.off()

svg("tools_to_cmp/comparison_highly_lowly_methylated_sites2.boxplot.rice2samples.raw2.svg", 
    width = 39/2.54, 
    height = 39/2.54)
grid.arrange(arrangeGrob(p_cg_o, grid.rect(gp=gpar(col="white")), p_cg_o2,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chg_o, grid.rect(gp=gpar(col="white")), p_chg_o2,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh_o, grid.rect(gp=gpar(col="white")), p_chh_o2,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12, 1, 12))
dev.off()





