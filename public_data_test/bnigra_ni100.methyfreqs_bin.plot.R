library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(extrafont)
library(plyr)

# 254.6 bnigra
# === cg
df_cg <- read.table("public_data_test/bnigra_ni100.methyfreq_bins.CG.comb_dp2p0.8_meg_CNN.tsv", 
                 header = F, sep = "\t", stringsAsFactors = F)
colnames(df_cg) <- c("key", "rmet", "bins", "method")
df_cg <- df_cg[, c("rmet", "bins", "method")]
df_cg$binname <- "low frequency"
df_cg[df_cg$bins=="mediate",]$binname <- "intermediate frequency"
df_cg[df_cg$bins=="high",]$binname <- "high frequency"
df_cg <- df_cg[, c("rmet", "binname", "method")]


df_cg$binname <- factor(df_cg$binname, levels=c("low frequency", "intermediate frequency", 
                                                "high frequency"))
df_cg$method <- factor(df_cg$method, levels=c("bisulfite",  
                                              "deepsignal2.arab", 
                                              "deepsignal2.rice", "deepsignal2", 
                                              "megalodon"))

df_cg <- df_cg[df_cg$method %in% c("bisulfite",  
                                   "deepsignal2", 
                                   "megalodon"), ]

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
# cbPalette <- c("#bababa", "#ffffbf", "#abdda4", "#2b83ba", "#fdae61")
# cbPalette <- c("#bababa", "#fdae61", "#ffffbf", "#2b83ba")
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_cg <- ggplot(df_cg, aes(x=binname, y=rmet, fill=method)) + 
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
        # plot.title = element_text(size=33, hjust = -0.05, vjust = -0.01, 
        #                          family="Arial", face="bold")) + 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="Arial")) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) + 
  scale_x_discrete(breaks=c("low frequency", "intermediate frequency", 
                            "high frequency"), 
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("bisulfite",  
                             "deepsignal2", 
                             "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", 
                             "Megalodon")) +
  xlab("Bisulfite methylation frequency") + 
  ylab("Methylation frequency") + 
  # ggtitle("a")
  ggtitle(expression(paste(italic("B. nigra"), "     CpG", sep = "")))

# p_cg



# === chg
df_chg <- read.table("public_data_test/bnigra_ni100.methyfreq_bins.CHG.comb_dp2p0.0_meg_CNN.tsv",
                    header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chg) <- c("key", "rmet", "bins", "method")
df_chg <- df_chg[, c("rmet", "bins", "method")]
df_chg$binname <- "low frequency"
df_chg[df_chg$bins=="mediate",]$binname <- "intermediate frequency"
df_chg[df_chg$bins=="high",]$binname <- "high frequency"
df_chg <- df_chg[, c("rmet", "binname", "method")]


df_chg$binname <- factor(df_chg$binname, levels=c("low frequency", "intermediate frequency",
                                                "high frequency"))
df_chg$method <- factor(df_chg$method, levels=c("bisulfite",  
                                                "deepsignal2.arab", 
                                                "deepsignal2.rice", "deepsignal2", 
                                                "megalodon"))

df_chg <- df_chg[df_chg$method %in% c("bisulfite",  
                                   "deepsignal2", 
                                   "megalodon"), ]

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
# cbPalette <- c("#bababa", "#ffffbf", "#abdda4", "#2b83ba", "#fdae61")
# cbPalette <- c("#bababa", "#fdae61", "#ffffbf", "#2b83ba")
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chg <- ggplot(df_chg, aes(x=binname, y=rmet, fill=method)) +
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
        # plot.title = element_text(size=33, hjust = -0.05, vjust = -0.01, 
        #                          family="Arial", face="bold")) + 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="Arial")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(breaks=c("low frequency", "intermediate frequency",
                            "high frequency"),
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette,
                    breaks=c("bisulfite",  
                             "deepsignal2", 
                             "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", 
                             "Megalodon")) +
  labs(x="Bisulfite methylation frequency",
       y="Methylation frequency", 
       title = expression(paste(italic("B. nigra"), "     CHG", sep = "")))

# p_chg




# === chh
df_chh <- read.table("public_data_test/bnigra_ni100.methyfreq_bins.CHH.comb_dp2p0.8_meg_CNN.tsv",
                     header = F, sep = "\t", stringsAsFactors = F)
colnames(df_chh) <- c("key", "rmet", "bins", "method")
df_chh <- df_chh[, c("rmet", "bins", "method")]
df_chh$binname <- "low frequency"
df_chh[df_chh$bins=="mediate",]$binname <- "intermediate frequency"
df_chh[df_chh$bins=="high",]$binname <- "high frequency"
df_chh <- df_chh[, c("rmet", "binname", "method")]


df_chh$binname <- factor(df_chh$binname, levels=c("low frequency", "intermediate frequency",
                                                  "high frequency"))
df_chh$method <- factor(df_chh$method, levels=c("bisulfite",  
                                                "deepsignal2.arab", 
                                                "deepsignal2.rice", "deepsignal2", 
                                                "megalodon"))

df_chh <- df_chh[df_chh$method %in% c("bisulfite",  
                                      "deepsignal2", 
                                      "megalodon"), ]

# ggplot(df_high, aes(x=method, y=rmet, fill=method)) + geom_violin()
# cbPalette <- c("#bababa", "#ffffbf", "#abdda4", "#2b83ba", "#fdae61")
# cbPalette <- c("#bababa", "#fdae61", "#ffffbf", "#2b83ba")
cbPalette <- c("#bababa", "#2b83ba", "#fdae61")
p_chh <- ggplot(df_chh, aes(x=binname, y=rmet, fill=method)) +
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
        # plot.title = element_text(size=33, hjust = -0.05, vjust = -0.01, 
        #                          family="Arial", face="bold")) + 
        plot.title = element_text(size=17, hjust = 0.5, vjust = 0, 
                                  family="Arial")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(breaks=c("low frequency", "intermediate frequency",
                            "high frequency"),
                   labels=c("0.0-0.3", "0.3-0.7", "0.7-1.0")) +
  scale_fill_manual(values=cbPalette,
                    breaks=c("bisulfite", "deepsignal2", 
                             "megalodon"), 
                    labels=c("bisulfite   ", "DeepSignal-plant   ", 
                             "Megalodon")) +
  labs(x="Bisulfite methylation frequency",
       y="Methylation frequency", 
       title = expression(paste(italic("B. nigra"), "     CHH", sep = "")))

# p_chh



ppi= 300
png("public_data_test/bnigra_ni100.methyfreqs_bin.tools_cmp.raw.png", 
    width = 39, 
    height = 24, units = "cm", res=ppi) # 37, 38
# grid.arrange(arrangeGrob(p_cg,
#                          grid.rect(gp=gpar(col="white")),
#                          p_chg,
#                          grid.rect(gp=gpar(col="white")),
#                          p_chh,
#                          nrow=5,
#                          ncol=1,
#                          heights=c(12, 1, 12, 1, 12)))
grid.arrange(arrangeGrob(p_cg, grid.rect(gp=gpar(col="white")), p_chg,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh, grid.rect(gp=gpar(col="white")), 
                         grid.rect(gp=gpar(col="white")),
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12))
dev.off()

svg("public_data_test/bnigra_ni100.methyfreqs_bin.tools_cmp.raw.svg",
    width = 39/2.54,
    height = 24/2.54)
# grid.arrange(arrangeGrob(p_cg,
#                          grid.rect(gp=gpar(col="white")),
#                          p_chg,
#                          grid.rect(gp=gpar(col="white")),
#                          p_chh,
#                          nrow=5,
#                          ncol=1,
#                          heights=c(12, 1, 12, 1, 12)))
grid.arrange(arrangeGrob(p_cg, grid.rect(gp=gpar(col="white")), p_chg,
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_chh, grid.rect(gp=gpar(col="white")), 
                         grid.rect(gp=gpar(col="white")),
                         nrow=1,
                         ncol=3,
                         widths=c(19, 1, 19)),
             heights=c(12, 1, 12))
dev.off()



