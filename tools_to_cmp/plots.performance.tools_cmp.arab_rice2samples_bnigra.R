library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(gridExtra)
library(grid)

# font_import()
# loadfonts(device = "win")


df.vmeg <- read.table("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice2-1_rep2-1.dp0.8.txt", 
                      header = T, 
                      sep = "\t", stringsAsFactors = F)
df.vmeg$coverage <- paste(substr(df.vmeg$coverage, 1, nchar(df.vmeg$coverage)-1), 
                          "\U00D7", sep = "")
df.vmeg$coverage <- factor(df.vmeg$coverage, levels = c("20\U00D7", "40\U00D7", 
                                                        "60\U00D7", "80\U00D7", "100\U00D7"))


# facet with motif, species single plot ========================================================
df.vmeg.arab <- df.vmeg[df.vmeg$species=="A. thaliana", ]
df.vmeg.arab[df.vmeg.arab$motif=="CG", ]$motif = "CpG"
df.vmeg.arab$motif <- factor(df.vmeg.arab$motif, levels = c("CpG", "CHG", "CHH"))
df.vmeg.arab$method <- factor(df.vmeg.arab$method, levels = c("DeepSignal2", "Megalodon"))
# cbPalettec <- c("#abdda4", "#fdae61")
cbPalettec <- c("#2b83ba", "#fdae61")
p_arab <- ggplot(data = df.vmeg.arab, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", size=0.3, stat="identity", position="dodge") + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                              round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=3, angle=90, hjust=1.2,
            family="Arial") + 
  facet_grid(. ~ motif, scales = "free_y") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("DeepSignal2", "Megalodon"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon (retrain)")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        # legend.text = element_text(size=12, family = "Arial"), 
        legend.key.size = unit(0.25, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(), 
        text=element_text(size=12,  family="Arial"), 
        axis.title = element_text(size = 12, family = "Arial"), 
        axis.text.x=element_text(size=9, family = "Arial"), 
        axis.text.y=element_text(size=12, family = "Arial"), 
        plot.title = element_text(hjust = 0.5, size = 15, family="Arial", face = "italic")) + 
  labs(x="Coverage", y="Pearson correlation", title = "A. thaliana")
#p_arab


df.vmeg.rice <- df.vmeg[df.vmeg$species=="O. sativa", ]
df.vmeg.rice[df.vmeg.rice$motif=="CG", ]$motif = "CpG"
df.vmeg.rice$motif <- factor(df.vmeg.rice$motif, levels = c("CpG", "CHG", "CHH"))
df.vmeg.rice$method <- factor(df.vmeg.rice$method, levels = c("DeepSignal2", "Megalodon"))
# cbPalettec <- c("#abdda4", "#fdae61")
cbPalettec <- c("#2b83ba", "#fdae61")
p_rice <- ggplot(data = df.vmeg.rice, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", size=0.3, stat="identity", position="dodge") + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                              round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=3, angle=90, hjust=1.2,
            family="Arial") + 
  facet_grid(. ~ motif, scales = "free_y") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("DeepSignal2", "Megalodon"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon (retrain)")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        # legend.text = element_text(size=12, family = "Arial"), 
        legend.key.size = unit(0.25, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(), 
        text=element_text(size=12,  family="Arial"), 
        axis.title = element_text(size = 12, family = "Arial"), 
        axis.text.x=element_text(size=9, family = "Arial"), 
        axis.text.y=element_text(size=12, family = "Arial"), 
        plot.title = element_text(hjust = 0.5, size = 15, family="Arial", face = "italic")) + 
  labs(x="Coverage", y="Pearson correlation", title = expression(paste(italic("O. sativa"), 
                                                                       " (sample1)", 
                                                                       sep = "")))
#p_rice


# rice samples =======
df.rice2 <- read.table("tools_to_cmp/eval.tools_cmp.10-50x.rice1-1_rep1-2.dp0.8.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)

df.rice2$coverage <- paste(substr(df.rice2$coverage, 1, nchar(df.rice2$coverage)-1), 
                           "\U00D7", sep = "")
df.rice2$coverage <- factor(df.rice2$coverage, levels = c("20\U00D7", "40\U00D7", 
                                                          "60\U00D7", "80\U00D7", "100\U00D7"))
df.rice2$method <- factor(df.rice2$method, levels = c("DeepSignal2", "Megalodon"))
df.rice2[df.rice2$motif=="CG", ]$motif = "CpG"
df.rice2$motif <- factor(df.rice2$motif, levels = c("CpG", "CHG", "CHH"))

cbPalettec <- c("#2b83ba", "#fdae61")
p_rice2 <- ggplot(data = df.rice2, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", size=0.3, stat="identity", position="dodge") + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                              round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=3, angle=90, hjust=1.2,
            family="Arial") + 
  facet_grid(. ~ motif, scales = "free_y") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("DeepSignal2", "Megalodon"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon (retrain)")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        # legend.text = element_text(size=12, family = "Arial"), 
        legend.key.size = unit(0.25, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(), 
        text=element_text(size=12,  family="Arial"), 
        axis.title = element_text(size = 12, family = "Arial"), 
        axis.text.x=element_text(size=9, family = "Arial"), 
        axis.text.y=element_text(size=12, family = "Arial"), 
        plot.title = element_text(hjust = 0.5, size = 15, family="Arial", face = "italic")) + 
  labs(x="Coverage", y="Pearson correlation", title = expression(paste(italic("O. sativa"), 
                                                                       " (sample2)", 
                                                                       sep = "")))
# p_rice2



# bnigra =========
df.bnigra <- read.table("public_data_test/bnigra_ni100.correlation.coverage_effect.cmp_tools.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)

df.bnigra$coverage <- paste(substr(df.bnigra$coverage, 1, 2), "\U00D7", sep = "")
df.bnigra$method <- factor(df.bnigra$method, levels = c("DeepSignal2_arab", "DeepSignal2_rice", 
                                                        "DeepSignal2_comb", "Megalodon"))
df.bnigra[df.bnigra$motif=="CG", ]$motif = "CpG"
df.bnigra$motif <- factor(df.bnigra$motif, levels = c("CpG", "CHG", "CHH"))

df.bnigra <- df.bnigra[df.bnigra$method %in% c("DeepSignal2_comb", "Megalodon"), ]

# cbPalette <- c("#ffffbf", "#abdda4", "#2b83ba", "#fdae61")
cbPalette <- c("#2b83ba", "#fdae61")
p_bnigra <- ggplot(data = df.bnigra, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", size=0.3, stat="identity", position=position_dodge(0.9)) + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                              round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=3, angle=90, hjust=1.2,
            family="Arial") + 
  facet_grid(. ~ motif, scales = "free") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  # scale_fill_manual(values=cbPalette, 
  #                   breaks=c("DeepSignal2_arab", "DeepSignal2_rice", 
  #                            "DeepSignal2_comb", "Megalodon"), 
  #                   labels=c("DeepSignal2_arab   ", "DeepSignal2_rice   ", 
  #                            "DeepSignal2_comb   ", "Megalodon")) +
  scale_fill_manual(values=cbPalette, 
                    breaks=c("DeepSignal2_comb", "Megalodon"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon (retrain)")) +
  theme_bw() + 
  theme(panel.spacing.x = unit(2, "lines"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        # legend.text = element_text(size=12, family = "Arial"), 
        legend.key.size = unit(0.25, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(), 
        text=element_text(size=12,  family="Arial"), 
        axis.title = element_text(size = 12, family = "Arial"), 
        axis.text.x=element_text(size=9, family = "Arial"), 
        axis.text.y=element_text(size=12, family = "Arial"), 
        plot.title = element_text(hjust = 0.5, size = 15, family="Arial", face = "italic")) + 
  labs(x="Coverage", y="Pearson correlation", title = "B. nigra")
#p_bnigra



# ============
ppi= 300
png("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice_2samples.bnigra.dp0.8.raw.png", 
    width = 36, 
    height = 23, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice,
                         nrow=1, 
                         ncol = 3,
                         widths=c(17.5, 1, 17.5)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_bnigra, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice2,
                         nrow=1, 
                         ncol = 3,
                         widths=c(17.5, 1, 17.5)), 
             nrow=3, 
             heights=c(11, 1, 11))
dev.off()

svg("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice_2samples.bnigra.dp0.8.raw.svg", 
    width = 36/2.54, 
    height = 23/2.54)
grid.arrange(arrangeGrob(p_arab, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice,
                         nrow=1, 
                         ncol = 3,
                         widths=c(17.5, 1, 17.5)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(p_bnigra, 
                         grid.rect(gp=gpar(col="white")),
                         p_rice2,
                         nrow=1, 
                         ncol = 3,
                         widths=c(17.5, 1, 17.5)), 
             nrow=3, 
             heights=c(11, 1, 11))
dev.off()


