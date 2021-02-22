library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(gridExtra)
library(grid)

font_import()
# loadfonts(device = "win")


# cmp other tools ===
# df.cot <- read.table("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_1-1_rep1-2.txt",
#                      header = T,
#                      sep = "\t", stringsAsFactors = F)
df.cot <- read.table("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.txt",
                     header = T,
                     sep = "\t", stringsAsFactors = F)
df.cot[df.cot$method=="DeepSignal2", ]$method = "DeepSignal-plant"
df.cot[df.cot$motif == "CG", ]$motif = "CpG"
df.cot$motif <- factor(df.cot$motif, levels=c("CpG", "CHG", "CHH"))
df.cot <- ddply(df.cot, .(motif, method),
                transform,
                pos = cumsum(Pearson_Correlation))
df.cot$data <- factor(df.cot$data, levels = c("arab.20x", "rice.20x"))
df.cot$method <- factor(df.cot$method, levels = c("Tombo", "DeepSignal",
                                                  "Megalodon",
                                                  "DeepSignal-plant"))

# charts.data <- read.csv("data/copper-data-for-tutorial.csv")
# cbPalettec <- c("#66c2a5", "#fc8d62")
# cbPalettec <- c("#fc8d62", "#66c2a5")
cbPalettec <- c("#1f78b4", "#33a02c")
# p1 <- ggplot() +
#   geom_bar(data = df.cot, aes(y=Pearson_Correlation, x=method, fill = data),
#            colour="black", stat="identity", width = 1) +
#   geom_text(data=df.cot, aes(x = method, y = pos,
#                              label=sprintf("%.4f",
#                                            round(Pearson_Correlation, digits = 4))),
#             size=5,
#             family="Arial", vjust=-0.2) +
#   facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
#   scale_y_continuous(breaks = seq(0, 2, 0.2)) +
#   scale_fill_manual(values=cbPalettec,
#                     breaks=c("arab.20x", "rice.20x"),
#                     labels=c("A. thaliana     ", "O. sativa")) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_text(size=18, face="italic", family = "Arial"),
#         legend.key.size = unit(0.35, "cm"),
#         legend.margin=margin(-5, 0, 0, 0),
#         plot.title = element_text(size=30, hjust = -0.06, face = "bold", family="Arial", 
#                                   vjust = -0.1),
#         strip.background = element_rect(colour="white", fill="white",
#                                         size=1, linetype="solid"),
#         strip.text.x = element_text(size = 18),
#         panel.border = element_blank(),
#         axis.text.x  = element_text(size=12, angle=15, vjust=1, hjust = 1),
#         axis.text.y  = element_text(size=15),
#         text=element_text(size=16,  family="Arial")) +
#   labs(x="", y="Pearson correlation", title = "a")
p1 <- ggplot(data = df.cot, aes(y=Pearson_Correlation, x=method, fill = data, 
                                label=sprintf("%.4f",
                                              round(Pearson_Correlation, digits = 4)))) +
  geom_bar(colour="black", stat="identity", position = position_dodge()) +
  geom_text(size=4, position=position_dodge(width=0.9),
            family="Arial", angle = 90, hjust=1.2) +
  facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values=cbPalettec,
                    breaks=c("arab.20x", "rice.20x"),
                    labels=c(" A. thaliana      ", " O. sativa")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, face="italic", family = "Arial"),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=18, hjust = -0.06, face = "bold", family="Arial", 
                                  vjust = -0.1),
        strip.background = element_rect(colour="white", fill="white",
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(),
        axis.text.x  = element_text(size=12, angle=15, vjust=1, hjust = 1),
        axis.text.y  = element_text(size=12),
        text=element_text(size=12,  family="Arial")) +
  labs(x="", y="Pearson correlation", title = "a")
# p1


# cmp other tools CG_fb_comb ===
# df.cotfb <- read.table("tools_to_cmp/eval.tools_cmp.10x_fb_comb.CG.arab_rep123.rice_1-1_rep1-2.txt",
#                        header = T,
#                        sep = "\t", stringsAsFactors = F)
df.cotfb <- read.table("tools_to_cmp/eval.tools_cmp.10x_fb_comb.CG.arab_rep123.rice_2-1_rep2-1.txt",
                       header = T,
                       sep = "\t", stringsAsFactors = F)
df.cotfb[df.cotfb$method=="DeepSignal2", ]$method = "DeepSignal-plant"
df.cotfb[df.cotfb$motif == "CG", ]$motif = "CpG"
df.cotfb$motif <- factor(df.cotfb$motif, levels=c("CpG", "CHG", "CHH"))
df.cotfb <- ddply(df.cotfb, .(motif, method),
                transform,
                pos = cumsum(Pearson_Correlation))
df.cotfb$data <- factor(df.cotfb$data, levels = c("arab.20x", "rice.20x"))
df.cotfb$method <- factor(df.cotfb$method, levels = c("Tombo", "DeepSignal",
                                                      "Megalodon", "nanopolish",
                                                      "DeepSignal-plant"))

# cbPalettec <- c("#fc8d62", "#66c2a5")
# cbPalettec <- c("#66c2a5", "#fc8d62")
cbPalettec <- c("#1f78b4", "#33a02c")
# p2 <- ggplot() +
#   geom_bar(data = df.cotfb, aes(y=Pearson_Correlation, x=method, fill = data),
#            colour="black", stat="identity", width = 1) +
#   geom_text(data=df.cotfb, aes(x = method, y = pos,
#                                label=sprintf("%.4f",
#                                              round(Pearson_Correlation, digits = 4))),
#             size=5,
#             family="Arial", vjust=-0.2) +
#   facet_grid(. ~ motif) +
#   scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.2)) +
#   scale_fill_manual(values=cbPalettec,
#                     breaks=c("arab.20x", "rice.20x"),
#                     labels=c("A. thaliana  ", "O. sativa")) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_text(size=16, face="italic", family = "Arial"),
#         legend.key.size = unit(0.35, "cm"),
#         legend.margin=margin(-5, 0, 0, 0),
#         plot.title = element_text(size=30, hjust = -0.13, face = "bold", family="Arial", 
#                                   vjust = -0.1),
#         strip.background = element_rect(colour="white", fill="white",
#                                         size=1, linetype="solid"),
#         strip.text.x = element_text(size = 18),
#         panel.border = element_blank(),
#         axis.text.x  = element_text(size=12, angle=15, vjust=1, hjust = 1), 
#         axis.text.y  = element_text(size=15),
#         text=element_text(size=16,  family="Arial")) +
#   labs(x="", y="Pearson correlation", title ="b")
p2 <- ggplot(data = df.cotfb, aes(y=Pearson_Correlation, x=method, fill = data, 
                                label=sprintf("%.4f",
                                              round(Pearson_Correlation, digits = 4)))) +
  geom_bar(colour="black", stat="identity", position = position_dodge()) +
  geom_text(size=4, position=position_dodge(width=0.9),
            family="Arial", angle = 90, hjust=1.2) +
  facet_grid(. ~ motif, scales = "free_x", space = "free_x") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values=cbPalettec,
                    breaks=c("arab.20x", "rice.20x"),
                    labels=c(" A. thaliana      ", " O. sativa")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, face="italic", family = "Arial"),
        legend.key.size = unit(0.35, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=18, hjust = -0.15, face = "bold", family="Arial", 
                                  vjust = -0.1),
        strip.background = element_rect(colour="white", fill="white",
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 12),
        panel.border = element_blank(),
        axis.text.x  = element_text(size=12, angle=15, vjust=1, hjust = 1),
        axis.text.y  = element_text(size=12),
        text=element_text(size=12,  family="Arial")) +
  labs(x="", y="Pearson correlation", title = "b")
# p2


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
p_legend = g_legend(p1)
ppi= 300
# 32, 20
png("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.png",
    width = 36,
    height = 16, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p2+ theme(legend.position = "none"),
                         nrow = 1,
                         ncol = 3,
                         widths=c(20, 1, 10.5)),
             p_legend,
             nrow = 2,
             heights = c(15, 1))
dev.off()

svg("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.svg",
    width = 36/2.54,
    height = 16/2.54)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p2+ theme(legend.position = "none"),
                         nrow = 1,
                         ncol = 3,
                         widths=c(20, 1, 10.5)),
             p_legend,
             nrow = 2,
             heights = c(15, 1))
dev.off()



# cmp megalodon under diff coverage ==
# df.vmeg <- read.table("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice_1-1_rep1-2.txt", 
#                       header = T, 
#                       sep = "\t", stringsAsFactors = F)
df.vmeg <- read.table("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice2-1_rep2-1.dp0.8.txt", 
                      header = T, 
                      sep = "\t", stringsAsFactors = F)
df.vmeg$coverage <- paste(substr(df.vmeg$coverage, 1, nchar(df.vmeg$coverage)-1), 
                          "\U00D7", sep = "")


# facet with species, motif single plot ====================================================
# CG
df.vmeg.cg <- df.vmeg[df.vmeg$motif=="CG", ]
df.vmeg.cg$method <- factor(df.vmeg.cg$method, levels = c("DeepSignal2", "Megalodon"))
# cbPalettec <- c("#abdda4", "#fdae61")
cbPalettec <- c("#2b83ba", "#fdae61")
p_cg <- ggplot(data = df.vmeg.cg, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", stat="identity", position="dodge") + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=5, angle=90, hjust=1.2,
            family="serif") + 
  facet_grid(. ~ species, scales = "free_y") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("DeepSignal2", "Megalodon"), 
                    labels=c("DeepSignal-plant         ", "Megalodon")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15, family = "serif"), 
        # legend.key.size = unit(0.4, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 17, face = "italic"),
        panel.border = element_blank(), 
        text=element_text(size=17,  family="serif"), 
        axis.title = element_text(size = 15, family = "serif"), 
        axis.text=element_text(size=15, family = "serif"), 
        plot.title = element_text(hjust = -0, family="serif")) + 
  labs(x="Coverage", y="Pearson correlation", title = "CG")

# CHG
df.vmeg.chg <- df.vmeg[df.vmeg$motif=="CHG", ]
df.vmeg.chg$method <- factor(df.vmeg.chg$method, levels = c("DeepSignal2", "Megalodon"))
# cbPalettec <- c("#abdda4", "#fdae61")
cbPalettec <- c("#2b83ba", "#fdae61")
p_chg <- ggplot(data = df.vmeg.chg, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", stat="identity", position="dodge") + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=5, angle=90, hjust=1.2,
            family="serif") + 
  facet_grid(. ~ species, scales = "free_y") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("DeepSignal2", "Megalodon"), 
                    labels=c("DeepSignal-plant         ", "Megalodon")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15, family = "serif"), 
        # legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 17, face = "italic"),
        panel.border = element_blank(), 
        text=element_text(size=17,  family="serif"), 
        axis.title = element_text(size = 15, family = "serif"), 
        axis.text=element_text(size=15, family = "serif"), 
        plot.title = element_text(hjust = -0, family="serif")) + 
  labs(x="Coverage", y="Pearson correlation", title = "CHG")


# CHH
df.vmeg.chh <- df.vmeg[df.vmeg$motif=="CHH", ]
df.vmeg.chh$method <- factor(df.vmeg.chh$method, levels = c("DeepSignal2", "Megalodon"))
# cbPalettec <- c("#abdda4", "#fdae61")
cbPalettec <- c("#2b83ba", "#fdae61")
p_chh <- ggplot(data = df.vmeg.chh, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", stat="identity", position="dodge") + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                              round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=5, angle=90, hjust=1.2,
            family="serif") + 
  facet_grid(. ~ species, scales = "free_y") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("DeepSignal2", "Megalodon"), 
                    labels=c("DeepSignal-plant         ", "Megalodon")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=15, family = "serif"), 
        # legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size = 17, face = "italic"),
        panel.border = element_blank(), 
        text=element_text(size=17,  family="serif"), 
        axis.title = element_text(size = 15, family = "serif"), 
        axis.text=element_text(size=15, family = "serif"), 
        plot.title = element_text(hjust = -0, family="serif")) + 
  labs(x="Coverage", y="Pearson correlation", title = "CHH")


# ============
ppi= 300
png("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice_2-1_rep2-1.raw.png", 
    width = 44, 
    height = 13, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_cg, 
                         grid.rect(gp=gpar(col="white")),
                         p_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_chh,
                         nrow=1, 
                         ncol = 5,
                         widths=c(12, 1, 12, 1, 12)))
dev.off()

svg("tools_to_cmp/eval.tools_cmp.10-50x.arab_rep123.rice_2-1_rep2-1.raw.svg", 
    width = 44/2.54, 
    height = 13/2.54)
grid.arrange(arrangeGrob(p_cg, 
                         grid.rect(gp=gpar(col="white")),
                         p_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_chh,
                         nrow=1, 
                         ncol = 5,
                         widths=c(12, 1, 12, 1, 12)))
dev.off()















