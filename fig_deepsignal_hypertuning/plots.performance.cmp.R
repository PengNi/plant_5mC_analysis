library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)

font_import()
loadfonts(device = "win")

# 1 inch = 2.54 cm

# diff kmers =
df.klen <- read.table("fig_deepsignal_hypertuning/eval.diff_kmer.arab.10x.vs_rep1.txt", 
                     header = T, 
                     sep = "\t", stringsAsFactors = F)
df.klen$klen <- factor(df.klen$klen, levels = c(9, 11, 13, 15, 17))

cbPalette <- c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177")
labelx_face=c(rep("italic", 1), rep("plain", 11))
p <- ggplot() + 
  geom_bar(data = df.klen, aes(y=Pearson_Correlation, x=klen, fill = klen), 
           stat="identity", position=position_dodge(), colour="black", 
           width=1) + 
  geom_text(data = df.klen, aes(x=klen, y=Pearson_Correlation, 
                                label=sprintf("%.4f", 
                                              round(Pearson_Correlation, digits = 4))), 
            size = 2.1, vjust=-0.5, 
            family="serif") +
  facet_grid(. ~ motif) + 
  theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_rect(colour="white", fill="white", 
                                        size=2.0, linetype="solid"), 
        panel.border = element_blank(), 
        legend.title = element_blank(), 
        text=element_text(size=12,  family="serif")) + 
  scale_fill_manual(values=cbPalette) +
  coord_cartesian(ylim=c(0.6, 1)) +
  scale_y_continuous(breaks = seq(0.6, 1, 0.05)) + 
  labs(x=expression(paste(italic("k"), "-mer length")), y="Pearson correlation")

ppi= 300
png("fig_deepsignal_hypertuning/eval.diff_kmer.arab.10x.vs_rep1.plot.png", 
     width = 12, 
     height = 9, units = "cm", res=ppi)
p
dev.off()

# svg 
svg("fig_deepsignal_hypertuning/eval.diff_kmer.arab.10x.vs_rep1.plot.svg", 
     width = 12/2.54, 
     height = 9/2.54)
p
dev.off()



# denoise selection ===
df.denoise <- read.table("fig_deepsignal_hypertuning/eval.diff_feature.denoise.arab.10x.vs_bsrep123.txt", 
                      header = T, 
                      sep = "\t", stringsAsFactors = F)
df.denoise$feature <- factor(df.denoise$feature, levels = c("signal", "sequence", "signal+sequence"))

cbPalette <- c("#fc8d59", "#ffffbf", "#91cf60")
p1 <- ggplot(data = df.denoise, aes(y=Pearson_Correlation, x=motif, fill = feature, 
                                  label=sprintf("%.4f", 
                                                round(Pearson_Correlation, digits = 4)))) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  geom_text(size = 3.3, position=position_dodge(width=0.9), vjust=-0.5, 
            family="serif") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        text=element_text(size=12,  family="serif"), 
        plot.title = element_text(size=22, hjust = -0.18, face = "bold", 
                                  vjust=-0.1),
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-8, 0, 0, 0)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks=c("signal", "sequence", "signal+sequence"), 
                    labels=c("signal  ", "sequence  ", "signal+sequence")) +
  coord_cartesian(ylim=c(0.65, 0.97)) +
  scale_y_continuous(breaks = seq(0.65, 0.97, 0.05)) + 
  labs(x="", y="Pearson correlation", 
       title = "a")

# call_methyl selection
df.call <- read.table("fig_deepsignal_hypertuning/eval.diff_feature.call_methyl.arab.10x.vs_bsrep123.txt", 
                         header = T, 
                         sep = "\t", stringsAsFactors = F)
df.call$feature <- factor(df.call$feature, levels = c("signal", "sequence", "signal+sequence"))

cbPalette2 <- c("#fc8d59", "#ffffbf", "#91cf60")

p2 <- ggplot(data = df.call, aes(y=Pearson_Correlation, x=motif, fill = feature, 
                                    label=sprintf("%.4f", 
                                                  round(Pearson_Correlation, digits = 4)))) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  geom_text(size = 3.3, position=position_dodge(width=0.9), vjust=-0.5, 
            family="serif") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        text=element_text(size=12, family="serif"), 
        plot.title = element_text(size=22, hjust = -0.13, face = "bold", family="serif",
                                  vjust=-0.1),
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-8, 0, 0, 0)) + 
  scale_fill_manual(values=cbPalette2, 
                    breaks=c("signal", "sequence", "signal+sequence"), 
                    labels=c("signal  ", "sequence  ", "signal+sequence")) +
  coord_cartesian(ylim=c(0.65, 0.97)) +
  scale_y_continuous(breaks = seq(0.65, 0.97, 0.05)) + 
  labs(x="", y="Pearson correlation", 
       title = "b")



ppi= 300
png("fig_deepsignal_hypertuning/eval.diff_feature.arab.10x.vs_bsrep123.plot.png", 
     width = 22, 
     height = 13, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p1,
                         grid.rect(gp=gpar(col="white")),
                         p2, 
                         nrow=1, 
                         widths = c(9, 1, 12)),
             heights=c(13))
dev.off()

# svg
svg("fig_deepsignal_hypertuning/eval.diff_feature.arab.10x.vs_bsrep123.plot.svg", 
     width = 22/2.54, 
     height = 13/2.54)
grid.arrange(arrangeGrob(p1,
                         grid.rect(gp=gpar(col="white")),
                         p2, 
                         nrow=1, 
                         widths = c(9, 1, 12)),
             heights=c(13))
dev.off()







