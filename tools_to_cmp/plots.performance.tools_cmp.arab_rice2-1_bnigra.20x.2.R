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


# cmp megalodon megalodon_retrain ===
df.cot <- read.table("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.bnigra.txt",
                     header = T,
                     sep = "\t", stringsAsFactors = F)
df.cot[df.cot$method=="DeepSignal2", ]$method = "DeepSignal-plant"
df.cot$method = factor(df.cot$method, levels = c("DeepSignal-plant", 
                                                 "Megalodon", 
                                                 "Megalodon (retrain)"))
df.cot[df.cot$motif == "CG", ]$motif = "CpG"
df.cot$motif <- factor(df.cot$motif, levels=c("CpG", "CHG", "CHH"))

# cbPalettec <- c("#66c2a5", "#fc8d62")
# cbPalettec <- c("#fc8d62", "#66c2a5")
# cbPalettec <- c("#1f78b4", "#33a02c")
cbPalettec <- c("#2b83ba", "#abdda4", "#fdae61")


df.cot_arab <- df.cot[df.cot$data == "arab.20x", ]
p_a <- ggplot(data = df.cot_arab, aes(y=Pearson_Correlation, x=motif, fill = method, 
                                      label=sprintf("%.4f", 
                                                    round(Pearson_Correlation, digits = 4)))) +
  geom_bar(colour="black", size=0.3, stat="identity", position = position_dodge()) +
  geom_text(size=2.2, position=position_dodge(width=0.9),
            family="Arial", angle = 90, hjust=1.1) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  scale_fill_manual(values = cbPalettec, 
                    breaks=c("DeepSignal-plant", "Megalodon", "Megalodon (retrain)"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon        ", 
                             " Megalodon (retrain)")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=10, family = "Arial"),
        legend.key.size = unit(0.22, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=12, hjust = 0.5, face = "italic", family="Arial"),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=10),
        text=element_text(size=10,  family="Arial")) +
  labs(x="", y="Pearson correlation", title = "A. thaliana")
# p_a

df.cot_rice <- df.cot[df.cot$data == "rice.20x", ]
p_o <- ggplot(data = df.cot_rice, aes(y=Pearson_Correlation, x=motif, fill = method, 
                                      label=sprintf("%.4f", 
                                                    round(Pearson_Correlation, digits = 4)))) +
  geom_bar(colour="black", size=0.3, stat="identity", position = position_dodge()) +
  geom_text(size=2.2, position=position_dodge(width=0.9),
            family="Arial", angle = 90, hjust=1.1) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  scale_fill_manual(values = cbPalettec, 
                    breaks=c("DeepSignal-plant", "Megalodon", "Megalodon (retrain)"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon        ", 
                             " Megalodon (retrain)")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=10, family = "Arial"),
        legend.key.size = unit(0.22, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=12, hjust = 0.5, face = "italic", family="Arial"),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=10),
        text=element_text(size=10,  family="Arial")) +
  labs(x="", y="Pearson correlation", title = "O. sativa")
# p_o

df.cot_bnigra <- df.cot[df.cot$data == "bnigra.20x", ]
p_b <- ggplot(data = df.cot_bnigra, aes(y=Pearson_Correlation, x=motif, fill = method, 
                                        label=sprintf("%.4f", 
                                                      round(Pearson_Correlation, digits = 4)))) +
  geom_bar(colour="black", size=0.3, stat="identity", position = position_dodge()) +
  geom_text(size=2.2, position=position_dodge(width=0.9),
            family="Arial", angle = 90, hjust=1.1) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  scale_fill_manual(values = cbPalettec, 
                    breaks=c("DeepSignal-plant", "Megalodon", "Megalodon (retrain)"), 
                    labels=c(" DeepSignal-plant         ", " Megalodon        ", 
                             " Megalodon (retrain)")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=10, family = "Arial"),
        legend.key.size = unit(0.22, "cm"),
        legend.margin=margin(-5, 0, 0, 0),
        plot.title = element_text(size=12, hjust = 0.5, face = "italic", family="Arial"),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=10),
        text=element_text(size=10,  family="Arial")) +
  labs(x="", y="Pearson correlation", title = "B. nigra")
# p_b



# =====
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ppi= 300
png("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.bnigra.raw.2.png", 
    width = 15, 
    height = 7.5, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_a + theme(legend.position = "none"), 
                         grid.rect(gp=gpar(col="white")),
                         p_o + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_b + theme(legend.position = "none"),
                         nrow = 1, 
                         ncol = 5,
                         widths=c(5.66, 0.5, 5.66, 0.5, 5.66)), 
             g_legend(p_a), 
             nrow=2, 
             heights = c(16, 1))
dev.off()

svg("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.bnigra.raw.2.svg", 
    width = 15/2.54, 
    height = 7.5/2.54)
grid.arrange(arrangeGrob(p_a + theme(legend.position = "none"), 
                         grid.rect(gp=gpar(col="white")),
                         p_o + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_b + theme(legend.position = "none"),
                         nrow = 1, 
                         ncol = 5,
                         widths=c(5.66, 0.5, 5.66, 0.5, 5.66)), 
             g_legend(p_a), 
             nrow=2, 
             heights = c(16, 1))
dev.off()

pdf("tools_to_cmp/eval.tools_cmp.10x.arab_rep123.rice_2-1_rep2-1.bnigra.raw.2.pdf", 
    width = 15/2.54, 
    height = 7.5/2.54)
grid.arrange(arrangeGrob(p_a + theme(legend.position = "none"), 
                         grid.rect(gp=gpar(col="white")),
                         p_o + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_b + theme(legend.position = "none"),
                         nrow = 1, 
                         ncol = 5,
                         widths=c(5.66, 0.5, 5.66, 0.5, 5.66)), 
             g_legend(p_a), 
             nrow=2, 
             heights = c(16, 1))
dev.off()

