library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)


library(eulerr)
library(extrafont)
library(dplyr)

# prop ==
# arab high=================================
m_arab_highly <- euler(c("bisulfite"= 180774, 
                        "nanopore" = 408232,
                        "bisulfite&nanopore"=1359931))
p_arab_highly_venn <- plot(m_arab_highly, 
                          quantities = list(fontsize=12, fontfamily="serif", 
                                            labels=c(180774, 
                                                     408232, 
                                                     1359931)), 
                          fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                          edges =TRUE,
                          labels = list(font=2, fontsize=12, fontfamily="serif"))
# p_arab_highly_venn
m_arab_highly_ratio <- data.frame(group=c("CG", "CHG", "CHH"), 
                                  value=c(1187348, 133892,38691))
m_arab_highly_ratio <- m_arab_highly_ratio %>% 
  arrange(desc(group)) %>%
  mutate(prop = round(value / sum(m_arab_highly_ratio$value) *100, 1)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# cbpalette = c("#e31a1c", "#fd8d3c", "#fecc5c")
p_arab_highly_pie <- ggplot(m_arab_highly_ratio, aes(x="", 
                                                     y=prop, 
                                                     fill=group)) +
  geom_bar(stat="identity", width=2, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="bottom", 
        legend.title = element_blank(),
        text=element_text(family="serif", size=13)) +
  geom_text(aes(y = ypos, label = paste(prop, "%", sep="")), 
            color = "white", 
            size=5,
            angle=75, 
            fontface="bold") +
  scale_fill_brewer(palette="Set2")
# p_arab_highly_pie
# arab low ===============================================
m_arab_lowly <- euler(c("bisulfite"= 294339, 
                         "nanopore" = 1176548,
                         "bisulfite&nanopore"=37626441))
p_arab_lowly_venn <- plot(m_arab_lowly, 
                           quantities = list(fontsize=12, fontfamily="serif", 
                                             labels=c(294339, 
                                                      1176548, 
                                                      37626441)), 
                           fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                           edges =TRUE,
                           labels = list(font=2, fontsize=12, fontfamily="serif"))
# p_arab_lowly_venn
m_arab_lowly_ratio <- data.frame(group=c("CG", "CHG", "CHH"), 
                                  value=c(3725710,
                                          5071373,
                                          28829358))
m_arab_lowly_ratio <- m_arab_lowly_ratio %>% 
  arrange(desc(group)) %>%
  mutate(prop = round(value / sum(m_arab_lowly_ratio$value) *100, 1)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# cbpalette = c("#66c2a5", "#fc8d62", "#8da0cb")
p_arab_lowly_pie <- ggplot(m_arab_lowly_ratio, aes(x="", 
                                                     y=prop, 
                                                     fill=group)) +
  geom_bar(stat="identity", width=2, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="bottom", 
        legend.title = element_blank(),
        text=element_text(family="serif", size=13)) +
  geom_text(aes(y = ypos, label = paste(prop, "%", sep="")), 
            color = "white", 
            size=5,
            angle=70, 
            fontface="bold") +
  scale_fill_brewer(palette="Set2")
# p_arab_lowly_pie

# ppi= 300
# # 36, 30
# png("stats_correlation_with_bisulfite/comparison_highly_lowly_methylated_sites.venn.arab.raw.png",
#      width = 24,
#      height = 20, units = "cm", res=ppi)
# grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_arab_lowly_venn,
#                          p_arab_lowly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 16)),
#              arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_arab_highly_venn,
#                          p_arab_highly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 16)),
#              heights = c(12, 12))
# dev.off()
# 
# svg("stats_correlation_with_bisulfite/comparison_highly_lowly_methylated_sites.venn.arab.raw.svg",
#      width = 24/2.54,
#      height = 20/2.54)
# grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_arab_lowly_venn,
#                          p_arab_lowly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 16)), # 1, 9, 16
#              arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_arab_highly_venn,
#                          p_arab_highly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 16)),
#              heights = c(12, 12))
# dev.off()



# rice high=================================
m_rice_highly <- euler(c("bisulfite"= 1173015, 
                         "nanopore" = 3956150,
                         "bisulfite&nanopore"=20254672))
p_rice_highly_venn <- plot(m_rice_highly, 
                           quantities = list(fontsize=12, fontfamily="serif", 
                                             labels=c(1173015, 
                                                      3956150, 
                                                      20254672)), 
                           fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                           edges =TRUE,
                           labels = list(font=2, fontsize=12, fontfamily="serif"))
# p_rice_highly_venn
m_rice_highly_ratio <- data.frame(group=c("CG", "CHG", "CHH"), 
                                  value=c(14907937, 4897659,449076))
m_rice_highly_ratio <- m_rice_highly_ratio %>% 
  arrange(desc(group)) %>%
  mutate(prop = round(value / sum(m_rice_highly_ratio$value) *100, 1)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# cbpalette = c("#e31a1c", "#fd8d3c", "#fecc5c")
p_rice_highly_pie <- ggplot(m_rice_highly_ratio, aes(x="", 
                                                     y=prop, 
                                                     fill=group)) +
  geom_bar(stat="identity", width=2, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="bottom", 
        legend.title = element_blank(),
        text=element_text(family="serif", size=13)) +
  geom_text(aes(y = ypos, label = paste(prop, "%", sep="")), 
            color = "white", 
            size=5,
            angle=75, 
            fontface="bold") +
  scale_fill_brewer(palette="Set2")
# p_rice_highly_pie
# rice low ===============================================
m_rice_lowly <- euler(c("bisulfite"= 2970271, 
                        "nanopore" = 6861540,
                        "bisulfite&nanopore"=121840663))
p_rice_lowly_venn <- plot(m_rice_lowly, 
                          quantities = list(fontsize=12, fontfamily="serif", 
                                            labels=c(2970271, 
                                                     6861540, 
                                                     121840663)), 
                          fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                          edges =TRUE,
                          labels = list(font=2, fontsize=12, fontfamily="serif"))
# p_rice_lowly_venn
m_rice_lowly_ratio <- data.frame(group=c("CG", "CHG", "CHH"), 
                                 value=c(12388200,
                                         16272073,
                                         93180390))
m_rice_lowly_ratio <- m_rice_lowly_ratio %>% 
  arrange(desc(group)) %>%
  mutate(prop = round(value / sum(m_rice_lowly_ratio$value) *100, 1)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# cbpalette = c("#66c2a5", "#fc8d62", "#8da0cb")
p_rice_lowly_pie <- ggplot(m_rice_lowly_ratio, aes(x="", 
                                                   y=prop, 
                                                   fill=group)) +
  geom_bar(stat="identity", width=2, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="bottom", 
        legend.title = element_blank(),
        text=element_text(family="serif", size=13)) +
  geom_text(aes(y = ypos, label = paste(prop, "%", sep="")), 
            color = "white", 
            size=5,
            angle=70, 
            fontface="bold") +
  scale_fill_brewer(palette="Set2")
# p_rice_lowly_pie

# ppi= 300
# png("stats_correlation_with_bisulfite/comparison_highly_lowly_methylated_sites.venn.rice.raw.png",
#      width = 24,
#      height = 20, units = "cm", res=ppi)
# grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_rice_lowly_venn,
#                          p_rice_lowly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 16)),
#              arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_rice_highly_venn,
#                          p_rice_highly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 16)),
#              heights = c(12, 12))
# dev.off()
# 
# svg("stats_correlation_with_bisulfite/comparison_highly_lowly_methylated_sites.venn.rice.raw.svg",
#      width = 24/2.54,
#      height = 20/2.54)
# grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_rice_lowly_venn,
#                          p_rice_lowly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 17)),
#              arrangeGrob(grid.rect(gp=gpar(col="white")),
#                          p_rice_highly_venn,
#                          p_rice_highly_pie,
#                          nrow=1,
#                          widths = c(1, 9, 17)),
#              heights = c(12, 12))
# dev.off()












# ===
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}



cbPalette <- c("#d7301f", "#fc8d59", "#2b8cbe", "#7bccc4")


# read data ===
f_arab_cg <- read.table("stats_correlation_with_bisulfite/athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.high_low_in_bs.tsv", 
                        sep = "\t", header = F, stringsAsFactors = F)
colnames(f_arab_cg) <- c("rmet", "methylevel")
f_arab_chg <- read.table("stats_correlation_with_bisulfite/athaliana.guppy.pass.part2.CHG.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.freq.high_low_in_bs.tsv", 
                         sep = "\t", header = F, stringsAsFactors = F)
colnames(f_arab_chg) <- c("rmet", "methylevel")
f_arab_chh <- read.table("stats_correlation_with_bisulfite/athaliana.guppy.pass.part2.CHH.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.freq.high_low_in_bs.tsv", 
                         sep = "\t", header = F, stringsAsFactors = F)
colnames(f_arab_chh) <- c("rmet", "methylevel")

f_rice_cg <- read.table("stats_correlation_with_bisulfite/shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.high_low_in_bs.tsv", 
                        sep = "\t", header = F, stringsAsFactors = F)
colnames(f_rice_cg) <- c("rmet", "methylevel")
f_rice_chg <- read.table("stats_correlation_with_bisulfite/shuidao1-1.guppy.pass.part2.CHG.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.freq.high_low_in_bs.tsv", 
                        sep = "\t", header = F, stringsAsFactors = F)
colnames(f_rice_chg) <- c("rmet", "methylevel")
f_rice_chh <- read.table("stats_correlation_with_bisulfite/shuidao1-1.guppy.pass.part2.CHH.bn13_sn16.arabnrice2-1.denoise_signal_bilstm.both_bilstm.50x_12345.freq.high_low_in_bs.tsv", 
                         sep = "\t", header = F, stringsAsFactors = F)
colnames(f_rice_chh) <- c("rmet", "methylevel")



p_arab_cg <- ggplot(f_arab_cg, aes(x=rmet, fill=methylevel)) + 
  # geom_histogram(aes(y=..density..),  
  #                binwidth=0.02, position="dodge") +
  geom_density(colour="black", alpha=0.7, size=0.3) + 
  annotate("text",x=0.4, y=15,
           hjust=0,vjust=0,label="CG", 
           size=9, family="serif") + 
  geom_vline(xintercept=0.3, 
             colour=cbPalette[3], 
             linetype="dashed", size=1) +
  geom_vline(xintercept=0.7, 
             colour=cbPalette[1], 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=14, family = "serif"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=14, family = "serif"), 
        legend.key.size = unit(0.4, "cm"), 
        plot.title = element_text(size=32, face = "bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks = c("lowly_bs", "lowly_nano",  
                               "highly_bs", "highly_nano"), 
                    labels=c(" lowly methylated (bisulfite)    ",
                             " lowly methylated (nanopore)    ", 
                             " highly methylated (bisulfite)    ", 
                             " highly methylated (nanopore)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 30)) +
  labs(x="Methylation Frequency", y="Density", title = "a")
# p_arab_cg

p_arab_chg <- ggplot(f_arab_chg, aes(x=rmet, fill=methylevel)) + 
  # geom_histogram(aes(y=..density..),  
  #                binwidth=0.02, position="dodge") +
  geom_density(colour="black", alpha=0.7, size=0.3) + 
  annotate("text",x=0.4, y=15,
           hjust=0,vjust=0,label="CHG", 
           size=9, family="serif") + 
  geom_vline(xintercept=0.3, 
             colour=cbPalette[3], 
             linetype="dashed", size=1) +
  geom_vline(xintercept=0.7, 
             colour=cbPalette[1], 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=14, family = "serif"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        plot.title = element_text(size=32, face = "bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks = c("lowly_bs", "lowly_nano",  
                               "highly_bs", "highly_nano"), 
                    labels=c("lowly methylated (bisulfite)",
                             "lowly methylated (nanopore)", 
                             "highly methylated (bisulfite)", 
                             "highly methylated (nanopore)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 30)) +
  labs(x="Methylation Frequency", y="Density", title = "b")
# p_arab_chg

p_arab_chh <- ggplot(f_arab_chh, aes(x=rmet, fill=methylevel)) + 
  # geom_histogram(aes(y=..density..),  
  #                binwidth=0.02, position="dodge") +
  geom_density(colour="black", alpha=0.7, size=0.3) + 
  annotate("text",x=0.4, y=15,
           hjust=0,vjust=0,label="CHH", 
           size=9, family="serif") + 
  geom_vline(xintercept=0.3, 
             colour=cbPalette[3], 
             linetype="dashed", size=1) +
  geom_vline(xintercept=0.7, 
             colour=cbPalette[1], 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=14, family = "serif"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        plot.title = element_text(size=32, face = "bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks = c("lowly_bs", "lowly_nano",  
                               "highly_bs", "highly_nano"), 
                    labels=c("lowly methylated (bisulfite)",
                             "lowly methylated (nanopore)", 
                             "highly methylated (bisulfite)", 
                             "highly methylated (nanopore)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 30)) +
  labs(x="Methylation Frequency", y="Density", title = "c")
# p_arab_chh

p_rice_cg <- ggplot(f_rice_cg, aes(x=rmet, fill=methylevel)) + 
  # geom_histogram(aes(y=..density..),  
  #                binwidth=0.02, position="dodge") +
  geom_density(colour="black", alpha=0.7, size=0.3) + 
  annotate("text",x=0.4, y=15,
           hjust=0,vjust=0,label="CG", 
           size=9, family="serif") + 
  geom_vline(xintercept=0.3, 
             colour=cbPalette[3], 
             linetype="dashed", size=1) +
  geom_vline(xintercept=0.7, 
             colour=cbPalette[1], 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=14, family = "serif"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        plot.title = element_text(size=32, face = "bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks = c("lowly_bs", "lowly_nano",  
                               "highly_bs", "highly_nano"), 
                    labels=c("lowly methylated (bisulfite)",
                             "lowly methylated (nanopore)", 
                             "highly methylated (bisulfite)", 
                             "highly methylated (nanopore)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 30)) +
  labs(x="Methylation Frequency", y="Density", title = "d")
# p_rice_cg

p_rice_chg <- ggplot(f_rice_chg, aes(x=rmet, fill=methylevel)) + 
  # geom_histogram(aes(y=..density..),  
  #                binwidth=0.02, position="dodge") +
  geom_density(colour="black", alpha=0.7, size=0.3) + 
  annotate("text",x=0.4, y=15,
           hjust=0,vjust=0,label="CHG", 
           size=9, family="serif") + 
  geom_vline(xintercept=0.3, 
             colour=cbPalette[3], 
             linetype="dashed", size=1) +
  geom_vline(xintercept=0.7, 
             colour=cbPalette[1], 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=14, family = "serif"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        plot.title = element_text(size=32, face = "bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks = c("lowly_bs", "lowly_nano",  
                               "highly_bs", "highly_nano"), 
                    labels=c("lowly methylated (bisulfite)",
                             "lowly methylated (nanopore)", 
                             "highly methylated (bisulfite)", 
                             "highly methylated (nanopore)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 30)) +
  labs(x="Methylation Frequency", y="Density", title = "e")
# p_rice_chg

p_rice_chh <- ggplot(f_rice_chh, aes(x=rmet, fill=methylevel)) + 
  # geom_histogram(aes(y=..density..),  
  #                binwidth=0.02, position="dodge") +
  geom_density(colour="black", alpha=0.7, size=0.3) + 
  annotate("text",x=0.4, y=15,
           hjust=0,vjust=0,label="CHH", 
           size=9, family="serif") + 
  geom_vline(xintercept=0.3, 
             colour=cbPalette[3], 
             linetype="dashed", size=1) +
  geom_vline(xintercept=0.7, 
             colour=cbPalette[1], 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=14, family = "serif"), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        plot.title = element_text(size=32, face = "bold", hjust = -0.1)) + 
  scale_fill_manual(values=cbPalette, 
                    breaks = c("lowly_bs", "lowly_nano",  
                               "highly_bs", "highly_nano"), 
                    labels=c("lowly methylated (bisulfite)",
                             "lowly methylated (nanopore)", 
                             "highly methylated (bisulfite)", 
                             "highly methylated (nanopore)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(limits = c(0, 30)) +
  labs(x="Methylation Frequency", y="Density", title = "f")
# p_rice_chh


# multiple 
p_legend = g_legend(p_arab_cg)
ppi= 300
# 55, 30
png("stats_correlation_with_bisulfite/comparison_highly_lowly_methylated_sites.png",
     width = 33,
     height = 18, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab_cg+ theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chg+ theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chh+ theme(legend.position = "none"),
                         nrow=1,
                         widths = c(12, 1, 12, 1, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")),
                         nrow = 1),
             arrangeGrob(p_rice_cg + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg+ theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh+ theme(legend.position = "none"),
                         nrow=1,
                         widths = c(12, 1, 12, 1, 12)),
             p_legend,
             heights = c(12, 1, 12, 2))
dev.off()


svg("stats_correlation_with_bisulfite/comparison_highly_lowly_methylated_sites.svg", 
     width = 33/2.54, 
     height = 18/2.54)
grid.arrange(arrangeGrob(p_arab_cg+ theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chg+ theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chh+ theme(legend.position = "none"),
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg + theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg+ theme(legend.position = "none"),
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh+ theme(legend.position = "none"),
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             p_legend,
             heights = c(12, 1, 12, 2))
dev.off()
