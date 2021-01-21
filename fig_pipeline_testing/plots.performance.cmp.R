library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)

font_import()
loadfonts(device = "win")


# preprocess (neg_aspos, denoise) ===
df.pre <- read.table("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.txt", 
                     header = T, 
                     sep = "\t", stringsAsFactors = F)

df.pre$preprocess <- factor(df.pre$preprocess, levels = c("random", "balance", 
                                                          "balance+denoise"))
# cbPalette <- c("#fdcc8a", "#fc8d59", "#d7301f")
cbPalettea <- c("#bae4bc", "#7bccc4", "#2b8cbe")
pa <- ggplot(data = df.pre, aes(y=Pearson_Correlation, x=motif, fill = preprocess, 
                               label=sprintf("%.4f", 
                                             round(Pearson_Correlation, digits = 4)))) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  geom_text(size = 4, position=position_dodge(width=0.9), # vjust=-0.5, # 2
            angle=90, hjust=1.2,
            family="serif") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 9, family="serif"),
        text=element_text(size=14, family="serif"), # 12
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-10, 0, 0, 0)) + 
  scale_fill_manual(values=cbPalettea) +
  coord_cartesian(ylim=c(0.2, 0.95)) +
  scale_y_continuous(breaks = seq(0.2, 0.95, 0.1)) + 
  labs(x="", y="Pearson correlation")

ppi= 300
png("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.plot.png", 
     width = 8, 
     height = 14.5, units = "cm", res=ppi)  # 8/14.5
pa
dev.off()


# cross species validation ===
# df.csv <- read.table("fig_pipeline_testing/eval.cross_species_validation.10x.arab_rep123.rice1-1_rep1-2.txt", 
#                      header = T, 
#                      sep = "\t", stringsAsFactors = F)
df.csv <- read.table("fig_pipeline_testing/eval.cross_species_validation.10x.arab_rep123.rice2-1_rep2-1.txt", 
                     header = T, 
                     sep = "\t", stringsAsFactors = F)
df.csv <- ddply(df.csv, .(motif, model),
                     transform, 
                pos = cumsum(Pearson_Correlation) - (0.5 * Pearson_Correlation))
df.csv$data <- factor(df.csv$data, levels = c("rice.20x", "arab.20x"))

# charts.data <- read.csv("data/copper-data-for-tutorial.csv")
cbPalettec <- c("#66c2a5", "#fc8d62")
# cbPalettec <- c("#fc8d62", "#66c2a5")
p1 <- ggplot() + 
  geom_bar(data = df.csv, aes(y=Pearson_Correlation, x=model, fill = data),
           colour="black", stat="identity", width = 1) + 
  geom_text(data=df.csv, aes(x = model, y = pos, 
                             label=sprintf("%.4f", 
                                           round(Pearson_Correlation, digits = 4))), 
            size=7, angle=90, 
            family="serif") +
  facet_grid(. ~ motif) + 
  scale_x_discrete(limits=c("m_arab","m_rice","m_comb")) + 
  scale_y_continuous(breaks = seq(0, 2, 0.2)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("arab.20x", "rice.20x"), 
                    labels=c("A. thaliana", "O. sativa")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=30, face="italic", family = "serif"), 
        strip.background = element_rect(colour="white", fill="white", 
                                        size=2, linetype="solid"), 
        panel.border = element_blank(), 
        axis.text.x  = element_text(angle=45, vjust=1, hjust = 1),
        text=element_text(size=30,  family="serif")) + 
  labs(x="Model", y="Pearson correlation")

pc <- ggplot() + 
  geom_bar(data = df.csv, aes(y=Pearson_Correlation, x=model, fill = data),
           colour="black", stat="identity", width = 1) + 
  geom_text(data=df.csv, aes(x = model, y = pos, 
                             label=sprintf("%.4f", 
                                           round(Pearson_Correlation, digits = 4))), size=4, # 3.2
            family="serif") +
  facet_grid(motif ~ .) + 
  scale_x_discrete(limits=c("m_comb", "m_rice", "m_arab"), 
                   labels=c("m_comb", "m_rice", "m_arab")) + 
  scale_y_continuous(breaks = seq(0, 2, 0.2)) + 
  scale_fill_manual(values=cbPalettec, 
                    breaks=c("arab.20x", "rice.20x"), 
                    labels=c("A. thaliana   ", "O. sativa")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=12, face="italic", family = "serif"), # 12
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-5, 0, 0, 0),
        strip.background = element_rect(colour="white", fill="white", 
                                        size=1, linetype="solid"), # 1
        panel.border = element_blank(), 
        text=element_text(size=12,  family="serif"), # 12
        axis.text.y = element_text(size=12,  family="serif")) + 
  labs(x="Model", y="Pearson correlation") + 
  coord_flip()

ppi= 300
png("fig_pipeline_testing/eval.cross_species_validation.10x.arab_rep123.rice2-1_rep2-1.plot.png", 
     width = 16, 
     height = 7.5, units = "cm", res=ppi)
pc
dev.off()



# combine pa, pb, pc, svg/jpg
svg("fig_pipeline_testing/fig_pipeline_testing.abc.raw.svg", 
    width = 24/2.5, 
    height = 15/2.5)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), pa,
                         ncol=1, nrow = 2,
                         heights = c(0.5, 14.5)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(grid.rect(gp=gpar(col="white")), pb, 
                         grid.rect(gp=gpar(col="white")), pc,
                         ncol=1,
                         nrow=4,
                         heights = c(0.5, 7, 0.5, 7)), 
             nrow = 1,
             ncol = 3,
             widths=c(8, 1, 15))
dev.off()
ppi= 300
png("fig_pipeline_testing/fig_pipeline_testing.abc.raw.png", 
     width = 24, 
     height = 15, units = "cm", res=ppi)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), pa,
                         ncol=1, nrow = 2,
                         heights = c(0.5, 14.5)), 
             grid.rect(gp=gpar(col="white")),
             arrangeGrob(grid.rect(gp=gpar(col="white")), pb, 
                         grid.rect(gp=gpar(col="white")), pc,
                         ncol=1,
                         nrow=4,
                         heights = c(0.5, 7, 0.5, 7)), 
             nrow = 1,
             ncol = 3,
             widths=c(8, 1, 15))
dev.off()







