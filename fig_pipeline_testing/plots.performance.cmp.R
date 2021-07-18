library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)

# font_import()
# loadfonts(device = "win")


# preprocess (neg_aspos, denoise) ===
df.pre <- read.table("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.txt", 
                     header = T, 
                     sep = "\t", stringsAsFactors = F)

df.pre$preprocess <- factor(df.pre$preprocess, levels = c("random", "balance", 
                                                          "balance+denoise"))
df.pre[df.pre$motif=="CG", ]$motif <- "CpG"
df.pre$motif <- factor(df.pre$motif, levels = c("CpG", "CHG", "CHH"))
# cbPalette <- c("#fdcc8a", "#fc8d59", "#d7301f")
cbPalettea <- c("#bae4bc", "#7bccc4", "#2b8cbe")
pa <- ggplot(data = df.pre, aes(y=Pearson_Correlation, x=motif, fill = preprocess, 
                               label=sprintf("%.4f", 
                                             round(Pearson_Correlation, digits = 4)))) + 
  geom_bar(stat="identity", size=0.3, position=position_dodge(), colour="black") + 
  geom_text(size = 3, position=position_dodge(width=0.9), # vjust=-0.5, # 2
            angle=90, hjust=1.2,
            family="Arial") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 6, family="Arial"),
        text=element_text(size=9, family="Arial"), # 12
        legend.key.size = unit(0.25, "cm"), 
        legend.margin=margin(-20, 0, 0, 0)) + 
  scale_fill_manual(values=cbPalettea, 
                    labels=c(" random  ", " balance  ", " balance+denoise")) +
  coord_cartesian(ylim=c(0.2, 0.95)) +
  scale_y_continuous(breaks = seq(0.2, 0.95, 0.1)) + 
  labs(x="", y="Pearson correlation")

ppi= 300
png("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.plot.png", 
     width = 5.5, 
     height = 9, units = "cm", res=ppi)  # 8/14.5
pa
dev.off()
svg("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.plot.svg", 
    width = 5.5/2.54, 
    height = 9/2.54)  # 8/14.5
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
df.csv[df.csv$motif=="CG", ]$motif = "CpG"
df.csv$motif <- factor(df.csv$motif, levels = c("CpG", "CHG", "CHH"))
df.csv$model <- factor(df.csv$model, levels = c("m_arab", "m_rice", "m_comb"))


cbPalettec <- c("#1f78b4", "#33a02c")
df.csv_a <- df.csv[df.csv$data=="arab.20x", ]
pc_a <- ggplot() + 
  geom_bar(data = df.csv_a, 
           aes(y=Pearson_Correlation, x=model),
           fill = cbPalettec[1],
           colour="black", stat="identity", position=position_dodge(), 
           size=0.3) + 
  geom_text(data=df.csv_a, aes(x = model, y = Pearson_Correlation, 
                             label=sprintf("%.4f", 
                                           round(Pearson_Correlation, digits = 4))), size=5, # 3.2
            position=position_dodge(width=0.9), angle=90, hjust=1.2, family="Arial") +
  facet_grid(. ~ motif) + 
  scale_x_discrete(limits=c("m_arab", "m_rice", "m_comb"), 
                   labels=c("m_arab", "m_rice", "m_comb")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  coord_cartesian(ylim = c(0.3, 1)) +
  theme_bw() + 
  theme(# strip.background = element_rect(colour="white", fill="white", 
        #                                 size=1, linetype="solid"), # 1
        # panel.border = element_blank(), 
        text=element_text(size=12,  family="Arial"), # 12
        axis.text.y = element_text(size=12,  family="Arial"),
        axis.text.x = element_text(size=12,  family="Arial"),
        plot.title = element_text(size=15, face = "italic", family = "Arial", 
                                  hjust = 0.5)) + 
  labs(x="Model", y="Pearson correlation",title = "A. thaliana")
pc_a

df.csv_o <- df.csv[df.csv$data=="rice.20x", ]
pc_o <- ggplot() + 
  geom_bar(data = df.csv_o, 
           aes(y=Pearson_Correlation, x=model),
           fill = cbPalettec[2],
           colour="black", stat="identity", position=position_dodge(), 
           size=0.3) + 
  geom_text(data=df.csv_o, aes(x = model, y = Pearson_Correlation, 
                               label=sprintf("%.4f", 
                                             round(Pearson_Correlation, digits = 4))), size=5, # 3.2
            position=position_dodge(width=0.9), angle=90, hjust=1.2, family="Arial") +
  facet_grid(. ~ motif) + 
  scale_x_discrete(limits=c("m_arab", "m_rice", "m_comb"), 
                   labels=c("m_arab", "m_rice", "m_comb")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  coord_cartesian(ylim = c(0.3, 1)) +
  theme_bw() + 
  theme(# strip.background = element_rect(colour="white", fill="white", 
        #                                 size=1, linetype="solid"), # 1
        # panel.border = element_blank(), 
        text=element_text(size=12,  family="Arial"), # 12
        axis.text.y = element_text(size=12,  family="Arial"), 
        axis.text.x = element_text(size=12,  family="Arial"), 
        plot.title = element_text(size=15, face = "italic", family = "Arial", 
                                  hjust = 0.5)) + 
  labs(x="Model", y="Pearson correlation", title = "O. sativa")
pc_o

ppi= 300
png("fig_pipeline_testing/eval.cross_species_validation.10x.arab_rep123.rice2-1_rep2-1.plot.raw.png", 
    width = 36, 
    height = 13, units = "cm", res=ppi)
grid.arrange(pc_a, 
             grid.rect(gp=gpar(col="white")),
             pc_o, 
             nrow = 1,
             ncol = 3,
             widths=c(17.5, 1, 17.5))
dev.off()

svg("fig_pipeline_testing/eval.cross_species_validation.10x.arab_rep123.rice2-1_rep2-1.plot.raw.svg", 
    width = 36/2.54, 
    height = 13/2.54)
grid.arrange(pc_a, 
             grid.rect(gp=gpar(col="white")),
             pc_o, 
             nrow = 1,
             ncol = 3,
             widths=c(17.5, 1, 17.5))
dev.off()



# ==============
cbPalettec <- c("#66c2a5", "#fc8d62")
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







