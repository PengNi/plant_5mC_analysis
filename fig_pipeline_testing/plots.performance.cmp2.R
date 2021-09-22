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
  geom_bar(stat="identity", size=0.25, position=position_dodge(), colour="black") + 
  geom_text(size = 2.5, position=position_dodge(width=0.9), # vjust=-0.5, # 2
            angle=90, hjust=1.2,
            family="Arial") +
  theme_bw() + 
  theme(legend.position = c(0.4, -0.20),
        legend.direction = "horizontal",
        legend.title = element_blank(), 
        legend.text = element_text(size = 5, family="Arial"),
        text=element_text(size=6, family="Arial"), # 12
        axis.text = element_text(size=6, family = "Arial"),
        legend.key.size = unit(0.15, "cm"), 
        legend.spacing.x = unit(0.02, 'cm'),
        legend.margin=margin(-17, 0, 0, 0)) + 
  scale_fill_manual(values=cbPalettea, 
                    labels=c("random ", "balance  ", "balance+denoise")) +
  coord_cartesian(ylim=c(0.2, 0.95)) +
  scale_y_continuous(breaks = seq(0.2, 0.95, 0.1)) + 
  labs(x="", y="Pearson correlation") + 
  guides(fill = guide_legend(label.hjust = 10))

ppi= 500
png("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.plot2.png", 
    width = 3.61, 
    height = 5.9, units = "cm", res=ppi)  # 8/14.5
pa
dev.off()
svg("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.plot2.svg", 
    width = 3.61/2.54, 
    height = 5.9/2.54)  # 8/14.5
pa
dev.off()
pdf("fig_pipeline_testing/eval.preprocess.arab.10x_1.vs_bsrep123.plot2.pdf", 
    width = 3.61/2.54, 
    height = 5.9/2.54)  # 8/14.5
pa
dev.off()







