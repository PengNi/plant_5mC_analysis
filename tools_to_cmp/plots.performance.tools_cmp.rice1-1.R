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

# cmp deepsignal2 and megalodon in rice1-1

df <- read.table("tools_to_cmp/eval.tools_cmp.10-50x.rice1-1_rep1-2.dp0.8.txt", 
                 header = T, sep = "\t", stringsAsFactors = F)

df$coverage <- paste(substr(df$coverage, 1, 2), "\U00D7", sep = "")
df$method <- factor(df$method, levels = c("DeepSignal2", "Megalodon"))

cbPalette <- c("#2b83ba", "#fdae61")
p <- ggplot(data = df, aes(x=coverage)) + 
  geom_bar(aes(y=Pearson_Correlation, fill=method),
           colour="black", stat="identity", position=position_dodge(0.9)) + 
  geom_text(aes(y = Pearson_Correlation, 
                label=sprintf("%.4f", 
                              round(Pearson_Correlation, digits = 4)), 
                group=method),
            position = position_dodge(0.9),
            size=8, angle=90, hjust=1.2,
            family="serif") + 
  facet_grid(. ~ motif, scales = "free") + 
  coord_cartesian(ylim=c(0.3, 1)) +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) + 
  scale_fill_manual(values=cbPalette, 
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
        strip.text.x = element_text(size = 20),
        panel.border = element_blank(), 
        text=element_text(size=17,  family="serif"), 
        axis.title = element_text(size = 15, family = "serif"), 
        axis.text=element_text(size=15, family = "serif"), 
        plot.title = element_text(hjust = -0, family="serif")) + 
  labs(x="Coverage", y="Pearson correlation", title = "")
p

ppi= 300
png("tools_to_cmp/eval.tools_cmp.10-50x.rice1-1_rep1-2.dp0.8.raw.png", 
    width = 44, 
    height = 13, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p, 
                         nrow=1, 
                         ncol=1))
dev.off()

svg("tools_to_cmp/eval.tools_cmp.10-50x.rice1-1_rep1-2.dp0.8.raw.svg", 
    width = 44/2.54, 
    height = 13/2.54)
grid.arrange(arrangeGrob(p, 
                         nrow=1, 
                         ncol=1))
dev.off()







