library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)


# read data
df_signalinfo <- read.table("fig_pipeline_testing/athaliana.guppy.pass.part1.CHH.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_6m.denoised_kmer_stats.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)

# kmer = "GGTCCTCAAATGT"
# kmer = "ACATGGCAAACAT"
# kmer = "ATCCAACAACTGG"
kmer = "TTTGTTCAATCAA"
df_signalinfo_one <- df_signalinfo[df_signalinfo$kmer == kmer, ]
df_signalinfo_one$label <- factor(df_signalinfo_one$label, levels = c(1, -1, 0))
positions = seq(-6, 6, 1)
colnames(df_signalinfo_one)[4:16] = positions
cnt_pos_kept = sum(df_signalinfo_one$label==1)
cnt_pos_removed = sum(df_signalinfo_one$label==-1)
cnt_neg = sum(df_signalinfo_one$label==0)


# plot 
df_signalinfo_one_m <- melt(df_signalinfo_one, id.vars=c("kmer", "label", "sid"))
x_breaks = positions

signal_text = data.frame(xaxis=as.factor(positions), 
                         kmer=unlist(strsplit(kmer,split = "")), 
                         face=c(rep("plain", 6), rep("bold", 3), rep("plain", 4)))
text_face=c(rep("plain", 6), rep("bold", 3), rep("plain", 4))

cbPaletteb <- c("#d7191c", "#fdae61", "#abdda4")
min_y = -2.2
pb <- ggplot(df_signalinfo_one_m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=label), position="dodge", lwd=0.5, outlier.size=0.5, fatten=1) + # 0.5
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-5, 0, 0, 0),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size=10, family = "serif"), # 10
        text = element_text(size=12, family="serif"), # 12
        plot.title = element_text(size=10, family="serif")) + # 10
  scale_fill_manual(values=cbPaletteb, breaks=c(1, -1, 0), 
                    labels=c("positive_kept  ", "positive_removed  ", "negative")) +
  scale_x_discrete(labels=x_breaks, 
                   breaks=positions, 
                   limits=as.factor(positions)) +
  scale_y_continuous(limits = c(min_y, 2.2), 
                     breaks = seq(-2, 2, 1)) + 
  geom_text(data=signal_text, 
            aes(x=xaxis, y=min_y, label=kmer, colour=kmer),
            hjust=0.5, size=3.5, # 3
            show.legend=FALSE, vjust=0, angle=0,
            fontface=text_face, 
            family="serif") +
  ylab("signal") + 
  xlab("position") + 
  ggtitle(sprintf("denoise training samples (positive_kept: %d, positive_removed: %d, negative: %d)", 
                  cnt_pos_kept, cnt_pos_removed, cnt_neg))
pb

ppi= 300
jpeg("fig_pipeline_testing/effect_of_denoise.plot.jpg", 
     width = 15, 
     height = 7, units = "cm", res=ppi)
pb
dev.off()






