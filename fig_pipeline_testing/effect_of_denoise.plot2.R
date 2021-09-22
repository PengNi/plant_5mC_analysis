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
colnames(df_signalinfo_one)[2] <- "tlabel"
df_signalinfo_one$tlabel <- factor(df_signalinfo_one$tlabel, levels = c(1, -1, 0))
positions = seq(-6, 6, 1)
colnames(df_signalinfo_one)[4:16] = positions
cnt_pos_kept = sum(df_signalinfo_one$tlabel==1)
cnt_pos_removed = sum(df_signalinfo_one$tlabel==-1)
cnt_neg = sum(df_signalinfo_one$tlabel==0)


# plot 
df_signalinfo_one_m <- melt(df_signalinfo_one, id.vars=c("kmer", "tlabel", "sid"))
x_breaks = positions

signal_text = data.frame(xaxis=as.factor(positions), 
                         kmer=unlist(strsplit(kmer,split = "")), 
                         face=c(rep("plain", 6), rep("bold", 3), rep("plain", 4)))
text_face=c(rep("plain", 6), rep("bold", 3), rep("plain", 4))
bcolor = c(rep("red", 3), "yellow", rep("red", 2), 
           "blue", rep("green",2), "red", "blue", 
           rep("green", 2))

cbPaletteb <- c("#d7191c", "#fdae61", "#abdda4")
min_y = -2.2
pb <- ggplot(df_signalinfo_one_m, aes(x=variable, y=value)) + 
  stat_boxplot(aes(group=interaction(variable, tlabel)), geom = "errorbar", 
               size=0.3, 
               position="dodge") +
  geom_boxplot(aes(fill=tlabel), position="dodge", lwd=0.3, outlier.size=0.15, fatten=1) + # 0.5
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin=margin(-5, 0, 0, 0),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size=6, family = "Arial"), # 10
        axis.text = element_text(size=6, family = "Arial"),
        text = element_text(size=6, family="Arial"), # 12
        plot.title = element_text(size=6, family="Arial")) + # 10
  scale_fill_manual(values=cbPaletteb, breaks=c(1, -1, 0), 
                    labels=c(bquote(paste(" positive_kept (", italic("n"), " = ", .(cnt_pos_kept), ")   ")),
                             bquote(paste(" positive_removed (", italic("n"), " = ", .(cnt_pos_removed), ")   ")), 
                             bquote(paste(" negative (", italic("n"), " = ", .(cnt_neg), ")")))) +
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
            family="mono") +
  scale_color_manual(values = c("green", "blue", "yellow", "red")) +
  ylab("Signal") + 
  xlab("Position") + 
  # ggtitle(sprintf("Denoise training samples (positive_kept: %d, positive_removed: %d, negative: %d)", 
  #                 cnt_pos_kept, cnt_pos_removed, cnt_neg))
  ggtitle(sprintf("Denoise training samples"))
pb

ppi= 500
png("fig_pipeline_testing/effect_of_denoise.plot2.png", 
    width = 9.5, 
    height = 5.9, units = "cm", res=ppi) # 16/7.5
pb
dev.off()

svg("fig_pipeline_testing/effect_of_denoise.plot2.svg", 
    width = 9.5/2.54, 
    height = 5.9/2.54)
pb
dev.off()

pdf("fig_pipeline_testing/effect_of_denoise.plot2.pdf", 
    width = 9.5/2.54, 
    height = 5.9/2.54)
pb
dev.off()





