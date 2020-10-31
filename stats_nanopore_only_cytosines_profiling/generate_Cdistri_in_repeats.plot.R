library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)


# arab repeats
arab_rep <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.c_count.bs_3reps_vs_nano_guppy_part2_50x.only_nanopore_detected_Cs.pos.cnum_in_repeats.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)
arab_rep$motif_ratio <- arab_rep$moitf_num / arab_rep$motif_total
arab_rep <- arab_rep[arab_rep$region!="Repeats(WindowMasker)", ]
arab_rep$region <- factor(arab_rep$region, 
                           levels = c("Repeats(RepeatMasker)", 
                                      "Tandem repeats", 
                                      "Inverted repeats"))

p_arab_rep <- ggplot(data=arab_rep, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text=element_text(size=16, family = "serif"),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=32, face="bold", hjust = -0.15), 
        axis.text.x  = element_text(size=10)) + 
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(limits=c(0, 0.5), breaks = seq(0, 0.5, 0.1), 
                     labels = seq(0, 0.5, 0.1) * 100) +
  labs(x="", y="ratio (%)", title="a")

# ppi= 300
# jpeg("stats_nanopore_only_cytosines_profiling/generate_Cdistri_in_repeats.arab.raw.jpg", 
#      width = 32, 
#      height = 20, units = "cm", res=ppi)
# p_arab_rep
# dev.off()


# rice repeats
rice_rep <- read.table("stats_nanopore_only_cytosines_profiling/shuidao.c_count.bs_2-1n1-2_vs_nano_guppy_shuidao1-1_50x.only_nanopore_detected_Cs.pos.cnum_in_repeats.txt", 
                       header = T, sep = "\t", stringsAsFactors = F)
rice_rep$motif_ratio <- rice_rep$moitf_num / rice_rep$motif_total
rice_rep <- rice_rep[rice_rep$region!="Repeats(WindowMasker)", ]
rice_rep$region <- factor(rice_rep$region, 
                          levels = c("Repeats(RepeatMasker)", 
                                     "Tandem repeats", 
                                     "Inverted repeats"))

p_rice_rep <- ggplot(data=rice_rep, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text = element_text(size=16, family = "serif"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=32, face="bold", hjust = -0.15), 
        axis.text.x = element_text(size=10)) + 
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(limits=c(0, 0.85), breaks = seq(0, 0.8, 0.1), 
                     labels = seq(0, 0.8, 0.1) * 100) +
  labs(x="", y="ratio (%)", title="d")
# p
# ppi= 300
# jpeg("stats_nanopore_only_cytosines_profiling/generate_Cdistri_in_repeats.rice.raw.jpg", 
#      width = 32, 
#      height = 20, units = "cm", res=ppi)
# p
# dev.off()





