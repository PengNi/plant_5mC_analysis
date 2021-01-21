library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)


# arab repeats
arab_rep <- read.table("stats_nanopore_only_cytosines_profiling/ninanjie.c_count.bs_3reps_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_repeats.txt", 
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
        plot.title = element_text(size=32, face="bold", hjust = -0.15, vjust = -0.01), 
        axis.text.x  = element_text(size=10)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c("CG   ", "CHG   ", "CHH")) +
  scale_y_continuous(limits=c(0, 0.25), breaks = seq(0, 0.25, 0.05), 
                     labels = seq(0, 0.25, 0.05) * 100) +
  labs(x="", y="Ratio (%)", title="a")

# ppi= 300
# jpeg("stats_nanopore_only_cytosines_profiling/generate_Cdistri_in_repeats.arab.raw.jpg", 
#      width = 32, 
#      height = 20, units = "cm", res=ppi)
# p_arab_rep
# dev.off()


# rice repeats
rice_rep <- read.table("stats_nanopore_only_cytosines_profiling/shuidao2-1.c_count.bs_rep2-1_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_repeats.txt", 
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
        plot.title = element_text(size=32, face="bold", hjust = -0.15, vjust = -0.01), 
        axis.text.x = element_text(size=10)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c("CG   ", "CHG   ", "CHH")) +
  scale_y_continuous(limits=c(0, 0.4), breaks = seq(0, 0.4, 0.1), 
                     labels = seq(0, 0.4, 0.1) * 100) +
  labs(x="", y="Ratio (%)", title="d")


# rice repeats rep2
rice_rep2 <- read.table("stats_nanopore_only_cytosines_profiling/shuidao1-1.c_count.bs_rep1-2_vs_nano50x_dp2_p0.8.only_nanopore_detected_Cs.pos.cnum_in_repeats.txt", 
                        header = T, sep = "\t", stringsAsFactors = F)
rice_rep2$motif_ratio <- rice_rep2$moitf_num / rice_rep2$motif_total
rice_rep2 <- rice_rep2[rice_rep2$region!="Repeats(WindowMasker)", ]
rice_rep2$region <- factor(rice_rep2$region, 
                          levels = c("Repeats(RepeatMasker)", 
                                     "Tandem repeats", 
                                     "Inverted repeats"))

p_rice_rep2 <- ggplot(data=rice_rep2, aes(x=region, y=motif_ratio, fill=motif)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw() + 
  theme(text = element_text(size=16, family = "serif"),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.margin = margin(-5, 0, 0, 0),
        plot.title = element_text(size=32, face="bold", hjust = -0.15, vjust = -0.01), 
        axis.text.x = element_text(size=10)) + 
  scale_fill_brewer(palette="Set2", 
                    labels = c("CG   ", "CHG   ", "CHH")) +
  scale_y_continuous(limits=c(0, 0.48), breaks = seq(0, 0.5, 0.1), 
                     labels = seq(0, 0.5, 0.1) * 100) +
  labs(x="", y="Ratio (%)", title="g")




