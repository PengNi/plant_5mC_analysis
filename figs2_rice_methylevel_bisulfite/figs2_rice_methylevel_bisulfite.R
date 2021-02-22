library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

# function 
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}





cbPalette <- c("#1b9e77", "#d95f02", "#7570b3")
# fig a
s1a <- read.table("figs2_rice_methylevel_bisulfite/bisulfite_methy_genome_level_2reps.txt", 
                  header = T, sep = "\t", stringsAsFactors = F, row.names = NULL)
s1a[s1a$motif=="CG", ]$motif <- "CpG"
s1a$motif <- factor(s1a$motif, levels = c("CpG", "CHG", "CHH"))
pa <- ggplot(data=s1a, aes(x=motif, y=methylevel, fill=replicate)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-3, 0, 0, 0),
        text = element_text(size = 12, family="Arial"), 
        plot.title = element_text(size = 18, hjust = -0.13, vjust=-0.13, face = "bold"), 
        axis.title.y = element_text(size=10)) + 
  scale_y_continuous(limits = c(0, 64), breaks = seq(0, 64, 10)) +
  scale_fill_manual(values=cbPalette, 
                    breaks = c("replicate1", "replicate2"), 
                    labels = c(" sample1    ", " sample2")) + 
  labs(x="", y="Average methylation level (%)", title="a")
pa

# fig b-d
s1b_1 <- read.table("figs2_rice_methylevel_bisulfite/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt.distr_plot.tsv", 
                    header = T, sep = "\t", stringsAsFactors = F)
s1b_2 <- read.table("figs2_rice_methylevel_bisulfite/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt.distr_plot.tsv", 
                    header = T, sep = "\t", stringsAsFactors = F)
s1b <- rbind(s1b_1, s1b_2)
s1b$ratio <- s1b$ratio * 100
pb <- ggplot(data=s1b, aes(x=ranges, y=ratio, fill=replicate)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  annotate("text",x="0.4-0.5",y=26,hjust=0,vjust=0,label="CpG", size = 6, family="Arial") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-3, 0, 0, 0),
        text = element_text(size = 12, family="Arial"), 
        axis.text.x = element_text(size=7),
        plot.title = element_text(size = 18, hjust = -0.13, vjust=-0.13, face = "bold")) + 
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  scale_fill_manual(values=cbPalette, 
                    breaks = c("replicate1", "replicate2"), 
                    labels = c(" sample1    ", " sample2")) + 
  labs(x="Methylation frequency", y="Percent of counts (%)", title="b")
pb


s1c_1 <- read.table("figs2_rice_methylevel_bisulfite/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt.distr_plot.tsv", 
                    header = T, sep = "\t", stringsAsFactors = F)
s1c_2 <- read.table("figs2_rice_methylevel_bisulfite/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt.distr_plot.tsv", 
                    header = T, sep = "\t", stringsAsFactors = F)
s1c <- rbind(s1c_1, s1c_2)
s1c$ratio <- s1c$ratio * 100
pc <- ggplot(data=s1c, aes(x=ranges, y=ratio, fill=replicate)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  annotate("text",x="0.4-0.5",y=31,hjust=0,vjust=0,label="CHG", size = 6, family="Arial") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-3, 0, 0, 0),
        text = element_text(size = 12, family="Arial"), 
        axis.text.x = element_text(size = 7),
        plot.title = element_text(size = 18, hjust = -0.13, vjust=-0.13, face = "bold")) + 
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  scale_fill_manual(values=cbPalette, 
                    breaks = c("replicate1", "replicate2"), 
                    labels = c(" sample1    ", " sample2")) + 
  labs(x="Methylation frequency", y="Percent of counts (%)", title="c")
pc


s1d_1 <- read.table("figs2_rice_methylevel_bisulfite/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt.distr_plot.tsv", 
                    header = T, sep = "\t", stringsAsFactors = F)
s1d_2 <- read.table("figs2_rice_methylevel_bisulfite/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt.distr_plot.tsv", 
                    header = T, sep = "\t", stringsAsFactors = F)
s1d <- rbind(s1d_1, s1d_2)
s1d$ratio <- s1d$ratio * 100
pd <- ggplot(data=s1d, aes(x=ranges, y=ratio, fill=replicate)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  annotate("text",x="0.4-0.5",y=46,hjust=0,vjust=0,label="CHH", size = 6, family="Arial") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(0.35, "cm"), 
        legend.margin=margin(-3, 0, 0, 0),
        text = element_text(size = 12, family="Arial"), 
        axis.text.x = element_text(size = 7),
        plot.title = element_text(size = 18, hjust = -0.13, vjust=-0.13, face = "bold")) + 
  scale_y_continuous(limits = c(0, 91), breaks = seq(0, 91, 10)) +
  scale_fill_manual(values=cbPalette, 
                    breaks = c("replicate1", "replicate2"), 
                    labels = c(" sample1    ", " sample2")) + 
  labs(x="Methylation frequency", y="Percent of counts (%)", title="d")
pd




# multiple jpg
ppi= 300
png("figs2_rice_methylevel_bisulfite/Figs2_rice_methylevel_bisulfite.png", 
     width = 20, 
     height = 16, units = "cm", res=ppi)
grid.arrange(arrangeGrob(pa,
                         grid.rect(gp=gpar(col="white")),
                         pb, 
                         nrow=1, 
                         widths = c(24, 1, 24)),
             arrangeGrob(pc,
                         grid.rect(gp=gpar(col="white")),
                         pd,
                         nrow=1, 
                         widths = c(24, 1, 24)), 
             heights=c(24, 24))
dev.off()


# svg
svg("figs2_rice_methylevel_bisulfite/Figs2_rice_methylevel_bisulfite.svg", 
     width = 20/2.54, 
     height = 16/2.54)
grid.arrange(arrangeGrob(pa,
                         grid.rect(gp=gpar(col="white")),
                         pb, 
                         nrow=1, 
                         widths = c(24, 1, 24)),
             arrangeGrob(pc,
                         grid.rect(gp=gpar(col="white")),
                         pd,
                         nrow=1, 
                         widths = c(24, 1, 24)), 
             heights=c(24, 24))
dev.off()



