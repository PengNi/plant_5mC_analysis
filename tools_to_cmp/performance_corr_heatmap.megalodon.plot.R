## geom_bin2d() may be a better way ##
# https://stackoverflow.com/questions/50331320/how-to-move-tick-marks-and-labels-at-right-left-end-of-tiles-in-geom-tile-ggplot
# https://stackoverflow.com/questions/10710463/how-to-force-the-x-axis-tick-marks-to-appear-at-the-end-of-bar-in-heatmap-graph
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

############
# arab, rice2-1; megalodon bs_cov_15
############


classify_one_type <- function(rmet, sranges){
  for(i in 2:(length(sranges))){
    if(rmet<=sranges[i]){return(i-1)}
  }
}

classify_rtype <- function(Rmet, sranges){
  return(sapply(Rmet, classify_one_type, sranges))
}

generate_rtype_matrix <- function(Rmet1, Rmet2, rescale_num=10){
  scale_range <- 1 / rescale_num
  sranges <- seq(0, 1, scale_range)
  
  res_mat <- matrix(data=0, nrow=length(sranges)-1, ncol=length(sranges)-1, 
                    dimnames = list(sranges[2:length(sranges)], 
                                    sranges[2:length(sranges)]))
  # res_mat <- matrix(data=0, nrow=length(sranges)+1, ncol=length(sranges)+1, 
  #                   dimnames = list(c(sranges, '1+'), c(sranges, '1+')))
  rtype1 = classify_rtype(Rmet1, sranges)
  rtype2 = classify_rtype(Rmet2, sranges)
  for(i in 1:length(rtype1)){
    tmp = res_mat[rtype1[i], rtype2[i]]
    res_mat[rtype1[i], rtype2[i]] = tmp + 1
  }
  res_mat
}

plot_rmet_heatmap <- function(rmet_df_file, hunit=40, 
                              data_label="", 
                              name_bs="",
                              name_nano="", 
                              titlehjust=-0.3, 
                              pearcor=0.0){
  rmet_df <- read.table(rmet_df_file, 
                        header = T, sep = '\t', stringsAsFactors = F)
  rmet_bis <- rmet_df$rmet_bis
  rmet_nan <- rmet_df$rmet_nan
  # hunit=40
  rmet_mat_10 <- as.data.frame(generate_rtype_matrix(rmet_nan, rmet_bis, hunit))
  rmet_mat_10[rmet_mat_10==0] = 1
  sim.df <- cbind(ID=as.numeric(rownames(rmet_mat_10)), rmet_mat_10)
  sim.df <- melt(sim.df, id.vars = 'ID')
  sim.df$variable <- as.numeric(as.character(sim.df$variable))
  sim.df$logvalue <- log10(sim.df$value)
  sim.df$label <- as.character(round(sim.df$value, 2))
  
  maxvalue = max(sim.df$value)
  message(sprintf("maxvalue: %d", maxvalue))
  if(maxvalue < 1000000){
    maxvalue = 1000000
  }
  
  # data_label = "A   HX1"
  # name_bs = 'Bisulfite Methylation Frequency'
  # name_nano = "DeepSignal Methylation Frequency"
  
  # RdYlGn
  # RdYlBu
  rf <- colorRampPalette(rev(brewer.pal(11,'RdYlBu')))
  r <- rf(32)
  
  if(pearcor<=0){
    pearcor = cor(rmet_bis, rmet_nan, method='pearson')
  }
  pearcor_fmt <- format(round(pearcor, 4), nsmall = 4)
  
  p <- ggplot(rmet_df, aes(rmet_bis, rmet_nan)) + 
    geom_bin2d(bins=hunit) +
    theme_bw() + 
    scale_fill_gradientn(colors = r, 
                         limit = c(1, maxvalue), 
                         space = "Lab", name="count", trans="log10", 
                         labels=comma) + 
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme(text=element_text(size=17, family = "serif"), 
          axis.title = element_text(size = 15, family = "serif"), 
          axis.text=element_text(size=15, family = "serif"), 
          legend.text = element_text(size=12, family = "serif"), 
          legend.title = element_text(size=15, family = "serif"), 
          plot.title = element_text(hjust = titlehjust)) +
    xlab(name_bs) + ylab(name_nano) + 
    ggtitle(bquote(paste(.(data_label), ",  ", italic("r"), " = ", .(pearcor_fmt))))
  
  return(p)
}

p_arab_cg <- plot_rmet_heatmap("tools_to_cmp/megalodon_results.arab.pass.part2_guppy.50x_12345.CG_retrain_comb.bs_cov_15.vs_bsrep123.rmet.tsv", 
                               hunit = 50,
                               data_label = "CG", 
                               name_bs = "Bisulfite methylation frequency", 
                               name_nano = "Megalodon methylation frequency",  
                               titlehjust = -0.0, 
                               pearcor = 0.965437)
p_arab_cg

p_arab_chg <- plot_rmet_heatmap("tools_to_cmp/megalodon_results.arab.pass.part2_guppy.50x_12345.CHG_retrain_comb.bs_cov_15.vs_bsrep123.rmet.tsv", 
                                hunit = 50,
                                data_label = "CHG", 
                                name_bs = "Bisulfite methylation frequency", 
                                name_nano = "Megalodon methylation frequency", 
                                titlehjust = -0.0, 
                                pearcor = 0.935341)
p_arab_chg

p_arab_chh <- plot_rmet_heatmap("tools_to_cmp/megalodon_results.arab.pass.part2_guppy.50x_12345.CHH_retrain_comb.bs_cov_15.vs_bsrep123.rmet.tsv", 
                                hunit = 50,
                                data_label = "CHH", 
                                name_bs = "Bisulfite methylation frequency", 
                                name_nano = "Megalodon methylation frequency", 
                                titlehjust = -0.0, 
                                pearcor = 0.5782865)
p_arab_chh


# in methylbed, two digits of rmet are kept, which will cause imprecise, 
# then leads to lower performance,
# case1: rmet=0.4652, case2: round(rmet*100)/100 turns to 0.47
# e.g. rice1-1 CHH pearson cors with BS: case1 0.72196 vs case2 0.7219
# so here we use manually filled pearson correlations

p_rice_cg <- plot_rmet_heatmap("tools_to_cmp/megalodon_results.rice2-1.pass2_guppy.50x_12345.CG_retrain_comb.bs_cov_15.vs_bs2-1.rmet.tsv", 
                               hunit = 50,
                               data_label = "CG", 
                               name_bs = "Bisulfite methylation frequency", 
                               name_nano = "Megalodon methylation frequency", 
                               titlehjust = -0.0, 
                               pearcor = 0.987095)
p_rice_cg

p_rice_chg <- plot_rmet_heatmap("tools_to_cmp/megalodon_results.rice2-1.pass2_guppy.50x_12345.CHG_retrain_comb.bs_cov_15.vs_bs2-1.rmet.tsv", 
                                hunit = 50,
                                data_label = "CHG", 
                                name_bs = "Bisulfite methylation frequency", 
                                name_nano = "Megalodon methylation frequency", 
                                titlehjust = -0.0, 
                                pearcor = 0.9550643)
p_rice_chg

p_rice_chh <- plot_rmet_heatmap("tools_to_cmp/megalodon_results.rice2-1.pass2_guppy.50x_12345.CHH_retrain_comb.bs_cov_15.vs_bs2-1.rmet.tsv", 
                                hunit = 50,
                                data_label = "CHH", 
                                name_bs = "Bisulfite methylation frequency", 
                                name_nano = "Megalodon methylation frequency", 
                                titlehjust = -0.0, 
                                pearcor = 0.715651)
p_rice_chh


svg("tools_to_cmp/performance_corr_heatmap.megalodon.raw.svg", 
    width = 44/2.54, 
    height = 23/2.54)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_arab_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chh,
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh,
                         nrow=1, 
                         widths = c(12, 1, 12, 1, 12)), 
             heights = c(1, 12, 2, 12))
dev.off()
ppi= 300
png("tools_to_cmp/performance_corr_heatmap.megalodon.raw.png",
    width = 44,
    height = 23, units = "cm", res=ppi)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")),
                         nrow = 1),
             arrangeGrob(p_arab_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_arab_chh,
                         nrow=1,
                         widths = c(12, 1, 12, 1, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")),
                         nrow = 1),
             arrangeGrob(p_rice_cg,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_rice_chh,
                         nrow=1,
                         widths = c(12, 1, 12, 1, 12)),
             heights = c(1, 12, 2, 12))
dev.off()






