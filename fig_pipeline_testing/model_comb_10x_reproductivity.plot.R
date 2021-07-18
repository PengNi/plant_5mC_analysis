library(ggplot2)
library(reshape2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)

library(RColorBrewer)
library(grid)
library(gridExtra)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

make_matrix <- function(heatdf){
  idx1 <- vector(mode = "character", length = 10)
  idx2 <- vector(mode = "character", length = 10)
  count <- 1
  for(i in c(1:4)){
    for(j in c((i+1):5)){
      idx1[count]=paste("repeat", i, sep = " ")
      idx2[count]=paste("repeat", j, sep = " ")
      count = count + 1
    }
  }
  heatdf$idx1 = idx1
  heatdf$idx2 = idx2
  names <- paste("repeat", c(1:5), sep = " ")
  heatdf.mat <- matrix(data=0, nrow=5, ncol=5, 
                             dimnames = list(names, names))
  for(i in 1:nrow(heatdf)){
    id1 = heatdf[i, "idx1"]
    id2 = heatdf[i, "idx2"]
    value = heatdf[i, "pearson"]
    heatdf.mat[id1, id2] = value
    heatdf.mat[id2, id1] = value
  }
  diag(heatdf.mat) <- 1.0000
  heatdf.mat <- apply(heatdf.mat,2,rev)
  return(heatdf.mat)
}

heatplot <- function(heat_mat, xlab="", ylab="", title="", 
                     labs=TRUE, digits=4, labs.size=5) {
  sim.df <- as.data.frame(heat_mat)
  rn <- row.names(sim.df)
  
  sim.df <- cbind(ID=rownames(sim.df), sim.df)
  sim.df <- melt(sim.df)
  
  sim.df[,1] <- factor(sim.df[,1], levels=rev(rn))
  if (labs == TRUE) {
    ## lbs <- c(apply(round(sim, digits), 2, as.character))
    sim.df$label <- as.character(format(round(sim.df$value, digits), nsmall = digits))
  }
  variable <- ID <- value <- label <- NULL ## to satisfy codetools
  if (labs == TRUE)
    p <- ggplot(sim.df, aes(variable, ID, fill=value, label=label))
  else
    p <- ggplot(sim.df, aes(variable, ID, fill=value))
  
  p <- p + geom_tile(color="black")+
    # scale_fill_gradient(low=color.low, high=color.high, limits=c(0, 1)) +
    scale_fill_distiller(palette = "RdYlBu", 
                         limit = c(0, 1)) + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))+
    theme(axis.ticks=element_blank())
  if (labs == TRUE)
    p <- p+geom_text(size=labs.size, family = "Arial")
  p <- p + theme(axis.text.x=element_text( hjust=0.5, angle=0)) +
    theme(axis.text.y=element_text(vjust=0.5, hjust=0.5, angle=90))
  p <- p+theme(legend.title=element_blank())
  p <- p+theme(text = element_text(size = 17, family = "Arial"))
  ##geom_point(aes(size=value))
  p <- p+xlab(xlab)+ylab(ylab) + ggtitle(" ")
  
  return(p)
}


# ============ read ===========
df_comb <- read.table("fig_pipeline_testing/model_comb_10x_reproductivity.txt", 
                      header = T, sep = "\t", stringsAsFactors = F)
df_comb_meanstd <- summarySE(df_comb, "Pearson_correlation", c("motif", "species"))



# plot ==============================
cbPalettec <- c("#1f78b4", "#33a02c")

# ====== arab
species = "A. thaliana"

motif = "CG"
pearson_corrs <- data.frame(pearson_corr = df_comb[df_comb$motif==motif & 
                                                    df_comb$species==species,
                                                  ]$Pearson_correlation, 
                            xaxis = paste(species, motif, sep = " "))
p_a_cg <- ggplot(pearson_corrs, aes(y=pearson_corr, x=xaxis, fill=xaxis)) +
  geom_boxplot(position=position_dodge(),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=17, family="Arial"),
        axis.title = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title=element_text(size=17, family="Arial", hjust=0.5)) +
  scale_fill_manual(values=cbPalettec[1]) +
  scale_y_continuous(labels = fmt_dcimals(4)) +
  labs(title=bquote(paste(italic(.(species)), "   CpG", sep = "")),
       y="Pearson correlation", 
       x="") + 
  annotation_custom(grobTree(textGrob(paste("mean: ", 
                                            round(mean(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.92, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial")))) + 
  annotation_custom(grobTree(textGrob(paste("   s.d.: ", 
                                            round(sd(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.87, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial"))))
p_a_cg


df_a_cg_heat <- read.table("fig_pipeline_testing/model_cnn_comb_p0.8_10x_self_corr.arab.part2.CG.txt", 
                           header = T, sep = "\t", stringsAsFactors = F)
df_a_cg_heat.mat <- make_matrix(df_a_cg_heat)
p_a_cg_heat <- heatplot(df_a_cg_heat.mat, digits = 4)
p_a_cg_heat


motif = "CHG"
pearson_corrs <- data.frame(pearson_corr = df_comb[df_comb$motif==motif & 
                                                     df_comb$species==species,
                                                   ]$Pearson_correlation, 
                            xaxis = paste(species, motif, sep = " "))
p_a_chg <- ggplot(pearson_corrs, aes(y=pearson_corr, x=xaxis, fill=xaxis)) +
  geom_boxplot(position=position_dodge(),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=17, family="Arial"),
        axis.title = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title=element_text(size=17, family="Arial", hjust=0.5)) +
  scale_fill_manual(values=cbPalettec[1]) +
  scale_y_continuous(labels = fmt_dcimals(4)) +
  labs(title=bquote(paste(italic(.(species)), "   ", .(motif), sep = "")),
       y="Pearson correlation", 
       x="") + 
  annotation_custom(grobTree(textGrob(paste("mean: ", 
                                            format(round(mean(pearson_corrs$pearson_corr), 4),
                                                   nsmall = 4)), 
                                      x=0.51,  y=0.92, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial")))) + 
  annotation_custom(grobTree(textGrob(paste("   s.d.: ", 
                                            round(sd(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.87, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial"))))
p_a_chg


df_a_chg_heat <- read.table("fig_pipeline_testing/model_cnn_comb_p0.8_10x_self_corr.arab.part2.CHG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
df_a_chg_heat.mat <- make_matrix(df_a_chg_heat)
p_a_chg_heat <- heatplot(df_a_chg_heat.mat, digits = 4)
p_a_chg_heat


motif = "CHH"
pearson_corrs <- data.frame(pearson_corr = df_comb[df_comb$motif==motif & 
                                                     df_comb$species==species,
                                                   ]$Pearson_correlation, 
                            xaxis = paste(species, motif, sep = " "))
p_a_chh <- ggplot(pearson_corrs, aes(y=pearson_corr, x=xaxis, fill=xaxis)) +
  geom_boxplot(position=position_dodge(),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=17, family="Arial"),
        axis.title = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title=element_text(size=17, family="Arial", hjust=0.5)) +
  scale_fill_manual(values=cbPalettec[1]) +
  scale_y_continuous(labels = fmt_dcimals(4)) +
  labs(title=bquote(paste(italic(.(species)), "   ", .(motif), sep = "")),
       y="Pearson correlation", 
       x="") + 
  annotation_custom(grobTree(textGrob(paste("mean: ", 
                                            round(mean(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.92, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial")))) + 
  annotation_custom(grobTree(textGrob(paste("   s.d.: ", 
                                            round(sd(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.87, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial"))))
p_a_chh

df_a_chh_heat <- read.table("fig_pipeline_testing/model_cnn_comb_p0.8_10x_self_corr.arab.part2.CHH.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
df_a_chh_heat.mat <- make_matrix(df_a_chh_heat)
p_a_chh_heat <- heatplot(df_a_chh_heat.mat, digits = 4)
p_a_chh_heat


# ==== rice
species = "O. sativa"

motif = "CG"
pearson_corrs <- data.frame(pearson_corr = df_comb[df_comb$motif==motif & 
                                                     df_comb$species==species,
                                                   ]$Pearson_correlation, 
                            xaxis = paste(species, motif, sep = " "))
p_o_cg <- ggplot(pearson_corrs, aes(y=pearson_corr, x=xaxis, fill=xaxis)) +
  geom_boxplot(position=position_dodge(),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=17, family="Arial"),
        axis.title = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title=element_text(size=17, family="Arial", hjust=0.5)) +
  scale_fill_manual(values=cbPalettec[2]) +
  scale_y_continuous(labels = fmt_dcimals(4)) +
  labs(title=bquote(paste(italic(.(species)), "   CpG", sep = "")),
       y="Pearson correlation", 
       x="") + 
  annotation_custom(grobTree(textGrob(paste("mean: ", 
                                            round(mean(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.92, hjust=0, 
                                      gp=gpar(col="black", fontsize=9.5, fontfamily="Arial")))) + 
  annotation_custom(grobTree(textGrob(paste("   s.d.: ", 
                                            round(sd(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.87, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial"))))
p_o_cg

df_o_cg_heat <- read.table("fig_pipeline_testing/model_cnn_comb_p0.8_10x_self_corr.rice2-1.part2.CG.txt", 
                           header = T, sep = "\t", stringsAsFactors = F)
df_o_cg_heat.mat <- make_matrix(df_o_cg_heat)
p_o_cg_heat <- heatplot(df_o_cg_heat.mat, digits = 4)
p_o_cg_heat


motif = "CHG"
pearson_corrs <- data.frame(pearson_corr = df_comb[df_comb$motif==motif & 
                                                     df_comb$species==species,
                                                   ]$Pearson_correlation, 
                            xaxis = paste(species, motif, sep = " "))
p_o_chg <- ggplot(pearson_corrs, aes(y=pearson_corr, x=xaxis, fill=xaxis)) +
  geom_boxplot(position=position_dodge(),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=17, family="Arial"),
        axis.title = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title=element_text(size=17, family="Arial", hjust=0.5)) +
  scale_fill_manual(values=cbPalettec[2]) +
  scale_y_continuous(labels = fmt_dcimals(4)) +
  labs(title=bquote(paste(italic(.(species)), "   ", .(motif), sep = "")),
       y="Pearson correlation", 
       x="") + 
  annotation_custom(grobTree(textGrob(paste("mean: ", 
                                            round(mean(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.92, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial")))) + 
  annotation_custom(grobTree(textGrob(paste("   s.d.: ", 
                                            round(sd(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.87, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial"))))
p_o_chg


df_o_chg_heat <- read.table("fig_pipeline_testing/model_cnn_comb_p0.8_10x_self_corr.rice2-1.part2.CHG.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
df_o_chg_heat.mat <- make_matrix(df_o_chg_heat)
p_o_chg_heat <- heatplot(df_o_chg_heat.mat, digits = 4)
p_o_chg_heat


motif = "CHH"
pearson_corrs <- data.frame(pearson_corr = df_comb[df_comb$motif==motif & 
                                                     df_comb$species==species,
                                                   ]$Pearson_correlation, 
                            xaxis = paste(species, motif, sep = " "))
p_o_chh <- ggplot(pearson_corrs, aes(y=pearson_corr, x=xaxis, fill=xaxis)) +
  geom_boxplot(position=position_dodge(),
               lwd=0.5, outlier.size=-1, fatten=2) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=17, family="Arial"),
        axis.title = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title=element_text(size=17, family="Arial", hjust=0.5)) +
  scale_fill_manual(values=cbPalettec[2]) +
  scale_y_continuous(labels = fmt_dcimals(4)) +
  labs(title=bquote(paste(italic(.(species)), "   ", .(motif), sep = "")),
       y="Pearson correlation", 
       x="") + 
  annotation_custom(grobTree(textGrob(paste("mean: ", 
                                            round(mean(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.92, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial")))) + 
  annotation_custom(grobTree(textGrob(paste("   s.d.: ", 
                                            round(sd(pearson_corrs$pearson_corr), 4)), 
                                      x=0.51,  y=0.87, hjust=0, 
                                      gp=gpar(col="black", fontsize=10, fontfamily="Arial"))))
p_o_chh


df_o_chh_heat <- read.table("fig_pipeline_testing/model_cnn_comb_p0.8_10x_self_corr.rice2-1.part2.CHH.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
df_o_chh_heat.mat <- make_matrix(df_o_chh_heat)
p_o_chh_heat <- heatplot(df_o_chh_heat.mat, digits = 4)
p_o_chh_heat






# save ===
ppi= 300
png("fig_pipeline_testing/model_comb_10x_reproductivity.raw.png",
    width = 65,
    height = 26, units = "cm", res=ppi)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_a_cg,
                         p_a_cg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_a_chg,
                         p_a_chg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_a_chh,
                         p_a_chh_heat,
                         nrow=1, 
                         widths = c(7.2, 14, 0.9, 7.2, 14, 0.9, 7.2, 14)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_o_cg,
                         p_o_cg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_o_chg,
                         p_o_chg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_o_chh,
                         p_o_chh_heat,
                         nrow=1, 
                         widths = c(7.2, 14, 0.9, 7.2, 14, 0.9, 7.2, 14)), 
             heights = c(1, 12, 1, 12))
dev.off()

svg("fig_pipeline_testing/model_comb_10x_reproductivity.raw.svg",
    width = 65/2.54,
    height = 26/2.54)
grid.arrange(arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_a_cg,
                         p_a_cg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_a_chg,
                         p_a_chg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_a_chh,
                         p_a_chh_heat,
                         nrow=1, 
                         widths = c(7.2, 14, 0.9, 7.2, 14, 0.9, 7.2, 14)), 
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_o_cg,
                         p_o_cg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_o_chg,
                         p_o_chg_heat,
                         grid.rect(gp=gpar(col="white")),
                         p_o_chh,
                         p_o_chh_heat,
                         nrow=1, 
                         widths = c(7.2, 14, 0.9, 7.2, 14, 0.9, 7.2, 14)), 
             heights = c(1, 12, 1, 12))
dev.off()





