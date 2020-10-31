library(eulerr)
library(extrafont)
library(gridExtra)
library(grid)

font_import()
loadfonts(device = "win")


kmer_diver_chg <- euler(c("union"=25, 
                          "union+denoise"=0,
                          "intersection"=0,
                          "union&union+denoise"=22670, 
                          "union&intersection"=0,
                          "union+denoise&intersection"=0,
                          "union&union+denoise&intersection"=8847))
p_kmer_diver_chg <- plot(kmer_diver_chg, 
                            main = list(label="", 
                                        fontsize=18, 
                                        hjust=7,
                                        vjust=1), 
                            quantities = list(fontsize=17, fontfamily="serif"), 
                            fills = c('#7bccc4', '#bae4bc', '#f0f9e8'),
                            edges =TRUE,
                            labels = list(font=2, fontsize=17, fontfamily="serif"))
p_kmer_diver_chg

kmer_diver_chh <- euler(c("union"=15, 
                          "union+denoise"=0,
                          "intersection"=0,
                          "union&union+denoise"=10733, 
                          "union&intersection"=0,
                          "union+denoise&intersection"=0,
                          "union&union+denoise&intersection"=1088))
p_kmer_diver_chh <- plot(kmer_diver_chh, 
                         main = list(label="", 
                                     fontsize=18, 
                                     hjust=7,
                                     vjust=1), 
                         quantities = list(fontsize=17, fontfamily="serif"), 
                         fills = c('#7bccc4', '#bae4bc', '#f0f9e8'),
                         edges =TRUE,
                         labels = list(font=2, fontsize=17, fontfamily="serif"))
p_kmer_diver_chh

ppi= 300
jpeg("stats_high_confidence_kmers/stats_kmer_diversity_comparison.plot.raw.jpg", 
     width = 30, 
     height = 12, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_kmer_diver_chg,
                         grid.rect(gp=gpar(col="white")),
                         p_kmer_diver_chh, 
                         nrow=1, 
                         widths = c(12, 1, 12)))
dev.off()



