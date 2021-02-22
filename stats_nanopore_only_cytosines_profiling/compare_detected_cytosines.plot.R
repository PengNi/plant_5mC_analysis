library(ggplot2)
library(reshape2)
library(eulerr)
library(extrafont)
library(dplyr)
library(grid)
library(gridExtra)


# arab cytosine ==
m_arab_c <- euler(c("bisulfite"= 167689, 
                    "Nanopore" = 644273,
                    "bisulfite&Nanopore"=41951037))
p_arab_c <- plot(m_arab_c, 
                 quantities = list(fontsize=15, fontfamily="Arial", 
                                   labels=c("167,689", 
                                            "644,273", 
                                            "41,951,037")), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                 main = list(label="A. thaliana", 
                             fontsize=15,
                             font=3,
                             vjust=0.8, 
                             fontfamily="Arial"))
# p_arab_c

# rice cytosine =
m_rice_c <- euler(c("bisulfite"= 796695, 
                    "Nanopore" = 10446128,
                    "bisulfite&Nanopore"=150936284))
p_rice_c <- plot(m_rice_c, 
                 quantities = list(fontsize=15, fontfamily="Arial", 
                                   labels=c("796,695", 
                                            "10,446,128", 
                                            "150,936,284")), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                 main = list(label=expression(paste(italic("O. sativa"),  " (sample1)", sep = "")), 
                             fontsize=15,
                             font=3,
                             vjust=0.8, 
                             fontfamily="Arial"))
# p_rice_c


# rice cytosine rep2 =
m_rice_c2 <- euler(c("bisulfite"= 810454, 
                    "Nanopore" = 9432531,
                    "bisulfite&Nanopore"=151935813))
p_rice_c2 <- plot(m_rice_c2, 
                 quantities = list(fontsize=15, fontfamily="Arial", 
                                   labels=c("810,454", 
                                            "9,432,531", 
                                            "151,935,813")), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                 main = list(label=expression(paste(italic("O. sativa"),  " (sample2)", sep = "")), 
                             fontsize=15,
                             # font=3,
                             vjust=0.8, 
                             fontfamily="Arial"))
# p_rice_c2


ppi= 300
png("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice1n2.c.raw.png",
     width = 36,
     height = 10, units = "cm", res=ppi)
grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c,
             grid.rect(gp=gpar(col="white")), p_rice_c2, 
             widths = c(11.3, 1, 11.3, 1, 11.3))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice1n2.c.raw.svg",
    width = 36/2.54,
    height = 10/2.54)
grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c,
             grid.rect(gp=gpar(col="white")), p_rice_c2, 
             widths = c(11.3, 1, 11.3, 1, 11.3))
dev.off()



# arab cg, chg, chh
m_arab_cg <- euler(c("bisulfite"= 32636, 
                     "Nanopore" = 84684,
                     "bisulfite&Nanopore"=5436360))
p_arab_cg <- plot(m_arab_cg, 
                  quantities = list(fontsize=15, fontfamily="Arial", 
                                    labels=trimws(format(c(32636, 
                                             84684, 
                                             5436360), 
                                             big.mark = ",",
                                             scientific = F))),
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                  main = list(label="CpG", 
                              fontsize=15,
                              font=1,
                              vjust=0.8, 
                              fontfamily="Arial"))

m_arab_chg <- euler(c("bisulfite"= 18891, 
                      "Nanopore" = 85575,
                      "bisulfite&Nanopore"=5977439))
p_arab_chg <- plot(m_arab_chg, 
                   quantities = list(fontsize=15, fontfamily="Arial", 
                                     labels=trimws(format(c(18891, 
                                              85575, 
                                              5977439), 
                                              big.mark = ",",
                                              scientific = F))),
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                   main = list(label="CHG", 
                               fontsize=15,
                               font=1,
                               vjust=0.8, 
                               fontfamily="Arial"))


m_arab_chh <- euler(c("bisulfite"= 116162, 
                      "Nanopore" = 474014,
                      "bisulfite&Nanopore"=30537238))
p_arab_chh <- plot(m_arab_chh, 
                   quantities = list(fontsize=15, fontfamily="Arial", 
                                     labels=trimws(format(c(116162, 
                                              474014, 
                                              30537238), 
                                              big.mark = ",",
                                              scientific = F))), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                   main = list(label="CHH", 
                               fontsize=15,
                               font=1,
                               vjust=0.8, 
                               fontfamily="Arial"))


# rice cg, chg, chh
m_rice_cg <- euler(c("bisulfite"= 224432, 
                     "Nanopore" = 2010752,
                     "bisulfite&Nanopore"=28488226))
p_rice_cg <- plot(m_rice_cg, 
                  quantities = list(fontsize=15, fontfamily="Arial", 
                                    labels=trimws(format(c(224432, 
                                             2010752, 
                                             28488226), 
                                             big.mark = ",",
                                             scientific = F))), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                  main = list(label="CpG", 
                              fontsize=15,
                              font=1,
                              vjust=0.8, 
                              fontfamily="Arial"))

m_rice_chg <- euler(c("bisulfite"= 122935, 
                      "Nanopore" = 1552407,
                      "bisulfite&Nanopore"=25644528))
p_rice_chg <- plot(m_rice_chg, 
                   quantities = list(fontsize=15, fontfamily="Arial", 
                                     labels=trimws(format(c(122935, 
                                              1552407, 
                                              25644528), 
                                              big.mark = ",",
                                              scientific = F))),
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                   main = list(label="CHG", 
                               fontsize=15,
                               font=1,
                               vjust=0.8, 
                               fontfamily="Arial"))

m_rice_chh <- euler(c("bisulfite"= 449328, 
                      "Nanopore" = 6882969,
                      "bisulfite&Nanopore"=96803530))
p_rice_chh <- plot(m_rice_chh, 
                   quantities = list(fontsize=15, fontfamily="Arial", 
                                     labels=trimws(format(c(449328, 
                                              6882969, 
                                              96803530), 
                                              big.mark = ",",
                                              scientific = F))),
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                   main = list(label="CHH", 
                               fontsize=15,
                               font=1,
                               vjust=0.8, 
                               fontfamily="Arial"))


# rice rep2, cg, chg, chh
m_rice_cg2 <- euler(c("bisulfite"= 233818, 
                     "Nanopore" = 2177632,
                     "bisulfite&Nanopore"=28308883))
p_rice_cg2 <- plot(m_rice_cg2, 
                  quantities = list(fontsize=15, fontfamily="Arial", 
                                    labels=trimws(format(c(233818, 
                                             2177632, 
                                             28308883), 
                                             big.mark = ",",
                                             scientific = F))),
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                  main = list(label="CpG", 
                              fontsize=15,
                              font=1,
                              vjust=0.8, 
                              fontfamily="Arial"))

m_rice_chg2 <- euler(c("bisulfite"= 121819, 
                      "Nanopore" = 1595249,
                      "bisulfite&Nanopore"=25600473))
p_rice_chg2 <- plot(m_rice_chg2, 
                   quantities = list(fontsize=15, fontfamily="Arial", 
                                     labels=trimws(format(c(121819, 
                                              1595249, 
                                              25600473), 
                                              big.mark = ",",
                                              scientific = F))), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                   main = list(label="CHG", 
                               fontsize=15,
                               font=1,
                               vjust=0.8, 
                               fontfamily="Arial"))

m_rice_chh2 <- euler(c("bisulfite"= 454817, 
                      "Nanopore" = 5659650,
                      "bisulfite&Nanopore"=98026457))
p_rice_chh2 <- plot(m_rice_chh2, 
                   quantities = list(fontsize=15, fontfamily="Arial", 
                                     labels=trimws(format(c(454817, 
                                              5659650, 
                                              98026457), 
                                              big.mark = ",",
                                              scientific = F))),
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=15, fontfamily="Arial"), 
                   main = list(label="CHH", 
                               fontsize=15,
                               font=1,
                               vjust=0.8, 
                               fontfamily="Arial"))




ppi= 300
png("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.cg_chg_chh.png", 
    width = 36, 
    height = 32, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab_cg,
                         p_arab_chg,
                         p_arab_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg,
                         p_rice_chg,
                         p_rice_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg2,
                         p_rice_chg2,
                         p_rice_chh2,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             heights = c(10, 1, 10, 1, 10))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.cg_chg_chh.svg", 
    width = 36/2.54, 
    height = 32/2.54)
grid.arrange(arrangeGrob(p_arab_cg,
                         p_arab_chg,
                         p_arab_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg,
                         p_rice_chg,
                         p_rice_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(grid.rect(gp=gpar(col="white")), 
                         nrow = 1), 
             arrangeGrob(p_rice_cg2,
                         p_rice_chg2,
                         p_rice_chh2,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             heights = c(10, 1, 10, 1, 10))
dev.off()








