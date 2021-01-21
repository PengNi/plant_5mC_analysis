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
                 quantities = list(fontsize=16, fontfamily="serif", 
                                   labels=c(167689, 
                                            644273, 
                                            41951037)), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=16, fontfamily="serif"), 
                 main = list(label="", 
                             fontsize=24,
                             font=2,
                             hjust=8,
                             vjust=1, 
                             fontfamily="serif"))
# p_arab_c

# rice cytosine =
m_rice_c <- euler(c("bisulfite"= 796695, 
                    "Nanopore" = 10446128,
                    "bisulfite&Nanopore"=150936284))
p_rice_c <- plot(m_rice_c, 
                 quantities = list(fontsize=16, fontfamily="serif", 
                                   labels=c(796695, 
                                            10446128, 
                                            150936284)), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=16, fontfamily="serif"), 
                 main = list(label="", 
                             fontsize=24,
                             font=2,
                             hjust=8,
                             vjust=1, 
                             fontfamily="serif"))
# p_rice_c


# rice cytosine rep2 =
m_rice_c2 <- euler(c("bisulfite"= 810454, 
                    "Nanopore" = 9432531,
                    "bisulfite&Nanopore"=151935813))
p_rice_c2 <- plot(m_rice_c2, 
                 quantities = list(fontsize=16, fontfamily="serif", 
                                   labels=c(810454, 
                                            9432531, 
                                            151935813)), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=16, fontfamily="serif"), 
                 main = list(label="", 
                             fontsize=24,
                             font=2,
                             hjust=8,
                             vjust=1, 
                             fontfamily="serif"))
# p_rice_c2


# ppi= 300
# png("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.c.png", 
#      width = 23, 
#      height = 10, units = "cm", res=ppi)
# grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c,
#              widths = c(12, 1, 12))
# dev.off()



# arab cg, chg, chh
m_arab_cg <- euler(c("bisulfite"= 32636, 
                     "Nanopore" = 84684,
                     "bisulfite&Nanopore"=5436360))
p_arab_cg <- plot(m_arab_cg, 
                  quantities = list(fontsize=16, fontfamily="serif", 
                                    labels=c(32636, 
                                             84684, 
                                             5436360)), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=16, fontfamily="serif"), 
                  main = list(label="a", 
                              fontsize=24,
                              font=2,
                              hjust=8,
                              vjust=1, 
                              fontfamily="serif"))

m_arab_chg <- euler(c("bisulfite"= 18891, 
                      "Nanopore" = 85575,
                      "bisulfite&Nanopore"=5977439))
p_arab_chg <- plot(m_arab_chg, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(18891, 
                                              85575, 
                                              5977439)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="b", 
                               fontsize=24,
                               font=2,
                               hjust=8,
                               vjust=1, 
                               fontfamily="serif"))


m_arab_chh <- euler(c("bisulfite"= 116162, 
                      "Nanopore" = 474014,
                      "bisulfite&Nanopore"=30537238))
p_arab_chh <- plot(m_arab_chh, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(116162, 
                                              474014, 
                                              30537238)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="c", 
                               fontsize=24,
                               font=2,
                               hjust=8,
                               vjust=1, 
                               fontfamily="serif"))


# rice cg, chg, chh
m_rice_cg <- euler(c("bisulfite"= 224432, 
                     "Nanopore" = 2010752,
                     "bisulfite&Nanopore"=28488226))
p_rice_cg <- plot(m_rice_cg, 
                  quantities = list(fontsize=16, fontfamily="serif", 
                                    labels=c(224432, 
                                             2010752, 
                                             28488226)), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=16, fontfamily="serif"), 
                  main = list(label="d", 
                              fontsize=24,
                              font=2,
                              hjust=7.2,
                              vjust=1, 
                              fontfamily="serif"))

m_rice_chg <- euler(c("bisulfite"= 122935, 
                      "Nanopore" = 1552407,
                      "bisulfite&Nanopore"=25644528))
p_rice_chg <- plot(m_rice_chg, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(122935, 
                                              1552407, 
                                              25644528)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="e", 
                               fontsize=24,
                               font=2,
                               hjust=9.5,
                               vjust=1, 
                               fontfamily="serif"))

m_rice_chh <- euler(c("bisulfite"= 449328, 
                      "Nanopore" = 6882969,
                      "bisulfite&Nanopore"=96803530))
p_rice_chh <- plot(m_rice_chh, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(449328, 
                                              6882969, 
                                              96803530)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="f", 
                               fontsize=24,
                               font=2,
                               hjust=10.3,
                               vjust=1, 
                               fontfamily="serif"))


# rice rep2, cg, chg, chh
m_rice_cg2 <- euler(c("bisulfite"= 233818, 
                     "Nanopore" = 2177632,
                     "bisulfite&Nanopore"=28308883))
p_rice_cg2 <- plot(m_rice_cg2, 
                  quantities = list(fontsize=16, fontfamily="serif", 
                                    labels=c(233818, 
                                             2177632, 
                                             28308883)), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=16, fontfamily="serif"), 
                  main = list(label="g", 
                              fontsize=24,
                              font=2,
                              hjust=7.6,
                              vjust=1, 
                              fontfamily="serif"))

m_rice_chg2 <- euler(c("bisulfite"= 121819, 
                      "Nanopore" = 1595249,
                      "bisulfite&Nanopore"=25600473))
p_rice_chg2 <- plot(m_rice_chg2, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(121819, 
                                              1595249, 
                                              25600473)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="h", 
                               fontsize=24,
                               font=2,
                               hjust=7.5,
                               vjust=1, 
                               fontfamily="serif"))

m_rice_chh2 <- euler(c("bisulfite"= 454817, 
                      "Nanopore" = 5659650,
                      "bisulfite&Nanopore"=98026457))
p_rice_chh2 <- plot(m_rice_chh2, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(454817, 
                                              5659650, 
                                              98026457)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="i", 
                               fontsize=24,
                               font=2,
                               hjust=11.5,
                               vjust=1, 
                               fontfamily="serif"))




ppi= 300
png("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.cg_chg_chh.png", 
    width = 36, 
    height = 30, units = "cm", res=ppi)
grid.arrange(arrangeGrob(p_arab_cg,
                         p_arab_chg,
                         p_arab_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(p_rice_cg,
                         p_rice_chg,
                         p_rice_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(p_rice_cg2,
                         p_rice_chg2,
                         p_rice_chh2,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             heights = c(10, 10, 10))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.cg_chg_chh.svg", 
    width = 36/2.54, 
    height = 30/2.54)
grid.arrange(arrangeGrob(p_arab_cg,
                         p_arab_chg,
                         p_arab_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(p_rice_cg,
                         p_rice_chg,
                         p_rice_chh,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             arrangeGrob(p_rice_cg2,
                         p_rice_chg2,
                         p_rice_chh2,
                         nrow=1, 
                         widths = c(12, 12, 12)),
             heights = c(10, 10, 10))
dev.off()








