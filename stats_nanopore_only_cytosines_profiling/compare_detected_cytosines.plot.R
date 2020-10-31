library(ggplot2)
library(reshape2)
library(eulerr)
library(extrafont)
library(dplyr)
library(grid)
library(gridExtra)


# arab cytosine ==
m_arab_c <- euler(c("bisulfite"= 134141, 
                         "nanopore" = 654737,
                         "bisulfite&nanopore"=41984585))
p_arab_c <- plot(m_arab_c, 
                 quantities = list(fontsize=16, fontfamily="serif", 
                                   labels=c(134141, 
                                            654737, 
                                            41984585)), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=16, fontfamily="serif"), 
                 main = list(label="a", 
                             fontsize=24,
                             font=2,
                             hjust=8,
                             vjust=1, 
                             fontfamily="serif"))
# p_arab_c

# rice cytosine =
m_rice_c <- euler(c("bisulfite"= 652667, 
                    "nanopore" = 6920145,
                    "bisulfite&nanopore"=154696947))
p_rice_c <- plot(m_rice_c, 
                 quantities = list(fontsize=16, fontfamily="serif", 
                                   labels=c(652667, 
                                            6920145, 
                                            154696947)), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=16, fontfamily="serif"), 
                 main = list(label="b", 
                             fontsize=24,
                             font=2,
                             hjust=8,
                             vjust=1, 
                             fontfamily="serif"))
# p_rice_c

# ppi= 300
# png("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.c.png", 
#      width = 23, 
#      height = 10, units = "cm", res=ppi)
# grid.arrange(p_arab_c, grid.rect(gp=gpar(col="white")), p_rice_c,
#              widths = c(12, 1, 12))
# dev.off()



# arab cg, chg, chh
m_arab_cg <- euler(c("bisulfite"= 17464, 
                    "nanopore" = 87720,
                    "bisulfite&nanopore"=5451532))
p_arab_cg <- plot(m_arab_cg, 
                 quantities = list(fontsize=16, fontfamily="serif", 
                                   labels=c(17464, 
                                            87720, 
                                            5451532)), 
                 fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                 edges =TRUE,
                 labels = list(font=2, fontsize=16, fontfamily="serif"), 
                 main = list(label="a", 
                             fontsize=24,
                             font=2,
                             hjust=8,
                             vjust=1, 
                             fontfamily="serif"))

m_arab_chg <- euler(c("bisulfite"= 15548, 
                     "nanopore" = 86732,
                     "bisulfite&nanopore"=5980782))
p_arab_chg <- plot(m_arab_chg, 
                  quantities = list(fontsize=16, fontfamily="serif", 
                                    labels=c(15548, 
                                             86732, 
                                             5980782)), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=16, fontfamily="serif"), 
                  main = list(label="b", 
                              fontsize=24,
                              font=2,
                              hjust=8,
                              vjust=1, 
                              fontfamily="serif"))


m_arab_chh <- euler(c("bisulfite"= 101129, 
                      "nanopore" = 480285,
                      "bisulfite&nanopore"=30552271))
p_arab_chh <- plot(m_arab_chh, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(101129, 
                                              480285, 
                                              30552271)), 
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
m_rice_cg <- euler(c("bisulfite"= 141448, 
                     "nanopore" = 1515588,
                     "bisulfite&nanopore"=29098689))
p_rice_cg <- plot(m_rice_cg, 
                  quantities = list(fontsize=16, fontfamily="serif", 
                                    labels=c(141448, 
                                             1515588, 
                                             29098689)), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=16, fontfamily="serif"), 
                  main = list(label="d", 
                              fontsize=24,
                              font=2,
                              hjust=7.2,
                              vjust=1, 
                              fontfamily="serif"))

m_rice_chg <- euler(c("bisulfite"= 105322, 
                     "nanopore" = 1158737,
                     "bisulfite&nanopore"=26068103))
p_rice_chg <- plot(m_rice_chg, 
                  quantities = list(fontsize=16, fontfamily="serif", 
                                    labels=c(105322, 
                                             1158737, 
                                             26068103)), 
                  fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                  edges =TRUE,
                  labels = list(font=2, fontsize=16, fontfamily="serif"), 
                  main = list(label="e", 
                              fontsize=24,
                              font=2,
                              hjust=10,
                              vjust=1, 
                              fontfamily="serif"))

m_rice_chh <- euler(c("bisulfite"= 405897, 
                      "nanopore" = 4245820,
                      "bisulfite&nanopore"=99530155))
p_rice_chh <- plot(m_rice_chh, 
                   quantities = list(fontsize=16, fontfamily="serif", 
                                     labels=c(405897, 
                                              4245820, 
                                              99530155)), 
                   fills = list(fill=c('#fc8d62', '#66c2a5'), alpha=0.7),
                   edges =TRUE,
                   labels = list(font=2, fontsize=16, fontfamily="serif"), 
                   main = list(label="f", 
                               fontsize=24,
                               font=2,
                               hjust=10.5,
                               vjust=1, 
                               fontfamily="serif"))

ppi= 300
png("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.cg_chg_chh.png", 
     width = 36, 
     height = 20, units = "cm", res=ppi)
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
             heights = c(10, 10))
dev.off()
svg("stats_nanopore_only_cytosines_profiling/compare_detected_cytosines.venn.arab_rice.cg_chg_chh.svg", 
    width = 36/2.54, 
    height = 20/2.54)
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
             heights = c(10, 10))
dev.off()


