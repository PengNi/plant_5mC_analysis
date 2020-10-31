library(ggplot2)
library(scales)
library(plyr)

signalnumdf <- read.table("fig_deepsignal_hypertuning/athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.1mer.signalnum.txt", 
                          header = F, sep = '\t', stringsAsFactors = F)
colnames(signalnumdf) <- c("signalnum")
signalnumdf$count <- c(1:nrow(signalnumdf))
signalnum.mean <- mean(signalnumdf$signalnum)
signalnum.std <- sd(signalnumdf$signalnum)
signalnumdf$norm <- (signalnumdf$signalnum - signalnum.mean) / signalnum.std

signalnum.mean
signalnum.std
signalnum.mean + signalnum.std
signalnum.mean + 2*signalnum.std

signalnumdf_30 <- signalnumdf[signalnumdf$signalnum<=30,]
p <- ggplot(signalnumdf_30, 
            aes(x=signalnum)) + 
  geom_histogram(binwidth = 1, position = "dodge", colour="black", size=0.5, fill="red") + 
  geom_vline(xintercept = ceiling(signalnum.mean) + ceiling(signalnum.std), 
             linetype="dashed", size=1) +
  theme_bw() + 
  theme(text=element_text(size=12, family="serif"), 
        axis.text.x=element_text(size=7, family="serif")) + 
  scale_x_continuous(breaks=seq(0, 30, 1)) + 
  ggtitle('') + 
  labs(x="number of signals of one base")


ppi=300
png("fig_deepsignal_hypertuning/plot.number_of_signals_of_one_base.png", 
     width = 12, height = 9, units = "cm", res=ppi)
p
dev.off()


# svg
svg("fig_deepsignal_hypertuning/plot.number_of_signals_of_one_base.svg", 
     width = 12/2.54, height = 9/2.54)
p
dev.off()

