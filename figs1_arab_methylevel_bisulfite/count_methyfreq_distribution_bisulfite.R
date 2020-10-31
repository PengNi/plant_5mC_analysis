## #!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# arg1: bs rmet filepath
# arg2: replicate name, e.g.: replicate1
# arg3: motif, e.g.: CG

# log cmd:
# Rscript ~/tools/plant_5mC_analysis/figs1_arab_methylevel_bisulfite/count_methyfreq_distribution_bisulfite.R ninanjie-2-1_D1902826A-ZJ/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt replicate1 CG
# Rscript ~/tools/plant_5mC_analysis/figs1_arab_methylevel_bisulfite/count_methyfreq_distribution_bisulfite.R ninanjie-2-1_D1902826A-ZJ/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt replicate1 CHG &
# Rscript ~/tools/plant_5mC_analysis/figs1_arab_methylevel_bisulfite/count_methyfreq_distribution_bisulfite.R ninanjie-2-1_D1902826A-ZJ/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt replicate1 CHH &
# Rscript ~/tools/plant_5mC_analysis/figs1_arab_methylevel_bisulfite/count_methyfreq_distribution_bisulfite.R ninanjie-2-2_D1902827A-ZJ/D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt replicate2 CG &
# Rscript ~/tools/plant_5mC_analysis/figs1_arab_methylevel_bisulfite/count_methyfreq_distribution_bisulfite.R ninanjie-2-3_D1902828A-ZJ/D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt replicate3 CG &



# functions ==================
ratio_stats <- function(values, lb=0, ub=1, xbin=0.1){
  sranges <- seq(lb, ub, xbin)
  bincounts <- vector('numeric', length=length(sranges)-1)
  bincounts[1] <- sum(sranges[1] <= values & values <= sranges[2])
  for(i in 2:(length(sranges)-1)){
    bincounts[i] <- sum(sranges[i] < values & values <= sranges[i+1])
  }
  x <- data.frame(idx=c(1:length(bincounts)), 
                  ranges=paste(sranges[1:(length(sranges)-1)], 
                               sranges[2:length(sranges)], 
                               sep = "-"), 
                  counts=bincounts, 
                  ratio=bincounts/length(values))
}

# pipeline ===================
covcf = 5
width = 0.1

bs_rmet_file <- args[1]
be_rmet <- read.table(bs_rmet_file, 
                      header = T, sep = "\t", stringsAsFactors = F)
# coverage_bis <- bisulfite_data[bisulfite_data$coverage<=cov_cf,]$coverage
be_rmet$key <- paste(be_rmet$chromosome, be_rmet$pos, sep=" ")

be_rmet <- be_rmet[be_rmet$coverage>=covcf, 
                   c("key", "Rmet")]
rmet_bis <- be_rmet$Rmet
rmet_bis[rmet_bis<0] = 0
rmet_bis <- ratio_stats(rmet_bis, 0, 1, width)
rmet_bis$replicate <- args[2]
rmet_bis$motif <- args[3]

wfile <- paste(bs_rmet_file, ".distr_plot.tsv", sep="")
write.table(rmet_bis, file=wfile, quote = F, row.names = F, col.names = T, sep = "\t")






