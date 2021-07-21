#!/usr/bin/env bash



# arab, 254.2 ===
# =tsv2bed
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.freq.tsv --covcf 5 --wfile athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed --conv_chrom --sort &
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHG.freq.tsv --covcf 5 --wfile athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHG.cov5.rmet.bed --conv_chrom --sort &
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.freq.tsv --covcf 5 --wfile athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed --conv_chrom --sort &
# rmetbed2bedgraph, strand specific
awk '{if($6=="+") print}' athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.bedgraph &
awk '{if($6=="+") print}' athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHG.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHG.cov5.rmet.minus.bedgraph &
awk '{if($6=="+") print}' athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > athaliana.pass2.C.bn13_sn16_arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.bedgraph &






# rice2-1, 254.3 ===
# =tsv2bed
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.freq.tsv --covcf 5 --wfile shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed --sort &
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.freq.tsv --covcf 5 --wfile shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.bed --sort &
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.freq.tsv --covcf 5 --wfile shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed --sort &
# rmetbed2bedgraph, strand specific
awk '{if($6=="+") print}' shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.bedgraph &
awk '{if($6=="+") print}' shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.minus.bedgraph &
awk '{if($6=="+") print}' shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.bedgraph &
# bedgraph2tdf
igvtools toTDF -z 7 shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.bedgraph shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.bedgraph shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.plus.bedgraph shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.plus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.minus.bedgraph shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.minus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.bedgraph shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.bedgraph shuidao2-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &







# rice1-1, 254.3 ===
# =tsv2bed
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.freq.tsv --covcf 5 --wfile shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed --sort &
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.freq.tsv --covcf 5 --wfile shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.bed --sort &
python3 /homeb/nipeng/tools/plant_5mC_analysis/stats_correlation_with_bisulfite/convert_dp2_rmet2bedmethyl.py --freqfile shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.freq.tsv --covcf 5 --wfile shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed --sort &
# rmetbed2bedgraph, strand specific
awk '{if($6=="+") print}' shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.bedgraph &
awk '{if($6=="+") print}' shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.minus.bedgraph &
awk '{if($6=="+") print}' shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.bedgraph &
awk '{if($6=="-") print}' shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.bedgraph &
# bedgraph2tdf
igvtools toTDF -z 7 shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.bedgraph shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.plus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.bedgraph shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CG.cov5.rmet.minus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.plus.bedgraph shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.plus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.minus.bedgraph shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.0.CHG.cov5.rmet.minus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.bedgraph shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.plus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.bedgraph shuidao1-1.guppy.pass.part2.C.bn13_sn16.arabnrice2-1_120m.CNN.both_bilstm.50x_12345.call_mods.p0.8.CHH.cov5.rmet.minus.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &















