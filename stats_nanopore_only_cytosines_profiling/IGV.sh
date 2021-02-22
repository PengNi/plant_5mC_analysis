#!/usr/bin/env bash

# bed2bedgraph
awk '{print $1"\t"$2"\t"$3"\t"$5 }' foo.bed > foo.cov.bedgraph
awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' foo.bed > foo.rmet.bedgraph
# strand specific
awk '{if($6=="+") print}' ${foo}.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > ${foo}.rmet.plus.bedgraph &
awk '{if($6=="-") print}' ${foo}.bed | awk '{print $1"\t"$2"\t"$3"\t"$11/100 }' - > ${foo}.rmet.minus.bedgraph &


# bisulfite coverage
samtools merge -@ 40 merged.bam 1.bam 2.bam 3.bam
samtools sort -@ 40 -o bismark.sorted.mark_dup.sorted.sorted.bam bismark.sorted.mark_dup.sorted.bam &
bedtools genomecov -bga -ibam bismark.sorted.mark_dup.sorted.sorted.bam > bismark.sorted.mark_dup.sorted.sorted.bam.bedgraph

# nanopore coverage
minimap2 GCF_000001735.4_TAIR10.1_genomic.fna arab.part2_50x_12345.guppy.fastq.fq -ax map-ont -t 40 > arab.part2_50x_12345.guppy.fastq.fq.minimap2.sam &
samtools view -bS -@ 40 arab.part2_50x_12345.guppy.fastq.fq.minimap2.sam > arab.part2_50x_12345.guppy.fastq.fq.minimap2.bam
samtools sort -@ 40 -O bam -T arab.part2_50x_12345.guppy.fastq.fq.minimap2.bam.tmp arab.part2_50x_12345.guppy.fastq.fq.minimap2.bam > arab.part2_50x_12345.guppy.fastq.fq.minimap2.sorted.bam
bedtools genomecov -bga -ibam arab.part2_50x_12345.guppy.fastq.fq.minimap2.sorted.bam > arab.part2_50x_12345.guppy.fastq.fq.minimap2.sorted.bam.bedgraph

# to TDF
igvtools toTDF -z 7 shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.rmet.covcf_5.rmet.minus.bedgraph shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.rmet.covcf_5.rmet.minus.bedgraph.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &
igvtools toTDF -z 7 shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.rmet.covcf_5.rmet.plus.bedgraph shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.rmet.covcf_5.rmet.plus.bedgraph.tdf /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa &











