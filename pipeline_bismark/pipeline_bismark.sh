#!/usr/bin/env bash
# bismark_genome_preparation --bowtie2 genome.rice &

# get input params-----------------------------------
# pair-end reads, only for one group (2 fq) reads in .fastq format
read_fq1=$1
read_fq2=$2
ref_dir=$3
wk_dir=$4

mkdir $wk_dir

# 1. cutadapt
cutadapt_fq1="${wk_dir}/${read_fq1}.cut.fastq"
cutadapt_fq2="${wk_dir}/${read_fq2}.cut.fastq"
cutadapt_log="${wk_dir}/${read_fq1}.cut.log"
source ~/tools/cutadaptvenv/bin/activate
echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --max-n 0.1 -m 60 --trim-n -q 15,15 -j 20 -o $cutadapt_fq1 -p $cutadapt_fq2 $read_fq1 $read_fq2 > $cutadapt_log 2>&1"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --max-n 0.1 -m 60 --trim-n -q 15,15 -j 20 -o $cutadapt_fq1 -p $cutadapt_fq2 $read_fq1 $read_fq2 > $cutadapt_log 2>&1
deactivate

# 2. bismark
bismark_log="${wk_dir}/${read_fq1}.bismark.log"
# bismark_dir="${wk_dir}/bismark"
bismark_dir="${wk_dir}"
# --non_directional?
echo "bismark --parallel 20 -o $bismark_dir --temp_dir $bismark_dir --bowtie2 -N 1 -L 30 -most_valid_alignments 3 --genome $ref_dir -1 $cutadapt_fq1 -2 $cutadapt_fq2 > $bismark_log 2>&1"
bismark --parallel 20 -o $bismark_dir --temp_dir $bismark_dir --bowtie2 -N 1 -L 30 -most_valid_alignments 3 --genome $ref_dir -1 $cutadapt_fq1 -2 $cutadapt_fq2 > $bismark_log 2>&1

# 3. samtools sort / merge
readfq1_cut_prefix="${read_fq1}.cut"
bismark_bam="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.bam"
bismark_sbam="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.sorted.bam"
bismark_tbam="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.bam.tmp"
echo "samtools sort -o $bismark_sbam -T $bismark_tbam --threads 20 $bismark_bam"
samtools sort -o $bismark_sbam -T $bismark_tbam --threads 20 $bismark_bam
samtools index -@ 20 $bismark_sbam
rm $bismark_bam

# 4. picard / bam_sort
picard_log="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.sorted.bam.picard.log"
bismark_smbam="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.sorted.mark_dup.bam"
bismark_smmetrics="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.sorted.mark_dup.metrics.txt"
echo "java -Xmx200G -jar /home/nipeng/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$bismark_sbam O=$bismark_smbam M=$bismark_smmetrics >$picard_log 2>&1"
java -Xmx200G -jar /home/nipeng/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$bismark_sbam O=$bismark_smbam M=$bismark_smmetrics >$picard_log 2>&1

bismark_smsbam="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.sorted.mark_dup.sorted.bam"
bismark_smtbam="${bismark_dir}/${readfq1_cut_prefix}_bismark_bt2_pe.sorted.mark_dup.bam.tmp"
echo "samtools sort -n -o $bismark_smsbam -T $bismark_smtbam --threads 20 $bismark_smbam"
samtools sort -n -o $bismark_smsbam -T $bismark_smtbam --threads 20 $bismark_smbam
rm $bismark_smbam

# 5. bismark_methylation_extractor
# (-p --no_overlap --CX --ignore 4 --ignore_r2 4 --bedGraph --buffer_size 40G --cytosine_report)
ref_full_path=$(realpath $ref_dir)
bismark_full_smsbam=$(realpath $bismark_smsbam)
bismark_extract_log="${bismark_full_smsbam}.bme.log"
echo "cd $bismark_dir && bismark_methylation_extractor -p --no_overlap --ignore 4 --ignore_r2 4 --comprehensive --parallel 22 --bedGraph --cytosine_report --zero_based --CX --buffer_size 300G --genome_folder $ref_full_path $bismark_full_smsbam >$bismark_extract_log 2>&1"
cd $bismark_dir && bismark_methylation_extractor -p --no_overlap --ignore 4 --ignore_r2 4 --comprehensive --parallel 22 --bedGraph --cytosine_report --zero_based --CX --buffer_size 300G --genome_folder $ref_full_path $bismark_full_smsbam >$bismark_extract_log 2>&1
