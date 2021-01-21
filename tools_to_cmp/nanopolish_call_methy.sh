#!/usr/bin/env bash

read_f5=$1
read_fq_dir=$2
ref_fa=$3
nproc=$4
prefixname=$5

# 1. cat
read_fq_file="${prefixname}.fastq"
echo "cat ${read_fq_dir}/*.fastq > ${read_fq_file}"
cat ${read_fq_dir}/*.fastq > ${read_fq_file}

# 2. nanopolish index
echo "nanopolish index -d ${read_f5} ${read_fq_file}"
nanopolish index -d ${read_f5} ${read_fq_file}

# 3. minimap2
bam_tmp="${read_fq_file}.bam.tmp"
sorted_bam="${read_fq_file}.sorted.bam"
echo "minimap2 -a -x map-ont -t ${nproc} ${ref_fa} ${read_fq_file} | samtools view -@ ${nproc} -bS | samtools sort -@ ${nproc} -T ${bam_tmp} -o ${sorted_bam}"
minimap2 -a -x map-ont -t ${nproc} ${ref_fa} ${read_fq_file} | samtools view -@ ${nproc} -bS | samtools sort -@ ${nproc} -T ${bam_tmp} -o ${sorted_bam}
echo "samtools index -@ ${nproc} ${sorted_bam}"
samtools index -@ ${nproc} ${sorted_bam}


# 4. call methylation
methy_call_tsv="${read_fq_file}.methyl_calls.tsv"
echo "nanopolish call-methylation -t ${nproc} -r ${read_fq_file} -b ${sorted_bam} -g ${ref_fa} > ${methy_call_tsv}"
nanopolish call-methylation -t ${nproc} -r ${read_fq_file} -b ${sorted_bam} -g ${ref_fa} > ${methy_call_tsv}
echo "rm ${read_fq_file}"
echo "rm ${read_fq_file}.index*"
echo "rm ${sorted_bam}*"
rm ${read_fq_file}
rm ${read_fq_file}.index*
rm ${sorted_bam}
rm ${sorted_bam}.bai

# 5. methylation frequency
# scripts/calculate_methylation_frequency.py methylation_calls.tsv > methylation_frequency.tsv
