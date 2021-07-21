#!/usr/bin/env bash

inbam=$1

file_prefix=$(echo ${inbam%.bam})
sortedbam="${file_prefix}.sorted.bam"
samtools sort -@ 30 -o $sortedbam $inbam
echo $inbam
samtools depth -a $sortedbam | awk '{sum+=$3} END { print "Average = ",sum/NR}'

rm $sortedbam
