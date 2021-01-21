#!/usr/bin/env bash

python extract_contig_length.py ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.seq_len.txt
sort -k1,1V ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.seq_len.txt > ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.seq_len.sorted.txt



# =============================================================================
# 1. train CG (megalodon >=2.1.1)
# arab
export OMP_NUM_THREADS=1
megalodon_extras modified_bases create_ground_truth --bed-methyl-files athaliana.bsrep_123_comb.CG.methyl.bed \
--coverage-threshold 5 --pct-mod-thresholds 0 90 --out-csv athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.csv
megalodon ../reads_single_pass/part1 --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 60 \
    --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --mod-motif m CG 0 --outputs per_read_mods mappings --devices 0,1 --processes 25 \
    --output-directory megalodon_results.arab.pass.part1_guppy.CG --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CG.log 2>&1
awk '$2 > 90 && $7 > 90 && $3 - $6 > 5000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.arab.pass.part1_guppy.CG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.arab.pass.part1_guppy.CG/mappings.filtered.sorted.bed
cat megalodon_results.arab.pass.part1_guppy.CG/mappings.filtered.sorted.bed | \
    awk '{print $4}' > megalodon_results.arab.pass.part1_guppy.CG/train_read_ids.txt
megalodon_extras calibrate generate_modified_base_stats \
    --ground-truth-data athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.csv \
    --strand-specific-sites --out-filename athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    megalodon_results.arab.pass.part1_guppy.CG
megalodon_extras calibrate modified_bases \
    --ground-truth-llrs athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    --out-filename athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --out-pdf athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.megalodon_mod_calibration.pdf \
    --processes 30 \
    > megalodon_results.arab.pass.part1_guppy.CG.bsrep_123_comb.methyl.cov5.0_90.megalodon_mod_calibration.log
# megalodon_extras modified_bases estimate_threshold megalodon_results.arab.pass.part1_guppy.CG m
megalodon ../reads_single_pass/part1 --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 60 \
    --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --devices 0,1 --processes 25 \
    --mod-motif m CG 0 \
    --ref-include-mods \
    --mod-calibration-filename athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --read-ids-filename megalodon_results.arab.pass.part1_guppy.CG/train_read_ids.txt \
    --output-directory megalodon_results.arab.pass.part1_guppy.CG.retrain --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CG.retrain.log 2>&1
export OMP_NUM_THREADS=25
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.arab.pass.part1_guppy.CG.retrain/signal_mappings.hdf5 --device 0 \
    --outdir megalodon_results.arab.pass.part1_guppy.CG.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.arab.pass.part1_guppy.CG.retrain.taiyaki_train.log
dump_json.py megalodon_results.arab.pass.part1_guppy.CG.retrain/model_final.checkpoint \
    --output megalodon_results.arab.pass.part1_guppy.CG.retrain/model_final.jsn
# rice
export OMP_NUM_THREADS=1
megalodon_extras modified_bases create_ground_truth --bed-methyl-files shuidao2-1.bs.CG.methyl.bed \
--coverage-threshold 5 --pct-mod-thresholds 0 90 --out-csv shuidao2-1.bs.CG.methyl.cov5.0_90.csv
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain0/megalodon_results.arab.pass.part1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --mod-motif m CG 0 --outputs per_read_mods mappings --devices 0,1 --processes 25 \
    --mod-calibration-filename /home/nipeng/data/arab.nano/megalodon_retrain0/athaliana.bsrep_123_comb.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --output-directory megalodon_results.rice.pass1_guppy.CG --overwrite \
    > megalodon_results.rice.pass1_guppy.CG.log 2>&1
awk '$2 > 90 && $7 > 90 && $3 - $6 > 5000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.rice.pass1_guppy.CG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.rice.pass1_guppy.CG/mappings.filtered.sorted.bed
cat megalodon_results.rice.pass1_guppy.CG/mappings.filtered.sorted.bed | \
    awk '{print $4}' > megalodon_results.rice.pass1_guppy.CG/train_read_ids.txt
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .4) print $0}' megalodon_results.rice.pass1_guppy.CG/train_read_ids.txt > \
    megalodon_results.rice.pass1_guppy.CG/train_read_ids.f.txt
megalodon_extras calibrate generate_modified_base_stats \
    --ground-truth-data shuidao2-1.bs.CG.methyl.cov5.0_90.csv \
    --strand-specific-sites --out-filename shuidao2-1.bs.CG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    megalodon_results.rice.pass1_guppy.CG
megalodon_extras calibrate modified_bases \
    --ground-truth-llrs shuidao2-1.bs.CG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    --out-filename shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --out-pdf shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.pdf \
    --processes 30 \
    > megalodon_results.rice.pass1_guppy.CG.bs.methyl.cov5.0_90.megalodon_mod_calibration.log
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain0/megalodon_results.arab.pass.part1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --devices 0,1 --processes 25 \
    --mod-motif m CG 0 \
    --ref-include-mods \
    --mod-calibration-filename shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --read-ids-filename megalodon_results.rice.pass1_guppy.CG/train_read_ids.txt \
    --output-directory megalodon_results.rice.pass1_guppy.CG.retrain --overwrite \
    > megalodon_results.rice.pass1_guppy.CG.retrain.log 2>&1
export OMP_NUM_THREADS=30
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.rice.pass1_guppy.CG.retrain/signal_mappings.hdf5 --device 0 \
    --outdir megalodon_results.rice.pass1_guppy.CG.retrain --overwrite \
    --niteration 20000 --save_every 400 --sub_batches 1 \
    > megalodon_results.rice.pass1_guppy.CG.retrain.taiyaki_train.log
dump_json.py megalodon_results.rice.pass1_guppy.CG.retrain/model_final.checkpoint \
    --output megalodon_results.rice.pass1_guppy.CG.retrain/model_final.jsn
# test =
# arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_2.CG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_2.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_3.CG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_3.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_4.CG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_4.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_5.CG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_5.CG_retrain_comb.log 2>&1
# rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass2.single/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CG.retrain/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_1.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CG.retrain/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_2.CG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_2.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CG.retrain/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_3.CG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_3.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CG.retrain/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_4.CG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_4.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CG.retrain/shuidao2-1.bs.CG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_5.CG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_5.CG_retrain_comb.log 2>&1


# 2. train CHG (megalodon >=2.1.1)
export OMP_NUM_THREADS=1
megalodon_extras modified_bases create_ground_truth --bed-methyl-files athaliana.bsrep_123_comb.CHG.methyl.bed \
--coverage-threshold 5 --pct-mod-thresholds 0 90 --out-csv athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.csv
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 60 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --mod-motif Z CHG 0 --outputs per_read_mods mappings --devices 2,3 --processes 25 \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHG --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHG.log 2>&1
awk '$2 > 90 && $7 > 90 && $3 - $6 > 5000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.arab.pass.part1_guppy.CHG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.arab.pass.part1_guppy.CHG/mappings.filtered.sorted.bed
cat megalodon_results.arab.pass.part1_guppy.CHG/mappings.filtered.sorted.bed | \
    awk '{print $4}' > megalodon_results.arab.pass.part1_guppy.CHG/train_read_ids.txt
megalodon_extras calibrate generate_modified_base_stats \
    --ground-truth-data athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.csv \
    --strand-specific-sites --out-filename athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    megalodon_results.arab.pass.part1_guppy.CHG
megalodon_extras calibrate modified_bases \
    --ground-truth-llrs athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    --out-filename athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --processes 30
megalodon ../reads_single_pass/part1 --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 60 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --devices 2,3 --processes 25 \
    --mod-motif Z CHG 0 \
    --ref-include-mods \
    --mod-calibration-filename athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --read-ids-filename megalodon_results.arab.pass.part1_guppy.CHG/train_read_ids.txt \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHG.retrain --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHG.retrain.log 2>&1
export OMP_NUM_THREADS=25
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.arab.pass.part1_guppy.CHG.retrain/signal_mappings.hdf5 \
    --device 2 --outdir megalodon_results.arab.pass.part1_guppy.CHG.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.arab.pass.part1_guppy.CHG.retrain.taiyaki_train.log
dump_json.py megalodon_results.arab.pass.part1_guppy.CHG.retrain/model_final.checkpoint \
    --output megalodon_results.arab.pass.part1_guppy.CHG.retrain/model_final.jsn
# rice
export OMP_NUM_THREADS=1
megalodon_extras modified_bases create_ground_truth --bed-methyl-files shuidao2-1.bs.CHG.methyl.bed \
--coverage-threshold 5 --pct-mod-thresholds 0 90 --out-csv shuidao2-1.bs.CHG.methyl.cov5.0_90.csv
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain0/megalodon_results.arab.pass.part1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --mod-motif Z CHG 0 --outputs per_read_mods mappings --devices 2,3 --processes 25 \
    --mod-calibration-filename /home/nipeng/data/arab.nano/megalodon_retrain0/athaliana.bsrep_123_comb.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --output-directory megalodon_results.rice.pass1_guppy.CHG --overwrite \
    > megalodon_results.rice.pass1_guppy.CHG.log 2>&1
awk '$2 > 90 && $7 > 90 && $3 - $6 > 5000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.rice.pass1_guppy.CHG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.rice.pass1_guppy.CHG/mappings.filtered.sorted.bed
cat megalodon_results.rice.pass1_guppy.CHG/mappings.filtered.sorted.bed | \
    awk '{print $4}' > megalodon_results.rice.pass1_guppy.CHG/train_read_ids.txt
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .4) print $0}' megalodon_results.rice.pass1_guppy.CHG/train_read_ids.txt > \
    megalodon_results.rice.pass1_guppy.CHG/train_read_ids.f.txt
megalodon_extras calibrate generate_modified_base_stats \
    --ground-truth-data shuidao2-1.bs.CHG.methyl.cov5.0_90.csv \
    --strand-specific-sites --out-filename shuidao2-1.bs.CHG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    megalodon_results.rice.pass1_guppy.CHG
megalodon_extras calibrate modified_bases \
    --ground-truth-llrs shuidao2-1.bs.CHG.methyl.cov5.0_90.mod_calibration_statistics.npz \
    --out-filename shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --out-pdf shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.pdf \
    --processes 30 \
    > megalodon_results.rice.pass1_guppy.CHG.bs.methyl.cov5.0_90.megalodon_mod_calibration.log
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain0/megalodon_results.arab.pass.part1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --devices 2,3 --processes 25 \
    --mod-motif Z CHG 0 \
    --ref-include-mods \
    --mod-calibration-filename shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --read-ids-filename megalodon_results.rice.pass1_guppy.CHG/train_read_ids.f.txt \
    --output-directory megalodon_results.rice.pass1_guppy.CHG.retrain --overwrite \
    > megalodon_results.rice.pass1_guppy.CHG.retrain.log 2>&1
export OMP_NUM_THREADS=30
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.rice.pass1_guppy.CHG.retrain/signal_mappings.hdf5 --device 2 \
    --outdir megalodon_results.rice.pass1_guppy.CHG.retrain --overwrite \
    --niteration 20000 --save_every 400 --sub_batches 1 \
    > megalodon_results.rice.pass1_guppy.CHG.retrain.taiyaki_train.log
dump_json.py megalodon_results.rice.pass1_guppy.CHG.retrain/model_final.checkpoint \
    --output megalodon_results.rice.pass1_guppy.CHG.retrain/model_final.jsn
# test =
# arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_2.CHG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_2.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_3.CHG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_3.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_4.CHG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_4.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_5.CHG_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_5.CHG_retrain_comb.log 2>&1
# rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass2.single/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHG.retrain/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHG.retrain/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_2.CHG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_2.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHG.retrain/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_3.CHG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_3.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHG.retrain/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_4.CHG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_4.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHG.retrain/shuidao2-1.bs.CHG.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_5.CHG_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_5.CHG_retrain_comb.log 2>&1


# 3. train CHH (megalodon >=2.1.1)
# arab
export OMP_NUM_THREADS=1
megalodon_extras modified_bases create_ground_truth --bed-methyl-files athaliana.bsrep_123_comb.CHH.methyl.bed \
--coverage-threshold 5 --pct-mod-thresholds 0 90 --out-csv athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.csv
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 60 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --mod-motif Z CHH 0 --outputs per_read_mods mappings --devices 1,3 --processes 25 \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHH --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHH.log 2>&1
awk '$2 > 90 && $7 > 90 && $3 - $6 > 5000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.arab.pass.part1_guppy.CHH/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.arab.pass.part1_guppy.CHH/mappings.filtered.sorted.bed
cat megalodon_results.arab.pass.part1_guppy.CHH/mappings.filtered.sorted.bed | \
    awk '{print $4}' > megalodon_results.arab.pass.part1_guppy.CHH/train_read_ids.txt
megalodon_extras calibrate generate_modified_base_stats \
    --ground-truth-data athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.csv \
    --strand-specific-sites --out-filename athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.mod_calibration_statistics.npz \
    megalodon_results.arab.pass.part1_guppy.CHH
megalodon_extras calibrate modified_bases \
    --ground-truth-llrs athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.mod_calibration_statistics.npz \
    --out-filename athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --processes 30 \
     > megalodon_results.arab.pass.part1_guppy.CHH.bsrep_123_comb.methyl.cov5.0_90.megalodon_mod_calibration.log
megalodon ../reads_single_pass/part1 --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 60 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --devices 1,3 --processes 25 \
    --mod-motif Z CHH 0 \
    --ref-include-mods \
    --mod-calibration-filename athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --read-ids-filename megalodon_results.arab.pass.part1_guppy.CHH/train_read_ids.txt \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHH.retrain --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHH.retrain.log 2>&1
export OMP_NUM_THREADS=25
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.arab.pass.part1_guppy.CHH.retrain/signal_mappings.hdf5 \
    --device 3 --outdir megalodon_results.arab.pass.part1_guppy.CHH.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.arab.pass.part1_guppy.CHH.retrain.taiyaki_train.log
dump_json.py megalodon_results.arab.pass.part1_guppy.CHH.retrain/model_final.checkpoint \
    --output megalodon_results.arab.pass.part1_guppy.CHH.retrain/model_final.jsn
# rice
export OMP_NUM_THREADS=1
megalodon_extras modified_bases create_ground_truth --bed-methyl-files shuidao2-1.bs.CHH.methyl.bed \
--coverage-threshold 5 --pct-mod-thresholds 0 90 --out-csv shuidao2-1.bs.CHH.methyl.cov5.0_90.csv
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain0/megalodon_results.arab.pass.part1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --mod-motif Z CHH 0 --outputs per_read_mods mappings --devices 0,1 --processes 25 \
    --mod-calibration-filename /home/nipeng/data/arab.nano/megalodon_retrain0/athaliana.bsrep_123_comb.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --output-directory megalodon_results.rice.pass1_guppy.CHH --overwrite \
    > megalodon_results.rice.pass1_guppy.CHH.log 2>&1
awk '$2 > 90 && $7 > 90 && $3 - $6 > 5000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.rice.pass1_guppy.CHH/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.rice.pass1_guppy.CHH/mappings.filtered.sorted.bed
cat megalodon_results.rice.pass1_guppy.CHH/mappings.filtered.sorted.bed | \
    awk '{print $4}' > megalodon_results.rice.pass1_guppy.CHH/train_read_ids.txt
awk 'BEGIN {srand()} !/^$/ { if (rand() <= .4) print $0}' megalodon_results.rice.pass1_guppy.CHH/train_read_ids.txt > \
    megalodon_results.rice.pass1_guppy.CHH/train_read_ids.f.txt
megalodon_extras calibrate generate_modified_base_stats \
    --ground-truth-data shuidao2-1.bs.CHH.methyl.cov5.0_90.csv \
    --strand-specific-sites --out-filename shuidao2-1.bs.CHH.methyl.cov5.0_90.mod_calibration_statistics.npz \
    megalodon_results.rice.pass1_guppy.CHH
megalodon_extras calibrate modified_bases \
    --ground-truth-llrs shuidao2-1.bs.CHH.methyl.cov5.0_90.mod_calibration_statistics.npz \
    --out-filename shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --out-pdf shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.pdf \
    --processes 30 \
    > megalodon_results.rice.pass1_guppy.CHH.bs.methyl.cov5.0_90.megalodon_mod_calibration.log
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain0/megalodon_results.arab.pass.part1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --devices 0,1 --processes 25 \
    --mod-motif Z CHH 0 \
    --ref-include-mods \
    --mod-calibration-filename shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --read-ids-filename megalodon_results.rice.pass1_guppy.CHH/train_read_ids.f.txt \
    --output-directory megalodon_results.rice.pass1_guppy.CHH.retrain --overwrite \
    > megalodon_results.rice.pass1_guppy.CHH.retrain.log 2>&1
export OMP_NUM_THREADS=30
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.rice.pass1_guppy.CHH.retrain/signal_mappings.hdf5 --device 2 \
    --outdir megalodon_results.rice.pass1_guppy.CHH.retrain --overwrite \
    --niteration 20000 --save_every 400 --sub_batches 1 \
    > megalodon_results.rice.pass1_guppy.CHH.retrain.taiyaki_train.log
dump_json.py megalodon_results.rice.pass1_guppy.CHH.retrain/model_final.checkpoint \
    --output megalodon_results.rice.pass1_guppy.CHH.retrain/model_final.jsn
# test =
# arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_2.CHH_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_2.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_3.CHH_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_3.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_4.CHH_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_4.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --mod-calibration-filename /homeb/nipeng/data/shuidao2-1/megalodon_retrain0/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.arab.pass.part2_guppy.10x_5.CHH_retrain_comb \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_5.CHH_retrain_comb.log 2>&1
# rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass2.single/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHH.retrain/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CHH_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_1.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHH.retrain/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_2.CHH_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_2.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHH.retrain/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_3.CHH_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_3.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHH.retrain/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_4.CHH_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_4.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain0/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --mod-calibration-filename megalodon_results.rice.pass1_guppy.CHH.retrain/shuidao2-1.bs.CHH.methyl.cov5.0_90.megalodon_mod_calibration.npz \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_5.CHH_retrain_comb \
    --overwrite > megalodon_results.rice.pass2_guppy.10x_5.CHH_retrain_comb.log 2>&1












# =============================================================================
# 1. train CG (megalodon 2.2.3)
# arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --mod-motif m CG 0 --outputs per_read_mods mappings --devices 0,1 --processes 30 \
    --output-directory megalodon_results.arab.pass.part1_guppy.CG --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CG.log 2>&1
megalodon_extras \
    modified_bases per_site_thresholds \
    megalodon_results.arab.pass.part1_guppy.CG \
    athaliana.bsrep_123_comb.CG.methyl.bed \
    --mod-base m \
    --ground-truth-cov-min 25 \
    --nanopore-cov-min 30 \
    --ground-truth-coverage-pdf megalodon_results.arab.pass.part1_guppy.CG/gt_cov.CG.pdf \
    --out-blacklist-sites megalodon_results.arab.pass.part1_guppy.CG/low_coverage_sites.CG.bed \
    --out-per-site-mod-thresholds megalodon_results.arab.pass.part1_guppy.CG/site_mod_thresholds.CG.bed
sort -S 25% --parallel=30 -T /tmp/ \
    -k1,1V -k2,2n \
    -o megalodon_results.arab.pass.part1_guppy.CG/low_coverage_sites.CG.sorted.bed \
    megalodon_results.arab.pass.part1_guppy.CG/low_coverage_sites.CG.bed
awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.arab.pass.part1_guppy.CG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.arab.pass.part1_guppy.CG/mappings.filtered.sorted.bed
bedtools intersect \
    -a megalodon_results.arab.pass.part1_guppy.CG/mappings.filtered.sorted.bed \
    -b megalodon_results.arab.pass.part1_guppy.CG/low_coverage_sites.CG.sorted.bed -s -sorted -v | \
    awk '{print $4}' > megalodon_results.arab.pass.part1_guppy.CG/train_read_ids.txt
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --devices 0,1 --processes 30 \
    --mod-motif m CG 0 \
    --ref-include-mods \
    --mod-per-site-threshold megalodon_results.arab.pass.part1_guppy.CG/site_mod_thresholds.CG.bed \
    --read-ids-filename megalodon_results.arab.pass.part1_guppy.CG/train_read_ids.txt \
    --output-directory megalodon_results.arab.pass.part1_guppy.CG.retrain --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CG.retrain.log 2>&1

export OMP_NUM_THREADS=12
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.arab.pass.part1_guppy.CG.retrain/signal_mappings.hdf5 \
    --device 2 --outdir megalodon_results.arab.pass.part1_guppy.CG.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.arab.pass.part1_guppy.CG.retrain.taiyaki_train.log &
dump_json.py megalodon_results.arab.pass.part1_guppy.CG.retrain/model_final.checkpoint \
    --output megalodon_results.arab.pass.part1_guppy.CG.retrain/model_final.jsn
# test
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_arab \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_arab.log 2>&1 &
Rscript ../correlation_with_bs.cal_plot.general.R \
    ../bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt \
    megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_arab/modified_bases.5mC.bed \
    bisulfite.rep1 megalodon_retrain_arab.guppy.arab.10x.1 CG yes \
    megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_arab \
    CG_megalodon_retrain_arab.guppy.arab.10x.1_vs_bisulfite.rep1.tsv \
    > megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_arab/CG_megalodon_retrain_arab.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# rice =
export OMP_NUM_THREADS=1
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CG.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --mod-motif m CG 0 --outputs per_read_mods mappings --devices 2,3 --processes 20 \
    --output-directory megalodon_results.rice.pass1_guppy.CG --overwrite --disable-mod-calibration \
    > megalodon_results.rice.pass1_guppy.CG.log 2>&1
megalodon_extras \
    modified_bases per_site_thresholds \
    megalodon_results.rice.pass1_guppy.CG \
    shuidao2-1.bs.CG.methyl.bed \
    --mod-base m \
    --ground-truth-cov-min 25 \
    --nanopore-cov-min 30 \
    --ground-truth-coverage-pdf megalodon_results.rice.pass1_guppy.CG/gt_cov.CG.pdf \
    --out-blacklist-sites megalodon_results.rice.pass1_guppy.CG/low_coverage_sites.CG.bed \
    --out-per-site-mod-thresholds megalodon_results.rice.pass1_guppy.CG/site_mod_thresholds.CG.bed \
    > megalodon_results.rice.pass1_guppy.CG.extra.log 2>&1
sort -S 25% --parallel=30 -T /tmp/ \
    -k1,1V -k2,2n \
    -o megalodon_results.rice.pass1_guppy.CG/low_coverage_sites.CG.sorted.bed \
    megalodon_results.rice.pass1_guppy.CG/low_coverage_sites.CG.bed
awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.rice.pass1_guppy.CG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.rice.pass1_guppy.CG/mappings.filtered.sorted.bed
bedtools intersect \
    -a megalodon_results.rice.pass1_guppy.CG/mappings.filtered.sorted.bed \
    -b megalodon_results.rice.pass1_guppy.CG/low_coverage_sites.CG.sorted.bed -s -sorted -v \
    -g ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.seq_len.sorted.txt | \
    awk '{print $4}' > megalodon_results.rice.pass1_guppy.CG/train_read_ids.txt
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CG.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --devices 2,3 --processes 20 \
    --mod-motif m CG 0 \
    --ref-include-mods \
    --mod-per-site-threshold megalodon_results.rice.pass1_guppy.CG/site_mod_thresholds.CG.bed \
    --read-ids-filename megalodon_results.rice.pass1_guppy.CG/train_read_ids.txt \
    --output-directory megalodon_results.rice.pass1_guppy.CG.retrain --overwrite --disable-mod-calibration \
    > megalodon_results.rice.pass1_guppy.CG.retrain.log 2>&1
export OMP_NUM_THREADS=25
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.rice.pass1_guppy.CG.retrain/signal_mappings.hdf5 \
    --device 1 --outdir megalodon_results.rice.pass1_guppy.CG.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.rice.pass1_guppy.CG.retrain.taiyaki_train.log
dump_json.py megalodon_results.rice.pass1_guppy.CG.retrain/model_final.checkpoint \
    --output megalodon_results.rice.pass1_guppy.CG.retrain/model_final.jsn
# test =====
# arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_1.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_2.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_2.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_3.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_3.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_4.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_4.CG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_5.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_5.CG_retrain_comb.log 2>&1
# rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass2.single/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_1.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_2.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_2.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_3.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_3.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_4.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_4.CG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 1 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_5.CG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_5.CG_retrain_comb.log 2>&1






# 2. train CHG (megalodon 2.2.3)
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --mod-motif Z CHG 0 --outputs per_read_mods mappings --devices 0,1 --processes 30 \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHG --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHG.log 2>&1
megalodon_extras \
    modified_bases per_site_thresholds \
    megalodon_results.arab.pass.part1_guppy.CHG \
    athaliana.bsrep_123_comb.CHG.methyl.bed \
    --mod-base Z \
    --ground-truth-cov-min 25 \
    --nanopore-cov-min 30 \
    --ground-truth-coverage-pdf megalodon_results.arab.pass.part1_guppy.CHG/gt_cov.CHG.pdf \
    --out-blacklist-sites megalodon_results.arab.pass.part1_guppy.CHG/low_coverage_sites.CHG.bed \
    --out-per-site-mod-thresholds megalodon_results.arab.pass.part1_guppy.CHG/site_mod_thresholds.CHG.bed
sort -S 25% --parallel=30 -T /tmp/ \
    -k1,1V -k2,2n \
    -o megalodon_results.arab.pass.part1_guppy.CHG/low_coverage_sites.CHG.sorted.bed \
    megalodon_results.arab.pass.part1_guppy.CHG/low_coverage_sites.CHG.bed
awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.arab.pass.part1_guppy.CHG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.arab.pass.part1_guppy.CHG/mappings.filtered.sorted.bed
bedtools intersect \
    -a megalodon_results.arab.pass.part1_guppy.CHG/mappings.filtered.sorted.bed \
    -b megalodon_results.arab.pass.part1_guppy.CHG/low_coverage_sites.CHG.sorted.bed -s -sorted -v | \
    awk '{print $4}' > megalodon_results.arab.pass.part1_guppy.CHG/train_read_ids.txt
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --devices 0,1 --processes 30 \
    --mod-motif Z CHG 0 \
    --ref-include-mods \
    --mod-per-site-threshold megalodon_results.arab.pass.part1_guppy.CHG/site_mod_thresholds.CHG.bed \
    --read-ids-filename megalodon_results.arab.pass.part1_guppy.CHG/train_read_ids.txt \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHG.retrain --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHG.retrain.log 2>&1

export OMP_NUM_THREADS=12
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.arab.pass.part1_guppy.CHG.retrain/signal_mappings.hdf5 \
    --device 3 --outdir megalodon_results.arab.pass.part1_guppy.CHG.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.arab.pass.part1_guppy.CHG.retrain.taiyaki_train.log &
# test
dump_json.py megalodon_results.arab.pass.part1_guppy.CHG.retrain/model_final.checkpoint \
    --output megalodon_results.arab.pass.part1_guppy.CHG.retrain/model_final.jsn
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_00041.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 1 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_arab \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_arab.log 2>&1 &
Rscript ../correlation_with_bs.cal_plot.general.R \
    ../bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt \
    megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_arab/modified_bases.5mC.bed \
    bisulfite.rep1 megalodon_retrain_arab.guppy.arab.10x.1 CHG yes \
    megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_arab \
    CHG_megalodon_retrain_arab.guppy.arab.10x.1_vs_bisulfite.rep1.tsv \
    > megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_arab/CHG_megalodon_retrain_arab.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# rice =
export OMP_NUM_THREADS=1
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CHG.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --mod-motif Z CHG 0 --outputs per_read_mods mappings --devices 0,1 --processes 30 \
    --output-directory megalodon_results.rice.pass1_guppy.CHG --overwrite --disable-mod-calibration \
    > megalodon_results.rice.pass1_guppy.CHG.log 2>&1
megalodon_extras \
    modified_bases per_site_thresholds \
    megalodon_results.rice.pass1_guppy.CHG \
    shuidao2-1.bs.CHG.methyl.bed \
    --mod-base Z \
    --ground-truth-cov-min 25 \
    --nanopore-cov-min 30 \
    --ground-truth-coverage-pdf megalodon_results.rice.pass1_guppy.CHG/gt_cov.CHG.pdf \
    --out-blacklist-sites megalodon_results.rice.pass1_guppy.CHG/low_coverage_sites.CHG.bed \
    --out-per-site-mod-thresholds megalodon_results.rice.pass1_guppy.CHG/site_mod_thresholds.CHG.bed \
    > megalodon_results.rice.pass1_guppy.CHG.extra.log 2>&1
sort -S 25% --parallel=30 -T /tmp/ \
    -k1,1V -k2,2n \
    -o megalodon_results.rice.pass1_guppy.CHG/low_coverage_sites.CHG.sorted.bed \
    megalodon_results.rice.pass1_guppy.CHG/low_coverage_sites.CHG.bed
awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.rice.pass1_guppy.CHG/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.rice.pass1_guppy.CHG/mappings.filtered.sorted.bed
bedtools intersect \
    -a megalodon_results.rice.pass1_guppy.CHG/mappings.filtered.sorted.bed \
    -b megalodon_results.rice.pass1_guppy.CHG/low_coverage_sites.CHG.sorted.bed -s -sorted -v \
    -g ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.seq_len.sorted.txt | \
    awk '{print $4}' > megalodon_results.rice.pass1_guppy.CHG/train_read_ids.txt
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CHG.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --devices 0,1 --processes 30 \
    --mod-motif Z CHG 0 \
    --ref-include-mods \
    --mod-per-site-threshold megalodon_results.rice.pass1_guppy.CHG/site_mod_thresholds.CHG.bed \
    --read-ids-filename megalodon_results.rice.pass1_guppy.CHG/train_read_ids.txt \
    --output-directory megalodon_results.rice.pass1_guppy.CHG.retrain --overwrite --disable-mod-calibration \
    > megalodon_results.rice.pass1_guppy.CHG.retrain.log 2>&1
export OMP_NUM_THREADS=30
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.rice.pass1_guppy.CHG.retrain/signal_mappings.hdf5 \
    --device 0 --outdir megalodon_results.rice.pass1_guppy.CHG.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.rice.pass1_guppy.CHG.retrain.taiyaki_train.log &
dump_json.py megalodon_results.rice.pass1_guppy.CHG.retrain/model_final.checkpoint \
    --output megalodon_results.rice.pass1_guppy.CHG.retrain/model_final.jsn
# test
# test arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_1.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_2.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_2.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_3.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_3.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_4.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_4.CHG_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_54.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_5.CHG_retrain_comb.log 2>&1
## test rice
#megalodon ../fast5_pass2.single/10x.1/ \
#    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
#    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --gpu_runners_per_device 4 --num_callers 8" \
#    --guppy-timeout 60 \
#    --guppy-config model_final.cfg \
#    --outputs per_read_mods mod_mappings mods \
#    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
#    --write-mods-text \
#    --mod-output-formats bedmethyl \
#    --devices 0,1 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb \
#    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb.log 2>&1 &
#Rscript ../correlation_with_bs.cal_plot.general.R \
#    ../bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt \
#    megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb/modified_bases.5mC.bed bisulfite.rep1 \
#    megalodon_retrain_comb.guppy.rice1-1.10x.1 CHG yes \
#    megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb \
#    CHG_megalodon_retrain_comb.guppy.rice1-1.10x.1_vs_bisulfite.rep1.tsv \
#    > megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb/CHG_megalodon_retrain_comb.guppy.rice1-1.10x.1_vs_bisulfite.rep1.log \
#    2>&1 &
# test rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass2.single/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_1.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_2.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_2.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_3.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_3.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_4.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_4.CHG_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHG.retrain/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 0 --processes 20 --output-directory megalodon_results.rice.pass2_guppy.10x_5.CHG_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_5.CHG_retrain_comb.log 2>&1












# 3. train CHH (megalodon 2.2.3)
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --mod-motif Z CHH 0 --outputs per_read_mods mappings --devices 0,1 --processes 20 \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHH --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHH.log 2>&1
megalodon_extras \
    modified_bases per_site_thresholds \
    megalodon_results.arab.pass.part1_guppy.CHH \
    athaliana.bsrep_123_comb.CHH.methyl.bed \
    --mod-base Z \
    --ground-truth-cov-min 25 \
    --nanopore-cov-min 30 \
    --ground-truth-coverage-pdf megalodon_results.arab.pass.part1_guppy.CHH/gt_cov.CHH.pdf \
    --out-blacklist-sites megalodon_results.arab.pass.part1_guppy.CHH/low_coverage_sites.CHH.bed \
    --out-per-site-mod-thresholds megalodon_results.arab.pass.part1_guppy.CHH/site_mod_thresholds.CHH.bed
sort -S 25% --parallel=30 -T /tmp/ \
    -k1,1V -k2,2n \
    -o megalodon_results.arab.pass.part1_guppy.CHH/low_coverage_sites.CHH.sorted.bed \
    megalodon_results.arab.pass.part1_guppy.CHH/low_coverage_sites.CHH.bed
awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.arab.pass.part1_guppy.CHH/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.arab.pass.part1_guppy.CHH/mappings.filtered.sorted.bed
bedtools intersect \
    -a megalodon_results.arab.pass.part1_guppy.CHH/mappings.filtered.sorted.bed \
    -b megalodon_results.arab.pass.part1_guppy.CHH/low_coverage_sites.CHH.sorted.bed -s -sorted -v | \
    awk '{print $4}' > megalodon_results.arab.pass.part1_guppy.CHH/train_read_ids.txt
megalodon ../reads_single_pass/part1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --gpu_runners_per_device 4 --num_callers 8" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna \
    --devices 0,1 --processes 20 \
    --mod-motif Z CHH 0 \
    --ref-include-mods \
    --mod-per-site-threshold megalodon_results.arab.pass.part1_guppy.CHH/site_mod_thresholds.CHH.bed \
    --read-ids-filename megalodon_results.arab.pass.part1_guppy.CHH/train_read_ids.txt \
    --output-directory megalodon_results.arab.pass.part1_guppy.CHH.retrain --overwrite \
    > megalodon_results.arab.pass.part1_guppy.CHH.retrain.log 2>&1

export OMP_NUM_THREADS=12
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.arab.pass.part1_guppy.CHH.retrain/signal_mappings.hdf5 \
    --device 0 --outdir megalodon_results.arab.pass.part1_guppy.CHH.retrain --overwrite \
    --niteration 10000 --save_every 200 --sub_batches 1 \
    > megalodon_results.arab.pass.part1_guppy.CHH.retrain.taiyaki_train.log &
dump_json.py megalodon_results.arab.pass.part1_guppy.CHH.retrain/model_final.checkpoint \
    --output megalodon_results.arab.pass.part1_guppy.CHH.retrain/model_final.jsn
# test
megalodon ../reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CHH.retrain/ --num_callers 2" \
    --guppy-timeout 60 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_arab \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_arab.log 2>&1 &
Rscript ../correlation_with_bs.cal_plot.general.R \
    ../bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt \
    megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_arab/modified_bases.5mC.bed \
    bisulfite.rep1 megalodon_retrain_arab.guppy.arab.10x.1 CHH yes \
    megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_arab \
    CHH_megalodon_retrain_arab.guppy.arab.10x.1_vs_bisulfite.rep1.tsv \
    > megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_arab/CHH_megalodon_retrain_arab.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --mod-motif Z CHH 0 --outputs per_read_mods mappings --devices 2,3 --processes 30 \
    --output-directory megalodon_results.rice.pass1_guppy.CHH --overwrite --disable-mod-calibration \
    > megalodon_results.rice.pass1_guppy.CHH.log 2>&1
megalodon_extras \
    modified_bases per_site_thresholds \
    megalodon_results.rice.pass1_guppy.CHH \
    shuidao2-1.bs.CHH.methyl.bed \
    --mod-base Z \
    --ground-truth-cov-min 25 \
    --nanopore-cov-min 30 \
    --ground-truth-coverage-pdf megalodon_results.rice.pass1_guppy.CHH/gt_cov.CHH.pdf \
    --out-blacklist-sites megalodon_results.rice.pass1_guppy.CHH/low_coverage_sites.CHH.bed \
    --out-per-site-mod-thresholds megalodon_results.rice.pass1_guppy.CHH/site_mod_thresholds.CHH.bed \
    > megalodon_results.rice.pass1_guppy.CHH.extra.log 2>&1
sort -S 25% --parallel=30 -T /tmp/ \
    -k1,1V -k2,2n \
    -o megalodon_results.rice.pass1_guppy.CHH/low_coverage_sites.CHH.sorted.bed \
    megalodon_results.rice.pass1_guppy.CHH/low_coverage_sites.CHH.bed
awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
    megalodon_results.rice.pass1_guppy.CHH/mappings.summary.txt | \
    sort -S 25% --parallel=30 -T /tmp/ -k1,1V -k2,2n > \
    megalodon_results.rice.pass1_guppy.CHH/mappings.filtered.sorted.bed
bedtools intersect \
    -a megalodon_results.rice.pass1_guppy.CHH/mappings.filtered.sorted.bed \
    -b megalodon_results.rice.pass1_guppy.CHH/low_coverage_sites.CHH.sorted.bed -s -sorted -v \
    -g ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.seq_len.sorted.txt | \
    awk '{print $4}' > megalodon_results.rice.pass1_guppy.CHH/train_read_ids.txt
megalodon ../fast5_pass1.single/100x/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --outputs per_read_mods signal_mappings \
    --guppy-params "-d /home/nipeng/data/arab.nano/megalodon_retrain/megalodon_results.arab.pass.part1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
    --devices 2,3 --processes 30 \
    --mod-motif Z CHH 0 \
    --ref-include-mods \
    --mod-per-site-threshold megalodon_results.rice.pass1_guppy.CHH/site_mod_thresholds.CHH.bed \
    --read-ids-filename megalodon_results.rice.pass1_guppy.CHH/train_read_ids.txt \
    --output-directory megalodon_results.rice.pass1_guppy.CHH.retrain --overwrite --disable-mod-calibration \
    > megalodon_results.rice.pass1_guppy.CHH.retrain.log 2>&1
export OMP_NUM_THREADS=30
train_flipflop.py ~/tools/taiyaki/models/mLstm_cat_mod_flipflop.py \
    megalodon_results.rice.pass1_guppy.CHH.retrain/signal_mappings.hdf5 \
    --device 2 --outdir megalodon_results.rice.pass1_guppy.CHH.retrain --overwrite \
    --niteration 5000 --save_every 200 --sub_batches 1 \
    > megalodon_results.rice.pass1_guppy.CHH.retrain.taiyaki_train.log &
dump_json.py megalodon_results.rice.pass1_guppy.CHH.retrain/model_final.checkpoint \
    --output megalodon_results.rice.pass1_guppy.CHH.retrain/model_final.jsn
# test ==
# arab
export OMP_NUM_THREADS=1
megalodon ../reads_single_pass/part2/100x/10x.1/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --num_callers 2" \
    --guppy-timeout 60 --guppy-config model_final.cfg --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 --write-mods-text \
    --mod-output-formats bedmethyl --devices 2 --processes 20 \
    --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_1.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.2/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --num_callers 2" \
    --guppy-timeout 60 --guppy-config model_final.cfg --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 --write-mods-text \
    --mod-output-formats bedmethyl --devices 2 --processes 20 \
    --output-directory megalodon_results.arab.pass.part2_guppy.10x_2.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_2.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.3/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --num_callers 2" \
    --guppy-timeout 60 --guppy-config model_final.cfg --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 --write-mods-text \
    --mod-output-formats bedmethyl --devices 2 --processes 20 \
    --output-directory megalodon_results.arab.pass.part2_guppy.10x_3.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_3.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.4/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --num_callers 2" \
    --guppy-timeout 60 --guppy-config model_final.cfg --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 --write-mods-text \
    --mod-output-formats bedmethyl --devices 2 --processes 20 \
    --output-directory megalodon_results.arab.pass.part2_guppy.10x_4.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_4.CHH_retrain_comb.log 2>&1
megalodon ../reads_single_pass/part2/100x/10x.5/ --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homeb/nipeng/data/shuidao2-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --num_callers 2" \
    --guppy-timeout 60 --guppy-config model_final.cfg --outputs per_read_mods mod_mappings mods \
    --reference ../GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 --write-mods-text \
    --mod-output-formats bedmethyl --devices 2 --processes 20 \
    --output-directory megalodon_results.arab.pass.part2_guppy.10x_5.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.arab.pass.part2_guppy.10x_5.CHH_retrain_comb.log 2>&1
# test rice
export OMP_NUM_THREADS=1
megalodon ../fast5_pass2.single/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_1.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_1.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.2/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_2.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_2.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.3/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_3.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_3.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.4/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_4.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_4.CHH_retrain_comb.log 2>&1
megalodon ../fast5_pass2.single/10x.5/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /homec/nipeng/data/shuidao1-1/megalodon_retrain/megalodon_results.rice.pass1_guppy.CHH.retrain/ --gpu_runners_per_device 2 --num_callers 4" \
    --guppy-timeout 100 \
    --guppy-config model_final.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference ../Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl \
    --devices 2,3 --processes 30 --output-directory megalodon_results.rice.pass2_guppy.10x_5.CHH_retrain_comb \
    --overwrite --disable-mod-calibration > megalodon_results.rice.pass2_guppy.10x_5.CHH_retrain_comb.log 2>&1





















