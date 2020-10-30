# rice CG pipeline ====
# --
#(1) 选取BS数据标准集，主要是正样本
#(2) arab, rice各自denoise; resample负样本；训练模型
#(3) arab+rice的denoise和各自resample的负样本；训练模型
#(4) 测试训练好的模型

# CG pipeline (11/128) =========
# extract_features =====
# part1/0-337 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --methy_label 1 --kmer_len 11 --cent_signals_len 128 --motifs CG --mod_loc 0 --write_path rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.tsv --nproc 20 --positions bs.poses/shuidao.CG.2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.log &
# part1/0-337 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --methy_label 0 --kmer_len 11 --cent_signals_len 128 --motifs CG --mod_loc 0 --write_path rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.tsv --nproc 20 --positions bs.poses/shuidao.CG.2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.log &


# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.tsv --write_filepath rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
=python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.tsv --write_filepath rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/m5c_arab/get_kmer_dist_of_feafile.py --feafile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.27868, 0.2228643, 0.3009632, 0.1974925]
python /homeb/nipeng/tools/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.tsv --krfile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.asposkmer_10m.tsv &
# 261367 kmers from kmer2ratio file
# for 261244 common kmers, fill 9788186 samples, 222864 samples that can't filled
# totalline: 10000000, need to fill: 211814
# extract 155960 samples from 773 diff kmers


# get train data, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_pos.r10m.tsv --fp2 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_neg.asposkmer_10m.tsv --concated_fp rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv &
# 19944146 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv


# train, RNN+CNN, as_pos_kmer neg samples =====
# 19944146 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv
head -19934146 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv > rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv > rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 11 --signal_num 128 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 11 --signal_num 128 &
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.train.bin --valid_file rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.rice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn --log_dir model.CG.rice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn.vis_log --is_cnn yes --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 7 > model.CG.rice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn.log 2>&1 &
# train, denoise_rnnnobase, train, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal denoise --train_file rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv --is_cnn no --is_base no --is_rnn yes --seq_len 11 --cent_signals_len 128 --iterations 8 --rounds 3 --epoch_num 3 --score_cf 0.5 > rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.denoise_rnn_nobase.iter8_rd3_ep3.log 2>&1 &


# test RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.1/ --model_path cg.model.11mer/model.CG.rice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_6.ckpt --result_file shuidao3-1.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path shuidao3-1.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.tsv --result_file shuidao3-1.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1904165A-QJ_L3L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao3-1.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.freq.tsv bisulfite.rep3 deepsignal.11_0.9_rice_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_rice_negkmeraspos_brnncnn.10x.1.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_rice_negkmeraspos_brnncnn.10x.1.log 2>&1 &
# test arab 10x ==
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/10x.1/ --model_path rice.model.11mer/model.CG.rice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_6.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_rice_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_rice_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_rice_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_rice_negkmeraspos_brnncnn.10x.1.log 2>&1 &






# test rice3-1 40x ============================================
# CG 10x.1 ==
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/10x.1/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
# CG 10x.2 ==
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/10x.2/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.freq.tsv --prob_cf 0 &
# CG 10x.3 ==
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/10x.3/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.freq.tsv --prob_cf 0 &
# CG 10x.4 ==
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/10x.4/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.freq.tsv --prob_cf 0 &
# 40x_1234 call freq =
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --input_path shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --result_file shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.40x.1234_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1904165A-QJ_L3L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao3-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.40x.1234_call_mods.freq.tsv bisulfite.rep3 deepsignal.rice.11_0.9_comb_negkmeraspos_brnncnn.40x1234 CG yes analysis CG.bisulfite.rep3.vs.deepsignal.rice.11_0.9_comb_negkmeraspos_brnncnn.40x1234.tsv > analysis/CG.bisulfite.rep3.vs.deepsignal.rice.11_0.9_comb_negkmeraspos_brnncnn.40x1234.log 2>&1 &


# megalodon test ===
# rice3-1 10x.1 CG
megalodon reads_single_pass/part2/10x.1/ --load-default-model --outputs mods --write-mods-text --reference Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CG 0 --devices 2 3 --processes 16 --output-directory megalodon_results.rice3-1.CG.10x.1 --overwrite > megalodon_results.rice3-1.CG.10x.1.log 2>&1 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1904165A-QJ_L3L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt megalodon_results.rice3-1.CG.10x.1.tsv bisulfite.rep3 megalodon.taiyaki.rice3-1.10x.1 CG yes analysis CG.bisulfite.rep3.vs.megalodon.taiyaki.rice3-1.10x.1.tsv > analysis/CG.bisulfite.rep3.vs.megalodon.taiyaki.rice3-1.10x.1.log 2>&1 &
# rice3-1 10x.1 CHG
megalodon reads_single_pass/part2/10x.1/ --load-default-model --outputs mods --write-mods-text --reference Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHG 0 --devices 2 3 --processes 16 --output-directory megalodon_results.rice3-1.CHG.10x.1 --overwrite > megalodon_results.rice3-1.CHG.10x.1.log 2>&1 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1904165A-QJ_L3L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt megalodon_results.rice3-1.CHG.10x.1.tsv bisulfite.rep3 megalodon.taiyaki.rice3-1.10x.1 CHG yes analysis CHG.bisulfite.rep3.vs.megalodon.taiyaki.rice3-1.10x.1.tsv > analysis/CHG.bisulfite.rep3.vs.megalodon.taiyaki.rice3-1.10x.1.log 2>&1 &
# rice3-1 10x.1 CHH
megalodon reads_single_pass/part2/10x.1/ --load-default-model --outputs mods --write-mods-text --reference Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif Z CHH 0 --devices 2 3 --processes 16 --output-directory megalodon_results.rice3-1.CHH.10x.1 --overwrite > megalodon_results.rice3-1.CHH.10x.1.log 2>&1 &





# test rice2-1 100x ============================================
# CG 10x.1
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/10x.1 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/10x.1p --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1p_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1p_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --input_path shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1p_call_mods.tsv --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1a_call_mods.freq.tsv &
# CG 10x.2
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.2 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.log 2>&1 &
CUDA_VISIBLE_DEVICES=2 deepsignal call_mods --input_path reads_single_pass/10x.2p --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2p_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2p_call_mods.log 2>&1 &
# CG 10x.3
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.3 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.log 2>&1 &
# CG 10x.4
CUDA_VISIBLE_DEVICES=2 deepsignal call_mods --input_path reads_single_pass/10x.4 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.log 2>&1 &
# CG 10x.5
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/10x.5 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao2-1.pass.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.log 2>&1 &








# ======================================================================
# guppy basecalled deepsignal2 model ============================
# ======= CG pipeline (13/16) ================================
# extract_features shuidao2-1 =====
# part1/0-601 CG positive
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir fast5_pass1.single --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --methy_label 1 --seq_len 13 --signal_len 16 --motifs CG --mod_loc 0 --write_path shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.tsv --nproc 20 --positions bs.poses/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.rmet_0.9.cov_5_50000.hc_positions2_pos.txt > shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.log &
# part1/0-601 CG negative
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir fast5_pass1.single --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/rice/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --methy_label 0 --seq_len 13 --signal_len 16 --motifs CG --mod_loc 0 --write_path shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_neg.tsv --nproc 20 --positions bs.poses/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.rmet_0.0.cov_5_50000.hc_positions2_neg.txt > shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_neg.log &

# fliter samples, random shuidao2-1 =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.tsv --write_filepath shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos shuidao2-1=====
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.2569912, 0.2552761, 0.2919374, 0.1957953]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_neg.tsv --krfile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_neg.asposkmer_10m.tsv &
# 1949571 kmers from kmer2ratio file
# for 1443200 common kmers, fill 7468229 samples, 490362 samples that can't filled
# totalline: 10000000, need to fill: 2531771
# extract 2729174 samples from 1364757 diff kmers

# get train data, balance shuidao2-1 =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_pos.r10m.tsv --fp2 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_neg.asposkmer_10m.tsv --concated_fp shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv &
# 20197403 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv
head -20187403 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv > shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.train.tsv &
tail -10000 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv > shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.valid.tsv &
# train model, balance, seq+signal, shuidao2-1 =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/train.py --train_file shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.train.tsv --valid_file shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm.train.log 2>&1 &
# denoise signal, train seq+signal, shuidao2-1 =======
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/denoise.py --train_file shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv --is_filter_fn no --model_type signal_bilstm --is_base yes --epoch_num 3 --rounds 3 --iterations 8 --seq_len 13 > denoise.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.balance.signal_bilstm.iter8_rd3_ep3.log 2>&1 &
# train
# 17703610 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.tsv 
head -17693610 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.tsv > shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.train.tsv
tail -10000 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.tsv > shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.valid.tsv
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/train.py --train_file shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.train.tsv --valid_file shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp6.valid.tsv --model_dir model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.denoise_signal_bilstm.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.denoise_signal_bilstm.both_bilstm.train.log 2>&1 &


# test model, balance, seq+signal, shuidao1-1 =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.1 --model_path ../shuidao2-1/model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch5.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 15 --nproc_gpu 3 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.rice1-1.13_rice2-1_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.rice1-1.13_rice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.rice1-1.13_rice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test on arab.10x.1
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.rice2-1/model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch5.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.arab.13_rice2-1_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.arab.13_rice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.arab.13_rice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, denoise_signal_bilstm, seq+signal, shuidao1-1 =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.1 --model_path ../shuidao2-1/model.dp2.CG.rice2-1_R9.4plus_tem.bn13_sn16.denoise_signal_bilstm.both_bilstm/both_bilstm.b13_s16_epoch5.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.denoise_signal_bilstm.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.denoise_signal_bilstm.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.denoise_signal_bilstm.both_bilstm.10x_1.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao1-1.guppy.pass.part2.CG.bn13_sn16.rice2-1.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.rice1-1.13_rice2-1_denoise_signal_bilstm_both_bilstm.10x_1 CG yes analysis CG_deepsignal.rice1-1.13_rice2-1_denoise_signal_bilstm_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.rice1-1.13_rice2-1_denoise_signal_bilstm_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &






# guppy basecalled megalodon ============================
# rice 10x.1 CG
# example = 
# --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
# --write-mods-text \
# --mod-aggregate-method binary_threshold \
# --mod-binary-threshold 0.75 \
# --mod-output-formats bedmethyl wiggle \
megalodon fast5_pass2.single/10x.1/ \
    --guppy-server-path /usr/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
    --guppy-timeout 30 \
    --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 0 --processes 20 --output-directory megalodon_results.rice1-1.pass.part2_guppy.10x_1.CG \
    --overwrite > megalodon_results.rice1-1.pass.part2_guppy.10x_1.CG.log 2>&1 &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt megalodon_results.rice1-1.pass.part2_guppy.10x_1.CG/modified_bases.5mC.bed bisulfite.rep1 megalodon.guppy.rice1-1.10x.1 CG yes analysis CG_megalodon.guppy.rice1-1.10x.1_vs_bisulfite.rep1.tsv > analysis/CG_megalodon.guppy.rice1-1.10x.1_vs_bisulfite.rep1.log 2>&1 &















