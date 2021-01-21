# =====================================================================================================
# --
#(1) 选取BS数据标准集，主要是正样本
#(2) arab, rice各自denoise; resample负样本；训练模型
#(3) arab+rice的denoise和各自resample的负样本；训练模型
#(4) 测试训练好的模型

# CG pipeline (11/128) =========
# extract_features =====
# part1/0-601 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 11 --cent_signals_len 128 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --kmer_len 11 --cent_signals_len 128 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.log &



# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.r10m.tsv --num_lines 10000000 --header false &
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.r10m.tsv &
# (['CGA', 'CGC', 'CGG', 'CGT'], [0.449991, 0.1754438, 0.2665667, 0.1079985])
python /homeb/nipeng/tools/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 104809 kmers from kmer2ratio file
# using 101015 common kmers, fill 9463406 samples, 278983 samples that can't filled
# totalline: 10000000, need to fill: 536594
# extract 564613 samples from 141156 diff kmers



# get train data, random =====
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.r10m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.tsv &
head -19990000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.valid.tsv &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.train.tsv --base_num 11 --signal_num 128 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.valid.tsv --base_num 11 --signal_num 128 &
# get train data, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv &
# 20028019 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv
head -20018019 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 11 --signal_num 128 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 11 --signal_num 128 &



# train, RNN+CNN, random neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.r20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negrandom.brnncnn --log_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negrandom.brnncnn.vis_log --is_cnn yes --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 8 > model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negrandom.brnncnn.log 2>&1 &
# train, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn --log_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn.vis_log --is_cnn yes --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 7 > model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn.log 2>&1 &
# train, RNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnn --log_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnn.vis_log --is_cnn no --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 8 > model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnn.log 2>&1 &
# train, denoise_rnnnobase, train_RNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal denoise --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv --is_cnn no --is_base no --is_rnn yes --seq_len 11 --cent_signals_len 128 --iterations 8 --rounds 3 --epoch_num 3 --score_cf 0.5 > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.denoise_rnn_nobase.iter8_rd3_ep3.log 2>&1 &
# denoise denoise8 ==
python /media/cb201/disk1/nipeng/m5c_arab/filter_samples_by_label.py --sf_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.denoise8.tsv --midfix label1 --label 1 &
python /media/cb201/disk1/nipeng/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.denoise8.label1.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.4427166358220566, 0.17930791446128755, 0.27353358775642894, 0.10444186196022691]
# 7139551 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.denoise8.label1.tsv
python /media/cb201/disk1/nipeng/m5c_arab/filter_samples_by_label.py --sf_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv --midfix label0 --label 0 &
python /media/cb201/disk1/nipeng/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.denoise8.label1.kmer_distri.tsv --totalline 7139551 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.denoise8_kmeraspos.tsv &
# 104765 kmers from kmer2ratio file
# for 100972 common kmers, fill 6771917 samples, 182642 samples that can't filled
# totalline: 7139551, need to fill: 367634
# extract 423565 samples from 141199 diff kmers
python /media/cb201/disk1/nipeng/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.denoise8.label1.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.denoise8_kmeraspos.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.tsv &
# 14335033 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.tsv
head -14325033 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.valid.tsv &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.train.tsv --base_num 11 --signal_num 128 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.valid.tsv --base_num 11 --signal_num 128 &
# denoise8 train RNN
CUDA_VISIBLE_DEVICES=3 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnn --log_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnn.vis_log --is_cnn no --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 8 > model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnn.log 2>&1 &
# denoise8 train RNN+CNN
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs.hc_poses_all.negkmeraspos_20m.denoise8.negkmeraspos.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnncnn --log_dir model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnncnn.vis_log --is_cnn yes --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 8 > model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnncnn.log 2>&1 &




# test, RNN+CNN, random neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.1/ --model_path cg.model.11mer/model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negrandom.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negrandom_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negrandom_brnncnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negrandom_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negrandom_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negrandom_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_arab_negrandom_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negrandom_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negrandom_brnncnn.10x.1.log 2>&1 &
# test, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.1/ --model_path cg.model.11mer/model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_6.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_arab_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negkmeraspos_brnncnn.10x.1.log 2>&1 &
# test on rice.10x =
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/10x.1/ --model_path arab.model.11mer/model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_6.ckpt --result_file shuidao3-1.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path shuidao3-1.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.tsv --result_file shuidao3-1.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1904165A-QJ_L3L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao3-1.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.freq.tsv bisulfite.rep3 deepsignal.11_0.9_arab_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negkmeraspos_brnncnn.10x.1.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negkmeraspos_brnncnn.10x.1.log 2>&1 &
# test, RNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.1/ --model_path cg.model.11mer/model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnn/bn_11.sn_128.epoch_7.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnn.10x.1_call_mods.tsv --is_cnn no --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_negkmeraspos_brnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_arab_negkmeraspos_brnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negkmeraspos_brnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_negkmeraspos_brnn.10x.1.log 2>&1 &
# test rnnnobase_denoise8, train RNN =====
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/10x.1/ --model_path cg.model.11mer/model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnn/bn_11.sn_128.epoch_7.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnn.10x.1_call_mods.tsv --is_cnn no --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_arab_denoise8_negkmeraspos_brnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_denoise8_negkmeraspos_brnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_denoise8_negkmeraspos_brnn.10x.1.log 2>&1 &
# test rnnnobase_denoise8, train RNN+CNN =====
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/10x.1/ --model_path cg.model.11mer/model.CG.arab_R9.4plus_tem.bn11.sn128.0.9.denoise8_negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn11sn128_0.9_denoise8_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_arab_denoise8_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_denoise8_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_arab_denoise8_negkmeraspos_brnncnn.10x.1.log 2>&1 &







# model: CG negkmeraspos arab+rice =====================================================================
# random 10m negkmeraspos arab ===
# 20028019 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv
python /media/cb201/disk1/nipeng/m5c_arab/filter_samples_by_label.py --sf_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv --midfix label0 --label 0 &

python /media/cb201/disk1/nipeng/m5c_arab/filter_samples_by_label.py --sf_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv --midfix label1 --label 1 &
python /media/cb201/disk1/nipeng/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label1.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label1.r5m.tsv --num_lines 5000000 --header false &
python /media/cb201/disk1/nipeng/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label1.r5m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.4499446, 0.1754488, 0.2666818, 0.1079248]
python /media/cb201/disk1/nipeng/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label1.r5m.kmer_distri.tsv --totalline 5000000 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.asposkmer_5m.tsv &
# 104788 kmers from kmer2ratio file
# for 100995 common kmers, fill 4792049 samples, 83647 samples that can't filled
# totalline: 5000000, need to fill: 207951
# extract 282341 samples from 141176 diff kmers
python /media/cb201/disk1/nipeng/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label1.r5m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.asposkmer_5m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_10m.tsv &
# 5074390 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.label0.asposkmer_5m.tsv
# random 10m negkmeraspos rice ===
python /media/cb201/disk1/nipeng/m5c_arab/filter_samples_by_label.py --sf_path rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv --midfix label0 --label 0 &

python /media/cb201/disk1/nipeng/m5c_arab/filter_samples_by_label.py --sf_path rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.tsv --midfix label1 --label 1 &
python /media/cb201/disk1/nipeng/basemods_RNN_signal/randsel_file_rows.py --ori_filepath rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label1.tsv --write_filepath rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label1.r5m.tsv --num_lines 5000000 --header false &
python /media/cb201/disk1/nipeng/m5c_arab/get_kmer_dist_of_feafile.py --feafile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label1.r5m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.2785464, 0.2229132, 0.3011014, 0.197439]
python /media/cb201/disk1/nipeng/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label0.tsv --krfile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label1.r5m.kmer_distri.tsv --totalline 5000000 --wfile rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label0.asposkmer_5m.tsv &
# 258851 kmers from kmer2ratio file
# for 258733 common kmers, fill 4960505 samples, 44759 samples that can't filled
# totalline: 5000000, need to fill: 39495
# extract 16784 samples from 3284 diff kmers
python /media/cb201/disk1/nipeng/deepsignal/scripts/concat_two_files.py --fp1 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label1.r5m.tsv --fp2 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_20m.label0.asposkmer_5m.tsv --concated_fp rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_10m.tsv &
# combine arab+rice ===
python /media/cb201/disk1/nipeng/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_all.negkmeraspos_10m.tsv --fp2 rice.3-1.20190523-NPL0938-P5-PAD58112.pass.part1.CG.bn11.sn128.signal_feas.rep3irep2.hc_poses_all.negkmeraspos_10m.tsv --concated_fp arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.tsv &
# 20051679 arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.tsv
head -20041679 arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.tsv > arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.tsv > arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 11 --signal_num 128 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/generate_binary_feature_file.py --input_file arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 11 --signal_num 128 &
# train ===
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.train.bin --valid_file arab_and_rice_3-1.pass.part1.CG.bn11.sn128.signal_feas.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn --log_dir model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn.vis_log --is_cnn yes --is_base yes --is_rnn yes --kmer_len 11 --cent_signals_len 128 --is_binary yes --max_epoch_num 8 > model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn.log 2>&1 &
# test arab.10x ===
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/10x.1/ --model_path ../comb.arab_rice/cg.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /media/cb201/disk1/nipeng/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.11_0.9_comb_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_comb_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_comb_negkmeraspos_brnncnn.10x.1.log 2>&1 &
# test rice.10x ===
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/10x.1/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file shuidao3-1.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > shuidao3-1.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.log 2>&1 &
python ~/tools/deepsignal/scripts/call_modification_frequency.py --input_path shuidao3-1.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.tsv --result_file shuidao3-1.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1904165A-QJ_L3L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao3-1.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x1_call_mods.freq.tsv bisulfite.rep3 deepsignal.11_0.9_comb_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.11_0.9_comb_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.11_0.9_comb_negkmeraspos_brnncnn.10x.1.log 2>&1 &
# test arab.10x guppy (10.53 /home)===
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.freq.tsv bisulfite.rep1 deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.10x.1_a2g CG yes analysis CG.bisulfite.rep1.vs.deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.10x.1_a2g.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.10x.1_a2g.log 2>&1 &




# test arab.100x ===
# 10x.1 =
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
# 10x.2 =
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.2 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.freq.tsv --prob_cf 0 &
# 10x.3 =
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.3 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.freq.tsv --prob_cf 0 &
# 10x.4 =
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.4 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.freq.tsv --prob_cf 0 &
# 10x.5 =
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.5 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.freq.tsv --prob_cf 0 &
# 50x_12345 call freq =
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency2.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods.tsv --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods.tsv --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods.tsv --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods.freq.tsv bisulfite.rep1 deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.50x12345 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.50x12345.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.50x12345.log 2>&1 &


# test arab.100x guppy ===
# test arab.10x guppy (10.53 /home)===
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1/ --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.freq.tsv --prob_cf 0 &
# 10x.2 =
CUDA_VISIBLE_DEVICES=2 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.2 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods_a2g.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods_a2g.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods_a2g.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods_a2g.freq.tsv --prob_cf 0 &
# 10x.3 =
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.3 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods_a2g.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods_a2g.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods_a2g.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods_a2g.freq.tsv --prob_cf 0 &
# 10x.4 =
CUDA_VISIBLE_DEVICES=2 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.4 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods_a2g.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods_a2g.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods_a2g.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods_a2g.freq.tsv --prob_cf 0 &
# 10x.5 =
CUDA_VISIBLE_DEVICES=3 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.5 --model_path comb.model.11mer/model.CG.arabnrice_R9.4plus_tem.bn11.sn128.0.9.negkmeraspos.brnncnn/bn_11.sn_128.epoch_5.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods_a2g.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 11 --cent_signals_len 128 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods_a2g.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods_a2g.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods_a2g.freq.tsv --prob_cf 0 &
# combine 50x_12345
python ~/tools/deepsignal/scripts/combine_call_mods_freq_files.py --modsfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.1_call_mods_a2g.freq.tsv --modsfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.2_call_mods_a2g.freq.tsv --modsfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.3_call_mods_a2g.freq.tsv --modsfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.4_call_mods_a2g.freq.tsv --modsfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.10x.5_call_mods_a2g.freq.tsv --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods_a2g.freq.tsv &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_comb_bn11sn128_0.9_negkmeraspos_brnncnn.50x.12345_call_mods_a2g.freq.tsv bisulfite.rep1 deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.50x12345.a2g CG yes analysis CG.bisulfite.rep1.vs.deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.50x12345.a2g.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.arab.11_0.9_comb_negkmeraspos_brnncnn.50x12345.a2g.log 2>&1 &







# megalodon test ==================
# arab 10x.1 CG
megalodon reads_single_pass/part2/100x/10x.1/ --load-default-model --outputs mods --write-mods-text --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CG 0 --devices 0 1 --processes 16 --output-directory megalodon_results.arab.CG.10x.1 --overwrite > megalodon_results.arab.CG.10x.1.log 2>&1 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt modified_bases.5mC.megalodon.arab.10x.1.CG.bed bisulfite.rep1 megalodon.taiyaki.arab.10x.1 CG yes analysis CG.bisulfite.rep1.vs.megalodon.taiyaki.arab.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.megalodon.taiyaki.arab.10x.1.log 2>&1 &
# arab 10x.1 CHG
megalodon reads_single_pass/part2/100x/10x.1/ --load-default-model --outputs mods --write-mods-text --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 --devices 1 --processes 10 --output-directory megalodon_results.arab.CHG.10x.1 --overwrite > megalodon_results.arab.CHG.10x.1.log 2>&1 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt modified_bases.5mC.megalodon.arab.10x.1.CHG.bed bisulfite.rep1 megalodon.taiyaki.arab.10x.1 CHG yes analysis CHG.bisulfite.rep1.vs.megalodon.taiyaki.arab.10x.1.tsv > analysis/CHG.bisulfite.rep1.vs.megalodon.taiyaki.arab.10x.1.log 2>&1 &
# arab 10x.1 CHH
megalodon reads_single_pass/part2/100x/10x.1/ --load-default-model --outputs mods --write-mods-text --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 --devices 0 1 --processes 16 --output-directory megalodon_results.arab.CHH.10x.1 --overwrite > megalodon_results.arab.CHH.10x.1.log 2>&1 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt modified_bases.5mC.megalodon.arab.10x.1.CHH.bed bisulfite.rep1 megalodon.taiyaki.arab.10x.1 CHH yes analysis CHH.bisulfite.rep1.vs.megalodon.taiyaki.arab.10x.1.tsv > analysis/CHH.bisulfite.rep1.vs.megalodon.taiyaki.arab.10x.1.log 2>&1 &


# tombo test (guppy) =================
# dp2, 192.16.10.53
tombo resquiggle reads_single_pass/part2/100x/10x.1 GCF_000001735.4_TAIR10.1_genomic.fna --processes 20 --corrected-group RawGenomeCorrected_002 --basecall-group Basecall_1D_001 --overwrite --ignore-read-locks > tombo_arab_pass.part2_10x.1.guppy.log 2>&1 &
# arab 10x.1 CG
tombo detect_modifications alternative_model --fast5-basedirs reads_single_pass/part2/100x/10x.1 --alternate-bases CpG --statistics-file-basename tombo_results.arab.pass.part2_guppy.10x.1.CG --dna --multiprocess-region-size 1000 --processes 20 --corrected-group RawGenomeCorrected_002 > tombo_results.arab.pass.part2_guppy.10x.1.CG.log 2>&1 &
tombo text_output browser_files \
   --fast5-basedirs reads_single_pass/part2/100x/10x.1 \
   --statistics-filename tombo_results.arab.pass.part2_guppy.10x.1.CG.CpG.tombo.stats \
   --file-types coverage dampened_fraction fraction\
   --browser-file-basename tombo_results.arab.pass.part2_guppy.10x.1.CG --corrected-group RawGenomeCorrected_002 \
   > tombo_results.arab.pass.part2_guppy.10x.1.CG.output.log 2>&1 &
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig.bed
# wig2bed --multisplit bar < tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig.bed
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.plus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.plus.wig.bed
# arab 10x.1 5mC
tombo detect_modifications alternative_model --fast5-basedirs reads_single_pass/part2/100x/10x.1 --alternate-bases 5mC --statistics-file-basename tombo_results.arab.pass.part2_guppy.10x.1.5mC --dna --multiprocess-region-size 1000 --processes 20 --corrected-group RawGenomeCorrected_002 > tombo_results.arab.pass.part2_guppy.10x.1.5mC.log 2>&1 &
tombo text_output browser_files \
   --fast5-basedirs reads_single_pass/part2/100x/10x.1 \
   --statistics-filename tombo_results.arab.pass.part2_guppy.10x.1.5mC.5mC.tombo.stats \
   --file-types coverage dampened_fraction fraction\
   --browser-file-basename tombo_results.arab.pass.part2_guppy.10x.1.5mC --corrected-group RawGenomeCorrected_002 \
   > tombo_results.arab.pass.part2_guppy.10x.1.5mC.output.log 2>&1 &
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.minus.wig.bed
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.plus.wig > tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.plus.wig.bed

   




# CG pipeline, different kmer lens test, 17, 15, 13, 11, 9 =====================
# test kmer2signals len =======
# deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 17 --cent_signals_len 1 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn1.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn1.signal_feas.bs_inter.hc_poses_pos.log &
# athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn1.signal_feas.bs_inter.hc_poses_pos.tsv
python /homeb/nipeng/tools/m5c_arab/stat_count_kmer_signals.py  athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn1.signal_feas.bs_inter.hc_poses_pos.tsv 7

# 17/192 =========
# extract_features =====
# part1/0-601 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 17 --cent_signals_len 192 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --kmer_len 17 --cent_signals_len 192 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_neg.log &
# get train data =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos
python /homeb/nipeng/tools/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.r10m.tsv &
# (['CGA', 'CGC', 'CGG', 'CGT'], [0.4499013, 0.1755593, 0.2664465, 0.1080929])
python /homeb/nipeng/tools/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 199808 kmers from kmer2ratio file
# using 3136 common kmers, fill 161811 samples, 9105 samples that can't filled
# totalline: 10000000, need to fill: 9838189
# extract 11042101 samples from 2208662 diff kmers
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv &
# 21203912 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv
head -21193912 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 17 --signal_num 192 &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 17 --signal_num 192 &
# train, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn17.sn192.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn17.sn192.0.9.negkmeraspos.brnncnn --log_dir model.CG.arab_R9.4plus_tem.bn17.sn192.0.9.negkmeraspos.brnncnn.vis_log --is_cnn yes --is_base yes --is_rnn yes --kmer_len 17 --cent_signals_len 192 --is_binary yes --max_epoch_num 7 > model.CG.arab_R9.4plus_tem.bn17.sn192.0.9.negkmeraspos.brnncnn.log 2>&1 &
# test, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1/ --model_path model.CG.arab_R9.4plus_tem.bn17.sn192.0.9.negkmeraspos.brnncnn/bn_17.sn_192.epoch_4.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn17sn192_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 17 --cent_signals_len 192 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn17sn192_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn17sn192_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn17sn192_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn17sn192_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.17_0.9_arab_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.17_0.9_arab_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.17_0.9_arab_negkmeraspos_brnncnn.10x.1.log 2>&1 &


# 15/172 ==========
# extract_features =====
# part1/0-601 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 15 --cent_signals_len 172 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --kmer_len 15 --cent_signals_len 172 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_neg.log &
# get train data =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos
python /homeb/nipeng/tools/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.r10m.tsv &
# (['CGA', 'CGC', 'CGG', 'CGT'], [0.449753, 0.1754735, 0.266567, 0.1082065])
python /homeb/nipeng/tools/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 194838 kmers from kmer2ratio file
# using 15819 common kmers, fill 803447 samples, 74843 samples that can't filled
# totalline: 10000000, need to fill: 9196553
# extract 10342219 samples from 2068661 diff kmers
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv &
# 21145666 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv
head -21135666 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 15 --signal_num 172 &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 15 --signal_num 172 &
# train, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=0 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn15.sn172.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn15.sn172.0.9.negkmeraspos.brnncnn --is_cnn yes --is_base yes --is_rnn yes --kmer_len 15 --cent_signals_len 172 --is_binary yes --max_epoch_num 7 > model.CG.arab_R9.4plus_tem.bn15.sn172.0.9.negkmeraspos.brnncnn.log 2>&1 &
# test, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1/ --model_path model.CG.arab_R9.4plus_tem.bn15.sn172.0.9.negkmeraspos.brnncnn/bn_15.sn_172.epoch_4.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn15sn172_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 15 --cent_signals_len 172 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn15sn172_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn15sn172_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn15sn172_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn15sn172_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.15_0.9_arab_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.15_0.9_arab_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.15_0.9_arab_negkmeraspos_brnncnn.10x.1.log 2>&1 &



# 13/150 ==========
# extract_features =====
# part1/0-601 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 13 --cent_signals_len 150 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --kmer_len 13 --cent_signals_len 150 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_neg.log &
# get train data =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos
python /homeb/nipeng/tools/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.r10m.tsv &
# (['CGA', 'CGC', 'CGG', 'CGT'], [0.4498681, 0.1753567, 0.2665934, 0.1081818])
python /homeb/nipeng/tools/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 179560 kmers from kmer2ratio file
# using 88528 common kmers, fill 4819619 samples, 267986 samples that can't filled
# totalline: 10000000, need to fill: 5180381
# extract 6276979 samples from 1255485 diff kmers
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv &
# 21096598 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv
head -21086598 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 13 --signal_num 150 &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 13 --signal_num 150 &
# train, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=0 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn13.sn150.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn13.sn150.0.9.negkmeraspos.brnncnn --is_cnn yes --is_base yes --is_rnn yes --kmer_len 13 --cent_signals_len 150 --is_binary yes --max_epoch_num 7 > model.CG.arab_R9.4plus_tem.bn13.sn150.0.9.negkmeraspos.brnncnn.log 2>&1 &
# test, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1/ --model_path model.CG.arab_R9.4plus_tem.bn13.sn150.0.9.negkmeraspos.brnncnn/bn_13.sn_150.epoch_6.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn13sn150_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 13 --cent_signals_len 150 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn13sn150_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn13sn150_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn13sn150_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn13sn150_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.13_0.9_arab_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.13_0.9_arab_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.13_0.9_arab_negkmeraspos_brnncnn.10x.1.log 2>&1 &


# 11/128 ==========


# 9/106 ==========
# extract_features =====
# part1/0-601 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 9 --cent_signals_len 106 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_001 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --kmer_len 9 --cent_signals_len 106 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_neg.log &
# get train data =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos
python /homeb/nipeng/tools/m5c_arab/get_kmer_dist_of_feafile.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.r10m.tsv &
# (['CGA', 'CGC', 'CGG', 'CGT'], [0.4497634, 0.1755318, 0.2665534, 0.1081514])
python /homeb/nipeng/tools/m5c_arab/select_neg_samples_by_kmer_distri.py --feafile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 15456 kmers from kmer2ratio file
# using 15456 common kmers, fill 9944522 samples, 56319 samples that can't filled
# totalline: 10000000, need to fill: 55478
# extract 55620 samples from 927 diff kmers
python /homeb/nipeng/tools/deepsignal/scripts/concat_two_files.py --fp1 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_pos.r10m.tsv --fp2 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv &
# 20000142 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv
head -19990142 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv &
tail -10000 athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.tsv --base_num 9 --signal_num 106 &
python /homeb/nipeng/tools/deepsignal/scripts/generate_binary_feature_file.py --input_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.tsv --base_num 9 --signal_num 106 &
# train, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=0 deepsignal train --train_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.train.bin --valid_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part1.CG.bn9.sn106.signal_feas.bs_inter.hc_poses_all.negkmeraspos_20m.valid.bin --model_dir model.CG.arab_R9.4plus_tem.bn9.sn106.0.9.negkmeraspos.brnncnn --is_cnn yes --is_base yes --is_rnn yes --kmer_len 9 --cent_signals_len 106 --is_binary yes --max_epoch_num 7 > model.CG.arab_R9.4plus_tem.bn9.sn106.0.9.negkmeraspos.brnncnn.log 2>&1 &
# test, RNN+CNN, as_pos_kmer neg samples =====
CUDA_VISIBLE_DEVICES=0 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1/ --model_path model.CG.arab_R9.4plus_tem.bn9.sn106.0.9.negkmeraspos.brnncnn/bn9.sn106.epoch_6.ckpt --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn9sn106_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --is_cnn yes --is_rnn yes --is_base yes --kmer_len 9 --cent_signals_len 106 --corrected_group RawGenomeCorrected_001 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --nproc 10 --is_gpu yes > athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn9sn106_0.9_negkmeraspos_brnncnn.10x.1_call_mods.log 2>&1 &
python /home/nipeng/tools/deepsignal/scripts/call_modification_frequency.py --input_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn9sn106_0.9_negkmeraspos_brnncnn.10x.1_call_mods.tsv --result_file athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn9sn106_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv --prob_cf 0 &
Rscript /homeb/nipeng/tools/m5c_arab/correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.20190421-NPL0867-P1-E11-H11-barcode.pass.part2.CG.m_arab_bn9sn106_0.9_negkmeraspos_brnncnn.10x.1_call_mods.freq.tsv bisulfite.rep1 deepsignal.9_0.9_arab_negkmeraspos_brnncnn.10x.1 CG yes analysis CG.bisulfite.rep1.vs.deepsignal.9_0.9_arab_negkmeraspos_brnncnn.10x.1.tsv > analysis/CG.bisulfite.rep1.vs.deepsignal.9_0.9_arab_negkmeraspos_brnncnn.10x.1.log 2>&1 &









# guppy basecalled, model train ============================
# CG pipeline (11/128) =========
# extract_features =====
# part1/0-601 CG positive
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --kmer_len 11 --cent_signals_len 128 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.guppy.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 5 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.guppy.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
deepsignal extract --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --kmer_len 11 --cent_signals_len 128 --motifs CG --mod_loc 0 --write_path athaliana.20190421-NPL0867-P1-E11-H11-barcode.guppy.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 10 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.20190421-NPL0867-P1-E11-H11-barcode.guppy.pass.part1.CG.bn11.sn128.signal_feas.bs_inter.hc_poses_neg.log &







# guppy basecalled, deepsignal2 model =============================================================================================================================
# --corrected_group RawGenomeCorrected_002
# ======= CG pipeline (11/16) ================================
# extract_features =====
# part1/0-601 CG positive
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --seq_len 11 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --seq_len 11 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.log &

# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --num_lines 10000000 --header false &
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.4496942, 0.1754998, 0.2665265, 0.1082795]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 104804 kmers from kmer2ratio file
# for 101010 common kmers, fill 9476141 samples, 266286 samples that can't filled
# totalline: 10000000, need to fill: 523859
# extract 564635 samples from 141161 diff kmers


# get train data, random =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.r10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.tsv &
head -19990000 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.train.tsv &
tail -10000 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.valid.tsv &
# get train data, balance =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv &
head -19990000 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv



# train model, random, seq+signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.r20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.random.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 > model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.random.both_bilstm.train.log 2>&1 &
# train model, balance, seq+signal =====
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 > model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.both_bilstm.train.log 2>&1 &
# train model, balance, seq =====
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.seq_bilstm --model_type seq_bilstm --layernum1 3 --layernum2 1 > model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.seq_bilstm.train.log 2>&1 &
# train model, balance, signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.signal_bilstm --model_type signal_bilstm --layernum1 3 --layernum2 1 > model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.signal_bilstm.train.log 2>&1 &
# train model, balance+denoise, seq+signal =====
# denoise seq_nobase, train seq+signal
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/denoise.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv --is_filter_fn no --model_type seq_bilstm --is_base no --epoch_num 3 --rounds 3 --iterations 8 > denoise.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.seq_bilstm_nobase.iter8_rd3_ep3.log 2>&1 &
# denoise signal, train seq+signal
CUDA_VISIBLE_DEVICES=3 python ~/tools/deepsignal2/deepsignal2/denoise.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv --is_filter_fn no --model_type signal_bilstm --epoch_num 3 --rounds 3 --iterations 8 > denoise.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.signal_bilstm.iter8_rd3_ep3.log 2>&1 &
# 15383722 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.tsv
head -15373722 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.tsv > athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.valid.tsv
CUDA_VISIBLE_DEVICES=3 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn11_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.denoise_signal_bilstm.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 > model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.denoise_signal_bilstm.both_bilstm.train.log 2>&1 &




# test model, random, seq+signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.random.both_bilstm/both_bilstm.b11_s16_epoch5.ckpt --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.random.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 11 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn11_sn16.random.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn11_sn16.random.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.random.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn11_sn16.random.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.11_arab_random_both_bilstm.10x_1 CG yes analysis CG_deepsignal.11_arab_random_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.11_arab_random_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, balance, seq+signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.both_bilstm/both_bilstm.b11_s16_epoch4.ckpt --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 11 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 15 --nproc_gpu 3 > athaliana.guppy.pass.part2.CG.bn11_sn16.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn11_sn16.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn11_sn16.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.11_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.11_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.11_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# - rice shuidao1-1 10x.1
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.1 --model_path model.arab/model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.both_bilstm/both_bilstm.b11_s16_epoch4.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn11_sn16.arab.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 11 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > shuidao1-1.guppy.pass.part2.CG.bn11_sn16.arab.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn11_sn16.arab.balance.both_bilstm.10x_1.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn11_sn16.arab.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao1-1.guppy.pass.part2.CG.bn11_sn16.arab.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.rice1-1.11_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.rice1-1.11_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.rice1-1.11_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, balance, signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.signal_bilstm/signal_bilstm.b11_s16_epoch4.ckpt --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.signal_bilstm.10x_1.tsv --model_type signal_bilstm --seq_len 11 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn11_sn16.balance.signal_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn11_sn16.balance.signal_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.signal_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn11_sn16.balance.signal_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.11_arab_balance_signal_bilstm.10x_1 CG yes analysis CG_deepsignal.11_arab_balance_signal_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.11_arab_balance_signal_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, balance, seq =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.balance.seq_bilstm/seq_bilstm.b11_s16_epoch4.ckpt --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.seq_bilstm.10x_1.tsv --model_type seq_bilstm --seq_len 11 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn11_sn16.balance.seq_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn11_sn16.balance.seq_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.seq_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn11_sn16.balance.seq_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.11_arab_balance_seq_bilstm.10x_1 CG yes analysis CG_deepsignal.11_arab_balance_seq_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.11_arab_balance_seq_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, denoise signal_bilstm, seq+signal =====
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn11_sn16.denoise_signal_bilstm.both_bilstm/both_bilstm.b11_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.denoise_signal_bilstm.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 11 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 15 --nproc_gpu 3 > athaliana.guppy.pass.part2.CG.bn11_sn16.denoise_signal_bilstm.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn11_sn16.denoise_signal_bilstm.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn11_sn16.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn11_sn16.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.11_arab_denoise_signal_bilstm_both_bilstm.10x_1 CG yes analysis CG_deepsignal.11_arab_denoise_signal_bilstm_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.11_arab_denoise_signal_bilstm_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &





# ======= CG pipeline (17/16) ================================
# extract_features =====
# part1/0-601 CG positive
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --seq_len 17 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --seq_len 17 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_neg.log &

# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.4499883, 0.1753386, 0.2664226, 0.1082505]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 199781 kmers from kmer2ratio file
# for 3136 common kmers, fill 162468 samples, 8972 samples that can't filled
# totalline: 10000000, need to fill: 9837532
# extract 11042225 samples from 2208672 diff kmers

# get train data, balance =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv &
# wc -l athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv
head -21194693 athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv
# train model, balance, seq+signal
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn17_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn17_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 17 > model.dp2.CG.arab_R9.4plus_tem.bn17_sn16.balance.both_bilstm.train.log 2>&1 &

# test model, balance, seq+signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn17_sn16.balance.both_bilstm/both_bilstm.b17_s16_epoch7.ckpt --result_file athaliana.guppy.pass.part2.CG.bn17_sn16.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 17 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn17_sn16.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn17_sn16.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn17_sn16.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn17_sn16.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.17_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.17_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.17_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &




# ======= CG pipeline (15/16) ================================
# extract_features =====
# part1/0-601 CG positive
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --seq_len 15 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --seq_len 15 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_neg.log &

# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.449606, 0.1754845, 0.2665972, 0.1083123]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 194811 kmers from kmer2ratio file
# for 15818 common kmers, fill 806922 samples, 74252 samples that can't filled
# totalline: 10000000, need to fill: 9193078
# extract 10342339 samples from 2068670 diff kmers

# get train data, balance =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv &
# wc -l athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv
# 21149261 athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv
head -21139261 athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv
# train model, balance, seq+signal
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn15_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn15_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 15 > model.dp2.CG.arab_R9.4plus_tem.bn15_sn16.balance.both_bilstm.train.log 2>&1 &

# test model, balance, seq+signal =====
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn15_sn16.balance.both_bilstm/both_bilstm.b15_s16_epoch5.ckpt --result_file athaliana.guppy.pass.part2.CG.bn15_sn16.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 15 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn15_sn16.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn15_sn16.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn15_sn16.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn15_sn16.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.15_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.15_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.15_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &



# ======= CG pipeline (13/16) ================================
# extract_features =====
# part1/0-601 CG positive
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --seq_len 13 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --seq_len 13 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.log &

# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --num_lines 10000000 --header false &
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.449687, 0.1753789, 0.2667101, 0.108224]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 179532 kmers from kmer2ratio file
# for 88507 common kmers, fill 4827716 samples, 262584 samples that can't filled
# totalline: 10000000, need to fill: 5172284
# extract 6277131 samples from 1255509 diff kmers

# get train data, random =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.r10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.tsv &
head -19990000 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.valid.tsv
# train seq+siganl
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.r20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.random.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.random.both_bilstm.train.log 2>&1 &
# get train data, balance =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv &
# 21104847 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv
head -21094847 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv
# train model, balance, seq+signal
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.both_bilstm.train.log 2>&1 &
# train model, balance, seq
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.seq_bilstm --model_type seq_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.seq_bilstm.train.log 2>&1 &
# train model, balance, signal
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.signal_bilstm --model_type signal_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.signal_bilstm.train.log 2>&1 &
# denoise signal, train seq+signal =======
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/denoise.py --train_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv --is_filter_fn no --model_type signal_bilstm --is_base yes --epoch_num 3 --rounds 3 --iterations 8 --seq_len 13 > denoise.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.signal_bilstm.iter8_rd3_ep3.log 2>&1 &
# 16512362 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.tsv
head -16502362 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.tsv > athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.valid.tsv
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.signal_bilstm.denoise_fp7.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.denoise_signal_bilstm.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.denoise_signal_bilstm.both_bilstm.train.log 2>&1 &



# test model, random, seq+signal =====
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.random.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.random.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.random.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.random.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.random.both_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.random.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.13_arab_random_both_bilstm.10x_1 CG yes analysis CG_deepsignal.13_arab_random_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.13_arab_random_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, balance, seq+signal =====
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch4.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.13_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.13_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.13_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# combine fwd+rev
python ~/tools/deepsignal2/scripts/combine_two_strands_frequency.py --frequency_fp athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.freq.tsv --ref_fp GCF_000001735.4_TAIR10.1_genomic.fna &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.fb_combined.txt athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.freq.fb_combined.tsv bisulfite.rep1.fb deepsignal.13_arab_balance_both_bilstm.10x_1_fb CG yes analysis CG_deepsignal.13_arab_balance_both_bilstm.10x_1_fb_vs_bisulfite.rep1.fb.tsv > analysis/CG_deepsignal.13_arab_balance_both_bilstm.10x_1_fb_vs_bisulfite.rep1.fb.log 2>&1 &
# - rice shuidao1-1 10x.1
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.1 --model_path model.arab/model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch4.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arab.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arab.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arab.balance.both_bilstm.10x_1.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arab.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arab.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.rice1-1.13_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.rice1-1.13_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.rice1-1.13_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, balance, seq =====
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.seq_bilstm/seq_bilstm.b13_s16_epoch4.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.seq_bilstm.10x_1.tsv --model_type seq_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.balance.seq_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.balance.seq_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.seq_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.balance.seq_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.13_arab_balance_seq_bilstm.10x_1 CG yes analysis CG_deepsignal.13_arab_balance_seq_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.13_arab_balance_seq_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, balance, signal =====
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.balance.signal_bilstm/signal_bilstm.b13_s16_epoch7.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.signal_bilstm.10x_1.tsv --model_type signal_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.balance.signal_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.balance.signal_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.signal_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.balance.signal_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.13_arab_balance_signal_bilstm.10x_1 CG yes analysis CG_deepsignal.13_arab_balance_signal_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.13_arab_balance_signal_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test model, denoise signal_bilstm, signal =====
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn13_sn16.denoise_signal_bilstm.both_bilstm/both_bilstm.b13_s16_epoch4.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.denoise_signal_bilstm.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.denoise_signal_bilstm.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.denoise_signal_bilstm.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.13_arab_denoise_signal_bilstm_both_bilstm.10x_1 CG yes analysis CG_deepsignal.13_arab_denoise_signal_bilstm_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.13_arab_denoise_signal_bilstm_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &


# combine arab+rice2-1 to train ====
# arab 10m pos+neg ==
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.r5m.tsv --num_lines 5000000 --header false &
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.r5m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.4498146, 0.1751442, 0.2667586, 0.1082826]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --krfile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.r5m.kmer_distri.tsv --totalline 5000000 --wfile athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_5m.tsv &
# 179478 kmers from kmer2ratio file
# for 88483 common kmers, fill 2469504 samples, 79060 samples that can't filled
# totalline: 5000000, need to fill: 2530496
# extract 3766449 samples from 1255533 diff kmers
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_pos.r5m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_5m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_10m.tsv &
# rice2-1 10m pos+neg ==
python /home/nipeng/tools/m5c_arab/filter_samples_by_label.py --sf_path shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv --midfix label0 --label 0 &
python /home/nipeng/tools/m5c_arab/filter_samples_by_label.py --sf_path shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.tsv --midfix label1 --label 1 &
python /home/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label1.tsv --write_filepath shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label1.r5m.tsv --num_lines 5000000 --header false &
python /home/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label1.r5m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.2568006, 0.2553064, 0.2919978, 0.1958952]
python /home/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label0.tsv --krfile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label1.r5m.kmer_distri.tsv --totalline 5000000 --wfile shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label0.asposkmer_5m.tsv &
# 1449984 kmers from kmer2ratio file
# for 1091307 common kmers, fill 3855944 samples, 122875 samples that can't filled
# totalline: 5000000, need to fill: 1144056
# extract 1716650 samples from 1716650 diff kmers
python /home/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label1.r5m.tsv --fp2 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_20m.label0.asposkmer_5m.tsv --concated_fp shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_10m.tsv &
# = comb
python /home/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs_inter.hc_poses_all.balanced_10m.tsv --fp2 shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balanced_10m.tsv --concated_fp athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.tsv &
# 21808547 athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.tsv
head -21798547 athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.tsv > athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.train.tsv
tail -10000 athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.tsv > athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.valid.tsv
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.train.tsv --valid_file athaliana_shuidao2-1.guppy.pass.part1.CG.bn13_sn16.signal_feas.bs.hc_poses_all.balance.20m.valid.tsv --model_dir model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 13 > model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm.train.log 2>&1 &

# test arab.10x =
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path /homeb/nipeng/data/comb_arab_n_rice2-1/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.arab.13_arabnrice2-1_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.arab.13_arabnrice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.arab.13_arabnrice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test arab.10x.2 =
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.2 --model_path /homeb/nipeng/data/comb_arab_n_rice2-1/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.freq.tsv &
# test arab.10x.3 =
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.3 --model_path /homeb/nipeng/data/comb_arab_n_rice2-1/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.freq.tsv &
# test arab.10x.4 =
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.4 --model_path /homeb/nipeng/data/comb_arab_n_rice2-1/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.freq.tsv &
# test arab.10x.5 =
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.5 --model_path /homeb/nipeng/data/comb_arab_n_rice2-1/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.tsv --result_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.freq.tsv &
# test arab.50x_12345
python ~/tools/deepsignal2/scripts/combine_call_mods_freq_files.py --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.freq.tsv --modsfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.freq.tsv --wfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.tsv bisulfite.rep1 deepsignal.arab.13_arabnrice2-1_balance_both_bilstm.50x_12345 CG yes analysis CG_deepsignal.arab.13_arabnrice2-1_balance_both_bilstm.50x_12345_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.arab.13_arabnrice2-1_balance_both_bilstm.50x_12345_vs_bisulfite.rep1.log 2>&1 &

# test rice1-1.10x =
CUDA_VISIBLE_DEVICES=2 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.1 --model_path model.comb/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.rice1-1.13_arabnrice2-1_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.rice1-1.13_arabnrice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.rice1-1.13_arabnrice2-1_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &
# test rice1-1.10x.2 =
CUDA_VISIBLE_DEVICES=0 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.2 --model_path model.comb/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.freq.tsv &
# test rice1-1.10x.3 =
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.3 --model_path model.comb/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.freq.tsv &
# test rice1-1.10x.4 =
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.4 --model_path model.comb/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 25 --nproc_gpu 5 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.freq.tsv &
# test rice1-1.10x.5 =
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path fast5_pass2.single/10x.5 --model_path model.comb/model.dp2.CG.arabnrice2-1_R9.4plus_tem.bn13_sn16.balance.both_bilstm/both_bilstm.b13_s16_epoch6.ckpt --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.tsv --model_type both_bilstm --seq_len 13 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_001 --reference_path Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.tsv --result_file shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.freq.tsv &
# test rice1-1.50x_12345
python ~/tools/deepsignal2/scripts/combine_call_mods_freq_files.py --modsfile shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv --modsfile shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_2.freq.tsv --modsfile shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_3.freq.tsv --modsfile shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_4.freq.tsv --modsfile shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_5.freq.tsv --wfile shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt shuidao1-1.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.tsv bisulfite.rep1 deepsignal.rice1-1.13_arabnrice2-1_balance_both_bilstm.50x_12345 CG yes analysis CG_deepsignal.rice1-1.13_arabnrice2-1_balance_both_bilstm.50x_12345_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.rice1-1.13_arabnrice2-1_balance_both_bilstm.50x_12345_vs_bisulfite.rep1.log 2>&1 &





# ======= CG pipeline (9/16) ================================
# extract_features =====
# part1/0-601 CG positive
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 1 --seq_len 9 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.9.cov_5_50000.hc_poses_positive.inter.tsv > athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.log &
# part1/0-601 CG negative
python /homeb/nipeng/tools/deepsignal2/deepsignal2/extract_features.py --fast5_dir reads_single_pass/part1 --corrected_group RawGenomeCorrected_002 --basecall_subgroup BaseCalled_template --reference_path /homeb/nipeng/data/genome/athaliana/GCF_000001735.4_TAIR10.1_genomic.fna --methy_label 0 --seq_len 9 --signal_len 16 --motifs CG --mod_loc 0 --write_path athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --nproc 20 --positions bs.poses/ninanjie-2.CG.1n2n3.rmet_0.0.cov_5_50000.hc_poses_negative.inter.tsv > athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_neg.log &

# fliter samples, random =====
python /homeb/nipeng/tools/basemods_RNN_signal/randsel_file_rows.py --ori_filepath athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.tsv --write_filepath athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --num_lines 10000000 --header false &
# fliter samples, negkmeraspos =====
python /homeb/nipeng/tools/deepsignal2/scripts/get_kmer_dist_of_feafile.py --feafile athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv &
# ['CGA', 'CGC', 'CGG', 'CGT'] [0.4497627, 0.1754766, 0.2664961, 0.1082646]
python /homeb/nipeng/tools/deepsignal2/scripts/select_neg_samples_by_kmer_distri.py --feafile athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_neg.tsv --krfile athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.kmer_distri.tsv --totalline 10000000 --wfile athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv &
# 15455 kmers from kmer2ratio file
# for 15455 common kmers, fill 9946780 samples, 54105 samples that can't filled
# totalline: 10000000, need to fill: 53220
# extract 53824 samples from 928 diff kmers

# get train data, balance =====
python /homeb/nipeng/tools/deepsignal2/scripts/concat_two_files.py --fp1 athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_pos.tsv.r10m.tsv --fp2 athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_neg.asposkmer_10m.tsv --concated_fp athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv &
# 20000604 athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv
head -19990604 athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv
tail -10000 athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.tsv > athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv
# train model, balance, seq+signal
CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/train.py --train_file athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.train.tsv --valid_file athaliana.guppy.pass.part1.CG.bn9_sn16.signal_feas.bs_inter.hc_poses_all.balanced_20m.valid.tsv --model_dir model.dp2.CG.arab_R9.4plus_tem.bn9_sn16.balance.both_bilstm --model_type both_bilstm --layernum1 3 --layernum2 1 --seq_len 9 > model.dp2.CG.arab_R9.4plus_tem.bn9_sn16.balance.both_bilstm.train.log 2>&1 &

CUDA_VISIBLE_DEVICES=1 python ~/tools/deepsignal2/deepsignal2/call_modifications.py --input_path reads_single_pass/part2/100x/10x.1 --model_path model.dp2.CG.arab_R9.4plus_tem.bn9_sn16.balance.both_bilstm/both_bilstm.b9_s16_epoch5.ckpt --result_file athaliana.guppy.pass.part2.CG.bn9_sn16.balance.both_bilstm.10x_1.tsv --model_type both_bilstm --seq_len 9 --signal_len 16 --layernum1 3 --layernum2 1 --corrected_group RawGenomeCorrected_002 --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --motifs CG --mod_loc 0 --f5_batch_size 20 --nproc 20 --nproc_gpu 4 > athaliana.guppy.pass.part2.CG.bn9_sn16.balance.both_bilstm.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.bn9_sn16.balance.both_bilstm.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.bn9_sn16.balance.both_bilstm.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.bn9_sn16.balance.both_bilstm.10x_1.freq.tsv bisulfite.rep1 deepsignal.9_arab_balance_both_bilstm.10x_1 CG yes analysis CG_deepsignal.9_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.9_arab_balance_both_bilstm.10x_1_vs_bisulfite.rep1.log 2>&1 &







# guppy basecalled ================================================
# megalodon test ==================
# arab 10x.1 CG
# example = 
# --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
# --write-mods-text \
# --mod-aggregate-method binary_threshold \
# --mod-binary-threshold 0.75 \
# --mod-output-formats bedmethyl wiggle \
megalodon reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
    --guppy-timeout 30 \
    --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 2 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CG \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_1.CG.log 2>&1 &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt megalodon_results.arab.pass.part2_guppy.10x_1.CG/modified_bases.5mC.bed bisulfite.rep1 megalodon.guppy.arab.10x.1 CG yes analysis CG_megalodon.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CG_megalodon.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# =
python ~/tools/plant_5mC_analysis/tools_to_cmp/combine_cpg_two_strands_methy_freqs.py --report_fp megalodon_results.arab.pass.part2_guppy.10x_1.CG/modified_bases.5mC.bed -t bedmethyl --ref_fp GCF_000001735.4_TAIR10.1_genomic.fna &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.fb_combined.txt megalodon_results.arab.pass.part2_guppy.10x_1.CG/modified_bases.5mC.fb_combined.bed bisulfite.rep1.fb megalodon.guppy.arab.10x.1_fb CG yes analysis CG_megalodon.guppy.arab.10x.1_fb_vs_bisulfite.rep1.fb.tsv yes yes > analysis/CG_megalodon.guppy.arab.10x.1_fb_vs_bisulfite.rep1.fb.log 2>&1 &
# arab 10x.1 CHG
megalodon reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
    --guppy-timeout 30 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 2 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHG \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_1.CHG.log 2>&1 &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt megalodon_results.arab.pass.part2_guppy.10x_1.CHG/modified_bases.5mC.bed bisulfite.rep1 megalodon.guppy.arab.10x.1 CHG yes analysis CHG_megalodon.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CHG_megalodon.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# arab 10x.1 CHH
megalodon reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
    --guppy-timeout 30 \
    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif Z CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 2 --processes 20 --output-directory megalodon_results.arab.pass.part2_guppy.10x_1.CHH \
    --overwrite > megalodon_results.arab.pass.part2_guppy.10x_1.CHH.log 2>&1 &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt megalodon_results.arab.pass.part2_guppy.10x_1.CHH/modified_bases.5mC.bed bisulfite.rep1 megalodon.guppy.arab.10x.1 CHH yes analysis CHH_megalodon.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CHH_megalodon.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# ---- arab 10x.1 CHG new model
megalodon reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_prom_modbases_5mC_v001.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CHG 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 2 --processes 20 --output-directory megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHG.new \
    --overwrite --disable-mod-calibration > megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHG.new.log 2>&1 &
python ~/tools/plant_5mC_analysis/tools_to_cmp/correlation_with_bs.py --nano_file megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHG.new/modified_bases.5mC.bed --bs_file bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt --contig_prefix NC_003 --cov_cf 5 > megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHG.new.vs_bsrep1.cov5.nc_003.log &
# ---- arab 10x.1 CHH new model
megalodon reads_single_pass/part2/100x/10x.1/ \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server \
    --guppy-params "-d /home/nipeng/tools/rerio/basecall_models/ --num_callers 2" \
    --guppy-timeout 100 \
    --guppy-config res_dna_r941_prom_modbases_5mC_v001.cfg \
    --outputs per_read_mods mod_mappings mods \
    --reference GCF_000001735.4_TAIR10.1_genomic.fna --mod-motif m CHH 0 \
    --write-mods-text \
    --mod-output-formats bedmethyl wiggle \
    --devices 2 --processes 20 --output-directory megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHH.new \
    --overwrite --disable-mod-calibration > megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHH.new.log 2>&1 &
python ~/tools/plant_5mC_analysis/tools_to_cmp/correlation_with_bs.py --nano_file megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHH.new/modified_bases.5mC.bed --bs_file bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt --contig_prefix NC_003 --cov_cf 5 > megalodon_results/megalodon_results.arab.pass.part2_guppy.10x_1.CHH.new.vs_bsrep1.cov5.nc_003.log &





# tombo test (guppy) =================
# dp2, 192.16.10.53
# already resquiggled guppy
# tombo resquiggle reads_single_pass/part2/100x/10x.1 GCF_000001735.4_TAIR10.1_genomic.fna --processes 20 --corrected-group RawGenomeCorrected_002 --basecall-group Basecall_1D_001 --overwrite --ignore-read-locks > tombo_arab_pass.part2_10x.1.guppy.log 2>&1 &
# arab 10x.1 CG
tombo detect_modifications alternative_model --fast5-basedirs reads_single_pass/part2/100x/10x.1 --alternate-bases CpG --statistics-file-basename tombo_results.arab.pass.part2_guppy.10x.1.CG --dna --multiprocess-region-size 1000 --processes 20 --corrected-group RawGenomeCorrected_002 > tombo_results.arab.pass.part2_guppy.10x.1.CG.log 2>&1 &
tombo text_output browser_files \
   --fast5-basedirs reads_single_pass/part2/100x/10x.1 \
   --statistics-filename tombo_results.arab.pass.part2_guppy.10x.1.CG.CpG.tombo.stats \
   --file-types coverage dampened_fraction fraction\
   --browser-file-basename tombo_results.arab.pass.part2_guppy.10x.1.CG --corrected-group RawGenomeCorrected_002 \
   > tombo_results.arab.pass.part2_guppy.10x.1.CG.output.log 2>&1 &
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig.bed
# wig2bed --multisplit bar < tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig.bed
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.plus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.plus.wig.bed
python ~/tools/plant_5mC_analysis/tools_to_cmp/tombo_reformat_bed.py --cov_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.coverage.plus.bedgraph --cov_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.coverage.minus.bedgraph --rmet_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.plus.wig.bed --rmet_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.minus.wig.bed --wfile tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.methyl.bed &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.methyl.bed bisulfite.rep1 tombo.guppy.arab.10x.1 CG yes analysis CG_tombo.guppy.arab.10x.1_vs_bisulfite.rep1.tsv yes yes > analysis/CG_tombo.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# = fb_combined
python ~/tools/plant_5mC_analysis/tools_to_cmp/combine_cpg_two_strands_methy_freqs.py --report_fp tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.methyl.bed -t bedmethyl --ref_fp GCF_000001735.4_TAIR10.1_genomic.fna &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.fb_combined.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.dampened_fraction_modified_reads.methyl.fb_combined.bed bisulfite.rep1.fb tombo.guppy.arab.10x.1_fb CG yes analysis CG_tombo.guppy.arab.10x.1_fb_vs_bisulfite.rep1.fb.tsv yes yes > analysis/CG_tombo.guppy.arab.10x.1_fb_vs_bisulfite.rep1.fb.log 2>&1 &
# = no dampened
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.minus.wig.bed
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.plus.wig > tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.plus.wig.bed
python ~/tools/plant_5mC_analysis/tools_to_cmp/tombo_reformat_bed.py --cov_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.coverage.plus.bedgraph --cov_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.coverage.minus.bedgraph --rmet_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.plus.wig.bed --rmet_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.minus.wig.bed --wfile tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.methyl.bed &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.CG.fraction_modified_reads.methyl.bed bisulfite.rep1 tombo_nodampened.guppy.arab.10x.1 CG yes analysis CG_tombo_nodampened.guppy.arab.10x.1_vs_bisulfite.rep1.tsv yes yes > analysis/CG_tombo_nodampened.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# arab 10x.1 5mC
tombo detect_modifications alternative_model --fast5-basedirs reads_single_pass/part2/100x/10x.1 --alternate-bases 5mC --statistics-file-basename tombo_results.arab.pass.part2_guppy.10x.1.5mC --dna --multiprocess-region-size 1000 --processes 20 --corrected-group RawGenomeCorrected_002 > tombo_results.arab.pass.part2_guppy.10x.1.5mC.log 2>&1 &
tombo text_output browser_files \
   --fast5-basedirs reads_single_pass/part2/100x/10x.1 \
   --statistics-filename tombo_results.arab.pass.part2_guppy.10x.1.5mC.5mC.tombo.stats \
   --file-types coverage dampened_fraction fraction\
   --browser-file-basename tombo_results.arab.pass.part2_guppy.10x.1.5mC --corrected-group RawGenomeCorrected_002 \
   > tombo_results.arab.pass.part2_guppy.10x.1.5mC.output.log 2>&1 &
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.minus.wig.bed
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.plus.wig > tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.plus.wig.bed
python ~/tools/plant_5mC_analysis/tools_to_cmp/tombo_reformat_bed.py --cov_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.coverage.plus.bedgraph --cov_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.coverage.minus.bedgraph --rmet_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.plus.wig.bed --rmet_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.minus.wig.bed --wfile tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.methyl.bed &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.methyl.bed bisulfite.rep1 tombo.guppy.arab.10x.1 CHG yes analysis CHG_tombo.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CHG_tombo.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.dampened_fraction_modified_reads.methyl.bed bisulfite.rep1 tombo.guppy.arab.10x.1 CHH yes analysis CHH_tombo.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CHH_tombo.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
# = no dampened
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.minus.wig > tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.minus.wig.bed
wig2bed --multisplit < tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.plus.wig > tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.plus.wig.bed
python ~/tools/plant_5mC_analysis/tools_to_cmp/tombo_reformat_bed.py --cov_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.coverage.plus.bedgraph --cov_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.coverage.minus.bedgraph --rmet_plus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.plus.wig.bed --rmet_minus tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.minus.wig.bed --wfile tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.methyl.bed &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHG.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.methyl.bed bisulfite.rep1 tombo_nodampened.guppy.arab.10x.1 CHG yes analysis CHG_tombo_nodampened.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CHG_tombo_nodampened.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CHH.txt tombo_results.arab.pass.part2_guppy.10x.1/tombo_results.arab.pass.part2_guppy.10x.1.5mC.fraction_modified_reads.methyl.bed bisulfite.rep1 tombo_nodampened.guppy.arab.10x.1 CHH yes analysis CHH_tombo_nodampened.guppy.arab.10x.1_vs_bisulfite.rep1.tsv > analysis/CHH_tombo_nodampened.guppy.arab.10x.1_vs_bisulfite.rep1.log 2>&1 &




# deepsignal guppy ==================
# arab 10x.1 CG
CUDA_VISIBLE_DEVICES=1 deepsignal call_mods --input_path reads_single_pass/part2/100x/10x.1 --model_path model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+/bn_17.sn_360.epoch_9.ckpt --reference_path GCF_000001735.4_TAIR10.1_genomic.fna --corrected_group RawGenomeCorrected_002 --nproc 20 --is_gpu yes --result_file athaliana.guppy.pass.part2.CG.deepsignal.10x_1.tsv > athaliana.guppy.pass.part2.CG.deepsignal.10x_1.log 2>&1 &
python ~/tools/deepsignal2/scripts/call_modification_frequency.py --input_path athaliana.guppy.pass.part2.CG.deepsignal.10x_1.tsv --result_file athaliana.guppy.pass.part2.CG.deepsignal.10x_1.freq.tsv &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.deepsignal.10x_1.freq.tsv bisulfite.rep1 deepsignal.arab.dp1.10x_1 CG yes analysis CG_deepsignal.arab.dp1.10x_1_vs_bisulfite.rep1.tsv > analysis/CG_deepsignal.arab.dp1.10x_1_vs_bisulfite.rep1.log 2>&1 &
# combine fwd+rev
python ~/tools/deepsignal2/scripts/combine_two_strands_frequency.py --frequency_fp athaliana.guppy.pass.part2.CG.deepsignal.10x_1.freq.tsv --ref_fp GCF_000001735.4_TAIR10.1_genomic.fna &
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.fb_combined.txt athaliana.guppy.pass.part2.CG.deepsignal.10x_1.freq.fb_combined.tsv bisulfite.rep1.fb deepsignal.arab.dp1.10x_1_fb CG yes analysis CG_deepsignal.arab.dp1.10x_1_fb_vs_bisulfite.rep1.fb.tsv > analysis/CG_deepsignal.arab.dp1.10x_1_fb_vs_bisulfite.rep1.fb.log 2>&1 &


# nanopolish guppy =====================
python ~/tools/NanoSigFeaGen/test/test_combine_fb_of_report.py --report_fp bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt -t bs --ref_fp GCF_000001735.4_TAIR10.1_genomic.fna &
# guppy basecall reads_single_pass/part2_10x_1.guppy_r.fastq/
# guppy_basecaller -i reads_single_pass/part2/100x/10x.1 -r -s reads_single_pass/part2_10x_1.guppy_r.fastq --config dna_r9.4.1_450bps_hac_prom.cfg --gpu_runners_per_device 2 --chunks_per_runner 2500 --device CUDA:1 > reads_single_pass/part2_10x_1.guppy_prom.log
/usr/bin/time -v bash nanopolish_call_methy.sh reads_single_pass/part2/100x/10x.1 reads_single_pass/part2_10x_1.guppy_r.fastq GCF_000001735.4_TAIR10.1_genomic.fna 20 athaliana.guppy.pass.part2.CG.nanopolish.10x_1 > athaliana.guppy.pass.part2.CG.nanopolish.10x_1.log 2>&1 &
python ~/tools/nanopolish/scripts/calculate_methylation_frequency.py -s athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.tsv > athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.fb_combined.txt athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.freq.tsv bisulfite.rep1.fb nanopolish.arab.split.10x_1 CG yes analysis CG_nanopolish.arab.split.10x_1_vs_bisulfite.rep1.fb.tsv yes yes > analysis/CG_nanopolish.arab.split.10x_1_vs_bisulfite.rep1.fb.log 2>&1 &
# =
python ~/tools/nanopolish/scripts/calculate_methylation_frequency.s.py -s athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.tsv > athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.freq.strand_specific.tsv &
# =
python nanopolish_add_rev_freq.py --freq_file athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.freq.tsv
Rscript correlation_with_bs.cal_plot.general.R bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt athaliana.guppy.pass.part2.CG.nanopolish.10x_1.fastq.methyl_calls.freq.add_rev.tsv bisulfite.rep1 nanopolish.arab.split_addrev.10x_1 CG yes analysis CG_nanopolish.arab.split_addrev.10x_1_vs_bisulfite.rep1.tsv yes yes > analysis/CG_nanopolish.arab.split_addrev.10x_1_vs_bisulfite.rep1.log 2>&1 &


# deepsignal2 arab+rice->arab.10x.1 fb_combined ===
python ~/tools/deepsignal2/scripts/combine_two_strands_frequency.py --frequency_fp athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv --ref_fp GCF_000001735.4_TAIR10.1_genomic.fna &
















