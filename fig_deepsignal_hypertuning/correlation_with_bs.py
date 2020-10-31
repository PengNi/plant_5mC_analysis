#! /usr/bin/python
import argparse
import os
import pandas as pd
import scipy.stats
import numpy as np


sep = "||"


def read_rmetfile_of_dp2(dp2_file, contig_prefix, contig_names, cov_cf):
    rmet_dp2 = pd.read_csv(dp2_file, sep="\t", header=None,
                           names=["chromosome", "pos", "strand", "pos_in_strand", "prob0", "prob1", "met", "unmet",
                                  "coverage", "Rmet", "kmer"])

    if contig_prefix is not None:
        rmet_dp2 = rmet_dp2[rmet_dp2.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_dp2 = rmet_dp2[rmet_dp2.chromosome.isin(contigset)]
    else:
        pass
    rmet_dp2['key'] = rmet_dp2.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_dp2 = rmet_dp2[["chromosome", "pos", "coverage", "Rmet", "key"]]

    rmet_dp2 = rmet_dp2[rmet_dp2["coverage"] >= cov_cf]

    return rmet_dp2.sort_values(by=['chromosome', 'pos'])


def read_rmetfile_of_bs(bs_file, contig_prefix, contig_names, cov_cf):
    rmet_bs = pd.read_csv(bs_file, sep="\t", header=0)

    if contig_prefix is not None:
        rmet_bs = rmet_bs[rmet_bs.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_bs = rmet_bs[rmet_bs.chromosome.isin(contigset)]
    else:
        pass

    rmet_bs['key'] = rmet_bs.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_bs = rmet_bs[["chromosome", "pos", "coverage", "Rmet", "key"]]

    rmet_bs = rmet_bs[rmet_bs["coverage"] >= cov_cf]

    return rmet_bs.sort_values(by=['chromosome', 'pos'])


def cal_pearson_corr_df1_vs_df2(df1, df2):
    df1_inter = df1[df1.key.isin(df2.key)].sort_values(by=['chromosome', 'pos'])
    df2_inter = df2[df2.key.isin(df1.key)].sort_values(by=['chromosome', 'pos'])
    # df1_inter["Rmet"].corr(df2_inter['Rmet'], method='pearson'), wrong? 0.2660 vs scipy 0.9459
    corr, pvalue = scipy.stats.pearsonr(np.array(list(df1_inter["Rmet"])),
                                        np.array(list(df2_inter["Rmet"])))
    return len(df1_inter), corr


def correlation_with_bs_rmets(args):
    dp2_files = args.dp2_file
    bs_files = args.bs_file

    bs_fname2rmetinfo = dict()
    for bs_file in bs_files:
        bs_fname2rmetinfo[os.path.basename(bs_file)] = read_rmetfile_of_bs(bs_file, args.contig_prefix,
                                                                           args.contig_names, args.cov_cf)
    for dp2_file in dp2_files:
        print("====== {}".format(os.path.basename(dp2_file)))
        dp2rmetinfo = read_rmetfile_of_dp2(dp2_file, args.contig_prefix,
                                           args.contig_names, args.cov_cf)
        print("{}\t{}\t{}".format("bs_file", "sitenum", "pear corr"))
        sitecorrs = []
        sitenums = []
        for bs_fname in sorted(list(bs_fname2rmetinfo.keys())):
            sitenum, sitecorr = cal_pearson_corr_df1_vs_df2(dp2rmetinfo, bs_fname2rmetinfo[bs_fname])
            sitecorrs.append(sitecorr)
            sitenums.append(sitenum)
            print("{}\t{}\t{:.4f}".format(bs_fname, sitenum, sitecorr))
        print("{}\t{:.2f}\t{:.4f}".format("average", sum(sitenums)/len(sitenums), sum(sitecorrs)/len(sitecorrs)))


def main():

    # python ~/tools/plant_5mC_analysis/fig_pipeline_testing/correlation_with_bs.py \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.random.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.denoise_signal_bilstm.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn17_sn16.balance.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn15_sn16.balance.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn11_sn16.balance.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn9_sn16.balance.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.seq_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.balance.signal_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.rice2-1.balance.both_bilstm.10x_1.freq.tsv \
    # --dp2_file athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.10x_1.freq.tsv \
    # --bs_file bs.poses/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.\
    # sorted.CX_report.CG.txt --bs_file bs.poses/D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.\
    # sorted.mark_dup.sorted.CX_report.CG.txt --bs_file bs.poses/D1902828A-ZJ_363363.L2L4_merged.cutadapt.\
    # R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --contig_prefix NC_003 \
    # > athaliana.guppy.pass.part2.CG.10x_1.vs_bs_rep123.stat.pearson_corr.log &
    parser = argparse.ArgumentParser()
    parser.add_argument("--dp2_file", type=str, action="append", required=True, help="dp2 freq file")
    parser.add_argument("--bs_file", type=str, action="append",
                        required=True, help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")

    args = parser.parse_args()
    correlation_with_bs_rmets(args)


if __name__ == '__main__':
    main()
