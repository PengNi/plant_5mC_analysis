#! /usr/bin/python
import argparse
import os
import pandas as pd
import scipy.stats
import numpy as np
from sklearn.metrics import mean_squared_error
import math


sep = "||"


def read_nanopolish_freqfile(nanopolish_file, contig_prefix, contig_names, cov_cf):
    rmet_nano = pd.read_csv(nanopolish_file, sep="\t", header=0,
                            dtype={"chromosome": str})
    rmet_nano = rmet_nano.rename(columns={"start": "pos", "called_sites": "coverage",
                                 "methylated_frequency": "Rmet"})
    if contig_prefix is not None:
        rmet_nano = rmet_nano[rmet_nano.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_nano = rmet_nano[rmet_nano.chromosome.isin(contigset)]
    else:
        pass
    rmet_nano['key'] = rmet_nano.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_nano = rmet_nano[["chromosome", "pos", "coverage", "Rmet", "key"]]

    meancov = rmet_nano["coverage"].mean()
    rmet_nano = rmet_nano[rmet_nano["coverage"] >= cov_cf]
    return meancov, rmet_nano.sort_values(by=['chromosome', 'pos'])


def read_methylbed(bed_file, contig_prefix, contig_names, cov_cf):
    rmet_bed = pd.read_csv(bed_file, sep="\t", header=None,
                           names=["chromosome", "pos", "end", "na1", "na2", "strand",
                                  "na3", "na4", "na5", "coverage", "rpercent"],
                           dtype={"chromosome": str})
    rmet_bed["Rmet"] = rmet_bed.apply(lambda row: row["rpercent"] / 100.0, axis=1)
    if contig_prefix is not None:
        rmet_bed = rmet_bed[rmet_bed.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_bed = rmet_bed[rmet_bed.chromosome.isin(contigset)]
    else:
        pass
    rmet_bed['key'] = rmet_bed.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_bed = rmet_bed[["chromosome", "pos", "coverage", "Rmet", "key"]]

    meancov = rmet_bed["coverage"].mean()
    rmet_bed = rmet_bed[rmet_bed["coverage"] >= cov_cf]
    return meancov, rmet_bed.sort_values(by=['chromosome', 'pos'])


def read_methylbed_fb(bed_file, contig_prefix, contig_names, cov_cf):
    rmet_bed = pd.read_csv(bed_file, sep="\t", header=0,
                           dtype={"chromosome": str})
    if contig_prefix is not None:
        rmet_bed = rmet_bed[rmet_bed.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_bed = rmet_bed[rmet_bed.chromosome.isin(contigset)]
    else:
        pass
    rmet_bed['key'] = rmet_bed.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_bed = rmet_bed[["chromosome", "pos", "coverage", "Rmet", "key"]]

    meancov = rmet_bed["coverage"].mean()
    rmet_bed = rmet_bed[rmet_bed["coverage"] >= cov_cf]
    return meancov, rmet_bed.sort_values(by=['chromosome', 'pos'])


def read_rmetfile_of_dp2(dp2_file, contig_prefix, contig_names, cov_cf):
    rmet_dp2 = pd.read_csv(dp2_file, sep="\t", header=None,
                           names=["chromosome", "pos", "strand", "pos_in_strand", "prob0", "prob1", "met", "unmet",
                                  "coverage", "Rmet", "kmer"],
                           dtype={"chromosome": str})

    if contig_prefix is not None:
        rmet_dp2 = rmet_dp2[rmet_dp2.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_dp2 = rmet_dp2[rmet_dp2.chromosome.isin(contigset)]
    else:
        pass
    rmet_dp2['key'] = rmet_dp2.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_dp2 = rmet_dp2[["chromosome", "pos", "coverage", "Rmet", "key"]]

    meancov = rmet_dp2["coverage"].mean()
    rmet_dp2 = rmet_dp2[rmet_dp2["coverage"] >= cov_cf]

    return meancov, rmet_dp2.sort_values(by=['chromosome', 'pos'])


def read_rmetfile_of_bs(bs_file, contig_prefix, contig_names, cov_cf):
    if bs_file.endswith(".bed"):
        rmet_bs = pd.read_csv(bs_file, sep="\t", header=None,
                              names=["chromosome", "pos", "end", "na1", "na2", "strand",
                                     "na3", "na4", "na5", "coverage", "rpercent"],
                              dtype={"chromosome": str})
        rmet_bs["Rmet"] = rmet_bs.apply(lambda row: row["rpercent"] / 100.0, axis=1)
    else:
        rmet_bs = pd.read_csv(bs_file, sep="\t", header=0, dtype={"chromosome": str})

    if contig_prefix is not None:
        rmet_bs = rmet_bs[rmet_bs.apply(lambda row: row["chromosome"].startswith(contig_prefix), axis=1)]
    elif contig_names is not None:
        contigset = pd.Series(contig_names.split(","))
        rmet_bs = rmet_bs[rmet_bs.chromosome.isin(contigset)]
    else:
        pass

    rmet_bs['key'] = rmet_bs.apply(lambda row: row["chromosome"] + sep + str(row["pos"]), axis=1)
    rmet_bs = rmet_bs[["chromosome", "pos", "coverage", "Rmet", "key"]]

    meancov = rmet_bs["coverage"].mean()
    rmet_bs = rmet_bs[rmet_bs["coverage"] >= cov_cf]

    return meancov, rmet_bs.sort_values(by=['chromosome', 'pos'])


def cal_corr_df1_vs_df2(df1, df2):
    df1_inter = df1[df1.key.isin(df2.key)].sort_values(by=['chromosome', 'pos'])
    df2_inter = df2[df2.key.isin(df1.key)].sort_values(by=['chromosome', 'pos'])
    # df1_inter["Rmet"].corr(df2_inter['Rmet'], method='pearson'), wrong? 0.2660 vs scipy 0.9459
    df1_array, df2_array = np.array(list(df1_inter["Rmet"])), np.array(list(df2_inter["Rmet"]))
    pcorr, _ = scipy.stats.pearsonr(df1_array, df2_array)  # pearson
    scorr, _ = scipy.stats.spearmanr(df1_array, df2_array)  # spearman
    _, _, r_value, _, _ = scipy.stats.linregress(df1_array, df2_array)
    r_square = r_value ** 2  # coefficient of determination
    rmse = math.sqrt(mean_squared_error(df2_array, df1_array))  # RMSE

    return len(df1.index), len(df2.index), len(df1_inter.index), pcorr, scorr, r_square, rmse


def _read_rmet_file(dp2_file, args):
    if str(dp2_file).endswith("fb_combined.bed"):
        mean_cov, dp2rmetinfo = read_methylbed_fb(dp2_file, args.contig_prefix, args.contig_names, args.cov_cf)
    elif str(dp2_file).endswith(".bed"):
        mean_cov, dp2rmetinfo = read_methylbed(dp2_file, args.contig_prefix, args.contig_names, args.cov_cf)
    elif "nanopolish" in str(dp2_file):
        mean_cov, dp2rmetinfo = read_nanopolish_freqfile(dp2_file, args.contig_prefix, args.contig_names,
                                                         args.cov_cf)
    elif str(dp2_file).endswith("freq.tsv") or str(dp2_file).endswith("freq.fb_combined.tsv"):
        mean_cov, dp2rmetinfo = read_rmetfile_of_dp2(dp2_file, args.contig_prefix,
                                                     args.contig_names, args.cov_cf)
    else:
        return None, None
    return mean_cov, dp2rmetinfo


def correlation_rmets(args):
    nano_files = args.nano_file

    # print("==coverage cutoff: {}\n".format(args.cov_cf))
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("rmet1", "rmet2", "nanonum1", "nanonum2", "internum", "pearson",
                                                      "spearman", "rsquare", "RMSE"))
    for i in range(len(nano_files)):
        # print("====== {}".format(dp2_file))
        mean_cov1, dp2rmetinfo1 = _read_rmet_file(nano_files[i], args)
        # print("mean_covarge: {}".format(mean_cov))
        for j in range(i+1, len(nano_files)):
            mean_cov1, dp2rmetinfo2 = _read_rmet_file(nano_files[j], args)
            nanonum1, nanonum2, internum, \
                pcorr, scorr, r_square, rmse = cal_corr_df1_vs_df2(dp2rmetinfo1, dp2rmetinfo2)

            print("{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(nano_files[i], nano_files[j],
                                                                              nanonum1, nanonum2,
                                                                              internum, pcorr, scorr,
                                                                              r_square, rmse))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nano_file", type=str, action="append", required=True, help="nano freq file, .freq.tsv/.bed")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")

    args = parser.parse_args()
    correlation_rmets(args)


if __name__ == '__main__':
    main()
