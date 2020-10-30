#! /usr/bin/python
import argparse
import numpy as np
import scipy.stats
from copy import deepcopy


def _read_one_mod_freq_file(freqfile, args):
    contigset = set(args.contig_names.strip().split(",")) if args.contig_names is not None else None
    freqinfo = {}
    with open(freqfile, "r") as rf:
        for line in rf:
            words = line.strip().split("\t")
            chrom = words[0]
            m_key = "\t".join([words[0], words[1], words[2]])
            # pos_in_strand = words[3]
            # methy_prob = float(words[4])
            # unmethy_prob = float(words[5])
            # methy_cov = int(words[6])
            # unmethy_cov = int(words[7])
            cov = int(words[8])
            rmet = float(words[9])
            # kmer = words[10]

            cnt_flag = 0
            if args.contig_prefix is not None:
                if str(chrom).startswith(args.contig_prefix):
                    cnt_flag = 1
            elif args.contig_names is not None:
                if chrom in contigset:
                    cnt_flag = 1
            else:
                cnt_flag = 1
            if cnt_flag == 1:
                if cov >= args.covcf:
                    freqinfo[m_key] = rmet
    return freqinfo


def _read_one_bs_rmet_file(bs_file, args):
    contigset = set(args.contig_names.strip().split(",")) if args.contig_names is not None else None
    freqinfo = {}
    with open(bs_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage, rmet = words[0], words[1], words[2], int(words[6]), float(words[7])
            mkey = "\t".join([chrom, str(pos), strand])
            cnt_flag = 0
            if args.contig_prefix is not None:
                if str(chrom).startswith(args.contig_prefix):
                    cnt_flag = 1
            elif args.contig_names is not None:
                if chrom in contigset:
                    cnt_flag = 1
            else:
                cnt_flag = 1
            if cnt_flag == 1:
                if coverage >= args.covcf:
                    freqinfo[mkey] = rmet
    return freqinfo


def _cal_pearson_corr_of_rmet1_and_rmet2(freqinfo1, freqinfo2):
    keys_inter = set(freqinfo1.keys()).intersection(set(freqinfo2.keys()))
    keys_inter = sorted(list(keys_inter))
    rmet1, rmet2 = [], []
    for ktmp in keys_inter:
        rmet1.append(freqinfo1[ktmp])
        rmet2.append(freqinfo2[ktmp])
    corr, pvalue = scipy.stats.pearsonr(np.array(rmet1), np.array(rmet2))
    return rmet1, rmet2, len(keys_inter), corr


def cmp_sitesrmet_of_nano_and_bs(args):
    print("==nano: {}".format(args.modsfile))
    freqinfo_nano = _read_one_mod_freq_file(args.modsfile, args)
    freqinfo_bs = []
    for bsfile in args.bsfile:
        freqtmp = _read_one_bs_rmet_file(bsfile, args)
        _, _, sitesnum, corr = _cal_pearson_corr_of_rmet1_and_rmet2(freqinfo_nano, freqtmp)
        print("==bs: {}, sites: {}, corr: {:.4f}".format(bsfile, sitesnum, corr))
        freqinfo_bs.append(freqtmp)
    freqinfo_bs_comb = dict()
    keyset = set()
    for freqinfotmp in freqinfo_bs:
        for fkey in freqinfotmp.keys():
            if fkey not in keyset:
                keyset.add(fkey)
                freqinfo_bs_comb[fkey] = []
            freqinfo_bs_comb[fkey].append(freqinfotmp[fkey])
    for fkey in freqinfo_bs_comb.keys():
        freqinfo_bs_comb[fkey] = sum(freqinfo_bs_comb[fkey]) / len(freqinfo_bs_comb[fkey])
    rmet_nano, rmet_bs, sitesnum, corr = _cal_pearson_corr_of_rmet1_and_rmet2(freqinfo_nano, freqinfo_bs_comb)
    print("==bs: combined, sitesnum: {}, corr: {:.4f}".format(sitesnum, corr))
    wf = open(args.wfile, "w")
    wf.write("\t".join(["rmet_bis", "rmet_nan"]) + "\n")
    for ridx in range(0, len(rmet_bs)):
        wf.write("\t".join([str(rmet_bs[ridx]), str(rmet_nano[ridx])]) + "\n")
    wf.close()


def main():
    # python performance_corr_heatmap.data.py --modsfile ~/tools/data/arab/athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.tsv --bsfile ~/tools/data/arab/bs.result/D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile ~/tools/data/arab/bs.result/D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --bsfile ~/tools/data/arab/bs.result/D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt --covcf 5 --contig_prefix NC_003 --wfile athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.vs_bsrep123.rmet.tsv > athaliana.guppy.pass.part2.CG.bn13_sn16.arabnrice2-1.balance.both_bilstm.50x_12345.freq.vs_bsrep123.rmet.log 2>&1 &
    parser = argparse.ArgumentParser()
    parser.add_argument("--modsfile", action="store", type=str, required=True,
                        help="50x call_mods_freq file")
    parser.add_argument("--bsfile", type=str, action="append", required=True,
                        help="bs rmet file")
    parser.add_argument("--covcf", type=int, required=False, default=5,
                        help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--wfile", type=str, required=True)

    args = parser.parse_args()
    cmp_sitesrmet_of_nano_and_bs(args)


if __name__ == '__main__':
    main()
