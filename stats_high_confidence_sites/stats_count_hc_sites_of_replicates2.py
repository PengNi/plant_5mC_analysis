#! /usr/bin/python
import argparse
import os

sep = "\t"


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def read_one_bs_file(bsfile, contig_prefix, contig_names, cov_cf=5):
    count_e0, count_e100, count_g90 = 0, 0, 0
    contigset = set(contig_names.split(",")) if contig_names is not None else None

    with open(bsfile, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage, rmet = words[0], words[1], words[2], int(words[6]), float(words[7])
            if coverage < cov_cf:
                continue
            if contig_prefix is not None and (not chrom.startswith(contig_prefix)):
                continue
            if contig_names is not None and chrom not in contigset:
                continue
            if rmet == 0.0:
                count_e0 += 1
            else:
                if rmet == 1.0:
                    count_e100 += 1
                if rmet >= 0.9:
                    count_g90 += 1
    return count_e0, count_e100, count_g90


def stats_hc_sites(args):
    bs_files = args.bs_file
    motifname = args.motif
    wf = open(args.wfile, "w")

    wf.write("\t".join(["motif", "replicate", "range", "count"]) + "\n")
    rep_cnt = 0
    for bs_file in bs_files:
        rep_cnt += 1
        count_e0, count_e100, count_g90 = read_one_bs_file(bs_file, args.contig_prefix,
                                                           args.contig_names, args.cov_cf)
        wf.write("\t".join([motifname, "rep" + str(rep_cnt), "=0.0", str(count_e0)]) + "\n")
        wf.write("\t".join([motifname, "rep" + str(rep_cnt), "=1.0", str(count_e100)]) + "\n")
        wf.write("\t".join([motifname, "rep" + str(rep_cnt), ">=0.9", str(count_g90)]) + "\n")
    wf.close()


def main():
    # python stats_count_hc_sites_of_replicates.py --bs_file ~/data/nextomics/athaliana/bismark.out/
    # ninanjie-2-1_D1902826A-ZJ/
    # D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/data/nextomics/athaliana/bismark.out/ninanjie-2-2_D1902827A-ZJ/
    # D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/data/nextomics/athaliana/bismark.out/ninanjie-2-3_D1902828A-ZJ/
    # D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --cov_cf 5 --contig_prefix NC_003 --wfile ninanjie-2.bs_3replicates.CG.hc_sites.main_genome.stat.txt &

    parser = argparse.ArgumentParser("count number of sites of which rmet =0.0/=1.0/>=0.9")
    parser.add_argument("--bs_file", type=str, action="append",
                        required=True, help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")
    parser.add_argument("--motif", type=str, required=False, default="CG")
    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    stats_hc_sites(args)


if __name__ == '__main__':
    main()
