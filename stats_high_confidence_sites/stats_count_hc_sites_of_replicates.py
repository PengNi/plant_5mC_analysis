#! /usr/bin/python
import argparse
import os

sep = "\t"


def str2bool(v):
    # susendberg's function
    return v.lower() in ("yes", "true", "t", "1")


def read_one_bs_file(bsfile, contig_prefix, contig_names, cov_cf=5, rmet_l=0.0, rmet_g=0.9):
    sites_p, sites_n = set(), set()
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
            if rmet <= rmet_l:
                sites_n.add(sep.join([chrom, str(pos), strand]))
            elif rmet >= rmet_g:
                sites_p.add(sep.join([chrom, str(pos), strand]))
    return sites_p, sites_n


def write_pos2file(wfile, poses):
    wf = open(wfile, "w")
    for pos in poses:
        wf.write(pos + "\n")
    wf.close()


def stats_hc_sites(args):
    bs_files = args.bs_file
    site_ps, site_ns = [], []
    wf = open(args.wfile, "w")

    wf.write("\t".join(["rep", "positive", "negative"]) + "\n")
    for bs_file in bs_files:
        sites_p, sites_n = read_one_bs_file(bs_file, args.contig_prefix, args.contig_names,
                                            args.cov_cf, args.rmet_l, args.rmet_g)
        site_ps.append(sites_p)
        site_ns.append(sites_n)
        wf.write("\t".join([bs_file, str(len(sites_p)), str(len(sites_n))]) + "\n")
        if str2bool(args.output_single):
            fname, fext = os.path.splitext(bs_file)
            wf_pos = fname + ".hc_sites.pos.txt"
            wf_neg = fname + ".hc_sites.neg.txt"
            write_pos2file(wf_pos, sites_p)
            write_pos2file(wf_neg, sites_n)
    posstr_union_p = set.union(*site_ps)
    posstr_inter_p = set.intersection(*site_ps)
    posstr_union_n = set.union(*site_ns)
    posstr_inter_n = set.intersection(*site_ns)
    wf.write("\t".join(["intersection", str(len(posstr_inter_p)), str(len(posstr_inter_n))]) + "\n")
    wf.write("\t".join(["union", str(len(posstr_union_p)), str(len(posstr_union_n))]) + "\n")
    wf.close()

    fname, fext = os.path.splitext(args.wfile)
    wf_pos_i = fname + ".pos.inter.txt"
    wf_pos_u = fname + ".pos.union.txt"
    wf_neg_i = fname + ".neg.inter.txt"
    wf_neg_u = fname + ".neg.union.txt"
    write_pos2file(wf_pos_i, posstr_inter_p)
    write_pos2file(wf_pos_u, posstr_union_p)
    write_pos2file(wf_neg_i, posstr_inter_n)


def main():
    # python stats_count_hc_sites_of_replicates.py --bs_file ~/data/nextomics/athaliana/bismark.out/
    # ninanjie-2-1_D1902826A-ZJ/
    # D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/data/nextomics/athaliana/bismark.out/ninanjie-2-2_D1902827A-ZJ/
    # D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/data/nextomics/athaliana/bismark.out/ninanjie-2-3_D1902828A-ZJ/
    # D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --cov_cf 5 --contig_prefix NC_003 --wfile ninanjie-2.bs_3replicates.CG.hc_sites.main_genome.stat.txt &

    # python stats_count_hc_sites_of_replicates.py --bs_file ~/tools/data/rice/bs.result/D1904161A-QJ_L2L4.
    # cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/tools/data/rice/bs.result/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --contig_names 1,2,3,4,5,6,7,8,9,10,11,12 --wfile shuidao.rep1-2nrep2-1.CG.hc_sites.main_genome.cmp.txt

    parser = argparse.ArgumentParser()
    parser.add_argument("--bs_file", type=str, action="append",
                        required=True, help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")
    parser.add_argument("--rmet_l", type=float, required=False, default=0.0, help="")
    parser.add_argument("--rmet_g", type=float, required=False, default=0.9, help="")
    parser.add_argument("--output_single", type=str, required=False, default="no",
                        help="output hc sites of single bs_file")
    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    stats_hc_sites(args)


if __name__ == '__main__':
    main()
