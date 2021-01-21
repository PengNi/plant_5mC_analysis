#! /usr/bin/python
import argparse

alphabeta = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def _get_one_bs_file_positions(bs_file, cov_cf):
    chrom2poses = {}
    with open(bs_file, "r") as rf:
        next(rf)
        for line in rf:
            words = line.strip().split("\t")
            chrom, pos, strand, coverage = words[0], words[1], words[2], int(words[6])
            if coverage >= cov_cf:
                if chrom not in chrom2poses.keys():
                    chrom2poses[chrom] = set()
                chrom2poses[chrom].add("\t".join([chrom, pos, strand]))
    return chrom2poses


def union_poses_of_bsfiles(bs_files, cov_cf):
    chrom2poses = {}
    for bs_file in bs_files:
        tmp_chrom2poses = _get_one_bs_file_positions(bs_file, cov_cf)
        for chrom in tmp_chrom2poses.keys():
            if chrom not in chrom2poses.keys():
                chrom2poses[chrom] = set()
            chrom2poses[chrom].update(tmp_chrom2poses[chrom])
    return chrom2poses


def stat_poses(chrom2poses, contig_prefix, contig_names):
    poses_genome = set()
    totalcount = 0

    contigset = set(contig_names.strip().split(",")) if contig_names is not None else None

    for chrom in chrom2poses.keys():
        poses = chrom2poses[chrom]
        if contig_prefix is not None:
            if str(chrom).startswith(contig_prefix):
                poses_genome.update(poses)
        elif contig_names is not None:
            # contigset = contig_names
            if chrom in contigset:
                poses_genome.update(poses)
        else:
            poses_genome.update(poses)
        print("{}\t{}".format(len(poses), chrom))
        totalcount += len(poses)
    print("{}\tgenome(main)".format(len(poses_genome)))
    print("{}\tgenome(all)".format(totalcount))
    return poses_genome


def write_union_poses(poses, wfile):
    poses = sorted(list(poses))
    wf = open(wfile, "w")
    for mpos in poses:
        wf.write(mpos + "\n")
    wf.close()


def main():
    # python stats_count_c_of_bisulfite.py --bs_file ~/tools/data/rice/bs.result/
    # D1904161A-QJ_L2L4.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/tools/data/rice/bs.result/shuidao2-1_1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --cov_cf 5 --contig_names 1,2,3,4,5,6,7,8,9,10,11,12 --wfile shuidao.bs_rep2-1n1-2.CG.main_genome.stat.txt
    # > shuidao.bs_rep2-1n1-2.CG.main_genome.stat.log &

    # python stats_count_c_of_bisulfite.py --bs_file ~/data/nextomics/athaliana/bismark.out/
    # ninanjie-2-1_D1902826A-ZJ/
    # D1902826A-ZJ_301301.L1L3_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/data/nextomics/athaliana/bismark.out/ninanjie-2-2_D1902827A-ZJ/
    # D1902827A-ZJ_364364.L1L2_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --bs_file ~/data/nextomics/athaliana/bismark.out/ninanjie-2-3_D1902828A-ZJ/
    # D1902828A-ZJ_363363.L2L4_merged.cutadapt.R1_bismark_bt2_pe.sorted.mark_dup.sorted.CX_report.CG.txt
    # --cov_cf 5 --contig_prefix NC_003 --wfile ninanjie-2.bs_3replicates.CG.main_genome.stat.txt
    # >ninanjie-2.bs_3replicates.CG.main_genome.stat.log &

    parser = argparse.ArgumentParser("count motif num of each contig and genome in bisulfite sequencing")
    parser.add_argument("--bs_file", type=str, action="append", required=True, help="")
    parser.add_argument("--cov_cf", type=int, required=False, default=5, help="")
    parser.add_argument("--contig_prefix", type=str, required=False, default=None)  # NC_003,
    parser.add_argument("--contig_names", type=str, required=False, default=None)  # 1,2,3,4,5,6,7,8,9,10,11,12
    parser.add_argument("--wfile", type=str, required=True, default=None)

    args = parser.parse_args()
    print("=============coverage cutoff: 1")
    chrom2poses = union_poses_of_bsfiles(args.bs_file, 1)
    poses = stat_poses(chrom2poses, args.contig_prefix, args.contig_names)
    print("=============coverage cutoff: {}".format(args.cov_cf))
    chrom2poses = union_poses_of_bsfiles(args.bs_file, args.cov_cf)
    poses = stat_poses(chrom2poses, args.contig_prefix, args.contig_names)
    write_union_poses(poses, args.wfile)


if __name__ == '__main__':
    main()
